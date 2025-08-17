#include "dds/dds.h"
#include <stdio.h>
#include <iostream>
#include <memory>
#include <vector>
#include <stdlib.h>
#include <zephyr/kernel.h>
#include <pthread.h>
#include <zephyr/net/net_if.h>
#include <zephyr/net/ethernet.h>
#include <Eigen/Core>
// #include "TrajectorySample.hpp"
// #include "CartesianSample.hpp"
// #include "CurvilinearSample.hpp"
#include "InputData.h"
// #include "CartesianSampleData.h"
// #include "CurvilinearSampleData.h"
#include "ResultData.h"
// #include "TrajectoryEvaluator.hpp"
#include <dds/ddsi/ddsi_config.h>
#include "zephyr_app.hpp"

#include "TrajectoryHandler.hpp"
#include "CoordinateSystemWrapper.hpp"
#include "FillCoordinates.hpp"


#define ACTUATION_STACK_SIZE 16 * 1024
#define DDS_DOMAIN_ACTUATION 2
static K_THREAD_STACK_DEFINE(main_stack_area, ACTUATION_STACK_SIZE);

#if defined(CONFIG_NET_DHCPV4)
/* Semaphore to indicate a lease has been acquired. */
static K_SEM_DEFINE(got_address, 0, 1);

static struct net_mgmt_event_callback mgmt_cb;
#endif  // CONFIG_NET_DHCPV4
using InputData = custom_trajectory_msgs_msg_InputData;

using SamplingMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, 13, Eigen::RowMajor>;
using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Global params
static constexpr int num_samples_d = 9;
double sampling_d[num_samples_d + 1]; // +1 in case d is not present and needs to be added
static int actual_d_samples = 0;

static constexpr int num_samples_time = 7;
double sampling_time[num_samples_time + 1]; // +1 in case t is not present and needs to be added
static int actual_t_samples = 0;
static double N = 30; //  self.N = int(config_plan.planning.planning_horizon / config_plan.planning.dt)
static double dT = 0.1;

std::shared_ptr<CoordinateSystemWrapper> coordinate_system;

SamplingMatrixXd generate_sampling_matrix_cstyle(
    double* t0_range, size_t t0_size,
    double* t1_range, size_t t1_size,
    double* s0_range, size_t s0_size,
    double* ss0_range, size_t ss0_size,
    double* sss0_range, size_t sss0_size,
    double* ss1_range, size_t ss1_size,
    double* sss1_range, size_t sss1_size,
    double* d0_range, size_t d0_size,
    double* dd0_range, size_t dd0_size,
    double* ddd0_range, size_t ddd0_size,
    double* d1_range, size_t d1_size,
    double* dd1_range, size_t dd1_size,
    double* ddd1_range, size_t ddd1_size,
    bool debug_mode = true
) {
    size_t total_combinations = t0_size * t1_size * s0_size * ss0_size * sss0_size *
                              ss1_size * sss1_size * d0_size * dd0_size * ddd0_size *
                              d1_size * dd1_size * ddd1_size;
    
    SamplingMatrixXd sampling_matrix(total_combinations, 13);
    
    size_t idx = 0;
    for (size_t i0 = 0; i0 < t0_size; ++i0)
    for (size_t i1 = 0; i1 < t1_size; ++i1)
    for (size_t i2 = 0; i2 < s0_size; ++i2)
    for (size_t i3 = 0; i3 < ss0_size; ++i3)
    for (size_t i4 = 0; i4 < sss0_size; ++i4)
    for (size_t i5 = 0; i5 < ss1_size; ++i5)
    for (size_t i6 = 0; i6 < sss1_size; ++i6)
    for (size_t i7 = 0; i7 < d0_size; ++i7)
    for (size_t i8 = 0; i8 < dd0_size; ++i8)
    for (size_t i9 = 0; i9 < ddd0_size; ++i9)
    for (size_t i10 = 0; i10 < d1_size; ++i10)
    for (size_t i11 = 0; i11 < dd1_size; ++i11)
    for (size_t i12 = 0; i12 < ddd1_size; ++i12) {
        sampling_matrix(idx, 0) = t0_range[i0];
        sampling_matrix(idx, 1) = t1_range[i1];
        sampling_matrix(idx, 2) = s0_range[i2];
        sampling_matrix(idx, 3) = ss0_range[i3];
        sampling_matrix(idx, 4) = sss0_range[i4];
        sampling_matrix(idx, 5) = ss1_range[i5];
        sampling_matrix(idx, 6) = sss1_range[i6];
        sampling_matrix(idx, 7) = d0_range[i7];
        sampling_matrix(idx, 8) = dd0_range[i8];
        sampling_matrix(idx, 9) = ddd0_range[i9];
        sampling_matrix(idx, 10) = d1_range[i10];
        sampling_matrix(idx, 11) = dd1_range[i11];
        sampling_matrix(idx, 12) = ddd1_range[i12];
        ++idx;
    }
    if (debug_mode) {
        std::cout << "<ReactivePlanner>: " << idx << " trajectories sampled" << std::endl;
    }
    
    return sampling_matrix;
}

static bool get_msg(dds_entity_t rd, void * sample)
{
  dds_sample_info_t info;
  dds_return_t rc = dds_take(rd, &sample, &info, 1, 1);
  if (rc < 0) {
    dds_log(DDS_LC_WARNING, __FILE__, __LINE__, DDS_FUNCTION, "can't take msg\n");
  }
  if (rc > 0 && info.valid_data) {
    return true;
  }
  return false;
}

template<typename T>
static void on_msg_dds(dds_entity_t rd, void * arg)
{
  if (!arg) {
    return;
  }
  auto process = *reinterpret_cast<std::function<void(const T &)> *>(arg);

  // Rely on CycloneDDS for the memory management of the message's sequence buffers.
  // A static message ensures the buffers allocated dynamically on the previous calls are not lost.
  static T msg;
  if (get_msg(rd, reinterpret_cast<void *>(&msg))) {
    process(msg);
  }
}

void create_reader(
  dds_entity_t m_participant,
  const dds_topic_descriptor_t * desc,
  const char * name, const dds_qos_t * qos, 
  dds_on_data_available_fn callback, void * arg)
{
  dds_entity_t topic = dds_create_topic(m_participant, desc, name, NULL, NULL);
  if (topic < 0) {
    printf("dds_create_topic (%s): %s\n", name, dds_strretcode(-topic));
    std::exit(EXIT_FAILURE);
  }

  dds_listener_t * listener = dds_create_listener(arg);
  dds_lset_data_available(listener, callback);

  dds_entity_t reader = dds_create_reader(m_participant, topic, qos, listener);
  if (reader < 0) {
    printf("dds_create_reader (%s): %s\n", name, dds_strretcode(-reader));
    std::exit(EXIT_FAILURE);
  }

  std::cout << "waiting for writer to be discovered" << std::endl;

  uint32_t status = 0;

  dds_return_t rc = dds_set_status_mask(reader, DDS_SUBSCRIPTION_MATCHED_STATUS);
  if (rc != DDS_RETCODE_OK)
    DDS_FATAL("dds_set_status_mask: %s\n", dds_strretcode(-rc));

  while(!(status & DDS_SUBSCRIPTION_MATCHED_STATUS))
  {
    rc = dds_get_status_changes (reader, &status);
    if (rc != DDS_RETCODE_OK)
      DDS_FATAL("dds_get_status_changes: %s\n", dds_strretcode(-rc));

    dds_sleepfor (DDS_MSECS (20));
  }

  dds_delete_listener(listener);
}

dds_entity_t create_writer(
  dds_entity_t m_participant,
  const dds_topic_descriptor_t * desc,
  const char * name, const dds_qos_t * qos)
{
  dds_entity_t topic = dds_create_topic(m_participant, desc, name, NULL, NULL);
  if (topic < 0) {
    printf("dds_create_topic (%s): %s\n", name, dds_strretcode(-topic));
    std::exit(EXIT_FAILURE);
  }

  dds_entity_t writer = dds_create_writer(m_participant, topic, qos, NULL);
  if (writer < 0) {
    printf("dds_create_writer (%s): %s\n", name, dds_strretcode(-writer));
    std::exit(EXIT_FAILURE);
  }

  return writer;
}

Eigen::VectorXd ddsSequenceToVector(const dds_sequence_double& sequence) {
  if (sequence._length == 0 || sequence._buffer == nullptr)
      return Eigen::VectorXd{};
  
  Eigen::VectorXd vec(sequence._length);
  for (size_t i = 0; i < sequence._length; ++i) {
      vec(i) = sequence._buffer[i];
  }
  return vec;
}

dds_entity_t result_writer;

void on_msg(const InputData& msg) {
  static double vehicle_a_max = 11.5;
  static double horizon = 3.0;
  static double vehicle_v_max = 50.8;
  static double v_limit = 36.0;

  double s = msg.s;                   
  double ss = msg.ss;                   
  double sss = msg.sss;               
  double d = msg.d;            
  double dd = msg.dd;        
  double ddd = msg.ddd;

  double velocity = msg.velocity;
  double timestep = msg.timestep;
  double desired_velocity = msg.desired_velocity;

  std::cout << "s: " << s << std::endl;

  // Generate 9 linearly spaced samples between min_v and max_v
  double min_v = std::max(0.001, velocity - vehicle_a_max * horizon);
  double max_v = std::min({velocity + (vehicle_a_max / 6.0) * horizon, v_limit, vehicle_v_max});
  static constexpr int num_samples_v = 9;
  double sampling_v[num_samples_v + 1]; // +1 in case ss is not present and needs to be added
  int actual_v_samples = 0;

  double step_v = (max_v - min_v) / (num_samples_v - 1);
  sampling_v[0] = ss; ++actual_v_samples;
  for (int i = 1; i < num_samples_v+1; ++i) {
    if(sampling_v[i] != ss) {
      ++actual_v_samples;
      sampling_v[i] = min_v + i * step_v;
    }
  }
  
  // Add d to sampling_d if not present
  bool add_d = true;
  for (int i = 0; i < num_samples_d; ++i) {
    if(sampling_d[i] == d) {
      add_d = false;
      break;  
    }
  }
  if (add_d) {
    actual_d_samples++;
    sampling_d[actual_d_samples] = d;
  }

  // Add N*dT to sampling_time if not present
  bool add_time = true;
  for (int i = 0; i < num_samples_time; ++i) {
    if(sampling_time[i] == N*dT) {
      add_time = false;
      break;
    }
  }
  if (add_time) {
    actual_t_samples++;
    sampling_time[actual_t_samples] = N*dT;
  }

  size_t t0_size = 1;
  size_t t1_size = actual_t_samples;
  size_t s0_size = 1;
  size_t ss0_size = 1;
  size_t sss0_size = 1;
  size_t ss1_size = actual_v_samples;
  size_t sss1_size = 1;
  size_t d0_size = 1;
  size_t dd0_size = 1;
  size_t ddd0_size = 1;
  size_t d1_size = actual_d_samples;
  size_t dd1_size = 1;
  size_t ddd1_size = 1;

  double t0_range[1] = {0.0};
  double s0_range[1] = {s};
  double ss0_range[1] = {ss};
  double sss0_range[1] = {sss};
  double sss1_range[1] = {0.0};
  double d0_range[1] = {d};
  double dd0_range[1] = {dd};
  double ddd0_range[1] = {ddd};
  double dd1_range[1] = {0.0};
  double ddd1_range[1] = {0.0};

  SamplingMatrixXd sampling_matrix = generate_sampling_matrix_cstyle(
      t0_range, t0_size,
      sampling_time, t1_size,
      s0_range, s0_size,
      ss0_range, ss0_size,
      sss0_range, sss0_size,
      sampling_v, ss1_size,
      sss1_range, sss1_size,
      d0_range, d0_size,
      dd0_range, dd0_size,
      ddd0_range, ddd0_size,
      sampling_d, d1_size,
      dd1_range, dd1_size,
      ddd1_range, ddd1_size
  );

  TrajectoryHandler trajectory_handler(0.1, desired_velocity);

  trajectory_handler.addFillCoordinates(std::make_unique<FillCoordinates>(false, orientation, coordinate_system, 3.0));

  trajectory_handler.generateTrajectories(sampling_matrix, false);
    
  // trajectory_handler.evaluateTrajectory(trajectory_handler.m_trajectories[0]);

  trajectory_handler.evaluateAllTrajectories();

  for (const auto& trajectory : trajectory_handler.m_trajectories) {
    printf("ID: %d, feasibility: %d, cost: %f\n", trajectory.m_uniqueId, trajectory.m_feasible, trajectory.m_cost);
    // std::cout << "Feasibility map:" << std::endl;
    // for (const auto& pair : trajectory.m_feasabilityMap) {
    //     std::cout << pair.first << ": " << pair.second << std::endl;
    // }
  }
}

static void * main_thread(void * arg)
{
    std::cout << "entered main_thread" << std::endl;
    (void)arg;
    dds_entity_t participant;
    dds_entity_t topic;
    dds_entity_t reader;
    dds_return_t rc;
    dds_qos_t *qos;

    struct ddsi_config dds_cfg;
    init_config(dds_cfg);

    dds_entity_t domain = dds_create_domain_with_rawconfig(DDS_DOMAIN_ACTUATION, &dds_cfg);
    // The domain could have been set prior to this point, don't fail in this case.
    if (domain < 0 && domain != DDS_RETCODE_PRECONDITION_NOT_MET) {
      printf("dds_create_domain_with_rawconfig: %s\n", dds_strretcode(-domain));
      std::exit(EXIT_FAILURE);
    }

    participant = dds_create_participant (DDS_DOMAIN_ACTUATION, NULL, NULL);
    if (participant < 0)
      DDS_FATAL("dds_create_participant: %s\n", dds_strretcode(-participant));

    qos = dds_create_qos ();
    dds_qset_reliability (qos, DDS_RELIABILITY_RELIABLE, DDS_MSECS(30));  

    std::function<void(const InputData&)> trajectory_callback = [&](const InputData & msg) {on_msg(msg);};

    std::cout << "creating reader" << std::endl;
    create_reader(
      participant,
      &custom_trajectory_msgs_msg_InputData_desc,
      "trajectory_sample_msg",
      qos,
      on_msg_dds<InputData>,
      reinterpret_cast<void*>(&trajectory_callback)
    );
    
    std::cout << "creating writer" << std::endl;
    result_writer = create_writer(
      participant,
      &custom_trajectory_msgs_msg_ResultData_desc,
      "result_data_msg",
      qos
    );
    dds_delete_qos(qos);
    
    
    return NULL;
}

static void net_if_cb(struct net_if * iface, void * user_data)
{
  auto ifs = reinterpret_cast<std::vector<struct net_if *> *>(user_data);
  ifs->push_back(iface);
}

static int setup_iface(
  struct net_if * iface, const char * addr, const char * gw, const char * netmask, uint16_t tag)
{
  struct in_addr inaddr;

  if (net_addr_pton(AF_INET, addr, &inaddr)) {
    std::printf("Invalid address: %s", addr);
    return 1;
  }
  if (!net_if_ipv4_addr_add(iface, &inaddr, NET_ADDR_MANUAL, 0)) {
    std::printf("Cannot add %s to interface %p", addr, iface);
    return 1;
  }

  if (net_addr_pton(AF_INET, gw, &inaddr)) {
    std::printf("Invalid address: %s", gw);
    return 1;
  }
  net_if_ipv4_set_gw(iface, &inaddr);

  if (net_addr_pton(AF_INET, netmask, &inaddr)) {
    std::printf("Invalid address: %s", netmask);
    return 1;
  }
  net_if_ipv4_set_netmask(iface, &inaddr);

#if defined(CONFIG_NET_VLAN)
  if (tag > 0) {
    int ret = net_eth_vlan_enable(iface, tag);
    if (ret < 0) {
      std::printf("Cannot set VLAN tag %d to interface %p", tag, iface);
      return 1;
    }
  }
#endif

  return 0;
}

#if defined(CONFIG_NET_DHCPV4)
static void handler(
  struct net_mgmt_event_callback *cb, uint32_t mgmt_event, struct net_if *iface)
{
  int i = 0;

  if (mgmt_event != NET_EVENT_IPV4_ADDR_ADD) {
    return;
  }

  for (i = 0; i < NET_IF_MAX_IPV4_ADDR; i++) {
    char buf[NET_IPV4_ADDR_LEN];

    if (iface->config.ip.ipv4->unicast[i].addr_type != NET_ADDR_DHCP) {
      continue;
    }

    printf("  IP address: %s\n",
      net_addr_ntop(AF_INET,
          &iface->config.ip.ipv4->unicast[i].address.in_addr,
              buf, sizeof(buf)));
    printf("  Lease time: %u seconds\n",
       iface->config.dhcpv4.lease_time);
    printf("  Netmask:    %s\n",
      net_addr_ntop(AF_INET,
               &iface->config.ip.ipv4->netmask,
               buf, sizeof(buf)));
    printf("  Gateway:    %s\n",
      net_addr_ntop(AF_INET,
             &iface->config.ip.ipv4->gw,
             buf, sizeof(buf)));

    k_sem_give(&got_address);
  }
}
#endif  // CONFIG_NET_DHCPV4

int main(void)
{
  // Initialize network interfaces
  std::vector<struct net_if *> ifs {};
  net_if_foreach(net_if_cb, &ifs);
  if (ifs.size() >= 1 && sizeof(CONFIG_NET_IFACE1_ADDR) > 1) {
    int ret = setup_iface(
      ifs[0],
      CONFIG_NET_IFACE1_ADDR,
      CONFIG_NET_IFACE1_GW,
      CONFIG_NET_IFACE1_NETMASK,
      CONFIG_NET_IFACE1_VLAN);
    if (ret) {
      return 1;
    }
  }
  if (ifs.size() >= 2 && sizeof(CONFIG_NET_IFACE2_ADDR) > 1) {
    int ret = setup_iface(
      ifs[1],
      CONFIG_NET_IFACE2_ADDR,
      CONFIG_NET_IFACE2_GW,
      CONFIG_NET_IFACE2_NETMASK,
      CONFIG_NET_IFACE2_VLAN);
    if (ret) {
      return 1;
    }
  }

#if defined(CONFIG_NET_DHCPV4)
  if(ifs.size() >= 1) {
    std::printf("Requesting a DHCP lease...\n");
    net_mgmt_init_event_callback(&mgmt_cb, handler, NET_EVENT_IPV4_ADDR_ADD);
    net_mgmt_add_event_callback(&mgmt_cb);
    net_dhcpv4_start(ifs[0]);

    /* Wait for a lease. */
    if (k_sem_take(&got_address, K_SECONDS(10)) != 0) {
      std::printf("Did not get a DHCP lease\n");
    }
  }
#endif  // CONFIG_NET_DHCPV4

  pthread_attr_t main_attr {};
  pthread_t pthread_main {};
  void * retval = NULL;
  size_t main_stacksize = K_THREAD_STACK_SIZEOF(main_stack_area);

  RowMatrixXd path(409, 2);  // 409 points (x,y)
  path << -157.98592136094146, -49.45232504653796,
          -157.08161898224327, -49.02543254498661,
          -156.17731660354505, -48.598540043435264,
          -155.27301422484683, -48.17164754188394,
          -154.3687118461486, -47.74475504033259,
          -153.46440946745042, -47.31786253878124,
          -152.56010708875223, -46.89097003722989,
          -151.65580471005399, -46.46407753567856,
          -150.7515023313558, -46.03718503412723,
          -149.8471999526576, -45.61029253257587,
          -148.94289757395939, -45.18340003102453,
          -148.03859519526117, -44.756507529473204,
          -147.13429281656295, -44.32961502792186,
          -146.22999043786476, -43.90272252637052,
          -145.32568805916657, -43.47583002481918,
          -144.42138568046835, -43.04893752326784,
          -143.51708330177013, -42.6220450217165,
          -142.61278092307194, -42.19515252016515,
          -141.70847854437372, -41.768260018613816,
          -140.8041761656755, -41.34136751706249,
          -139.89987378697728, -40.914475015511144,
          -138.99557140827906, -40.487582513959794,
          -138.09126902958084, -40.060690012408465,
          -137.18696665088268, -39.6337975108571,
          -136.28266427218443, -39.20690500930578,
          -135.37836189348624, -38.780012507754435,
          -134.47405951478802, -38.35312000620309,
          -133.56975713608983, -37.92622750465175,
          -132.6654547573917, -37.49933500310023,
          -131.76115237871178, -37.072442501510174,
          -130.8568500037566, -36.64554999202984,
          -129.9525484097171, -36.21865582831096,
          -129.04824676041957, -35.7917617816469,
          -128.14394304572292, -35.3648721102075,
          -127.23963615294599, -34.937989171132116,
          -126.33532917383376, -34.511106414948586,
          -125.43102649187777, -34.08421455582229,
          -124.52673249512797, -33.657304298989665,
          -123.6224487448657, -33.23037233838011,
          -122.71816366821758, -32.80344318721643,
          -121.81386164802375, -32.376549926558,
          -120.90952705614286, -31.949725670939117,
          -120.00515471219664, -31.522981411210722,
          -119.10078776212464, -31.096225721067164,
          -118.19648427929071, -30.66933556304334,
          -117.2923024459948, -30.242187811207774,
          -116.3881797637249, -29.814914878114692,
          -115.48349784689115, -29.388827786269676,
          -114.57747275807421, -28.96560547214828,
          -113.66933564398911, -28.54693616805126,
          -112.75856704579161, -28.134022941319895,
          -111.84568944434065, -27.725791420673303,
          -110.93152423086336, -27.320450163363766,
          -110.01687696317951, -26.916197316734262,
          -109.10241619362972, -26.511522818641332,
          -108.1882139720808, -26.10626453564237,
          -107.27415680146255, -25.7006791754061,
          -106.3601300449471, -25.295025273783455,
          -105.44603862063151, -24.889517117877617,
          -104.53187696940184, -24.48416730163432,
          -103.61766739221265, -24.078925585822105,
          -102.70343230559844, -23.673741422407556,
          -101.78919448438971, -23.268563429146365,
          -100.87497837514658, -22.863336449563207,
          -99.96080895524196, -22.45800415397181,
          -99.04671122257113, -22.052510218945425,
          -98.13262436521929, -21.646991782310085,
          -97.21809937693634, -21.242462714458508,
          -96.30256982853363, -20.840213106894026,
          -95.38547822992038, -20.441538901026643,
          -94.46645023599086, -20.04734916887184,
          -93.54587763005655, -19.6567788668968,
          -92.62437784163217, -19.26840031918989,
          -91.70255927593936, -18.880778691091706,
          -90.78082502551231, -18.492956612228014,
          -89.85867178405695, -18.10613235245799,
          -88.93531782117152, -17.722184529415326,
          -88.009997279318, -17.343003213177713,
          -87.08219881388084, -16.969926420589466,
          -86.15243689542504, -16.601767337263535,
          -85.22152141323667, -16.236533234254537,
          -84.29024499906723, -15.872219865859499,
          -83.35920344825296, -15.507306743745046,
          -82.42814109956453, -15.142446705658472,
          -81.49653459939827, -14.77897884648814,
          -80.5638650329459, -14.418248666418323,
          -79.62976695965554, -14.061234126433671,
          -78.69451333545538, -13.707256703026282,
          -77.75856853088955, -13.355110207175056,
          -76.82239093091336, -13.003582829187241,
          -75.88635985149534, -12.651665533987728,
          -74.95051449754203, -12.299254614937656,
          -74.0147861415161, -11.946533140661437,
          -73.07910523173764, -11.593685815543331,
          -72.14340793984452, -11.240881936879415,
          -71.20765617272849, -10.88822257325113,
          -70.27182005012028, -10.535787126309621,
          -69.33586981077568, -10.183654854317119,
          -68.39961853508584, -9.8323241854733,
          -67.46219018740717, -9.484149186596753,
          -66.52251655125747, -9.142085800777323,
          -65.57958495657432, -8.80911705103523,
          -64.63279950624583, -8.487270469463436,
          -63.683059621176284, -8.174238646804483,
          -62.731652700657605, -7.8663040799928385,
          -61.77980809648088, -7.559723272216695,
          -60.82837790564595, -7.251858911757546,
          -59.87663909611942, -6.944951129541261,
          -58.923380711877314, -6.642802424273961,
          -57.96744515468536, -6.349245337005392,
          -57.00809260113385, -6.067058008989304,
          -56.04601894191799, -5.794279598847353,
          -55.08228711613065, -5.527410520093862,
          -54.117897171955015, -5.262926200548911,
          -53.15365629878775, -4.9978990064388595,
          -52.18964582001975, -4.732034827155767,
          -51.22571479406539, -4.4658826594857475,
          -50.2617103492301, -4.199996593735598,
          -49.29730592179328, -3.935566737141908,
          -48.331433858945914, -3.676559341746219,
          -47.36286293411753, -3.42785530884073,
          -46.390514550333194, -3.1943887538660753,
          -45.413788772511666, -2.9799851456827984,
          -44.433307132605705, -2.783444543214273,
          -43.44992970131604, -2.601922922961515,
          -42.46438399724467, -2.4325474908606464,
          -41.4772389361513, -2.272744688931793,
          -40.48878125561711, -2.121268203093463,
          -39.49920242435043, -1.9772935702385084,
          -38.5086748413572, -1.8399939816618744,
          -37.51734658550135, -1.7085970657749012,
          -36.52532008361497, -1.5825781710030185,
          -35.532679221710644, -1.4614918346886379,
          -34.53950111284894, -1.3448924121842367,
          -33.54584032137447, -1.2324805323404706,
          -32.551678930755685, -1.1245877508918554,
          -31.55698225966554, -1.0217494548118213,
          -30.561723621463972, -0.9245032338349484,
          -29.56592766825388, -0.8329178637269854,
          -28.569795928691974, -0.7450489953390833,
          -27.573565505498554, -0.658302673808,
          -26.577465054375743, -0.5700806566837289,
          -25.58164594604858, -0.478737317784068,
          -24.585890741565144, -0.3866970501685077,
          -23.58986045302053, -0.2976955932463605,
          -22.593249598652832, -0.21548037481744595,
          -21.595892993044153, -0.14288301231438674,
          -20.597951150778258, -0.07879986459171862,
          -19.59963265794245, -0.020854214220689813,
          -18.601102553937423, 0.033339388964252134,
          -17.60248551038177, 0.08591229633874237,
          -16.603838346320067, 0.1379107095471627,
          -15.605197520461456, 0.19003055625415255,
          -14.606599599947938, 0.24296503192573554,
          -13.60801810068384, 0.29620593123895916,
          -12.609171563889321, 0.34413889219548144,
          -11.60980956010062, 0.37948283807288885,
          -10.609956436032437, 0.3949152162872468,
          -9.610030770246793, 0.38542036500970234,
          -8.61048224030956, 0.35581970644419714,
          -7.611356500975159, 0.3140967945968245,
          -6.612409402031936, 0.26822177300084543,
          -5.613357873066582, 0.22469232338157635,
          -4.614192243777667, 0.18385512776272156,
          -3.614984038348573, 0.1440688350461996,
          -2.615800124608489, 0.10367951734589546,
          -1.61668562370448, 0.061608097363945985,
          -0.6175853486166965, 0.01919845434314802,
          0.38158883473077304, -0.021423284870488835,
          1.380913810898639, -0.05812591910446627,
          2.3804297512840047, -0.0891794895377925,
          3.3801063175903767, -0.11455237472644297,
          4.379900816449027, -0.1347644360144542,
          5.379778541437215, -0.15033983843369766,
          6.379722625825588, -0.1606592188286473,
          7.379712362538882, -0.16028877189537935,
          8.379526422179218, -0.14222458778522415,
          9.378570909954018, -0.09945302291769756,
          10.3755948667622, -0.02341763935640963,
          11.367649591925542, 0.10105005381777696,
          12.349180052307421, 0.2908639833319503,
          13.311489671269815, 0.5612730805616128,
          14.242883227304453, 0.9237646419860255,
          15.127934651094549, 1.3876500909537561,
          15.943699236277169, 1.9641860532280797,
          16.656691154495544, 2.663155169970869,
          17.232922635319504, 3.478575560417812,
          17.674374479975796, 4.374765547822999,
          18.0009693529536, 5.319294894153323,
          18.23333849363567, 6.291548651961431,
          18.387844195713782, 7.2792916675610195,
          18.47606429830874, 8.27522478811551,
          18.510694255797283, 9.274525183919827,
          18.505595089519815, 10.274463454201044,
          18.471549909927326, 11.273848241279998,
          18.408635765265, 12.271825414387665,
          18.314491908585804, 13.267335178683176,
          18.18672180398544, 14.259082394240892,
          18.024751395809655, 15.245829112045623,
          17.833642688100362, 16.227365825427974,
          17.619616975451052, 17.204174961849507,
          17.38879007156617, 18.177160595644068,
          17.14625619161778, 19.147298347394273,
          16.894225051482, 20.115013756877303,
          16.634294090764133, 21.080638474196775,
          16.36805284369085, 22.04454338946017,
          16.096583692229295, 23.006989122735845,
          15.819432582019287, 23.967813610100997,
          15.53582352204198, 24.926751231506255,
          15.244987619143402, 25.883521257977975,
          14.946750024308612, 26.838010612984068,
          14.642742269425069, 27.79067902570799,
          14.33496730718316, 28.742137821243347,
          14.025417979319657, 29.693021174417172,
          13.715659016836325, 30.643836290432247,
          13.405951166325623, 31.594668056382964,
          13.09628153765342, 32.545512271268755,
          12.786637236784571, 33.49636473453826,
          12.477005009205321, 34.447221129395935,
          12.167370493751449, 35.398076779234515,
          11.857719098012407, 36.34892693202381,
          11.548036301930052, 37.299766858441046,
          11.238322116166506, 38.2505965609198,
          10.928621377589591, 39.20143064337384,
          10.618988269979365, 40.152286751205274,
          10.309476085328647, 41.103182226852745,
          9.999948491265453, 42.05407268161099,
          9.689683002609794, 43.00472256233971,
          9.377835960927568, 43.95485467920166,
          9.063566693526496, 44.90418816488286,
          8.746405298194416, 45.85255939227609,
          8.427024759470417, 46.800185786457114,
          8.106332710680554, 47.747369264870066,
          7.785234568052715, 48.694415188440445,
          7.464416986708212, 49.641556185276045,
          7.143888722486576, 50.588795134665254,
          6.823518714676702, 51.53608762239247,
          6.503175992752016, 52.48338933797362,
          6.182759158049603, 53.43066598780968,
          5.862258856131586, 54.37791440079531,
          5.541684639466187, 55.32513780167456,
          5.221046318293755, 56.27233950495156,
          4.900426475425506, 57.2195474620469,
          4.580134690030284, 58.166866388462225,
          4.26052711574156, 59.11441635035737,
          3.941959652521081, 60.0623164702119,
          3.6246273258989214, 61.01063080634124,
          3.308221169970453, 61.95925459817455,
          2.9923291374399392, 62.908049736933314,
          2.6765393153603116, 63.85687890336948,
          2.360563123192616, 64.80564602086031,
          2.0444985241439, 65.75438369228812,
          1.7285219424943512, 66.7031506800233,
          1.4128094931633044, 67.65200558907924,
          1.0974031701303921, 68.60096230517125,
          0.7819236624713031, 69.54989469033617,
          0.46590655527603103, 70.49864815440041,
          0.14888703095408476, 71.44706709479104,
          -0.16996582339181987, 72.39487098137113,
          -0.49263448874775206, 73.34138184529017,
          -0.8213260435432994, 74.28581682906359,
          -1.1582322249434838, 75.22735104858587,
          -1.5044195998481826, 76.16551314964633,
          -1.8574701003676446, 77.10111609830822,
          -2.2142904650901207, 78.03528888108642,
          -2.5718036415486, 78.9691969370169,
          -2.927720846848168, 79.90371428714985,
          -3.28223090038702, 80.83876643744502,
          -3.6360177764207355, 81.77409249478201,
          -3.989765995934723, 82.70943317542049,
          -4.343974642732974, 83.64459958588111,
          -4.698555565065491, 84.57962491349379,
          -5.0533041495955, 85.51458664606778,
          -5.408015905930706, 86.44956235069992,
          -5.76277201845423, 87.38452121744838,
          -6.1185576946278015, 88.31908867358932,
          -6.476535329090846, 89.25281848082119,
          -6.837864066477433, 90.18525625912156,
          -7.203181946499499, 91.11613846519225,
          -7.571497835346944, 92.04583891300902,
          -7.941501718224973, 92.97486906864165,
          -8.31188660540387, 93.90374743923415,
          -8.68168932379698, 94.83285772001517,
          -9.051024112649596, 95.76215412341082,
          -9.42022021695592, 96.69150563731311,
          -9.789606951885428, 97.62078139229591,
          -10.159357770334397, 98.54991234436515,
          -10.52915846135274, 99.47902344999757,
          -10.898597133463719, 100.40827854874144,
          -11.267261710110352, 101.33784099602873,
          -11.634845256134017, 102.26783144021829,
          -12.001369783497214, 103.19823977912733,
          -12.366923635789284, 104.12902992882482,
          -12.731595334214294, 105.06016605644922,
          -13.095508291710162, 105.99159899546544,
          -13.458894219394791, 106.92323768755126,
          -13.822006852556116, 107.85498293895,
          -14.185099984885293, 108.78673579038431,
          -14.548359387854289, 109.7184238300705,
          -14.911760027210105, 110.65005679117908,
          -15.275234074728198, 111.58166111524743,
          -15.638713680644583, 112.51326327067929,
          -16.002148897859776, 113.44488274378851,
          -16.365545068304726, 114.37651744860341,
          -16.728918810378474, 115.30816090162155,
          -17.09228674372056, 116.23980662023877,
          -17.45565164861547, 117.17145352001164,
          -17.818973695345363, 118.1031171341648,
          -18.182204324561145, 119.03481639290489,
          -18.54529500896855, 119.9665701974834,
          -18.908250837090176, 120.89837654340512,
          -19.271241491072765, 121.83016932258879,
          -19.634470494309035, 122.76186921000487,
          -19.998141163815035, 123.69339678124784,
          -20.362255644405558, 124.62475097797497,
          -20.72620071708012, 125.55617137503722,
          -21.08923588751826, 126.4879467276345,
          -21.45062047736485, 127.42036336494104,
          -21.810360482184493, 128.3534159123437,
          -22.170744122189248, 129.28621987637945,
          -22.534530249754795, 130.21770119066215,
          -22.90446625506652, 131.14675545124373,
          -23.28244986161545, 132.07256428190007,
          -23.667845713091406, 132.99531305709726,
          -24.05950225528367, 133.9154227470873,
          -24.456277536147894, 134.83333736010047,
          -24.856999372416972, 135.7495364972136,
          -25.260383816578138, 136.66456685683204,
          -25.665126762135895, 137.5789973292936,
          -26.069926948809368, 138.49340247562867,
          -26.473851254077644, 139.4081948396289,
          -26.877077637945263, 140.32329508379732,
          -27.28001712000883, 141.2385217074257,
          -27.68308101568014, 142.1536935408974,
          -28.086518102838603, 143.0687009214038,
          -28.490087542206574, 143.9836499403258,
          -28.893445620021687, 144.89869215017487,
          -29.296248343116137, 145.81397894187631,
          -29.698140419957966, 146.72966591435423,
          -30.098732921150976, 147.6459221080406,
          -30.49762917312013, 148.5629179654272,
          -30.894431903713762, 149.48082157697357,
          -31.28899354536273, 150.3996908126717,
          -31.681918693501842, 151.31926117755017,
          -32.0739726303657, 152.23920339442074,
          -32.465921470061204, 153.15919040430938,
          -32.8583252463731, 154.07898345180067,
          -33.25112703230814, 154.99860660655784,
          -33.644138957507224, 155.91813998008556,
          -34.037173079311856, 156.83766386695206,
          -34.43009046951373, 157.75723763924165,
          -34.822898454674174, 158.6768581511383,
          -35.215635532067694, 159.5965089477677,
          -35.608340452731795, 160.51617347642124,
          -36.00116135058306, 161.4357884671726,
          -36.394572009306565, 162.35515127606268,
          -36.78911540882493, 163.2740284680703,
          -37.18533331157019, 164.19218475385955,
          -37.58351806068499, 165.10948979403517,
          -37.98322417005819, 166.02613304995148,
          -38.38384909035031, 166.94237517335594,
          -38.784791073006794, 167.85847862031036,
          -39.1855729012629, 168.77465214048345,
          -39.58608584094304, 169.6909432384851,
          -39.98630039361394, 170.60736470320978,
          -40.38618680902646, 171.5239293990578,
          -40.78554210096583, 172.44072560026547,
          -41.18365008291612, 173.35806399591718,
          -41.579682483897784, 174.2763001612699,
          -41.97280932796481, 175.1957837513835,
          -42.362622776975314, 176.11667703782862,
          -42.749971421245434, 177.03861012343583,
          -43.13597882248796, 177.96110576678757,
          -43.521770234281696, 178.88369179074823,
          -43.908158883180676, 179.80602782779997,
          -44.295036413769985, 180.72815892647904,
          -44.6820947897738, 181.65021413805135,
          -45.06903119295449, 182.5723205379255,
          -45.45582235748665, 183.49448787352082,
          -45.84262200623279, 184.4166516505853,
          -46.229421619216076, 185.33881544265046,
          -46.616221232258425, 186.26097923469084,
          -47.00302084530068, 187.18314302673124,
          -47.389820458342925, 188.10530681877168,
          -47.776620071385146, 189.0274706108121,
          -48.163419684427424, 189.9496344028525,
          -48.550219297469674, 190.87179819489293,
          -48.937018910511945, 191.79396198693334,
          -49.32381852355417, 192.71612577897378,
          -49.71061813659641, 193.6382895710142,
          -50.09741774963865, 194.5604533630546,
          -50.48421736268091, 195.48261715509503,
          -50.871016975723165, 196.40478094713546,
          -51.25781658876542, 197.32694473917587,
          -51.644616201807665, 198.24910853121628,
          -52.03141581484989, 199.17127232325674,
          -52.41821542789216, 200.09343611529712,
          -52.80501504093441, 201.01559990733756,
          -53.191814653976664, 201.93776369937797,
          -53.57861426701893, 202.85992749141838,
          -53.96541388006115, 203.78209128345884,
          -54.35221349310342, 204.70425507549922,
          -54.73901310614567, 205.62641886753966,
          -55.12581271918792, 206.54858265958006,
          -55.51261233223017, 207.47074645162047,
          -55.8994119452724, 208.39291024366094,
          -56.28621155831466, 209.31507403570134,
          -56.67301117135692, 210.23723782774175,
          -56.84360341775189, 210.64394452424278;
          
  coordinate_system = std::make_shared<CoordinateSystemWrapper>(path);
    
  // Generate 9 linearly spaced samples between d_min and d_max
  static double d_min = -3.0;
  static double d_max = 3.0;

  static double d_step = (d_max - d_min) / (num_samples_d - 1);
  for (int i = 0; i < num_samples_d; ++i) {
      sampling_d[i] = d_min + i * d_step;
  }
  actual_d_samples = num_samples_d;

  // Generate 7 linearly spaced samples between t_min and t_max
  static double t_min = 1.1;
  static double step_size = 3;
  static double t_max = 3.0;

  sampling_time[actual_t_samples] = t_min;
  ++actual_t_samples;
  static double t_val = t_min + step_size * dT;
  while(t_val < t_max) {
      sampling_time[actual_t_samples] = t_val;
      ++actual_t_samples;
      t_val += step_size * dT;
  }

  pthread_attr_init(&main_attr);
  pthread_attr_setstack(&main_attr, &main_stack_area, main_stacksize);
  pthread_create(&pthread_main, &main_attr, &main_thread, NULL);
  pthread_join(pthread_main, &retval);

  return 0;
}
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
#include "TrajectorySample.hpp"
#include "CartesianSample.hpp"
#include "CurvilinearSample.hpp"
#include "TrajectorySampleData.h"
#include "CartesianSampleData.h"
#include "CurvilinearSampleData.h"
#include "ResultData.h"
#include "TrajectoryEvaluator.hpp"
#include <dds/ddsi/ddsi_config.h>
#include "zephyr_app.hpp"


#define ACTUATION_STACK_SIZE 16 * 1024
#define DDS_DOMAIN_ACTUATION 2
static K_THREAD_STACK_DEFINE(main_stack_area, ACTUATION_STACK_SIZE);

#if defined(CONFIG_NET_DHCPV4)
/* Semaphore to indicate a lease has been acquired. */
static K_SEM_DEFINE(got_address, 0, 1);

static struct net_mgmt_event_callback mgmt_cb;
#endif  // CONFIG_NET_DHCPV4

using TrajectorySampleData = custom_trajectory_msgs_msg_TrajectorySampleData;
using XClData = custom_trajectory_msgs_msg_XClData;

// Global params
static constexpr int num_samples_d = 9;
double sampling_d[num_samples_d + 1]; // +1 in case d is not present and needs to be added
static int actual_d_samples = 0;

static constexpr int num_samples_time = 7;
double sampling_time[num_samples_time + 1]; // +1 in case t is not present and needs to be added
static int actual_t_samples = 0;
static double N = 30; //  self.N = int(config_plan.planning.planning_horizon / config_plan.planning.dt)
static double dT = 0.1;

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

void on_msg(const XClData& msg) {
  // size_t m_size = msg.size;
  // size_t m_actualSize = msg.m_actual_size;
  // double m_dT = msg.m_d_t;

  // printf("size: %zu\n", m_size);
  // printf("d: %d\n", m_dT);

  float s = msg.s;                   
  float ss = msg.ss;                   
  float sss = msg.sss;               
  float d = msg.d;            
  float dd = msg.dd;        
  float ddd = msg.ddd;

  float velocity = msg.velocity;
  float timestep = msg.timestep;
  float desired_velocity = msg.desired_velocity;

  std::cout << "s: " << s << std::endl;

  // Eigen::VectorXd m_curvilinearSample_s = ddsSequenceToVector(msg.m_curvilinear_sample.s);
  // Eigen::VectorXd m_curvilinearSample_d = ddsSequenceToVector(msg.m_curvilinear_sample.d);
  // Eigen::VectorXd m_curvilinearSample_theta = ddsSequenceToVector(msg.m_curvilinear_sample.theta);
  // Eigen::VectorXd m_curvilinearSample_dd = ddsSequenceToVector(msg.m_curvilinear_sample.dd);
  // Eigen::VectorXd m_curvilinearSample_ddd = ddsSequenceToVector(msg.m_curvilinear_sample.ddd);
  // Eigen::VectorXd m_curvilinearSample_ss = ddsSequenceToVector(msg.m_curvilinear_sample.ss);
  // Eigen::VectorXd m_curvilinearSample_sss = ddsSequenceToVector(msg.m_curvilinear_sample.sss);

  // Eigen::VectorXd m_cartesianSample_x = ddsSequenceToVector(msg.m_cartesian_sample.x);
  // std::cout << "x: " << m_cartesianSample_x << std::endl;
  // Eigen::VectorXd m_cartesianSample_y = ddsSequenceToVector(msg.m_cartesian_sample.y);
  // Eigen::VectorXd m_cartesianSample_theta = ddsSequenceToVector(msg.m_cartesian_sample.theta);
  // Eigen::VectorXd m_cartesianSample_velocity = ddsSequenceToVector(msg.m_cartesian_sample.velocity);
  // Eigen::VectorXd m_cartesianSample_acceleration = ddsSequenceToVector(msg.m_cartesian_sample.acceleration);
  // Eigen::VectorXd m_cartesianSample_kappa = ddsSequenceToVector(msg.m_cartesian_sample.kappa);
  // Eigen::VectorXd m_cartesianSample_kappaDot = ddsSequenceToVector(msg.m_cartesian_sample.kappa_dot);

  static float vehicle_a_max = 0;
  static float horizon = 0;
  static float vehicle_v_max = 0;
  static float v_limit = 0;
  static float max_density = 0;

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

  // CartesianSample cartesianSample(
  //     m_cartesianSample_x,
  //     m_cartesianSample_y,
  //     m_cartesianSample_theta,
  //     m_cartesianSample_velocity,
  //     m_cartesianSample_acceleration,
  //     m_cartesianSample_kappa,
  //     m_cartesianSample_kappaDotä
  // );

  // CurviLinearSample curvilinearSample(
  //     m_curvilinearSample_s,
  //     m_curvilinearSample_d,
  //     m_curvilinearSample_theta,
  //     m_curvilinearSample_dd,
  //     m_curvilinearSample_ddd,
  //     m_curvilinearSample_ss,
  //     m_curvilinearSample_sss
  // );

  // TrajectorySample trajectory_sample(
  //     cartesianSample,
  //     curvilinearSample,
  //     m_size,
  //     m_actualSize,
  //     m_dT
  // );

  // TrajectoryEvaluator trajectory_evaluator;
  // trajectory_evaluator.evaluateTrajectory(trajectory_sample);
  // custom_trajectory_msgs_msg_ResultData result_msg;
  // result_msg.m_cost = trajectory_sample.m_cost;
  // result_msg.m_feasible = trajectory_sample.m_feasible;
  // dds_write(result_writer, &result_msg);

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

  size_t total_combinations = t0_size * t1_size * s0_size * ss0_size * sss0_size *
                            ss1_size * sss1_size * d0_size * dd0_size * ddd0_size *
                            d1_size * dd1_size * ddd1_size;

  double sampling_matrix[total_combinations][13] = {0}; // All elements initialized to 0

  generate_sampling_matrix_cstyle(
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
      ddd1_range, ddd1_size,
      sampling_matrix
  );

  // Print the sampling matrix
  std::cout << "Sampling Matrix:" << std::endl;
  for (size_t i = 0; i < total_combinations; ++i) {
      for (size_t j = 0; j < 13; ++j) {
          std::cout << sampling_matrix[i][j] << " ";
      }
      std::cout << std::endl;
  }
}

void generate_sampling_matrix_cstyle(
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
    double sampling_matrix[][13],
    bool debug_mode = true
) {
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
        sampling_matrix[idx][0] = t0_range[i0];
        sampling_matrix[idx][1] = t1_range[i1];
        sampling_matrix[idx][2] = s0_range[i2];
        sampling_matrix[idx][3] = ss0_range[i3];
        sampling_matrix[idx][4] = sss0_range[i4];
        sampling_matrix[idx][5] = ss1_range[i5];
        sampling_matrix[idx][6] = sss1_range[i6];
        sampling_matrix[idx][7] = d0_range[i7];
        sampling_matrix[idx][8] = dd0_range[i8];
        sampling_matrix[idx][9] = ddd0_range[i9];
        sampling_matrix[idx][10] = d1_range[i10];
        sampling_matrix[idx][11] = dd1_range[i11];
        sampling_matrix[idx][12] = ddd1_range[i12];
        ++idx;
    }
    if (debug_mode) {
        std::cout << "<ReactivePlanner>: " << idx << " trajectories sampled" << std::endl;
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

    std::function<void(const TrajectorySampleData&)> trajectory_callback = [&](const TrajectorySampleData & msg) {on_msg(msg);};

    std::cout << "creating reader" << std::endl;
    create_reader(
      participant,
      &custom_trajectory_msgs_msg_TrajectorySampleData_desc,
      "trajectory_sample_msg",
      qos,
      on_msg_dds<TrajectorySampleData>,
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
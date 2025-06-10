#pragma once

#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

// #include "trajectory/TrajectorySample.hpp"
#include "TrajectorySample.hpp"
#include "TrajectoryStrategy.hpp"
#include "polynomial.hpp"

#include "FillCoordinates.hpp"

using SamplingMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, 13, Eigen::RowMajor>;

class TrajectoryHandler
{
public:
    TrajectoryHandler(double dt, double desired_velocity);


    /**
     * @brief Generates a trajectory for every row in samplingMatrix
     * 
     * @param samplingMatrix The Matrix needs to be of the form: [ t0
     *                                                           , t1
     *                                                           , s0
     *                                                           , ss0
     *                                                           , sss0
     *                                                           , ss1
     *                                                           , sss1
     *                                                           , d0
     *                                                           , dd0
     *                                                           , ddd0
     *                                                           , d1
     *                                                           , dd1
     *                                                           , ddd1]
     * 
     * @param lowVelocityMode If the vehicle is currently in the lowVelocityMode
     * 
     * t0: starting time
     * t1: end time
     * s0/s1:       start/end longitudinal 
     * ss0/ss1:     start/end longitudinal velocity
     * sss0/sss1:   start/end longitudinal acceleration
     * d0/d1:       start/end lateral
     * dd0/dd1:     start/end lateral velocity
     * ddd0/ddd1:   start/end lateral acceleration
     */
    void generateTrajectories(const SamplingMatrixXd& samplingMatrix, bool lowVelocityMode);
    void evaluateTrajectory(TrajectorySample& trajectory);
    void evaluateAllTrajectories();
    void addFillCoordinates(std::unique_ptr<FillCoordinates> function);
    

    std::vector<TrajectorySample> m_trajectories;
private:

    double m_dt;

    std::vector<std::unique_ptr<TrajectoryStrategy>> functions_;

};



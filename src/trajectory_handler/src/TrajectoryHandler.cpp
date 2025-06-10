#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/QR>

#include "TrajectoryHandler.hpp"

#include "CheckAccelerationConstraint.hpp"
#include "CheckCurvatureConstraints.hpp"
#include "CheckCurvatureRateConstrains.hpp"
#include "CheckYawRateConstraint.hpp"

#include "CalculateAccelerationCost.hpp"
#include "CalculateDistanceToReferencePathCost.hpp"
#include "CalculateJerkCost.hpp"
#include "CalculateLateralAccelerationCost.hpp"
#include "CalculateLateralVelocityCost.hpp"
#include "CalculateLongitudinalAccelerationCost.hpp"
#include "CalculateLongitudinalVelocityCost.hpp"
#include "CalculateNegativeAccelerationCost.hpp"
#include "CalculateNegativeOrientationOffsetCost.hpp"
#include "CalculateNegativeVelocityOffsetCost.hpp"
#include "CalculateOrientationOffsetCost.hpp"
#include "CalculatePositiveAccelerationCost.hpp"
#include "CalculatePositiveOrientationOffsetCost.hpp"
#include "CalculatePositiveVelocityOffsetCost.hpp"
#include "CalculateVelocityOffsetCost.hpp"

#include "FillCoordinates.hpp"

#include "TrajectorySample.hpp"
#include "TrajectoryStrategy.hpp"

#include <map>
#include <string>

std::map<std::string, double> cost_weights = {
    {"acceleration", 0.0},
    {"jerk", 0.0},
    {"lateral_jerk", 0.2},
    {"longitudinal_jerk", 0.2},
    {"orientation_offset", 0.0},
    {"path_length", 0.0},
    {"lane_center_offset", 0.0},
    {"velocity_offset", 1.0},
    {"velocity", 0.0},
    {"distance_to_reference_path", 5.0},
    {"responsibility", 0}
};

std::map<std::string, double> vehicle_params = {
    {"length", 4.508},
    {"width", 1.61},
    {"mass", 1093.3},
    {"wheelbase", 2.5789},
    {"delta_max", 1.066},
    {"delta_min", -1.066},
    {"a_max", 11.5},
    {"v_max", 50.8},
    {"v_switch", 7.319},
    {"v_delta_max", 0.4},
    {"v_delta_min", -0.4},
    {"wb_front_axle", 1.156},
    {"wb_rear_axle", 1.422},
    {"cr_vehicle_id", 2} 
};

TrajectoryHandler::TrajectoryHandler(double dt, double desired_velocity)
    : m_dt(dt)
{
    functions_.reserve(9); // 1 for fill coordinates
    
    functions_.emplace_back(std::make_unique<CheckAccelerationConstraint>(vehicle_params["v_switch"], vehicle_params["a_max"], false));
    functions_.emplace_back(std::make_unique<CheckCurvatureConstraint>(vehicle_params["delta_max"], vehicle_params["wheelbase"], false));
    functions_.emplace_back(std::make_unique<CheckCurvatureRateConstraint>(vehicle_params["wheelbase"], vehicle_params["v_delta_max"], false));
    functions_.emplace_back(std::make_unique<CheckYawRateConstraint>(vehicle_params["delta_max"], vehicle_params["wheelbase"], false));

    // functions_.emplace_back(std::make_unique<CalculateAccelerationCost>("CalculateAccelerationCost", cost_weights["acceleration"]));
    functions_.emplace_back(std::make_unique<CalculateDistanceToReferencePathCost>("CalculateDistanceToReferencePathCost", cost_weights["distance_to_reference_path"]));
    functions_.emplace_back(std::make_unique<CalculateJerkCost>("CalculateJerkCost", cost_weights["jerk"]));
    // functions_.emplace_back(std::make_unique<CalculateLateralAccelerationCost>("CalculateLateralAccelerationCost", cost_weights["lateral_jerk"]));
    // functions_.emplace_back(std::make_unique<CalculateLateralVelocityCost>("CalculateLateralVelocityCost", cost_weights["velocity"]));
    // functions_.emplace_back(std::make_unique<CalculateLongitudinalAccelerationCost>("CalculateLongitudinalAccelerationCost", cost_weights["longitudinal_jerk"]));
    // functions_.emplace_back(std::make_unique<CalculateLongitudinalVelocityCost>("CalculateLongitudinalVelocityCost", cost_weights["velocity"]));
    // functions_.emplace_back(std::make_unique<CalculateNegativeAccelerationCost>("CalculateNegativeAccelerationCost", cost_weights["acceleration"]));
    // functions_.emplace_back(std::make_unique<CalculateNegativeOrientationOffsetCost>("CalculateNegativeOrientationOffsetCost", cost_weights["orientation_offset"]));
    // functions_.emplace_back(std::make_unique<CalculateNegativeVelocityOffsetCost>("CalculateNegativeVelocityOffsetCost", cost_weights["velocity_offset"], desired_velocity));
    functions_.emplace_back(std::make_unique<CalculateOrientationOffsetCost>("CalculateOrientationOffsetCost", cost_weights["orientation_offset"]));
    // functions_.emplace_back(std::make_unique<CalculatePositiveAccelerationCost>("CalculatePositiveAccelerationCost", cost_weights["acceleration"]));
    // functions_.emplace_back(std::make_unique<CalculatePositiveOrientationOffsetCost>("CalculatePositiveOrientationOffsetCost", cost_weights["orientation_offset"]));
    // functions_.emplace_back(std::make_unique<CalculatePositiveVelocityOffsetCost>("CalculatePositiveVelocityOffsetCost", cost_weights["velocity_offset"], desired_velocity));
    functions_.emplace_back(std::make_unique<CalculateVelocityOffsetCost>("CalculateVelocityOffsetCost", cost_weights["velocity_offset"], desired_velocity, 0.1, 1.1, false, 2));

// ('velocity_offset', <frenetix._frenetix.trajectory_functions.cost_functions.CalculateVelocityOffsetCost object at 0x7d7e703d8620>)
// ('distance_to_reference_path', <frenetix._frenetix.trajectory_functions.cost_functions.CalculateDistanceToReferencePathCost object at 0x7d7e7af80390>)
// ('longitudinal_jerk', <frenetix._frenetix.trajectory_functions.cost_functions.CalculateLongitudinalJerkCost object at 0x7d7e7af80180>)
// ('prediction', <frenetix._frenetix.trajectory_functions.cost_functions.CalculateCollisionProbabilityFast object at 0x7d7e7a588f30>)
// ('lateral_jerk', <frenetix._frenetix.trajectory_functions.cost_functions.CalculateLateralJerkCost object at 0x7d7e7af6af10>)
}

void TrajectoryHandler::generateTrajectories(const SamplingMatrixXd& samplingMatrix, bool lowVelocityMode)
{
    m_trajectories.clear();
    m_trajectories.reserve(samplingMatrix.rows());

    for(Eigen::Index iii = 0; iii < samplingMatrix.rows(); iii++)
    {
        Eigen::Vector3d x0_lon {samplingMatrix.row(iii)[2], samplingMatrix.row(iii)[3], samplingMatrix.row(iii)[4]};
        Eigen::Vector2d x1_lon {samplingMatrix.row(iii)[5], samplingMatrix.row(iii)[6]};

        TrajectorySample::LongitudinalTrajectory longitudinalTrajectory (
            samplingMatrix.row(iii)[0],
            samplingMatrix.row(iii)[1],
            x0_lon,
            x1_lon,
            TrajectorySample::LongitudinalX0Order,
            TrajectorySample::LongitudinalXDOrder
        );

        double t1 = 0.0;
        if (lowVelocityMode) {
            t1 = longitudinalTrajectory(samplingMatrix.row(iii)[1]) - x0_lon[0];
            if (t1 <= 0.0) {
                t1 = samplingMatrix.row(iii)[1];
            }
        } else {
            t1 = samplingMatrix.row(iii)[1];
        }

        Eigen::Vector3d x0_lat {samplingMatrix.row(iii)[7], samplingMatrix.row(iii)[8], samplingMatrix.row(iii)[9]};
        Eigen::Vector3d x1_lat {samplingMatrix.row(iii)[10], samplingMatrix.row(iii)[11], samplingMatrix.row(iii)[12]};

        TrajectorySample::LateralTrajectory lateralTrajectory(
            samplingMatrix.row(iii)[0],
            t1,
            x0_lat,
            x1_lat
            );

        m_trajectories.emplace_back(
            m_dt,
            longitudinalTrajectory,
            lateralTrajectory,
            iii,
            samplingMatrix.row(iii)
            );
    }
}

void TrajectoryHandler::addFillCoordinates(std::unique_ptr<FillCoordinates> function) {
    if (functions_.size() >= 9) {
        functions_.erase(functions_.begin()); // Remove the frist function if it exists
    }
    functions_.emplace(functions_.begin(), std::move(function));
}

void TrajectoryHandler::evaluateTrajectory(TrajectorySample& trajectory) {

    // Evaluate fillcoordinates_func first

    for (auto& func : functions_) {
        // std::string funName = func->getFunctionName();
        // std::cout << "Evaluating function: " << funName << std::endl;
        func->evaluateTrajectory(trajectory);
    }

}

void TrajectoryHandler::evaluateAllTrajectories() {
    for (auto& trajectory : m_trajectories) {
        evaluateTrajectory(trajectory);
    }
}
#include "TrajectoryEvaluator.hpp"

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
    {"wheelbase", 2.59},
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

double desired_velocity = 7.0;

TrajectoryEvaluator::TrajectoryEvaluator() {
    functions_.reserve(19);
    
    functions_.emplace_back(std::make_unique<CheckAccelerationConstraint>(vehicle_params["v_switch"], vehicle_params["a_max"], true));
    functions_.emplace_back(std::make_unique<CheckCurvatureConstraint>(vehicle_params["delta_max"], vehicle_params["wheelbase"], true));
    functions_.emplace_back(std::make_unique<CheckCurvatureRateConstraint>(vehicle_params["wheelbase"], vehicle_params["v_delta_max"], true));
    functions_.emplace_back(std::make_unique<CheckYawRateConstraint>(vehicle_params["delta_max"], vehicle_params["wheelbase"], true));

    functions_.emplace_back(std::make_unique<CalculateAccelerationCost>("CalculateAccelerationCost", cost_weights["acceleration"]));
    functions_.emplace_back(std::make_unique<CalculateDistanceToReferencePathCost>("CalculateDistanceToReferencePathCost", cost_weights["distance_to_reference_path"]));
    functions_.emplace_back(std::make_unique<CalculateJerkCost>("CalculateJerkCost", cost_weights["jerk"]));
    functions_.emplace_back(std::make_unique<CalculateLateralAccelerationCost>("CalculateLateralAccelerationCost", cost_weights["lateral_jerk"]));
    functions_.emplace_back(std::make_unique<CalculateLateralVelocityCost>("CalculateLateralVelocityCost", cost_weights["velocity"]));
    functions_.emplace_back(std::make_unique<CalculateLongitudinalAccelerationCost>("CalculateLongitudinalAccelerationCost", cost_weights["longitudinal_jerk"]));
    functions_.emplace_back(std::make_unique<CalculateLongitudinalVelocityCost>("CalculateLongitudinalVelocityCost", cost_weights["velocity"]));
    functions_.emplace_back(std::make_unique<CalculateNegativeAccelerationCost>("CalculateNegativeAccelerationCost", cost_weights["acceleration"]));
    functions_.emplace_back(std::make_unique<CalculateNegativeOrientationOffsetCost>("CalculateNegativeOrientationOffsetCost", cost_weights["orientation_offset"]));
    functions_.emplace_back(std::make_unique<CalculateNegativeVelocityOffsetCost>("CalculateNegativeVelocityOffsetCost", cost_weights["velocity_offset"], desired_velocity));
    functions_.emplace_back(std::make_unique<CalculateOrientationOffsetCost>("CalculateOrientationOffsetCost", cost_weights["orientation_offset"]));
    functions_.emplace_back(std::make_unique<CalculatePositiveAccelerationCost>("CalculatePositiveAccelerationCost", cost_weights["acceleration"]));
    functions_.emplace_back(std::make_unique<CalculatePositiveOrientationOffsetCost>("CalculatePositiveOrientationOffsetCost", cost_weights["orientation_offset"]));
    functions_.emplace_back(std::make_unique<CalculatePositiveVelocityOffsetCost>("CalculatePositiveVelocityOffsetCost", cost_weights["velocity_offset"], desired_velocity));
    functions_.emplace_back(std::make_unique<CalculateVelocityOffsetCost>("CalculateVelocityOffsetCost", cost_weights["velocity_offset"], desired_velocity, 0.1, 1.1, false, 2));
}


void TrajectoryEvaluator::evaluateTrajectory(TrajectorySample& trajectory) {
    for (auto& func : functions_) {
        func->evaluateTrajectory(trajectory);
    }

}

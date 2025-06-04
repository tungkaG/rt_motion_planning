#pragma once

#include "TrajectorySample.hpp"
#include "TrajectoryStrategy.hpp"
#include <memory>
#include <vector>

class TrajectoryEvaluator {
public:
    explicit TrajectoryEvaluator();
    
    void evaluateTrajectory(TrajectorySample& trajectory);

private:
    std::vector<std::unique_ptr<TrajectoryStrategy>> functions_;
};

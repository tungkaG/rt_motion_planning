#pragma once

#include <stddef.h>
#include <Eigen/Core>
#include <optional>
#include <string>
// #include <utility>
#include <unordered_map>
#include <memory>

#include "CartesianSample.hpp"
#include "CurvilinearSample.hpp"
#include "polynomial.hpp"


class TrajectorySample
{
public:
    using LongitudinalTrajectory = PolynomialTrajectory<4, 3, 2>;
    using LateralTrajectory = PolynomialTrajectory<5, 3, 3>;


    TrajectorySample(CartesianSample m_cartesianSample, CurviLinearSample m_curvilinearSample, size_t m_size, size_t m_actualSize, double m_dT);
    TrajectorySample(double dt,
                 LongitudinalTrajectory trajectoryLongitudinal,
                 LateralTrajectory trajectoryLateral,
                 int uniqueId,
                 Eigen::VectorXd samplingParameters);


    CartesianSample m_cartesianSample;
    CurviLinearSample m_curvilinearSample;

    static const LongitudinalTrajectory::OrderVectorX0 LongitudinalX0Order;
    static const LongitudinalTrajectory::OrderVectorXD LongitudinalXDOrder;

    std::shared_ptr<LinearTrajectory> m_trajectoryLongitudinal;
    std::shared_ptr<LinearTrajectory> m_trajectoryLateral;

    Eigen::Vector<double, 13> m_samplingParameters;


    std::unordered_map<std::string, std::pair<double,double>> m_costMap;
    std::unordered_map<std::string, double> m_feasabilityMap;

    size_t m_size;
    size_t m_actualSize;
    double m_dT; // = NaN
    double m_cost = 0.0;

    // These variables are only used by Python
    std::optional<double> m_harm_occ_module;

    std::optional<int> m_uniqueId;

    std::optional<int> m_currentTimeStep;

    bool m_valid = true;
    bool m_feasible = true;

    /**
     * @brief Initialize arrays of the curvilinear and cartesian samples with the specified size.
     *
     * @param size The size to resize the arrays to.
     */
    void initArraysWithSize(size_t size);

    /**
     * @brief Add a cost value to the list of cost values. The weight must therefore be speicified
     * in the costWeightMap. Add a costWeight to the map if it is not specified yet by
     * handler.addCostWeight(std::string functionName, double weight).
     *
     * @param costFunctionName
     * @param cost
     * @param costWeighted
     */
    void addCostValueToList(std::string costFunctionName, double cost, double costWeighted);

    void addFeasabilityValueToList(std::string costFunctionName, double value);


    size_t size();
};


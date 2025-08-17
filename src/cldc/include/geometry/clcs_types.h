#ifndef CLCS_TYPES_H
#define CLCS_TYPES_H

// #include <boost/geometry.hpp>
// #include <boost/geometry/geometries/point_xy.hpp>
// #include <boost/geometry/geometries/polygon.hpp>

#include "clipper2/clipper.h"

namespace geometry {

    // /// boost point type
    // typedef boost::geometry::model::d2::point_xy<double> point_type;

    // /// boost polygon type
    // typedef boost::geometry::model::polygon<point_type> polygon_type;

    // /// boost multi polygon type
    // typedef boost::geometry::model::multi_polygon<polygon_type> mpolygon_type;

    using point_type   = Clipper2Lib::PointD;  // 2D point with .x and .y (double)

    using polygon_type = Clipper2Lib::PathD;   // one closed ring = one polygon boundary

    using mpolygon_type = Clipper2Lib::PathsD;

    /// EigenPolyline
    typedef std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> EigenPolyline;

    /// EigenMatrix
    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

} // namespace geometry

#endif //CLCS_TYPES_H

#include "geometry/util_projection_domain.h"

namespace geometry {

namespace util_projection_domain {

// tiny epsilon check
static inline bool same_pt(const Eigen::Vector2d& a, const Eigen::Vector2d& b, double eps = 1e-12) {
    return std::abs(a.x() - b.x()) <= eps && std::abs(a.y() - b.y()) <= eps;
}

// void polylineToBoostPolygon(EigenPolyline& polyline, polygon_type& boost_polygon) {
//     boost_polygon.clear();
//     for (EigenPolyline::iterator it = polyline.begin(); it != polyline.end(); it++) {
//         boost::geometry::append(boost_polygon, point_type((*it)[0], (*it)[1]));
//     }
// }

void polylineToClipperPath(EigenPolyline& polyline, polygon_type& clipper_path) {
    clipper_path.clear();
    for (const auto& pt : polyline) {
        clipper_path.emplace_back(pt.x(), pt.y());
    }
}

// void polylineToBoostPolygon(const EigenPolyline& polyline, polygon_type& boost_polygon) {
//     boost_polygon.clear();
//     for (const auto &vert : polyline) {
//         double x = vert[0];
//         double y = vert[1];
//         boost::geometry::append(boost_polygon, point_type(vert[0], vert[1]));
//     }
// }
void polylineToClipperPath(const EigenPolyline& polyline, polygon_type& clipper_path) {
    clipper_path.clear();
    for (const auto& pt : polyline) {
        clipper_path.emplace_back(pt.x(), pt.y());
    }
}

void overapproximatePolygonAABB(const EigenPolyline& polyline, polygon_type& clipper_aabb) {
    // boost_poly_aabb.clear();
    clipper_aabb.clear();

    // init extremum points of polygon
    double x_min = std::numeric_limits<double>::infinity();
    double x_max = -std::numeric_limits<double>::infinity();
    double y_min = std::numeric_limits<double>::infinity();
    double y_max = -std::numeric_limits<double>::infinity();

    // iterate over vertices of polygon and store min/max
    for (const auto& vert: polyline) {
        // double x = vert[0];
        // double y = vert[1];
        double x = vert.x();
        double y = vert.y();
        x_min = x_min > x ? x : x_min;
        y_min = y_min > y ? y : y_min;
        x_max = x_max < x ? x : x_max;
        y_max = y_max < y ? y : y_max;
    }

    // // create polygon
    // boost::geometry::append(boost_poly_aabb, point_type(x_min, y_min));
    // boost::geometry::append(boost_poly_aabb, point_type(x_min, y_max));
    // boost::geometry::append(boost_poly_aabb, point_type(x_max, y_max));
    // boost::geometry::append(boost_poly_aabb, point_type(x_max, y_min));
    // boost::geometry::append(boost_poly_aabb, point_type(x_min, y_min));

    clipper_aabb.emplace_back(x_min, y_min);
    clipper_aabb.emplace_back(x_min, y_max);
    clipper_aabb.emplace_back(x_max, y_max);
    clipper_aabb.emplace_back(x_max, y_min);
    clipper_aabb.emplace_back(x_min, y_min);
}

// std::vector<EigenPolyline> polygonWithinPolygonBoost(const polygon_type& polygon_input,
//                                                      const polygon_type& polygon_other) {
//     // init output: vector of EigenPolylines
//     std::vector<EigenPolyline> polygons_within_polygon_other;

//     // create deque for parts of polygon_input in polygon_other
//     std::deque<polygon_type> parts_in_polygon;

//     // intersect with Boost
//     boost::geometry::intersection(polygon_input, polygon_other, parts_in_polygon);

//     for (polygon_type const &p : parts_in_polygon) {
//         EigenPolyline vertices;
//         for (auto it = boost::end(boost::geometry::exterior_ring(p)) - 1;
//              (it != boost::begin(boost::geometry::exterior_ring(p))); --it) {
//             double x = boost::geometry::get<0>(*it);
//             double y = boost::geometry::get<1>(*it);
//             vertices.emplace_back(x, y);
//         }
//         polygons_within_polygon_other.push_back(vertices);
//     }
//     return polygons_within_polygon_other;
// }

static inline EigenPolyline fromClipperPath(const polygon_type& path, bool close_ring = true) {
    EigenPolyline poly;
    poly.reserve(path.size() + 1);
    for (const auto& p : path) {
        poly.emplace_back(p.x, p.y);
    }
    if (close_ring && !poly.empty() && !same_pt(poly.front(), poly.back())) {
        poly.push_back(poly.front());
    }
    return poly;
}

std::vector<EigenPolyline> polygonWithinPolygonClipper(const polygon_type& polygon_input,
                                                       const polygon_type& polygon_other) {
    std::vector<EigenPolyline> polygons_within;

    mpolygon_type subject{ polygon_input };
    mpolygon_type clip{ polygon_other };
    mpolygon_type result = Clipper2Lib::Intersect(subject, clip, Clipper2Lib::FillRule::NonZero);

    for (const auto& path : result) {
        polygons_within.push_back(fromClipperPath(path));
    }
    return polygons_within;
}

// std::vector<EigenPolyline> polygonWithinPolygonBoost(const EigenPolyline& polygon_input,
//                                                      const polygon_type& polygon_other) {
//     // create Boost polygon from EigenPolyline
//     polygon_type poly_in;
//     polylineToBoostPolygon(polygon_input, poly_in);

//     // check intersection
//     return polygonWithinPolygonBoost(poly_in, polygon_other);
// }

std::vector<EigenPolyline> polygonWithinPolygonClipper(const EigenPolyline& polygon_input,
                                                       const polygon_type& polygon_other) {
    polygon_type input_path;
    polylineToClipperPath(polygon_input, input_path);
    return polygonWithinPolygonClipper(input_path, polygon_other);
}

bool pointInPolygon(const point_type& pt, const polygon_type& polygon) {
    bool inside = false;
    std::size_t n = polygon.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const auto& pi = polygon[i];
        const auto& pj = polygon[j];

        bool intersects = ((pi.y > pt.y) != (pj.y > pt.y)) &&
                          (pt.x < (pj.x - pi.x) * (pt.y - pi.y) / ((pj.y - pi.y) + 1e-12) + pi.x);
        if (intersects)
            inside = !inside;
    }
    return inside;
}

} // namespace util_projection_domain

} // namepsace geometry

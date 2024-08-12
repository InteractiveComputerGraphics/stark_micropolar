#pragma once

#include <vector>
#include <array>

#include <Eigen/Dense>

#include "fem_types.h"

namespace stark
{
namespace Fem
{
template<std::size_t DIM>
struct Quadrature
{
    std::vector<Eigen::Matrix<double, DIM, 1>> points;
    std::vector<double> weights;
    std::string name;

    /// Number of points of the quadrature rule
    int length() const {
        return (int)points.size();
    }

    /// Interpolates the quadrature points of this table to all elements of the given mesh
    template<std::size_t N, typename ElementInterpolator>
    std::vector<Eigen::Vector3d> interpolate_quad_points_generic(const std::vector<std::array<int, N>>& elements, const std::vector<Eigen::Vector3d>& vertices, ElementInterpolator&& interpolate) const {
        std::vector<Eigen::Vector3d> quad_points;
        quad_points.reserve(elements.size() * this->points.size());

        std::vector<Eigen::Vector3d> x_i;
        x_i.reserve(N);

        for (std::size_t i_element = 0; i_element < elements.size(); ++i_element) {
            const auto& element = elements[i_element];

            for (int j = 0; j < (int)element.size(); ++j) {
                x_i.push_back(vertices[element[j]]);
            }

            std::vector<Eigen::Vector3d> x_quad = interpolate(x_i);

            for (auto x : x_quad) {
                quad_points.emplace_back(x);
            }

            x_i.clear();
        }

        return quad_points;
    }
};

template<ElementShape SHAPE, std::size_t DIM>
struct ElementQuadrature : Quadrature<DIM> {};

struct TriQuadrature : ElementQuadrature<ElementShape::Tri2d, 2>
{
    static TriQuadrature tri_p1();
    static TriQuadrature tri_p2();
    static TriQuadrature tri_p2_edges();
    static TriQuadrature tri_p3();
    static TriQuadrature tri_p4();
    static TriQuadrature tri_p5();
    static TriQuadrature tri_p6();

    /// Converts  the quadrature rule into a vector of {weight, xi, eta} arrays, for use with symX add_for_each
    std::vector<std::array<double, 3>> to_symx_coefficients() const;
};

struct QuadQuadrature : ElementQuadrature<ElementShape::Quad2d, 2>
{
    static QuadQuadrature quad_p1();
    static QuadQuadrature quad_p3();
    static QuadQuadrature quad_p5();
};

struct TetQuadrature : ElementQuadrature<ElementShape::Tet3d, 3>
{
    static TetQuadrature tet_p1();
    static TetQuadrature tet_p2();
};
}
}

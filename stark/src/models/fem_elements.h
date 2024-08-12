#pragma once

#include <fmt/format.h>

#include "fem_types.h"

namespace stark
{
// TODO: Document reference element node coordinates
namespace Fem
{

namespace {
double _triangle_area(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
    return 0.5 * (p0 - p2).cross(p1 - p2).norm();
}

double _quad_area(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    return _triangle_area(p0, p1, p2) + _triangle_area(p0, p2, p3);
}
}

/// Common interface for all finite elements
template<typename BasisT>
struct Element : BasisT {
    /// Number of nodes of the element
    static constexpr std::size_t num_nodes = BasisT::num_nodes;
    /// Element type enum value
    static constexpr ElementType element_type = BasisT::element_type;
    /// Dimensionality of the reference element domain
    static constexpr std::size_t ref_dim = BasisT::ref_dim;
    /// Shape classification of the reference element
    static constexpr ElementShape ref_shape = BasisT::ref_shape;
    /// Area of the reference element
    static constexpr double ref_area = BasisT::ref_area;

    /// Returns the concrete element type from this interface
    BasisT& get_basis() const {
        return static_cast<BasisT&>(*this);
    }

    /// Computes or approximates the area of a physical element with the given nodal positions
    static double area(const std::array<Eigen::Vector3d, num_nodes>& vertices) {
        // TODO: Use Eigen Matrix instead of std::array?
        return BasisT::area(vertices);
    }

    /// Returns the nodal positions of the reference element
    static Eigen::Matrix<double, ref_dim, num_nodes> nodes() {
        return BasisT::nodes();
    }

    /// Computes the basis functions at the given reference coordinates
    static Eigen::Vector<double, num_nodes> basis(const Eigen::Vector2d& x_ref) {
        return BasisT::basis(x_ref);
    }

    /// Computes the gradient of the basis functions at the given reference coordinates
    static Eigen::Matrix<double, ref_dim, num_nodes> basis_grad(const Eigen::Vector2d& x_ref) {
        return BasisT::basis_grad(x_ref);
    }

    /// Interpolates nodal data to the given reference coordinates
    template<typename NodalDataT, typename PointT>
    static typename NodalDataT::value_type interpolate(const NodalDataT& nodal_data, const PointT& p) {
        return BasisT::interpolate(nodal_data, p);
    }

    /// Returns the symbolic gradient of the basis functions
    static symx::Matrix grad_symbolic(const symx::Vector& x_ref) {
        return BasisT::grad_symbolic(x_ref);
    }
};

/// Tri3: Degree 1 Lagrange basis on a triangle
/// Basis functions from https://defelement.com/elements/examples/triangle-lagrange-equispaced-1.html
/// Node numbering:
/// ```
///  2
///  | \
///  |   \
///  |     \
///  |       \
///  0---------1
/// ```
struct BasisTri3 {
    static constexpr std::size_t num_nodes = 3;
    static constexpr ElementType element_type = ElementType::Tri3;
    static constexpr std::size_t ref_dim = 2;
    static constexpr ElementShape ref_shape = ElementShape::Tri2d;
    static constexpr double ref_area = 0.5;

    static constexpr auto N0 = [](auto x, auto y) { return 1.0 - x - y; };
    static constexpr auto N1 = [](auto x, auto y) { return x; };
    static constexpr auto N2 = [](auto x, auto y) { return y; };

    static constexpr auto dN0dx = [](auto x, auto y) { return -1.0 + 0.0*x; };
    static constexpr auto dN0dy = [](auto x, auto y) { return -1.0 + 0.0*y; };
    static constexpr auto dN1dx = [](auto x, auto y) { return 1.0 + 0.0*x; };
    static constexpr auto dN1dy = [](auto x, auto y) { return 0.0*y; };
    static constexpr auto dN2dx = [](auto x, auto y) { return 0.0*x; };
    static constexpr auto dN2dy = [](auto x, auto y) { return 1.0 + 0.0*x; };

    static double area(const std::array<Eigen::Vector3d, 3>& vertices) {
        return _triangle_area(vertices[0], vertices[1], vertices[2]);
    }

    static Eigen::Matrix<double, 2, 3> nodes() {
        Eigen::Matrix<double, 2, 3> nodes = Eigen::Matrix<double, 2, 3>::Zero();
        nodes.col(0) = Eigen::Vector2d(0.0, 0.0);
        nodes.col(1) = Eigen::Vector2d(1.0, 0.0);
        nodes.col(2) = Eigen::Vector2d(0.0, 1.0);
        return nodes;
    }

    static Eigen::Vector<double, 3> basis(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        return {N0(x,y), N1(x,y), N2(x,y)};
    }

    static Eigen::Matrix<double, 2, 3> basis_grad(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        Eigen::Matrix<double, 2, 3> grad = Eigen::Matrix<double, 2, 3>::Zero();
        grad(0,0) = dN0dx(x,y);
        grad(1,0) = dN0dy(x,y);
        grad(0,1) = dN1dx(x,y);
        grad(1,1) = dN1dy(x,y);
        grad(0,2) = dN2dx(x,y);
        grad(1,2) = dN2dy(x,y);
        return grad;
    }

    template<typename NodalDataT, typename PointT>
    static typename NodalDataT::value_type interpolate(const NodalDataT& nodal_data, const PointT& p) {
        assert(nodal_data.size() == num_nodes);
        const auto x = p[0];
        const auto y = p[1];
        const auto& d = nodal_data;
        return (d[0] * N0(x,y) + d[1] * N1(x,y) + d[2] * N2(x,y));
    }

    static symx::Matrix grad_symbolic(const symx::Vector& x_ref) {
        symx::Scalar x = x_ref[0];
        symx::Scalar y = x_ref[1];
        return symx::Matrix({
            dN0dx(x,y), dN0dy(x,y),
            dN1dx(x,y), dN1dy(x,y),
            dN2dx(x,y), dN2dy(x,y),
        }, {3, 2}).transpose();
    }
};

/// Tri6: Degree 2 Lagrange basis on a triangle
/// Basis functions from https://defelement.com/elements/examples/triangle-lagrange-equispaced-2.html
/// Note: Re-ordered node order: 0-5-1-3-2-4
/// Node numbering:
/// ```
///  4
///  | \
///  |   \
///  5    3
///  |     \
///  |       \
///  0----1----2
/// ```
struct BasisTri6 {
    static constexpr std::size_t num_nodes = 6;
    static constexpr ElementType element_type = ElementType::Tri6;
    static constexpr std::size_t ref_dim = 2;
    static constexpr ElementShape ref_shape = ElementShape::Tri2d;
    static constexpr double ref_area = 0.5;

    static constexpr auto N0 = [](auto x, auto y) { return 2.0*x*x + 4.0*x*y - 3.0*x + 2.0*y*y - 3.0*y + 1.0; };
    static constexpr auto N1 = [](auto x, auto y) { return 4.0*x * (-x - y + 1.0); };
    static constexpr auto N2 = [](auto x, auto y) { return x * (2.0*x - 1.0); };
    static constexpr auto N3 = [](auto x, auto y) { return 4.0*x*y; };
    static constexpr auto N4 = [](auto x, auto y) { return y * (2.0*y - 1.0); };
    static constexpr auto N5 = [](auto x, auto y) { return 4.0*y * (-x - y + 1.0); };

    static constexpr auto dN0dx = [](auto x, auto y) { return 4.0*x + 4.0*y - 3.0; };
    static constexpr auto dN0dy = [](auto x, auto y) { return 4.0*x + 4.0*y - 3.0; };
    static constexpr auto dN1dx = [](auto x, auto y) { return -8.0*x - 4.0*y + 4.0; };
    static constexpr auto dN1dy = [](auto x, auto y) { return -4.0*x; };
    static constexpr auto dN2dx = [](auto x, auto y) { return 4.0*x - 1.0; };
    static constexpr auto dN2dy = [](auto x, auto y) { return 0.0*x; };
    static constexpr auto dN3dx = [](auto x, auto y) { return 4.0*y; };
    static constexpr auto dN3dy = [](auto x, auto y) { return 4.0*x; };
    static constexpr auto dN4dx = [](auto x, auto y) { return 0.0*x; };
    static constexpr auto dN4dy = [](auto x, auto y) { return 4.0*y - 1.0; };
    static constexpr auto dN5dx = [](auto x, auto y) { return -4.0*y; };
    static constexpr auto dN5dy = [](auto x, auto y) { return -4.0*x - 8.0*y + 4.0; };

    static double area(const std::array<Eigen::Vector3d, 6>& vertices) {
        return _triangle_area(vertices[0], vertices[2], vertices[4]);
    }

    static Eigen::Matrix<double, 2, 6> nodes() {
        Eigen::Matrix<double, 2, 6> nodes = Eigen::Matrix<double, 2, 6>::Zero();
        nodes.col(0) = Eigen::Vector2d(0.0, 0.0);
        nodes.col(1) = Eigen::Vector2d(0.5, 0.0);
        nodes.col(2) = Eigen::Vector2d(1.0, 0.0);
        nodes.col(3) = Eigen::Vector2d(0.5, 0.5);
        nodes.col(4) = Eigen::Vector2d(0.0, 1.0);
        nodes.col(5) = Eigen::Vector2d(0.0, 0.5);
        return nodes;
    }

    static Eigen::Vector<double, 6> basis(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        return {N0(x,y), N1(x,y), N2(x,y), N3(x,y), N4(x,y), N5(x,y)};
    }

    static Eigen::Matrix<double, 2, 6> basis_grad(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        Eigen::Matrix<double, 2, 6> grad = Eigen::Matrix<double, 2, 6>::Zero();
        grad(0,0) = dN0dx(x,y);
        grad(1,0) = dN0dy(x,y);
        grad(0,1) = dN1dx(x,y);
        grad(1,1) = dN1dy(x,y);
        grad(0,2) = dN2dx(x,y);
        grad(1,2) = dN2dy(x,y);
        grad(0,3) = dN3dx(x,y);
        grad(1,3) = dN3dy(x,y);
        grad(0,4) = dN4dx(x,y);
        grad(1,4) = dN4dy(x,y);
        grad(0,5) = dN5dx(x,y);
        grad(1,5) = dN5dy(x,y);
        return grad;
    }

    template<typename NodalDataT, typename PointT>
    static typename NodalDataT::value_type interpolate(const NodalDataT& nodal_data, const PointT& p) {
        assert(nodal_data.size() == num_nodes);
        const auto x = p[0];
        const auto y = p[1];
        const auto& d = nodal_data;
        return (
            d[0] * N0(x,y) + d[1] * N1(x,y) + d[2] * N2(x,y)+
            d[3] * N3(x,y) + d[4] * N4(x,y) + d[5] * N5(x,y)
        );
    }

    static symx::Matrix grad_symbolic(const symx::Vector& x_ref) {
        symx::Scalar x = x_ref[0];
        symx::Scalar y = x_ref[1];
        return symx::Matrix({
            dN0dx(x,y), dN0dy(x,y),
            dN1dx(x,y), dN1dy(x,y),
            dN2dx(x,y), dN2dy(x,y),
            dN3dx(x,y), dN3dy(x,y),
            dN4dx(x,y), dN4dy(x,y),
            dN5dx(x,y), dN5dy(x,y),
        }, {6, 2}).transpose();
    }
};

/// Tri10: Degree 3 Lagrange basis on a triangle
/// Basis functions from: https://defelement.com/elements/examples/triangle-lagrange-equispaced-3.html
/// Note: Re-ordered node order: 0-7-8-1-3-4-2-6-5-9
/// Node numbering:
/// ```
///  6
///  |  \
///  7    5
///  |      \
///  8   9    4
///  |          \
///  0---1---2---3
/// ```
struct BasisTri10 {
    static constexpr std::size_t num_nodes = 10;
    static constexpr ElementType element_type = ElementType::Tri10;
    static constexpr std::size_t ref_dim = 2;
    static constexpr ElementShape ref_shape = ElementShape::Tri2d;
    static constexpr double ref_area = 0.5;

    static constexpr auto N0 = [](auto x, auto y) { return -0.5 * (x*(9.0*x*x + 27.0*x*y + 11.0) + y*(9.0*y*y + 27.0*x*y + 11.0)) + 9.0*x*x + 18.0*x*y + 9.0*y*y + 1.0; };
    static constexpr auto N1 = [](auto x, auto y) { return 0.5 * 9.0*x * (3.0*x*x + 6.0*x*y - 5.0*x + 3.0*y*y - 5.0*y + 2.0); };
    static constexpr auto N2 = [](auto x, auto y) { return 0.5 * 9.0*x * (-3.0*x*x - 3.0*x*y + 4.0*x + y - 1.0); };
    static constexpr auto N3 = [](auto x, auto y) { return 0.5 * (x*(9.0*x*x - 9.0*x + 2.0)); };
    static constexpr auto N4 = [](auto x, auto y) { return 0.5 * (9.0*x*y*(3.0*x - 1.0)); };
    static constexpr auto N5 = [](auto x, auto y) { return 0.5 * (9.0*x*y*(3.0*y - 1.0)); };
    static constexpr auto N6 = [](auto x, auto y) { return 0.5 * (y*(9.0*y*y - 9.0*y + 2.0)); };
    static constexpr auto N7 = [](auto x, auto y) { return 0.5 * 9.0*y * (-3.0*x*y + x - 3.0*y*y + 4.0*y - 1.0); };
    static constexpr auto N8 = [](auto x, auto y) { return 0.5 * 9.0*y * (3.0*x*x + 6.0*x*y - 5.0*x + 3.0*y*y - 5.0*y + 2.0); };
    static constexpr auto N9 = [](auto x, auto y) { return 27.0*x*y*(1.0 - x - y);};

    static constexpr auto dN0dx = [](auto x, auto y) { return 0.5*(-27.0*x*x - 18.0*x*(3.0*y - 2.0) - 27.0*y*y + 36.0*y - 11.0); };
    static constexpr auto dN0dy = [](auto x, auto y) { return 0.5*(-27.0*x*x - 18.0*x*(3.0*y - 2.0) - 27.0*y*y + 36.0*y - 11.0); };
    static constexpr auto dN1dx = [](auto x, auto y) { return 40.5*x*x + x*(54.0*y - 45.0) + y*(13.5*y - 22.5) + 9.0; };
    static constexpr auto dN1dy = [](auto x, auto y) { return 27.0*x*x + x*(27.0*y - 22.5); };
    static constexpr auto dN2dx = [](auto x, auto y) { return -4.5*(9.0*x*x + 6.0*x*y - 8.0*x - y + 1.0); };
    static constexpr auto dN2dy = [](auto x, auto y) { return 4.5*(1.0 - 3.0*x)*x; };
    static constexpr auto dN3dx = [](auto x, auto y) { return 13.5*x*x - 9.0*x + 1.0; };
    static constexpr auto dN3dy = [](auto x, auto y) { return 0.0*x + 0.0*y; };
    static constexpr auto dN4dx = [](auto x, auto y) { return y*(27.0*x - 4.5); };
    static constexpr auto dN4dy = [](auto x, auto y) { return 4.5*x*(3.0*x - 1.0); };
    static constexpr auto dN5dx = [](auto x, auto y) { return 4.5*y*(3.0*y - 1.0); };
    static constexpr auto dN5dy = [](auto x, auto y) { return x*(27.0*y - 4.5); };
    static constexpr auto dN6dx = [](auto x, auto y) { return 0.0*x + 0.0*y; };
    static constexpr auto dN6dy = [](auto x, auto y) { return 13.5*y*y - 9.0*y + 1.0; };
    static constexpr auto dN7dx = [](auto x, auto y) { return 4.5*(1.0 - 3.0*y)*y; };
    static constexpr auto dN7dy = [](auto x, auto y) { return x*(4.5 - 27.0*y) - 40.5*y*y + 36*y - 4.5; };
    static constexpr auto dN8dx = [](auto x, auto y) { return 27.0*y*y + y*(27.0*x - 22.5); };
    static constexpr auto dN8dy = [](auto x, auto y) { return 13.5*x*x + x*(54.0*y - 22.5) + y*(40.5*y - 45.0) + 9.0; };
    static constexpr auto dN9dx = [](auto x, auto y) { return y*(-54.0*x - 27.0*y + 27.0); };
    static constexpr auto dN9dy = [](auto x, auto y) { return x*(-27.0*x - 54.0*y + 27.0); };

    static double area(const std::array<Eigen::Vector3d, 10>& vertices) {
        return _triangle_area(vertices[0], vertices[3], vertices[6]);
    }

    static Eigen::Matrix<double, 2, 10> nodes() {
        Eigen::Matrix<double, 2, 10> nodes = Eigen::Matrix<double, 2, 10>::Zero();
        nodes.col(0) = Eigen::Vector2d(0.0, 0.0);
        nodes.col(1) = Eigen::Vector2d(1.0/3.0, 0.0);
        nodes.col(2) = Eigen::Vector2d(2.0/3.0, 0.0);
        nodes.col(3) = Eigen::Vector2d(1.0, 0.0);
        nodes.col(4) = Eigen::Vector2d(2.0/3.0, 1.0/3.0);
        nodes.col(5) = Eigen::Vector2d(1.0/3.0, 2.0/3.0);
        nodes.col(6) = Eigen::Vector2d(0.0, 1.0);
        nodes.col(7) = Eigen::Vector2d(0.0, 2.0/3.0);
        nodes.col(8) = Eigen::Vector2d(0.0, 1.0/3.0);
        nodes.col(9) = Eigen::Vector2d(1.0/3.0, 1.0/3.0);
        return nodes;
    }

    static Eigen::Vector<double, 10> basis(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        return {N0(x,y), N1(x,y), N2(x,y), N3(x,y), N4(x,y), N5(x,y), N6(x,y), N7(x,y), N8(x,y), N9(x,y)};
    }

    static Eigen::Matrix<double, 2, 10> basis_grad(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        Eigen::Matrix<double, 2, 10> grad = Eigen::Matrix<double, 2, 10>::Zero();
        grad(0,0) = dN0dx(x,y);
        grad(1,0) = dN0dy(x,y);
        grad(0,1) = dN1dx(x,y);
        grad(1,1) = dN1dy(x,y);
        grad(0,2) = dN2dx(x,y);
        grad(1,2) = dN2dy(x,y);
        grad(0,3) = dN3dx(x,y);
        grad(1,3) = dN3dy(x,y);
        grad(0,4) = dN4dx(x,y);
        grad(1,4) = dN4dy(x,y);
        grad(0,5) = dN5dx(x,y);
        grad(1,5) = dN5dy(x,y);
        grad(0,6) = dN6dx(x,y);
        grad(1,6) = dN6dy(x,y);
        grad(0,7) = dN7dx(x,y);
        grad(1,7) = dN7dy(x,y);
        grad(0,8) = dN8dx(x,y);
        grad(1,8) = dN8dy(x,y);
        grad(0,9) = dN9dx(x,y);
        grad(1,9) = dN9dy(x,y);
        return grad;
    }

    template<typename NodalDataT, typename PointT>
    static typename NodalDataT::value_type interpolate(const NodalDataT& nodal_data, const PointT& p) {
        assert(nodal_data.size() == num_nodes);
        const auto x = p[0];
        const auto y = p[1];
        const auto& d = nodal_data;
        return (
            d[0] * N0(x,y) + d[1] * N1(x,y) + d[2] * N2(x,y) + d[3] * N3(x,y) + d[4] * N4(x,y) +
            d[5] * N5(x,y) + d[6] * N6(x,y) + d[7] * N7(x,y) + d[8] * N8(x,y) + d[9] * N9(x,y)
        );
    }

    static symx::Matrix grad_symbolic(const symx::Vector& x_ref) {
        symx::Scalar x = x_ref[0];
        symx::Scalar y = x_ref[1];
        return symx::Matrix({
            dN0dx(x,y), dN0dy(x,y),
            dN1dx(x,y), dN1dy(x,y),
            dN2dx(x,y), dN2dy(x,y),
            dN3dx(x,y), dN3dy(x,y),
            dN4dx(x,y), dN4dy(x,y),
            dN5dx(x,y), dN5dy(x,y),
            dN6dx(x,y), dN6dy(x,y),
            dN7dx(x,y), dN7dy(x,y),
            dN8dx(x,y), dN8dy(x,y),
            dN9dx(x,y), dN9dy(x,y)
        }, {10, 2}).transpose();
    }
};

/// Quad4: Degree 1 Lagrange basis on a quadrilateral
/// Basis functions from: https://defelement.com/elements/examples/quadrilateral-lagrange-equispaced-1.html
/// Note: Re-ordered node order: 0-1-3-2
/// Node numbering:
/// ```
///  3-----2
///  |     |
///  |     |
///  0-----1
/// ```
struct BasisQuad4 {
    static constexpr std::size_t num_nodes = 4;
    static constexpr ElementType element_type = ElementType::Quad4;
    static constexpr std::size_t ref_dim = 2;
    static constexpr ElementShape ref_shape = ElementShape::Quad2d;
    static constexpr double ref_area = 1.0;

    static constexpr auto N0 = [](auto x, auto y) { return x*y - x - y + 1.0; };
    static constexpr auto N1 = [](auto x, auto y) { return x*(1.0-y); };
    static constexpr auto N2 = [](auto x, auto y) { return x*y; };
    static constexpr auto N3 = [](auto x, auto y) { return y*(1.0-x); };

    static constexpr auto dN0dx = [](auto x, auto y) { return y - 1.0; };
    static constexpr auto dN0dy = [](auto x, auto y) { return x - 1.0; };
    static constexpr auto dN1dx = [](auto x, auto y) { return 1.0-y; };
    static constexpr auto dN1dy = [](auto x, auto y) { return -x; };
    static constexpr auto dN2dx = [](auto x, auto y) { return y; };
    static constexpr auto dN2dy = [](auto x, auto y) { return x; };
    static constexpr auto dN3dx = [](auto x, auto y) { return -y; };
    static constexpr auto dN3dy = [](auto x, auto y) { return 1.0 - x; };

    static double area(const std::array<Eigen::Vector3d, 4>& vertices) {
        return _quad_area(vertices[0], vertices[1], vertices[2], vertices[3]);
    }

    static Eigen::Matrix<double, 2, 4> nodes() {
        Eigen::Matrix<double, 2, 4> nodes = Eigen::Matrix<double, 2, 4>::Zero();
        nodes.col(0) = Eigen::Vector2d(0.0, 0.0);
        nodes.col(1) = Eigen::Vector2d(1.0, 0.0);
        nodes.col(2) = Eigen::Vector2d(1.0, 1.0);
        nodes.col(3) = Eigen::Vector2d(0.0, 1.0);
        return nodes;
    }

    static Eigen::Vector<double, 4> basis(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        return {N0(x,y), N1(x,y), N2(x,y), N3(x,y)};
    }

    static Eigen::Matrix<double, 2, 4> basis_grad(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        Eigen::Matrix<double, 2, 4> grad = Eigen::Matrix<double, 2, 4>::Zero();
        grad(0,0) = dN0dx(x,y);
        grad(1,0) = dN0dy(x,y);
        grad(0,1) = dN1dx(x,y);
        grad(1,1) = dN1dy(x,y);
        grad(0,2) = dN2dx(x,y);
        grad(1,2) = dN2dy(x,y);
        grad(0,3) = dN3dx(x,y);
        grad(1,3) = dN3dy(x,y);
        return grad;
    }

    template<typename NodalDataT, typename PointT>
    static typename NodalDataT::value_type interpolate(const NodalDataT& nodal_data, const PointT& p) {
        assert(nodal_data.size() == num_nodes);
        const auto x = p[0];
        const auto y = p[1];
        const auto& d = nodal_data;
        return (d[0] * N0(x,y) + d[1] * N1(x,y) + d[2] * N2(x,y) + d[3] * N3(x,y));
    }

    static symx::Matrix grad_symbolic(const symx::Vector& x_ref) {
        symx::Scalar x = x_ref[0];
        symx::Scalar y = x_ref[1];
        return symx::Matrix({
                dN0dx(x,y), dN0dy(x,y),
                dN1dx(x,y), dN1dy(x,y),
                dN2dx(x,y), dN2dy(x,y),
                dN3dx(x,y), dN3dy(x,y),
        }, {4, 2}).transpose();
    }
};

/// Quad4: Degree 2 Lagrange basis on a quadrilateral
/// Basis functions from: https://defelement.com/elements/examples/quadrilateral-lagrange-equispaced-2.html
/// Note: Re-ordered node order: 0-4-1-6-3-7-2-5-8
/// Node numbering:
/// ```
///  6---5---4
///  |       |
///  7   8   3
///  |       |
///  0---1---2
/// ```
struct BasisQuad9 {
    static constexpr std::size_t num_nodes = 9;
    static constexpr ElementType element_type = ElementType::Quad9;
    static constexpr std::size_t ref_dim = 2;
    static constexpr ElementShape ref_shape = ElementShape::Quad2d;
    static constexpr double ref_area = 1.0;

    static constexpr auto N0 = [](auto x, auto y) { return 4.0*x*x*y*y - 6.0*x*x*y + 2.0*x*x - 6.0*x*y*y + 9.0*x*y - 3.0*x + 2.0*y*y - 3.0*y + 1.0; };
    static constexpr auto N1 = [](auto x, auto y) { return 4.0*x*(-2.0*x*y*y + 3.0*x*y - x + 2.0*y*y - 3.0*y + 1.0); };
    static constexpr auto N2 = [](auto x, auto y) { return x*(4.0*x*y*y - 6.0*x*y + 2.0*x - 2.0*y*y + 3.0*y - 1.0); };
    static constexpr auto N3 = [](auto x, auto y) { return 4.0*x*y*(-2.0*x*y + 2.0*x + y - 1.0); };
    static constexpr auto N4 = [](auto x, auto y) { return x*y*(4.0*x*y - 2.0*x - 2.0*y + 1.0); };
    static constexpr auto N5 = [](auto x, auto y) { return 4.0*x*y*(-2.0*x*y + x + 2.0*y - 1.0); };
    static constexpr auto N6 = [](auto x, auto y) { return y*(4.0*x*x*y - 2.0*x*x - 6.0*x*y + 3.0*x + 2.0*y - 1.0); };
    static constexpr auto N7 = [](auto x, auto y) { return 4.0*y*(-2.0*x*x*y + 2.0*x*x + 3.0*x*y - 3.0*x - y + 1.0); };
    static constexpr auto N8 = [](auto x, auto y) { return 16.0*x*y*(x*y - x - y + 1.0); };

    static constexpr auto dN0dx = [](auto x, auto y) { return x*(8.0*y*y - 12.0*y + 4.0) - 6.0*y*y + 9.0*y - 3.0; };
    static constexpr auto dN0dy = [](auto x, auto y) { return x*x*(8.0*y - 6.0) + x*(9.0 - 12.0*y) + 4.0*y - 3.0; };
    static constexpr auto dN1dx = [](auto x, auto y) { return x*(-16.0*y*y + 24.0*y -8.0) + 8.0*y*y - 12.0*y + 4.0; };
    static constexpr auto dN1dy = [](auto x, auto y) { return -16.0*x*(x*(y - 0.75) - y + 0.75); };
    static constexpr auto dN2dx = [](auto x, auto y) { return x*(8.0*y*y - 12.0*y + 4.0) - 2.0*y*y + 3.0*y - 1.0; };
    static constexpr auto dN2dy = [](auto x, auto y) { return x*(x*(8.0*y - 6.0) - 4.0*y + 3.0); };
    static constexpr auto dN3dx = [](auto x, auto y) { return y*(x*(16.0 - 16.0*y) + 4.0*y - 4.0); };
    static constexpr auto dN3dy = [](auto x, auto y) { return x*(x*(8.0 - 16.0*y) + 8.0*y - 4.0); };
    static constexpr auto dN4dx = [](auto x, auto y) { return y*(x*(8.0*y-4.0)-2.0*y+1.0); };
    static constexpr auto dN4dy = [](auto x, auto y) { return x*(x*(8.0*y - 2.0) - 4.0*y + 1.0); };
    static constexpr auto dN5dx = [](auto x, auto y) { return y*(x*(8.0 - 16.0*y) + 8.0*y - 4.0); };
    static constexpr auto dN5dy = [](auto x, auto y) { return x*(x*(4.0 - 16.0*y) + 16.0*y - 4.0); };
    static constexpr auto dN6dx = [](auto x, auto y) { return y*(x*(8.0*y - 4.0) - 6.0*y + 3.0); };
    static constexpr auto dN6dy = [](auto x, auto y) { return x*x*(8.0*y - 2.0) + x*(3.0-12.0*y) + 4.0*y - 1.0; };
    static constexpr auto dN7dx = [](auto x, auto y) { return -16.0*y*(x*(y - 1.0) - 0.75*y + 0.75); };
    static constexpr auto dN7dy = [](auto x, auto y) { return x*x*(8.0 - 16.0*y) + x*(24.0*y-12.0) - 8.0*y + 4.0; };
    static constexpr auto dN8dx = [](auto x, auto y) { return y*(x*(32.0*y - 32.0) - 16.0*y + 16.0); };
    static constexpr auto dN8dy = [](auto x, auto y) { return x*(x*(32.0*y - 16.0) - 32.0*y + 16.0); };

    static double area(const std::array<Eigen::Vector3d, 9>& vertices) {
        return _quad_area(vertices[0], vertices[2], vertices[4], vertices[6]);
    }

    static Eigen::Matrix<double, 2, 9> nodes() {
        Eigen::Matrix<double, 2, 9> nodes = Eigen::Matrix<double, 2, 9>::Zero();
        nodes.col(0) = Eigen::Vector2d(0.0, 0.0);
        nodes.col(1) = Eigen::Vector2d(0.5, 0.0);
        nodes.col(2) = Eigen::Vector2d(1.0, 1.0);
        nodes.col(3) = Eigen::Vector2d(1.0, 0.5);
        nodes.col(4) = Eigen::Vector2d(1.0, 1.0);
        nodes.col(5) = Eigen::Vector2d(0.5, 1.0);
        nodes.col(6) = Eigen::Vector2d(0.0, 1.0);
        nodes.col(7) = Eigen::Vector2d(0.0, 0.5);
        nodes.col(8) = Eigen::Vector2d(0.5, 0.5);
        return nodes;
    }

    static Eigen::Vector<double, 9> basis(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        return {
            N0(x,y), N1(x,y), N2(x,y), N3(x,y), N4(x,y),
            N5(x,y), N6(x,y), N7(x,y), N8(x,y)
        };
    }

    static Eigen::Matrix<double, 2, 9> basis_grad(const Eigen::Vector2d& x_ref) {
        const auto x = x_ref[0];
        const auto y = x_ref[1];
        Eigen::Matrix<double, 2, 9> grad = Eigen::Matrix<double, 2, 9>::Zero();
        grad(0,0) = dN0dx(x,y);
        grad(1,0) = dN0dy(x,y);
        grad(0,1) = dN1dx(x,y);
        grad(1,1) = dN1dy(x,y);
        grad(0,2) = dN2dx(x,y);
        grad(1,2) = dN2dy(x,y);
        grad(0,3) = dN3dx(x,y);
        grad(1,3) = dN3dy(x,y);
        grad(0,4) = dN4dx(x,y);
        grad(1,4) = dN4dy(x,y);
        grad(0,5) = dN5dx(x,y);
        grad(1,5) = dN5dy(x,y);
        grad(0,6) = dN6dx(x,y);
        grad(1,6) = dN6dy(x,y);
        grad(0,7) = dN7dx(x,y);
        grad(1,7) = dN7dy(x,y);
        grad(0,8) = dN8dx(x,y);
        grad(1,8) = dN8dy(x,y);
        return grad;
    }

    template<typename NodalDataT, typename PointT>
    static typename NodalDataT::value_type interpolate(const NodalDataT& nodal_data, const PointT& p) {
        assert(nodal_data.size() == num_nodes);
        const auto x = p[0];
        const auto y = p[1];
        const auto& d = nodal_data;
        return (
            d[0] * N0(x,y) + d[1] * N1(x,y) + d[2] * N2(x,y) + d[3] * N3(x,y) + d[4] * N0(x,y) +
            d[5] * N1(x,y) + d[6] * N2(x,y) + d[7] * N3(x,y) + d[8] * N3(x,y)
        );
    }

    static symx::Matrix grad_symbolic(const symx::Vector& x_ref) {
        symx::Scalar x = x_ref[0];
        symx::Scalar y = x_ref[1];
        return symx::Matrix({
            dN0dx(x,y), dN0dy(x,y),
            dN1dx(x,y), dN1dy(x,y),
            dN2dx(x,y), dN2dy(x,y),
            dN3dx(x,y), dN3dy(x,y),
            dN4dx(x,y), dN4dy(x,y),
            dN5dx(x,y), dN5dy(x,y),
            dN6dx(x,y), dN6dy(x,y),
            dN7dx(x,y), dN7dy(x,y),
            dN8dx(x,y), dN8dy(x,y),
        }, {9, 2}).transpose();
    }
};
}
}

#include "fem_quadrature.h"

using namespace stark::Fem;

TriQuadrature TriQuadrature::tri_p1()
{
    return TriQuadrature {{{
        {{0.33333333333333333333, 0.33333333333333333333} },
        { 1.0 },
        "p1"
    }}};
}

TriQuadrature TriQuadrature::tri_p2()
{
    // See: https://mathsfromnothing.au/triangle-quadrature-rules
    return TriQuadrature {{{
        {
            {0.66666666666666666667, 0.16666666666666666667},
            {0.16666666666666666667, 0.66666666666666666667},
            {0.16666666666666666667, 0.16666666666666666667},
        },
        {
            0.33333333333333333333,
            0.33333333333333333333,
            0.33333333333333333333
        },
        "p2"
    }}};
}

TriQuadrature TriQuadrature::tri_p2_edges()
{
    return TriQuadrature {{{
        {
            {0.5, 0.0},
            {0.5, 0.5},
            {0.0, 0.5}
        },
        {
            0.33333333333333333333,
            0.33333333333333333333,
            0.33333333333333333333
        },
        "p2e"
    }}};
}

TriQuadrature TriQuadrature::tri_p4()
{
    return TriQuadrature {{{
        {
            {0.445948490915965, 0.108103018168070},
            {0.445948490915965, 0.445948490915965},
            {0.108103018168070, 0.445948490915965},
            {0.091576213509771, 0.816847572980459},
            {0.091576213509771, 0.091576213509771},
            {0.816847572980459, 0.091576213509771}
        },
        {
            0.223381589678011,
            0.223381589678011,
            0.223381589678011,
            0.109951743655322,
            0.109951743655322,
            0.109951743655322,
        },
        "p4"
    }}};
}

TriQuadrature TriQuadrature::tri_p3()
{
    return TriQuadrature {{{
        {
            {0.333333333333333 , 0.333333333333333},
            {0.2               , 0.6},
            {0.2               , 0.2},
            {0.6               , 0.2},
        },
        {
            -0.5625,
            0.520833333333333,
            0.520833333333333,
            0.520833333333333,
        },
        "p3"
    }}};
}

TriQuadrature TriQuadrature::tri_p5()
{
    return TriQuadrature {{{
        {
            {0.333333333333333, 0.333333333333333},
            {0.470142064105115, 0.059715871789770},
            {0.470142064105115, 0.470142064105115},
            {0.059715871789770, 0.470142064105115},
            {0.101286507323456, 0.797426985353087},
            {0.101286507323456, 0.101286507323456},
            {0.797426985353087, 0.101286507323456}
        },
        {
            0.225,
            0.132394152788506,
            0.132394152788506,
            0.132394152788506,
            0.125939180544827,
            0.125939180544827,
            0.125939180544827,
        },
        "p5"
    }}};
}

TriQuadrature TriQuadrature::tri_p6()
{
    return TriQuadrature {{{
        {
            {0.063089014491502, 0.873821971016996},
            {0.063089014491502, 0.063089014491502},
            {0.873821971016996, 0.063089014491502},
            {0.053145049844817, 0.636502499121399},
            {0.310352451033784, 0.053145049844817},
            {0.636502499121399, 0.310352451033784},
            {0.310352451033784, 0.636502499121399},
            {0.053145049844817, 0.310352451033784},
            {0.636502499121399, 0.053145049844817},
            {0.249286745170910, 0.501426509658179},
            {0.249286745170910, 0.249286745170910},
            {0.501426509658179, 0.249286745170910},
        },
        {
            0.050844906370207,
            0.050844906370207,
            0.050844906370207,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.116786275726379,
            0.116786275726379,
            0.116786275726379,
        },
        "p6"
    }}};
}

std::vector<std::array<double, 3>> TriQuadrature::to_symx_coefficients() const {
    std::vector<std::array<double, 3>> quadrature;
    for (int i = 0; i < this->length(); ++i) {
        quadrature.push_back({
            this->weights[i],
            this->points[i][0],
            this->points[i][1],
        });
    }
    return quadrature;
}

/// Source: Wriggers (2008): "Nonlinear finite element methods", Ch. 4.1, p. 117, Table 4.2
QuadQuadrature QuadQuadrature::quad_p1()
{
    return QuadQuadrature {{{
        {{0.5, 0.5} },
        { 1.0 },
        "p1"
    }}};
}

/// Source: Wriggers (2008): "Nonlinear finite element methods", Ch. 4.1, p. 117, Table 4.2
QuadQuadrature QuadQuadrature::quad_p3()
{
    const double a = 1.0 / std::sqrt(3.0);
    return QuadQuadrature {{{
        {
            0.5 * Eigen::Vector2d(-a + 1.0, -a + 1.0),
            0.5 * Eigen::Vector2d(a + 1.0, -a + 1.0),
            0.5 * Eigen::Vector2d(-a + 1.0, a + 1.0),
            0.5 * Eigen::Vector2d(a + 1.0, a + 1.0),
        },
        {
            0.25,
            0.25,
            0.25,
            0.25,
        },
        "p3"
    }}};
}

/// Source: Wriggers (2008): "Nonlinear finite element methods", Ch. 4.1, p. 117, Table 4.2
QuadQuadrature QuadQuadrature::quad_p5()
{
    const double a = std::sqrt(3.0 / 5.0);
    return QuadQuadrature {{{
        {
            0.5 * Eigen::Vector2d(-a + 1.0, -a + 1.0),
            0.5 * Eigen::Vector2d(0.0 + 1.0, -a + 1.0),
            0.5 * Eigen::Vector2d(a + 1.0, -a + 1.0),
            0.5 * Eigen::Vector2d(-a + 1.0, 0.0 + 1.0),
            0.5 * Eigen::Vector2d(0.0 + 1.0, 0.0 + 1.0),
            0.5 * Eigen::Vector2d(a + 1.0, 0.0 + 1.0),
            0.5 * Eigen::Vector2d(-a + 1.0, a + 1.0),
            0.5 * Eigen::Vector2d(0.0 + 1.0, a + 1.0),
            0.5 * Eigen::Vector2d(a + 1.0, a + 1.0),
        },
        {
            0.25 * (25.0/81.0),
            0.25 * (40.0/81.0),
            0.25 * (25.0/81.0),
            0.25 * (40.0/81.0),
            0.25 * (64.0/81.0),
            0.25 * (40.0/81.0),
            0.25 * (25.0/81.0),
            0.25 * (40.0/81.0),
            0.25 * (25.0/81.0),
        },
        "p5"
    }}};
}

TetQuadrature TetQuadrature::tet_p1()
{
    return TetQuadrature {{{
        { Eigen::Vector3d(0.25, 0.25, 0.25) },
        { 1.0/6.0 },
        "p1",
    }}};
}

TetQuadrature TetQuadrature::tet_p2()
{
    const double w = 1.0/24.0;
    const double a = 0.5854101966249685;    // 1/(3*sqrt(5) - 5)
    const double b = 0.1381966011250105;    // 1/(5 + sqrt(5))

    return TetQuadrature {{{
        {
            Eigen::Vector3d(a,b,b),
            Eigen::Vector3d(b,a,b),
            Eigen::Vector3d(b,b,a),
            Eigen::Vector3d(b,b,b)
        },
        { w, w, w, w },
        "p2",
    }}};
}

#pragma once
#include <vector>

#include "../../../core/Stark.h"
#include "../PointDynamics.h"
#include "../../types.h"

namespace stark
{
    /// Energy from "Kirchhoff-Love Shells with Arbitrary Hyperelastic Materials", Jiahao Wen, Jernej Barbič; 2023
    class EnergyTriangleStrainWen23
    {
    public:
        struct Params {
            STARK_PARAM_NON_NEGATIVE(double, thickness, 1e-3)
            STARK_PARAM_YOUNGS_MODULUS()
            STARK_PARAM_POISSONS_RATIO()
        };
        struct Handler { STARK_COMMON_HANDLER_CONTENTS(EnergyTriangleStrainWen23, Params) };

        /* Fields */
        const spPointDynamics dyn;
        symx::LabelledConnectivity<8> conn{ { "idx", "group", "i", "j", "k", "i_opposite", "j_opposite", "k_opposite" } };

        std::vector<int> triangle_offset;  // per group
        std::vector<std::vector<std::array<int, 3>>> triangles;  // per group
        std::vector<std::vector<std::array<int, 3>>> opposite_vertices;  // per group
        std::vector<PointSetHandler> point_sets;

        std::vector<double> thickness;  // per group
        std::vector<double> poissons_ratio;  // per group
        std::vector<double> youngs_modulus;  // per group

        std::vector<std::array<double, 3>> is_opposite_boundary_edge; // per triangle

        std::vector<double> triangle_area_rest;  // per triangle
        std::vector<double> h_rest;  // per triangle
        std::vector<double> k_rest;  // per triangle
        std::vector<std::array<double, 9>> t_rest_inv; // per triangle
        std::vector<std::array<double, 4>> l_rest; // per triangle
        std::vector<Eigen::Vector3d> n_rest;  // per triangle

        struct MidedgeNormal {
            Eigen::Vector3d n;  // normal on the edge
            int group;          // the mesh group of the edge
            int v0;             // vertex 1 of the edge
            int v1;             // vertex 2 of the edge
        };
        std::vector<MidedgeNormal> n_edge_rest;  // per edge per triangle

        /* Methods */
        EnergyTriangleStrainWen23(stark::core::Stark& stark, spPointDynamics dyn);
        Handler add(const PointSetHandler& set, const std::vector<std::array<int, 3>>& triangles, const Params& params);
        Params get_params(const Handler& handler) const;
        void set_params(const Handler& handler, const Params& params);

    private:
        void _before_simulation(stark::core::Stark& stark);
        void _write_frame(stark::core::Stark& stark);
    };
}

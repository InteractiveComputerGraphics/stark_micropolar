#pragma once
#include <vector>

#include "../../../core/Stark.h"
#include "../PointDynamics.h"
#include "../../types.h"

namespace stark
{
    /// Energy from "A Finite Element Formulation of Baraff-Witkin Cloth", Theodore Kim; 2020
    class EnergyTriangleStrainKim20
    {
    public:
        struct Params {
            STARK_PARAM_NON_NEGATIVE(double, thickness, 1e-3)
            STARK_PARAM_STIFFNESS()
            STARK_PARAM_STRAIN_LIMITING()
        };
        struct Handler { STARK_COMMON_HANDLER_CONTENTS(EnergyTriangleStrainKim20, Params) };

        /* Fields */
        const spPointDynamics dyn;
        symx::LabelledConnectivity<5> conn{ { "idx", "group", "i", "j", "k" } };

        std::vector<double> thickness;  // per group
        std::vector<double> stiffness;  // per group
        std::vector<double> strain_limit;  // per group
        std::vector<double> strain_limit_stiffness;  // per group

        std::vector<double> triangle_area_rest;  // per triangle
        std::vector<std::array<double, 4>> DXinv;  // per triangle

        /* Methods */
        EnergyTriangleStrainKim20(stark::core::Stark& stark, spPointDynamics dyn);
        Handler add(const PointSetHandler& set, const std::vector<std::array<int, 3>>& triangles, const Params& params);
        Params get_params(const Handler& handler) const;
        void set_params(const Handler& handler, const Params& params);
    };
}

#include "EnergyTriangleStrainKim20.h"

#include "../../time_integration.h"
#include "../../../utils/mesh_utils.h"


stark::EnergyTriangleStrainKim20::EnergyTriangleStrainKim20(stark::core::Stark& stark, spPointDynamics dyn)
	: dyn(dyn)
{
	// Energy from "A Finite Element Formulation of Baraff-Witkin Cloth", Theodore Kim; 2020
	stark.global_energy.add_energy("EnergyTriangleStrainKim20", this->conn,
		[&](symx::Energy& energy, symx::Element& conn)
		{
			// Unpack connectivity
			std::vector<symx::Index> triangle = conn.slice(2, 5);

			// Create symbols
			std::vector<symx::Vector> v1 = energy.make_dof_vectors(this->dyn->dof, this->dyn->v1.data, triangle);
			std::vector<symx::Vector> x0 = energy.make_vectors(this->dyn->x0.data, triangle);
			symx::Matrix DXinv = energy.make_matrix(this->DXinv, { 2, 2 }, conn["idx"]);
			symx::Scalar rest_area = energy.make_scalar(this->triangle_area_rest, conn["idx"]);
			symx::Scalar thickness = energy.make_scalar(this->thickness, conn["group"]);
			symx::Scalar stiffness = energy.make_scalar(this->stiffness, conn["group"]);
			symx::Scalar strain_limit = energy.make_scalar(this->strain_limit, conn["group"]);
			symx::Scalar strain_limit_stiffness = energy.make_scalar(this->strain_limit_stiffness, conn["group"]);
			symx::Scalar dt = energy.make_scalar(stark.dt);

			// Time integration
			std::vector<symx::Vector> x1 = time_integration(x0, v1, dt);

			// Kinematics
			symx::Matrix Dx_32 = symx::Matrix(symx::gather({ x1[1] - x1[0], x1[2] - x1[0] }), { 2, 3 }).transpose();
			symx::Matrix F_32 = Dx_32 * DXinv;  // 3x2

			// Anisotropic Baraff-Witkin energy

			// TODO: Make these parameters
			symx::Matrix au = energy.make_zero_matrix({2, 1});
			au(0, 0) = energy.make_one();
			symx::Matrix av = energy.make_zero_matrix({2, 1});
			av(1, 0) = energy.make_one();
			symx::Scalar bu = energy.make_one();
			symx::Scalar bv = energy.make_one();

			const auto I5 = [](const symx::Matrix& F, const symx::Matrix& a) {
				return (a.transpose() * F.transpose() * F * a)(0,0);
			};

			const auto I6 = [](const symx::Matrix& F, const symx::Matrix& a, const symx::Matrix& b) {
				return (a.transpose() * F.transpose() * F * b)(0,0);
			};

			symx::Scalar Psi_stretch = (I5(F_32, au).sqrt() - bu).powN(2) + (I5(F_32, av).sqrt() - bv).powN(2);
			symx::Scalar Psi_shear = I6(F_32, au, av).powN(2);
			symx::Scalar Psi = stiffness * (Psi_stretch + Psi_shear);

			symx::Scalar E_elasticity = thickness * rest_area * Psi;

			// Strain limiting
			symx::Matrix C = F_32.transpose() * F_32;
			symx::Scalar s1 = symx::sqrt(C.singular_values_2x2()[0]);
			symx::Scalar dl = s1 - (strain_limit + 1.0);
			symx::Scalar E_sl = symx::branch(dl > 0.0, thickness * rest_area * strain_limit_stiffness * dl.powN(3)/3.0, 0.0);

			energy.set(E_elasticity + E_sl);
		}
	);
}
stark::EnergyTriangleStrainKim20::Handler stark::EnergyTriangleStrainKim20::add(
        const PointSetHandler& set,
        const std::vector<std::array<int, 3>>& triangles,
        const Params& params)
{
    set.exit_if_not_valid("EnergyTriangleStrainKim20::add");
	const int group = (int)this->thickness.size();

	this->thickness.push_back(params.thickness);
	this->stiffness.push_back(params.stiffness);
	this->strain_limit.push_back(params.strain_limit);
	this->strain_limit_stiffness.push_back(params.strain_limit_stiffness);

	for (int tri_i = 0; tri_i < (int)triangles.size(); tri_i++) {
		// Connectivity
		const std::array<int, 3>& conn_loc = triangles[tri_i];
		const std::array<int, 3> conn_glob = set.get_global_indices(conn_loc);
		this->conn.numbered_push_back({ group, conn_glob[0], conn_glob[1], conn_glob[2] });

		// Fetch coordinates
		const Eigen::Vector3d& A = this->dyn->X[conn_glob[0]];
		const Eigen::Vector3d& B = this->dyn->X[conn_glob[1]];
		const Eigen::Vector3d& C = this->dyn->X[conn_glob[2]];

		// Area
		this->triangle_area_rest.push_back(triangle_area(A, B, C));

		// DXinv
		//// Projection matrix
		const Eigen::Vector3d u = (B - A).normalized();
		const Eigen::Vector3d n = u.cross(C - A);
		const Eigen::Vector3d v = u.cross(n).normalized();
		Eigen::Matrix<double, 2, 3> P;
		P.row(0) = u;
		P.row(1) = v;

		//// Projection
		const Eigen::Vector2d A_ = P*A;
		const Eigen::Vector2d B_ = P*B;
		const Eigen::Vector2d C_ = P*C;

		//// DX
		Eigen::Matrix2d DX;
		DX.col(0) = B_ - A_;
		DX.col(1) = C_ - A_;
		Eigen::Matrix2d DXinv = DX.inverse();
		this->DXinv.push_back({DXinv(0, 0), DXinv(0, 1), DXinv(1, 0), DXinv(1, 1)});
	}

    return Handler(this, group);
}
stark::EnergyTriangleStrainKim20::Params stark::EnergyTriangleStrainKim20::get_params(const Handler& handler) const
{
    handler.exit_if_not_valid("EnergyTriangleStrainKim20::get_params");
    const int group = handler.get_idx();

    Params params;
    params.thickness = this->thickness[group];
    params.stiffness = this->stiffness[group];
    params.strain_limit = this->strain_limit[group];
    params.strain_limit_stiffness = this->strain_limit_stiffness[group];
    return params;
}
void stark::EnergyTriangleStrainKim20::set_params(const Handler& handler, const Params& params)
{
    handler.exit_if_not_valid("EnergyTriangleStrainKim20::set_params");
    const int group = handler.get_idx();
    
    this->thickness[group] = params.thickness;
    this->stiffness[group] = params.stiffness;
    this->strain_limit[group] = params.strain_limit;
    this->strain_limit_stiffness[group] = params.strain_limit_stiffness;
}

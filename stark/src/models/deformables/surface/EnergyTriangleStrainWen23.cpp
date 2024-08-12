#include "EnergyTriangleStrainWen23.h"

#include "../../time_integration.h"
#include "../../../utils/mesh_utils.h"


stark::EnergyTriangleStrainWen23::EnergyTriangleStrainWen23(stark::core::Stark& stark, spPointDynamics dyn)
	: dyn(dyn)
{
	stark.callbacks.add_before_simulation([&]() { this->_before_simulation(stark); });
	stark.callbacks.add_write_frame([&]() { this->_write_frame(stark); });

	this->triangle_offset.push_back(0);

	// Energy from "Kirchhoff-Love Shells with Arbitrary Hyperelastic Materials", Jiahao Wen, Jernej Barbič; 2023
	stark.global_energy.add_energy("EnergyTriangleStrainWen23", this->conn,
		[&](symx::Energy& energy, symx::Element& conn)
		{
			energy.set_never_project_to_PD(stark.settings.models.never_project_tri_wen23);

			// Unpack connectivity
			std::vector<symx::Index> triangle = conn.slice(2, 5);
			std::vector<symx::Index> vertices_opposite = conn.slice(5, 8);

			// Create symbols
			std::vector<symx::Vector> v1 = energy.make_dof_vectors(this->dyn->dof, this->dyn->v1.data, triangle);
			std::vector<symx::Vector> x0 = energy.make_vectors(this->dyn->x0.data, triangle);
			symx::Scalar E = energy.make_scalar(this->youngs_modulus, conn["group"]);
			symx::Scalar nu = energy.make_scalar(this->poissons_ratio, conn["group"]);
			symx::Scalar h = energy.make_scalar(this->thickness, conn["group"]);
			symx::Scalar dt = energy.make_scalar(stark.dt);

			symx::Scalar area = energy.make_scalar(this->triangle_area_rest, conn["idx"]);
			symx::Scalar H = energy.make_scalar(this->h_rest, conn["idx"]);
			symx::Scalar K = energy.make_scalar(this->k_rest, conn["idx"]);
			symx::Matrix T_inv = energy.make_matrix(this->t_rest_inv, { 3, 3 }, conn["idx"]);
			symx::Matrix L = energy.make_matrix(this->l_rest, { 2, 2 }, conn["idx"]);

			symx::Vector is_opposite_boundary_edge = energy.make_vector(this->is_opposite_boundary_edge, conn["idx"]);
			symx::Scalar i_is_boundary = is_opposite_boundary_edge[0];
			symx::Scalar j_is_boundary = is_opposite_boundary_edge[1];
			symx::Scalar k_is_boundary = is_opposite_boundary_edge[2];

			// State of vertices opposite of triangle vertices
			std::vector<symx::Vector> v1_opp = energy.make_dof_vectors(this->dyn->dof, this->dyn->v1.data, vertices_opposite);
			std::vector<symx::Vector> x0_opp = energy.make_vectors(this->dyn->x0.data, vertices_opposite);

			// Time integration
			std::vector<symx::Vector> x1 = time_integration(x0, v1, dt);
			std::vector<symx::Vector> x1_opp = time_integration(x0_opp, v1_opp, dt);

			auto normal = [](symx::Vector a, symx::Vector b, symx::Vector c) {
				return (b - a).normalized().cross3((c - a).normalized()).normalized();
			};

			symx::Vector x_i = x1[0];
			symx::Vector x_j = x1[1];
			symx::Vector x_k = x1[2];

			// Compute triangle normal
			symx::Vector n_ijk = normal(x_i, x_j, x_k);
			// Compute neighboring triangles' normals
			symx::Vector n_i_opp_tri = normal(x1_opp[0], x_k, x_j);
			symx::Vector n_j_opp_tri = normal(x1_opp[1], x_i, x_k);
			symx::Vector n_k_opp_tri = normal(x1_opp[2], x_j, x_i);
			// Determine mid-edge normals
			symx::Vector n_i = i_is_boundary * n_ijk + (1.0 - i_is_boundary) * 0.5 * (n_ijk + n_i_opp_tri);
			symx::Vector n_j = j_is_boundary * n_ijk + (1.0 - j_is_boundary) * 0.5 * (n_ijk + n_j_opp_tri);
			symx::Vector n_k = k_is_boundary * n_ijk + (1.0 - k_is_boundary) * 0.5 * (n_ijk + n_k_opp_tri);

			//
			symx::Matrix t = energy.make_zero_matrix({3, 3});
			t.set_col(0, x_j - x_i);
			t.set_col(1, x_k - x_i);
			t.set_col(2, n_ijk);

			/*
			symx::Matrix q_half = energy.make_zero_matrix({3, 3});
			q_half.set_col(0, n_i - n_j);
			q_half.set_col(1, n_i - n_k);
			symx::Matrix q = 2.0 * q_half;
			*/

			symx::Matrix a = energy.make_zero_matrix({2, 2});
			a(0, 0) = (x_j - x_i).dot(x_j - x_i);
			a(0, 1) = (x_j - x_i).dot(x_k - x_i);
			a(1, 0) = (x_k - x_i).dot(x_j - x_i);
			a(1, 1) = (x_k - x_i).dot(x_k - x_i);
			symx::Matrix a_inv = a.inv();

			symx::Matrix b_half = energy.make_zero_matrix({2, 2});
			b_half(0, 0) = (n_i - n_j).dot(x_i - x_j);
			b_half(0, 1) = (n_i - n_j).dot(x_i - x_k);
			b_half(1, 0) = (n_i - n_k).dot(x_i - x_j);
			b_half(1, 1) = (n_i - n_k).dot(x_i - x_k);
			symx::Matrix b = 2.0 * b_half;

			symx::Matrix l = a_inv * b;
			symx::Matrix delta_l = L - l;
			symx::Matrix delta_l_L = delta_l * L;

			symx::Matrix F1_l = energy.make_zero_matrix({3, 3});
			F1_l(0,0) = delta_l(0,0);
			F1_l(0,1) = delta_l(0,1);
			F1_l(1,0) = delta_l(1,0);
			F1_l(1,1) = delta_l(1,1);

			symx::Matrix F2_l = energy.make_zero_matrix({3, 3});
			F2_l(0,0) = delta_l_L(0,0);
			F2_l(0,1) = delta_l_L(0,1);
			F2_l(1,0) = delta_l_L(1,0);
			F2_l(1,1) = delta_l_L(1,1);

			symx::Matrix F0 = t * T_inv;
			symx::Matrix F1 = t * (F1_l * T_inv);
			symx::Matrix F2 = t * (F2_l * T_inv);

			symx::Matrix I = energy.make_identity_matrix(3);
			auto Psi = [E, nu, I](const symx::Matrix& F) -> symx::Scalar {
				// StVK
				symx::Scalar mu = E / (2.0 * (1.0 + nu));
				symx::Scalar lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)); // 3D
				symx::Matrix E = (F.transpose() * F - I) * 0.5;

				symx::Scalar energy_density = mu * E.frobenius_norm_sq() + 0.5 * lambda * E.trace().powN(2);
				return energy_density;

				/*
				// Stable Neohookean
				symx::Scalar mu = E / (2.0 * (1.0 + nu));
				symx::Scalar lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)); // 3D
				symx::Scalar mu_ = (4.0 / 3.0) * mu;
				symx::Scalar lambda_ = lambda + (5.0 / 6.0) * mu;
				symx::Scalar detF = F.det();
				symx::Scalar Ic = F.frobenius_norm_sq();
				symx::Scalar alpha = 1.0 + mu_ / lambda_ - mu_ / (4.0 * lambda_);
				symx::Scalar energy_density = 0.5 * mu_ * (Ic - 3.0) + 0.5 * lambda_ * (detF - alpha).powN(2) - 0.5 * mu_ * symx::log(Ic + 1.0);
				return energy_density;
				*/
			};

			/*
			symx::Matrix factor_F1 = (std::sqrt(3.0)/6.0) * h * F1;
			symx::Matrix factor_F2 = (h*h / 24.0) * (F2 - 2.0 * H * F1);

			symx::Scalar energy_density = (h*h*K / 12.0) * Psi(F0) + 0.5 * (Psi(F0 + factor_F1) + Psi(F0 - factor_F1)) + Psi(F0 + factor_F2) - Psi(F0 - factor_F2);
			*/

			// Constants c_i for:
			// (c0 + c1 * h^2*K) * Psi(c2*F0 + (c3*h + c4*h^2*H)*F1 + c5*h^2*F2)
			std::vector<std::array<double, 6>> constants = {
				// c[0],		c[1],  c[2],					 c[3],		c[4],		c[5]
				{ 0.0,	1.0/12.0,	1.0,					  0.0,		 0.0,		0.0},
				{ 0.5,		 0.0,	1.0,  std::sqrt(3.0) / 6.0,		 0.0,		0.0},
				{ 0.5,		 0.0,	1.0, -std::sqrt(3.0) / 6.0,		 0.0,		0.0},
				{ 1.0,		 0.0,	1.0,					  0.0, -1.0/12.0,  1.0/24.0},
				{-1.0,		 0.0,	1.0,					  0.0,  1.0/12.0, -1.0/24.0},
			};

			symx::Scalar energy_density = energy.add_for_each(constants, [&](symx::Vector& c) {
				return (c[0] + c[1]*h*h*K) * Psi(c[2]*F0 + (c[3]*h + c[4]*h*h*H)*F1 + c[5]*h*h*F2);
			});

			symx::Scalar Energy = area * h * energy_density;
			energy.set(Energy);
		}
	);
}
stark::EnergyTriangleStrainWen23::Handler stark::EnergyTriangleStrainWen23::add(
    const PointSetHandler& set,
	const std::vector<std::array<int, 3>>& triangles,
	const Params& params)
{
    set.exit_if_not_valid("EnergyTriangleStrainWen23::add");
	const int group = (int)this->triangles.size();

	this->triangles.push_back(triangles);
	this->triangle_offset.push_back(triangles.size());
	this->point_sets.push_back(set);

	this->thickness.push_back(params.thickness);
	this->youngs_modulus.push_back(params.youngs_modulus);
	this->poissons_ratio.push_back(params.poissons_ratio);

	const int num_verts = set.size();
	const int num_tris = triangles.size();

	// Find internal_angles (dihedral) connectivity
	std::vector<std::array<int, 4>> internal_edges;
	find_internal_angles(internal_edges, triangles, num_verts);

	std::vector<std::array<int, 3>> opposite_vertices;
	opposite_vertices.resize(num_tris, {-1, -1, -1});

	std::vector<std::array<double, 3>> is_opposite_boundary_edge;
	is_opposite_boundary_edge.resize(num_tris, {true, true, true});

	// Collect "opposite vertex" per vertex per triangle
	//
	// For interior edges, opposite vertices can be determined using the information from `find_internal_angles`
	// If no corresponding interior edge was found, the vertex is marked to be opposite to itself
	{
		const auto hinge_cmp = [](const std::array<int, 4>& a, const std::array<int, 4>& b) {
			return (a[0] < b[0]) || (a[0] == b[0] && a[1] < b[1]);
		};

		std::vector<std::array<int, 4>> sorted_hinges = internal_edges;
		std::sort(sorted_hinges.begin(), sorted_hinges.end(), hinge_cmp);

		std::vector<std::array<int, 4>> hinge_buffer;

		for (int tri_idx = 0; tri_idx < num_tris; tri_idx++) {
			const auto& tri = triangles[tri_idx];
			for (int edge_idx = 0; edge_idx < 3; edge_idx++) {
				const int e0 = tri[edge_idx];
				const int e1 = tri[(edge_idx + 1) % 3];

				int vert_idx_loc = (edge_idx + 2) % 3;
				const int vert = tri[vert_idx_loc];

				std::array<int, 4> edge = {e0, e1, 0, 0};
				auto l = std::lower_bound(sorted_hinges.begin(), sorted_hinges.end(), edge, hinge_cmp);
				auto u = std::upper_bound(sorted_hinges.begin(), sorted_hinges.end(), edge, hinge_cmp);
				hinge_buffer.insert(hinge_buffer.end(), l, u);

				edge = {e1, e0, 0, 0};
				l = std::lower_bound(sorted_hinges.begin(), sorted_hinges.end(), edge, hinge_cmp);
				u = std::upper_bound(sorted_hinges.begin(), sorted_hinges.end(), edge, hinge_cmp);
				hinge_buffer.insert(hinge_buffer.end(), l, u);

				if (hinge_buffer.empty()) {
					opposite_vertices[tri_idx][vert_idx_loc] = vert;
				} else {
					assert(hinge_buffer.size() == 1);
					is_opposite_boundary_edge[tri_idx][vert_idx_loc] = false;

					edge = hinge_buffer[0];
					const int opposite_vert = edge[2] == vert ? edge[3] : edge[2];
					opposite_vertices[tri_idx][vert_idx_loc] = opposite_vert;
				}

				hinge_buffer.clear();
			}
		}

		/*
		for (int tri_idx = 0; tri_idx < num_tris; tri_idx++) {
			const auto& tri = triangles[tri_idx];
			for (int j = 0; j < 3; j++) {
				opposite_vertices[tri_idx][j] = tri[j];
			}
		}
		*/
	}

	this->opposite_vertices.push_back(opposite_vertices);
	this->is_opposite_boundary_edge.insert(this->is_opposite_boundary_edge.end(), is_opposite_boundary_edge.begin(), is_opposite_boundary_edge.end());

	// Initialize connectivity & rest state
	for (int tri_i = 0; tri_i < (int)triangles.size(); tri_i++) {
		const std::array<int, 3>& conn_loc = triangles[tri_i];
		const std::array<int, 3> conn_glob = set.get_global_indices(conn_loc);
		const std::array<int, 3> conn_opp_loc = opposite_vertices[tri_i];
		const std::array<int, 3> conn_opp_glob = set.get_global_indices(conn_opp_loc);

		// Initialize connectivity
		this->conn.numbered_push_back({
			group,
			conn_glob[0],
			conn_glob[1],
			conn_glob[2],
			conn_opp_glob[0],
			conn_opp_glob[1],
			conn_opp_glob[2]
		});
	}

    return Handler(this, group);
}
stark::EnergyTriangleStrainWen23::Params stark::EnergyTriangleStrainWen23::get_params(const Handler& handler) const
{
    handler.exit_if_not_valid("EnergyTriangleStrainWen23::get_params");
    const int group = handler.get_idx();

    Params params;
    params.thickness = this->thickness[group];
    params.youngs_modulus = this->youngs_modulus[group];
    params.poissons_ratio = this->poissons_ratio[group];
    return params;
}
void stark::EnergyTriangleStrainWen23::set_params(const Handler& handler, const Params& params)
{
    handler.exit_if_not_valid("EnergyTriangleStrainWen23::set_params");
    const int group = handler.get_idx();

    this->thickness[group] = params.thickness;
    this->youngs_modulus[group] = params.youngs_modulus;
    this->poissons_ratio[group] = params.poissons_ratio;
}

void stark::EnergyTriangleStrainWen23::_before_simulation(stark::core::Stark &stark)
{
	for(int group = 0; group < (int)this->triangles.size(); ++group) {
		// Initialize connectivity & rest state
		for(int tri_i = 0; tri_i < (int)this->triangles[group].size(); tri_i++) {
			const auto& set = this->point_sets[group];
			const auto& triangles = this->triangles[group];
			const int triangle_offset = this->triangle_offset[group];

			const std::array<int, 3>& conn_loc = triangles[tri_i];
			const std::array<int, 3> conn_glob = set.get_global_indices(conn_loc);
			const std::array<int, 3> conn_opp_loc = this->opposite_vertices[group][tri_i];
			const std::array<int, 3> conn_opp_glob = set.get_global_indices(conn_opp_loc);

			// Initialize rest state
			{
				const int i = conn_glob[0];
				const int j = conn_glob[1];
				const int k = conn_glob[2];

				// Fetch coordinates
				const Eigen::Vector3d& x_i = this->dyn->X[conn_glob[0]];
				const Eigen::Vector3d& x_j = this->dyn->X[conn_glob[1]];
				const Eigen::Vector3d& x_k = this->dyn->X[conn_glob[2]];

				// Area
				this->triangle_area_rest.push_back(triangle_area(x_i, x_j, x_k));

				// Fetch coordinates of neighboring vertices
				const Eigen::Vector3d& x_i_opp = this->dyn->X[conn_opp_glob[0]];
				const Eigen::Vector3d& x_j_opp = this->dyn->X[conn_opp_glob[1]];
				const Eigen::Vector3d& x_k_opp = this->dyn->X[conn_opp_glob[2]];

				auto normal = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) {
					return (b - a).normalized().cross((c - a).normalized()).normalized();
				};

				// Compute triangle normal
				const Eigen::Vector3d n_ijk = normal(x_i, x_j, x_k);
				// Compute neighboring triangles' normals
				const Eigen::Vector3d n_i_opp_tri = normal(x_i_opp, x_k, x_j);
				const Eigen::Vector3d n_j_opp_tri = normal(x_j_opp, x_i, x_k);
				const Eigen::Vector3d n_k_opp_tri = normal(x_k_opp, x_j, x_i);
				// Determine mid-edge normals
				const Eigen::Vector3d n_i = is_opposite_boundary_edge[triangle_offset + tri_i][0] == true ? n_ijk : 0.5 * (n_ijk + n_i_opp_tri);
				const Eigen::Vector3d n_j = is_opposite_boundary_edge[triangle_offset + tri_i][1] == true ? n_ijk : 0.5 * (n_ijk + n_j_opp_tri);
				const Eigen::Vector3d n_k = is_opposite_boundary_edge[triangle_offset + tri_i][2] == true ? n_ijk : 0.5 * (n_ijk + n_k_opp_tri);

				this->n_edge_rest.push_back({ n_i, group, j, k});
				this->n_edge_rest.push_back({ n_j, group, i, k});
				this->n_edge_rest.push_back({ n_k, group, j, i});

				Eigen::Matrix3d t = Eigen::Matrix3d::Zero();
				t.col(0) = x_j - x_i;
				t.col(1) = x_k - x_i;
				t.col(2) = n_ijk;
				Eigen::Matrix3d t_inv = t.inverse();

				/*
				Eigen::Matrix3d q_half = Eigen::Matrix3d::Zero();
				q_half.col(0) = n_i - n_j;
				q_half.col(1) = n_i - n_k;
				Eigen::Matrix3d q = 2.0 * q_half;
				*/

				Eigen::Matrix2d a = Eigen::Matrix2d::Zero();
				a(0, 0) = (x_j - x_i).dot(x_j - x_i);
				a(0, 1) = (x_j - x_i).dot(x_k - x_i);
				a(1, 0) = (x_k - x_i).dot(x_j - x_i);
				a(1, 1) = (x_k - x_i).dot(x_k - x_i);

				Eigen::Matrix2d b_half = Eigen::Matrix2d::Zero();
				b_half(0, 0) = (n_i - n_j).dot(x_i - x_j);
				b_half(0, 1) = (n_i - n_j).dot(x_i - x_k);
				b_half(1, 0) = (n_i - n_k).dot(x_i - x_j);
				b_half(1, 1) = (n_i - n_k).dot(x_i - x_k);
				Eigen::Matrix2d b = 2.0 * b_half;

				Eigen::Matrix2d a_inv = a.inverse();
				Eigen::Matrix2d l = a_inv * b;
				double H = l.trace() / 2.0;
				double K = l.determinant();

				this->t_rest_inv.push_back({
					t_inv(0,0), t_inv(0,1), t_inv(0,2),
					t_inv(1,0), t_inv(1,1), t_inv(1,2),
					t_inv(2,0), t_inv(2,1), t_inv(2,2)
				});
				this->l_rest.push_back({
					l(0,0), l(0,1),
					l(1,0), l(1,1)
				});
				this->h_rest.push_back(H);
				this->k_rest.push_back(K);
				this->n_rest.push_back(n_ijk);
			}
		}
	}
}

void stark::EnergyTriangleStrainWen23::_write_frame(stark::core::Stark &stark) {
	const std::string path = stark.get_frame_path("kl_") + ".vtk";
	const std::string path_edges = stark.get_frame_path("kl_edges_") + ".vtk";

	vtkio::VTKFile vtk_file;

	std::vector<Eigen::Vector3d> vertices;
	std::vector<std::array<int, 3>> triangles;

	std::vector<Eigen::Vector3d> v1;

	// Collect Tri3 vertices, triangles and v1
	for (int group = 0; group < this->triangles.size(); ++group) {
		const auto& set = this->point_sets[group];
		const int offset = (int)vertices.size();

		auto begin = set.get_global_index(0);
		auto end = begin + set.size();

		vertices.insert(vertices.end(), this->dyn->x1.data.begin() + begin, this->dyn->x1.data.begin() + end);
		v1.insert(v1.end(), this->dyn->v1.data.begin() + begin, this->dyn->v1.data.begin() + end);

		for (const std::array<int, 3>& tri : this->triangles[group]) {
			triangles.push_back({ tri[0] + offset, tri[1] + offset, tri[2] + offset });
		}
	}

	if (vertices.empty()) {
		vtk_file.write_empty(path);
	} else {
		vtk_file.set_points_from_twice_indexable(vertices);
		vtk_file.set_cells_from_twice_indexable(triangles, vtkio::CellType::Triangle);
		vtk_file.set_point_data_from_twice_indexable("v1", v1, vtkio::AttributeType::Vectors);
		vtk_file.set_cell_data_from_indexable("H", this->h_rest, vtkio::AttributeType::Scalars);
		vtk_file.set_cell_data_from_indexable("K", this->k_rest, vtkio::AttributeType::Scalars);
		vtk_file.set_cell_data_from_twice_indexable("n0", this->n_rest, vtkio::AttributeType::Vectors);

		vtk_file.write(path);
	}

	vtk_file = vtkio::VTKFile();

	std::vector<Eigen::Vector3d> midedge_points;
	std::vector<Eigen::Vector3d> midedge_normals;

	for (const auto& me : this->n_edge_rest) {
		const auto& set = this->point_sets[me.group];
		const auto& x_i = set.get_position(me.v0);
		const auto& x_j = set.get_position(me.v1);
		midedge_points.push_back(0.5*(x_i + x_j));
		midedge_normals.push_back(me.n);
	}

	if (vertices.empty()) {
		vtk_file.write_empty(path_edges);
	} else {
		vtk_file.set_points_from_twice_indexable(midedge_points);
		vtk_file.set_cells_as_particles(midedge_points.size());
		vtk_file.set_point_data_from_twice_indexable("n_edge", midedge_normals, vtkio::AttributeType::Vectors);

		vtk_file.write(path_edges);
	}
}

#pragma once

#include "../../../core/Stark.h"
#include "../PointDynamics.h"
#include "../../types.h"
#include "../../fem_types.h"
#include "../../fem_elements.h"
#include "../../fem_symx_types.h"
#include "../../../utils/include.h"
#include "../../matrix_wrapper.h"

namespace stark
{
class EnergyMicropolarShells {
public:
    /* Types */
    struct Params
    {
        STARK_PARAM_NON_NEGATIVE(double, thickness, 0.0)
        STARK_PARAM_YOUNGS_MODULUS()
        STARK_PARAM_POISSONS_RATIO()
        STARK_PARAM_NON_NEGATIVE(double, density, 0.0)
        STARK_PARAM_NON_NEGATIVE(double, angular_inertia_density, 0.0)
        STARK_PARAM_NON_NEGATIVE(double, inertia_damping, 0.0)
        STARK_PARAM_NON_NEGATIVE(double, mu_c_factor, 1.0)
        STARK_PARAM_NON_NEGATIVE(double, length_scale, 1e-2)
        STARK_PARAM_NON_NEGATIVE(double, shear_correction, 1.0)
        STARK_PARAM_NON_NEGATIVE(double, bending_correction, 1.0)
    };
    struct Handler {
        STARK_COMMON_HANDLER_CONTENTS(EnergyMicropolarShells, Params)

        /// Recomputes rest state quantities with the current positions of the shell
        void update_rest_configuration(bool reset_rotations = true);
        /// Interpolates the (curvature) quadrature points of all triangles in the mesh using the current nodal positions
        std::vector<Eigen::Vector3d> interpolate_quadrature_points_mesh() const;
        /// Get the rest deformation gradient (F0 or a) of the quadrature point with the given index (can only be called after the simulation was started)
        Eigen::Matrix3d get_rest_deformation_gradient(int quad_i) const;
        /// Sets the rest deformation gradient (F0 or a) of the quadrature point with the given index
        void set_rest_deformation_gradient(int quad_i, const Eigen::Matrix3d& F0);
        /// Sets the rest curvature of the quadrature point with the given index
        void set_rest_curvature(int quad_i, const Eigen::Matrix3d& curvature);
        /// Sets an angular velocity boundary condition for the vertices inside of the given AABB
        void prescribe_angular_velocity_inside_aabb(const Eigen::Vector3d& center, const Eigen::Vector3d& size, const Eigen::Vector3d& omega_bc, std::array<bool, 3> active = {true, true, true}, double stiffness = 1e3);
        /// Sets a target rotation boundary condition for the vertices inside of the given AABB
        int prescribe_rotation_inside_aabb(const Eigen::Vector3d& center, const Eigen::Vector3d& size, const Eigen::Matrix3d& rotation_bc, std::array<bool, 3> axis_active = {true, true, true}, double stiffness = 1e3);
        /// Updates the prescribed rotation for a boundary condition
        void set_prescribed_rotation(int bc_idx, const Eigen::Matrix3d& rotation_bc);
        /// Updates the prescribed rotation for a boundary condition
        void set_prescribed_rotation(int bc_idx, double angle_deg, const Eigen::Vector3d& axis);
        /// Returns the mesh of the displacement of the displacement degrees of freedom
        template<typename DisplacementElement>
        std::vector<std::array<int, DisplacementElement::num_nodes>> get_mesh() const {
            return this->model->displacement_meshes.get_mesh<DisplacementElement>(this->idx);
        }
    };

    Fem::MixedFemEnergyConnectivity<Fem::BasisTri3, Fem::BasisTri3> conn_disp_tri3_rot_tri3 = Fem::TriQuadrature::tri_p2();
    Fem::MixedFemEnergyConnectivity<Fem::BasisTri6, Fem::BasisTri3> conn_disp_tri6_rot_tri3 = Fem::TriQuadrature::tri_p2();
    Fem::MixedFemEnergyConnectivity<Fem::BasisTri6, Fem::BasisTri6> conn_disp_tri6_rot_tri6 = Fem::TriQuadrature::tri_p4();
    Fem::MixedFemEnergyConnectivity<Fem::BasisTri10, Fem::BasisTri3> conn_disp_tri10_rot_tri3 = Fem::TriQuadrature::tri_p6();
    Fem::MixedFemEnergyConnectivity<Fem::BasisTri10, Fem::BasisTri6> conn_disp_tri10_rot_tri6 = Fem::TriQuadrature::tri_p6();

    Fem::MixedFemEnergyConnectivity<Fem::BasisQuad4, Fem::BasisQuad4> conn_disp_quad4_rot_quad4 = Fem::QuadQuadrature::quad_p3();
    Fem::MixedFemEnergyConnectivity<Fem::BasisQuad9, Fem::BasisQuad4> conn_disp_quad9_rot_quad4 = Fem::QuadQuadrature::quad_p3();
    Fem::MixedFemEnergyConnectivity<Fem::BasisQuad9, Fem::BasisQuad9> conn_disp_quad9_rot_quad9 = Fem::QuadQuadrature::quad_p5();

    /// Flag to enable export of point cloud meshes of quadrature points with per quadrature point data
    bool WRITE_QUADRATURE_MESHES = false;
    /// Flag to enable export of subdivided meshes for visualization of curved meshes
    bool WRITE_SUBDIVIDED_MESHES = false;
    /// Number of subdivisions to apply for visualization meshes
    int MESH_SUBDIVISION_LEVEL = 3;

    struct ScalarFieldCallback {
        /// Name of the scalar field to be used in the VTK file
        std::string name;
        /// Callback mapping from rest state coordinates to scalar field values
        std::function<double(const Eigen::Vector3d&)> callback;
    };

    /// Callbacks for exporting scalar fields on the subdivided meshes
    std::vector<ScalarFieldCallback> scalar_field_callbacks;
private:
    /// Connectivity for lumped inertia energy for the displacement degrees of freedom
    symx::LabelledConnectivity<3> conn_lumped_disp_inertia {{"idx", "group", "dof" } };
    /// Connectivity for lumped inertia energy for the rotation degrees of freedom
    symx::LabelledConnectivity<3> conn_lumped_rot_inertia {{"idx", "group", "dof" } };

    /// Connectivity for the prescribed angular velocity penalty energy
    symx::LabelledConnectivity<3> conn_prescribed_ang_vel { { "idx", "group", "dof" } };
    /// Connectivity for the prescribed rotation matrix penalty energy
    symx::LabelledConnectivity<3> conn_prescribed_rot { { "idx", "group", "dof" } };

    const spPointDynamics dyn;

    std::vector<std::string> labels;
    std::vector<PointSetHandler> point_sets;

    symx::DoF micropolar_dof;
    std::vector<int> micropolar_dof_offsets;

    Fem::FemSurfaceMeshes displacement_meshes;
    Fem::FemSurfaceMeshes rotation_meshes;
    std::vector<std::vector<Eigen::Vector3d>> rotation_mesh_vertices;
    std::vector<stark::SubdividedMesh> subdivided_meshes;

    std::vector<QuaternionWrapperd> q0;     // Quaternions at time n
    std::vector<QuaternionWrapperd> q1;     // Quaternions at time n+1
    std::vector<Eigen::Vector3d> w0;        // Angular velocities at time n
    std::vector<Eigen::Vector3d> w1;        // Angular velocities at time n+1
    std::vector<Eigen::Vector3d> aa;        // Angular acceleration at time n+1
    std::vector<Eigen::Vector3d> t;         // Torque at time n+1

    std::vector<double> thickness;                  // per group
    std::vector<double> youngs_modulus;             // per group
    std::vector<double> poissons_ratio;             // per group
    std::vector<double> density;                    // per group
    std::vector<double> angular_inertia_density;    // per group
    std::vector<double> inertia_damping;            // per group
    std::vector<double> mu_c_factor;                // per group
    std::vector<double> length_scale;               // per group
    std::vector<double> shear_correction;           // per group
    std::vector<double> bending_correction;         // per group

    std::vector<double> nodal_lumped_area_disp;          // per displacement node
    std::vector<double> nodal_lumped_area_rot;           // per rotation node

    std::vector<Eigen::Vector3d> bc_target_ang_vel;         // Target angular velocity boundary condition
    std::vector<Eigen::Vector3d> bc_target_ang_vel_active;  // Active components of the boundary condition
    std::vector<double> bc_target_ang_vel_stiffness;        // Stiffness of the target angular velocity boundary condition

    std::vector<MatrixWrapper3d> bc_target_rot;             // Target rotation matrix boundary condition
    std::vector<Eigen::Vector3d> bc_target_rot_active;      // Active columns of the boundary condition
    std::vector<double> bc_target_rot_stiffness;            // Stiffness of the target rotation matrix boundary condition

    struct PerQuadraturePointData {
        /// Offset to the start of the quadrature data, per mesh
        std::vector<int> offsets;
        /// Orientation quaternions at time n, per quadrature point
        std::vector<QuaternionWrapperd> q0_quad;
        /// Rest state
        std::vector<MatrixWrapper3d> R_rest_quad;
        /// Rest state shape operator
        std::vector<MatrixWrapper2d> S_rest_quad;
        /// Rest state "First fundamental tensor" (not first fundamental form)
        std::vector<MatrixWrapper3d> a_rest_quad;
        /// Rest state "Second fundamental tensor" (not second fundamental form)
        std::vector<MatrixWrapper3d> b_rest_quad;
        /// Rest state curvature tensor, per quadrature point
        std::vector<MatrixWrapper3d> Gamma_rest_quad;
    };

    PerQuadraturePointData quadrature_data;

public:
    EnergyMicropolarShells(stark::core::Stark& stark, spPointDynamics dyn);

    /// Add a mesh for a micropolar shell. Only valid element combinations are explicitly instantiated in the source file.
    template<typename DisplacementElement, typename RotationElement>
    Handler add_mesh(const std::string& label, const PointSetHandler& set, const std::vector<std::array<int, DisplacementElement::num_nodes>>& primary_elements, const Mesh<RotationElement::num_nodes>& secondary_mesh, const Params& params) {
        static_assert(DisplacementElement::ref_dim == 2, "Displacement element type must be of reference dimension 2");
        static_assert(DisplacementElement::ref_dim == RotationElement::ref_dim, "Displacement and rotation elements must have the same reference dimension");
        static_assert(DisplacementElement::ref_shape == RotationElement::ref_shape, "Displacement and rotation elements must have the same reference shape");
        static_assert(DisplacementElement::num_nodes >= RotationElement::num_nodes, "Rotation discretization cannot have more DOF than the displacement discretization");
        return this->add_mesh_impl<DisplacementElement, RotationElement>(label, set, primary_elements, secondary_mesh, params);
    }

    Params get_params(const Handler& handler) const;
    void set_params(const Handler& handler, const Params& params);

private:
    void _initialize_fem_energies(stark::core::Stark& stark);
    void _initialize_bc_energies(stark::core::Stark& stark);
    void _initialize_rest_state_quantities(int group);

    // Stark callbacks
    void _before_simulation(stark::core::Stark& stark);
    void _before_time_step(stark::core::Stark& stark);
    void _after_time_step(stark::core::Stark& stark);
    void _write_frame(stark::core::Stark& stark);

    // Methods for handler
    void set_all_rotations(int group, const Eigen::Quaterniond& q);
    std::vector<Eigen::Vector3d> interpolate_quadrature_points_mesh(int group);
    Eigen::Matrix3d get_rest_deformation_gradient(int group, int quad_i) const;
    void set_rest_deformation_gradient(int group, int quad_i, const Eigen::Matrix3d& F0);
    void set_rest_curvature(int group, int quad_i, const Eigen::Matrix3d& curvature);
    void prescribe_angular_velocity_inside_aabb(int group, const Eigen::Vector3d& center, const Eigen::Vector3d& size, const Eigen::Vector3d& omega_bc, std::array<bool, 3> active, double stiffness);
    int prescribe_rotation_inside_aabb(int group, const Eigen::Vector3d& center, const Eigen::Vector3d& size, const Eigen::Matrix3d& rotation_bc, std::array<bool, 3> axis_active, double stiffness);
    void set_prescribed_rotation(int bc_idx, const Eigen::Matrix3d& rotation_bc);

    void add_bc_connectivity(symx::LabelledConnectivity<3>& bc_conn, int dof_group, int bc_group, const Eigen::AlignedBox3d& aabb);

    template<typename DisplacementElement, typename RotationElement>
    Handler add_mesh_impl(const std::string& label, const PointSetHandler& set, const std::vector<std::array<int, DisplacementElement::num_nodes>>& displacement_elements, const Mesh<RotationElement::num_nodes>& secondary_mesh, const Params& params);

    /// Returns whether the specified group is discretized with the given element types
    template<typename DisplacementElementT, typename RotationElementT>
    bool has_element_types(const int group) const {
        return this->displacement_meshes.element_type.at(group) == DisplacementElementT::element_type
            && this->rotation_meshes.element_type.at(group) == RotationElementT::element_type;
    }

    /// Run a function for each mixed FEM connectivity
    template<typename F>
    void for_each_connectivity(F&& f) {
        f(this->conn_disp_tri3_rot_tri3);
        f(this->conn_disp_tri6_rot_tri3);
        f(this->conn_disp_tri6_rot_tri6);
        f(this->conn_disp_tri10_rot_tri3);
        f(this->conn_disp_tri10_rot_tri6);
        f(this->conn_disp_quad4_rot_quad4);
        f(this->conn_disp_quad9_rot_quad4);
        f(this->conn_disp_quad9_rot_quad9);
    }

    /// Converts local micropolar DOF indices of a group to global micropolar DOF indices
    template<std::size_t N>
    std::array<int, N> micropolar_get_global_indices(const int group, const std::array<int, N>& local_indices)
    {
        std::array<int, N> global_indices;
        const int dof_start = this->micropolar_dof_offsets.at(group);
        for (int i = 0; i < N; ++i) {
            global_indices[i] = dof_start + local_indices[i];
        }
        return global_indices;
    }

};
}

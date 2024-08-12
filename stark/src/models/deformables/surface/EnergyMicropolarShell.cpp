#include "EnergyMicropolarShell.h"

#include <fmt/format.h>

#include "../../symx_quaternion.h"

/// Computes `(G^T*G).det().sqrt()`
const auto form = [](const symx::Matrix& G) -> symx::Scalar {
    return (G.transpose()*G).det().sqrt();
};

/// Computes `(G^T*G)^-1*G^T`, the pseudo-inverse of G
const auto pseudo_inv = [](const symx::Matrix& G) -> symx::Matrix {
    return (G.transpose()*G).inv() * G.transpose();
};

/// Returns the skew-symmetric part of the matrix
const auto sym = [](const symx::Matrix& m) -> symx::Matrix {
    return 0.5 * (m + m.transpose());
};

/// Returns the skew-symmetric part of the matrix
const auto skew = [](const symx::Matrix& m) -> symx::Matrix {
    return 0.5 * (m - m.transpose());
};

/// Returns the axial vector of a skew-symmetric matrix such that `cross(axl(A),x) = A*x` for all `x`
const auto axl = [](const symx::Matrix& m) -> symx::Vector {
    const symx::Scalar m_23 = m(1,2);
    const symx::Scalar m_31 = m(2,0);
    const symx::Scalar m_12 = m(0,1);
    // Source paper
    //return symx::Vector({m_23, m_31, m_12});
    // Corrected
    return symx::Vector({-m_23, -m_31, -m_12});
};

static Eigen::Matrix3d polar_decomposition_3d(Eigen::Matrix3d& F)
{
    // Eigen decomposition
    Eigen::Matrix3d FTF = F.transpose() * F; // U^2
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(FTF);
    Eigen::Matrix3d U_ = eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();
    Eigen::Matrix3d Q = eigen_solver.eigenvectors();

    // Check for reflection
    bool reflection = F.determinant() < 0.0;
    if (reflection) {
        // Find the smallest eigenvalue
        int min_eig_value_idx = -1;
        double min_eig_value = std::numeric_limits<double>::max();
        for (int i = 0; i < 3; i++) {
            if (U_(i, i) < min_eig_value) {
                min_eig_value = U_(i, i);
                min_eig_value_idx = i;
            }
        }

        // Swap
        Q.block(0, min_eig_value_idx, 3, 1) *= -1;
        U_(min_eig_value_idx, min_eig_value_idx) *= -1;
    }

    Eigen::Matrix3d U = Q * U_ * Q.transpose();
    Eigen::Matrix3d R = F * U.inverse();
    return R;
}

stark::EnergyMicropolarShells::EnergyMicropolarShells(stark::core::Stark& stark, spPointDynamics dyn)
        : dyn(dyn)
{
    this->micropolar_dof_offsets.push_back(0);
    this->quadrature_data.offsets.push_back(0);
    this->micropolar_dof = stark.global_energy.add_dof_array(this->w1, "EnergyTriangleMicropolarShells.w1");

    stark.callbacks.add_before_simulation([&]() { this->_before_simulation(stark); });
    stark.callbacks.add_before_time_step([&]() { this->_before_time_step(stark); });
    stark.callbacks.add_after_time_step([&]() { this->_after_time_step(stark); });
    stark.callbacks.add_write_frame([&]() { this->_write_frame(stark); });

    this->_initialize_fem_energies(stark);
    this->_initialize_bc_energies(stark);
}

void stark::EnergyMicropolarShells::_initialize_fem_energies(stark::core::Stark& stark)
{
    // TODO: Implement consistent mass matrix with separate quadrature

    stark.global_energy.add_energy("EnergyMpShellLumpedInertia", this->conn_lumped_disp_inertia, [&](symx::Energy& energy, symx::Element& NODE)
    {
        // Create symbols
        const symx::Vector v1 = energy.make_dof_vector(this->dyn->dof, this->dyn->v1.data, NODE["dof"]);
        const symx::Vector v0 = energy.make_vector(this->dyn->v0.data, NODE["dof"]);

        const symx::Vector a = energy.make_vector(this->dyn->a.data, NODE["dof"]);
        const symx::Vector f = energy.make_vector(this->dyn->f.data, NODE["dof"]);

        const symx::Scalar area = energy.make_scalar(this->nodal_lumped_area_disp, NODE["dof"]);
        const symx::Scalar h = energy.make_scalar(this->thickness, NODE["group"]);
        const symx::Scalar density = energy.make_scalar(this->density, NODE["group"]);
        const symx::Scalar damping = energy.make_scalar(this->inertia_damping, NODE["group"]);

        const symx::Scalar dt = energy.make_scalar(stark.dt);
        const symx::Vector gravity = energy.make_vector(stark.gravity);

        symx::Scalar e_tot = energy.make_zero();
        // External accelerations and forces
        {
            symx::Scalar mass = area * h * density;
            symx::Vector delta_x = dt * v1;
            e_tot += -mass * (gravity + a).dot(delta_x) - f.dot(delta_x);
        }
        // Classic kinetic energy
        {
            symx::Scalar mass = area * h * density;
            symx::Vector delta_v = v1 - v0;
            e_tot += 0.5 * mass * delta_v.squared_norm();
        }
        // Damping energy
        {
            symx::Scalar mass = area * h * density;
            e_tot += 0.5 * mass * dt * damping * v1.squared_norm();
        }
        energy.set(e_tot);
    });

    stark.global_energy.add_energy("EnergyMpShellLumpedAngularInertia", this->conn_lumped_rot_inertia, [&](symx::Energy& energy, symx::Element& NODE)
    {
        // Create symbols
        const symx::Vector w1 = energy.make_dof_vector(this->micropolar_dof, this->w1, NODE["dof"]);
        const symx::Vector w0 = energy.make_vector(this->w0, NODE["dof"]);

        const symx::Vector aa = energy.make_vector(this->aa, NODE["dof"]);
        const symx::Vector t = energy.make_vector(this->t, NODE["dof"]);

        const symx::Scalar area = energy.make_scalar(this->nodal_lumped_area_rot, NODE["dof"]);
        const symx::Scalar h = energy.make_scalar(this->thickness, NODE["group"]);
        const symx::Scalar angular_inertia_density = energy.make_scalar(this->angular_inertia_density, NODE["group"]);

        const symx::Scalar dt = energy.make_scalar(stark.dt);
        const symx::Vector gravity = energy.make_vector(stark.gravity);

        symx::Scalar e_tot = energy.make_zero();
        // External accelerations and torques
        {
            symx::Scalar mass = area * h * angular_inertia_density;
            symx::Vector delta_x = dt * w1;
            e_tot += -mass * aa.dot(delta_x) - t.dot(delta_x);
        }
        // Rotational inertia energy
        {
            symx::Scalar mass = area * h * angular_inertia_density;
            symx::Vector delta_w = w1 - w0;
            e_tot += 0.5 * mass * delta_w.squared_norm();
        }
        energy.set(e_tot);
    });

    // Split strain energy in individual terms to avoid ballooning compile times
    enum StrainEnergyTerm {
        FlatMembrane,
        FlatBending,
        FlatCurvature,
        FlatMembraneNoVol,
        FlatNonlinVol,
        CurvedMembrane1,
        CurvedMembrane2,
        CurvedMembrane3,
        CurvedMembrane4,
        CurvedBending1,
        CurvedBending2,
        CurvedBending3,
    };

    auto make_strain_energy = [&]<typename D, typename R>(StrainEnergyTerm term, Fem::MixedFemEnergyConnectivity<D, R>& CONNECTIVITY)
    {
        using DispEl = Fem::Element<D>;
        using RotEl = Fem::Element<R>;
        using Conn = Fem::MixedFemEnergyConnectivity<D, R>;

        return [&, term](symx::Energy& energy, symx::Element& ELE_WITH_QUAD)
        {
            energy.set_never_project_to_PD(stark.settings.models.never_project_mp_shell);

            // Unpack connectivity
            const std::vector<symx::Index> displacement_dof = Conn::element_primary_dof_slice(ELE_WITH_QUAD);
            const std::vector<symx::Index> rotation_dof = Conn::element_secondary_dof_slice(ELE_WITH_QUAD);

            // Create symbols
            const std::vector<symx::Vector> v1 = energy.make_dof_vectors(this->dyn->dof, this->dyn->v1.data, displacement_dof);
            const std::vector<symx::Vector> x0 = energy.make_vectors(this->dyn->x0.data, displacement_dof);
            const std::vector<symx::Vector> X = energy.make_vectors(this->dyn->X.data, displacement_dof);

            const std::vector<symx::Vector> w1 = energy.make_dof_vectors(this->micropolar_dof, this->w1, rotation_dof);
            const std::vector<Quaternion> q0 = make_quaternions(energy, this->q0, rotation_dof);

            const symx::Scalar dt = energy.make_scalar(stark.dt);

            const symx::Scalar h = energy.make_scalar(this->thickness, ELE_WITH_QUAD["group"]);
            const symx::Scalar youngs_modulus = energy.make_scalar(this->youngs_modulus, ELE_WITH_QUAD["group"]);
            const symx::Scalar poissons_ratio = energy.make_scalar(this->poissons_ratio, ELE_WITH_QUAD["group"]);
            const symx::Scalar mu_c_factor = energy.make_scalar(this->mu_c_factor, ELE_WITH_QUAD["group"]);
            const symx::Scalar L_c = energy.make_scalar(this->length_scale, ELE_WITH_QUAD["group"]);
            const symx::Scalar bending_correction = energy.make_scalar(this->bending_correction, ELE_WITH_QUAD["group"]);

            const auto E = youngs_modulus;
            const auto nu = poissons_ratio;
            const symx::Scalar mu = E / (2.0 * (1.0 + nu));
            const symx::Scalar lambda = (E * nu) / ( (1.0 + nu) * (1.0 - 2.0 * nu) );
            const symx::Scalar mu_c = mu * mu_c_factor;

            const symx::Vector q_point = energy.make_vector(CONNECTIVITY.quadrature_points, ELE_WITH_QUAD["q_idx_local"]);
            const symx::Scalar q_weight = energy.make_scalar(CONNECTIVITY.quadrature_weights, ELE_WITH_QUAD["q_idx_local"]);

            // Time integration of positions
            std::vector<symx::Vector> x1;
            for (int i = 0; i < x0.size(); i++) {
                x1.push_back(x0[i] + dt * v1[i]);
            }
            // Time integration of quaternions q_i (at nodes/midpoints)
            std::vector<Quaternion> q1_i;
            for (int i = 0; i < w1.size(); i++) {
                q1_i.push_back(q0[i].time_integration_with_global_w(w1[i], dt).normalized());
            }
            // Rotation matrices R_i (at nodes)
            std::vector<symx::Matrix> R1_i;
            for (int i = 0; i < w1.size(); i++) {
                R1_i.push_back(q1_i[i].to_rotation(energy));
            }

            {
                const symx::Matrix X_mat = symx::Matrix(
                        symx::gather(X), { (int)X.size(), 3 }).transpose();     // 3xN
                const symx::Matrix x1_mat = symx::Matrix(
                        symx::gather(x1), { (int)x1.size(), 3 }).transpose();   // 3xN

                const symx::Matrix dphidxi = DispEl::grad_symbolic(q_point).transpose();     // Nx2
                const symx::Matrix G = X_mat * dphidxi;										// 3x2
                const symx::Matrix g = x1_mat * dphidxi;									// 3x2
                const symx::Matrix Ginv = pseudo_inv(G);									// 2x3
                const symx::Matrix F = g * Ginv;											// 3x3

                symx::Matrix R1_quad = energy.make_identity_matrix(3);
                Quaternion q1_quad = Quaternion::from_vec4(energy.make_zero_vector(4));
                // Interpolation of rotation at quadrature point
                if (!stark.settings.models.enable_model_mp_shell_use_corot_rotation_interpolation) {
                    const Quaternion q0_quad = make_quaternion(energy, this->quadrature_data.q0_quad, ELE_WITH_QUAD["q_idx_global"]);
                    // Interpolate omega to quadrature point
                    const symx::Vector w1_quad = RotEl::interpolate(w1, q_point);
                    // Time integration at quadrature point
                    q1_quad = q0_quad.time_integration_with_global_w(w1_quad, dt).normalized();
                    R1_quad = q1_quad.to_rotation(energy);
                } else {
                    std::vector<Quaternion> q1;
                    for (int i = 0; i < w1.size(); ++i) {
                        q1.push_back(q0[i].time_integration_with_global_w(w1[i], dt).normalized());
                    }
                    // Make orientations relative to the first node
                    std::vector<Quaternion> q1_local;
                    for (int i = 0; i < w1.size(); ++i) {
                        q1_local.push_back(q1[0].inverse() * q1[i]);
                    }
                    // Convert to 4-component vectors
                    std::vector<symx::Vector> q1_local_vec;
                    for (int i = 0; i < w1.size(); ++i) {
                        q1_local_vec.push_back(q1_local[i].coeffs);
                    }
                    // Interpolate vectors to quadrature point
                    const symx::Vector q1_local_vec_quad = RotEl::interpolate(q1_local_vec, q_point);
                    const Quaternion q1_local_quad = Quaternion::from_vec4(q1_local_vec_quad).normalized();
                    // Recover absolute orientations
                    q1_quad = q1[0] * q1_local_quad;
                    R1_quad = q1_quad.to_rotation(energy);
                }

                const symx::Matrix dphidxi_rot = RotEl::grad_symbolic(q_point).transpose();   // Nx2
                symx::Matrix Gamma = energy.make_zero_matrix({3, 3});

                if (!stark.settings.models.enable_model_mp_shell_use_quaternion_gamma) {
                    symx::Matrix grad_R_alpha = energy.make_zero_matrix({3, 3});
                    symx::Matrix grad_R_beta = energy.make_zero_matrix({3, 3});
                    for (int i = 0; i < w1.size(); ++i) {
                        grad_R_alpha += R1_i[i] * dphidxi_rot(i, 0);
                        grad_R_beta += R1_i[i] * dphidxi_rot(i, 1);
                    }

                    const symx::Vector Gamma_alpha = axl(R1_quad.transpose() * grad_R_alpha);
                    const symx::Vector Gamma_beta = axl(R1_quad.transpose() * grad_R_beta);

                    symx::Matrix Gamma_omega = energy.make_zero_matrix({3, 2});
                    Gamma_omega.set_col(0, Gamma_alpha);
                    Gamma_omega.set_col(1, Gamma_beta);

                    Gamma = Gamma_omega * Ginv;
                } else {
                    Quaternion grad_q_alpha = Quaternion::from_vec4(energy.make_zero_vector(4));
                    Quaternion grad_q_beta = Quaternion::from_vec4(energy.make_zero_vector(4));

                    for (int i = 0; i < w1.size(); ++i) {
                        grad_q_alpha = grad_q_alpha + q1_i[i] * dphidxi_rot(i, 0);
                        grad_q_beta = grad_q_beta + q1_i[i] * dphidxi_rot(i, 1);
                    }

                    const symx::Vector Gamma_alpha = 2.0 * (q1_quad.inverse() * grad_q_alpha).to_vec3(energy);
                    const symx::Vector Gamma_beta = 2.0 * (q1_quad.inverse() * grad_q_beta).to_vec3(energy);

                    symx::Matrix Gamma_omega = energy.make_zero_matrix({3, 2});
                    Gamma_omega.set_col(0, Gamma_alpha);
                    Gamma_omega.set_col(1, Gamma_beta);

                    Gamma = Gamma_omega * Ginv;
                }

                const symx::Matrix Gamma_rest = energy.make_matrix(this->quadrature_data.Gamma_rest_quad, {3, 3}, ELE_WITH_QUAD["q_idx_global"]);
                Gamma = Gamma - Gamma_rest;

                //const symx::Matrix a = G * Ginv;
                const symx::Matrix a = energy.make_matrix(this->quadrature_data.a_rest_quad, {3, 3}, ELE_WITH_QUAD["q_idx_global"]);
                const symx::Matrix b = energy.make_matrix(this->quadrature_data.b_rest_quad, {3, 3}, ELE_WITH_QUAD["q_idx_global"]);

                const symx::Vector a1 = G.col(0);
                const symx::Vector a2 = G.col(1);
                const symx::Matrix c = form(G).inv() * (outer(a1, a2) - outer(a2, a1));

                const symx::Matrix E = R1_quad.transpose() * F - a;

                const symx::Scalar H = 0.5 * b.trace();
                const symx::Scalar K = 0.5 * (b.trace()*b.trace() - (b*b).trace());

                const auto W_m_bilinear = [&mu, &mu_c, &lambda](const symx::Matrix& S, const symx::Matrix& T) {
                    return mu * sym(S).double_dot(sym(T)) + mu_c * skew(S).double_dot(skew(T)) + ((lambda * mu) / (lambda + 2.0 * mu)) * S.trace() * T.trace();
                };
                const auto W_m_quadratic = [&mu, &mu_c, &lambda](const symx::Matrix& S) {
                    return mu * sym(S).frobenius_norm_sq() + mu_c * skew(S).frobenius_norm_sq() + ((lambda * mu) / (lambda + 2.0 * mu)) * S.trace().powN(2);
                };
                const auto W_mp = [&mu, &mu_c, &lambda](const symx::Matrix& S) {
                    return mu * sym(S).frobenius_norm_sq() + mu_c * skew(S).frobenius_norm_sq() + (lambda / 2.0) * S.trace().powN(2);
                };
                const auto W_curv = [&mu, &L_c](const symx::Matrix& S) {
                    return mu * L_c.powN(2) * S.frobenius_norm_sq();
                };
                const auto W_m_quadratic_no_vol = [&mu, &mu_c, &lambda](const symx::Matrix& S) {
                    return mu * sym(S).frobenius_norm_sq() + mu_c * skew(S).frobenius_norm_sq();
                };

                symx::Scalar e = energy.make_zero();

                // Terms of the material model, see:
                // "A geometrically nonlinear Cosserat shell model for orientable and non-orientable surfaces: Discretization with geometric finite elements"
                // by Nebel, Sander, Bîrsan & Neff, 2023
                {
                    // Terms for the reduced shell model for flat rest configurations
                    if (term == StrainEnergyTerm::FlatMembrane)
                        e += h * W_m_quadratic(E);
                    if (term == StrainEnergyTerm::FlatBending)
                        e += bending_correction * h.powN(3) / 12.0 * W_m_quadratic(c * Gamma);
                    if (term == StrainEnergyTerm::FlatCurvature)
                        e += h * mu * L_c.powN(2) * Gamma.frobenius_norm_sq();

                    // Nonlinear volume terms (measuring true volume change)
                    if (term == StrainEnergyTerm::FlatMembraneNoVol)
                        e += h * W_m_quadratic_no_vol(E);
                    if (term == StrainEnergyTerm::FlatNonlinVol) {
                        const symx::Matrix gtg = g.transpose() * g;
                        const symx::Matrix GtG = G.transpose() * G;
                        const symx::Scalar J = gtg.det().sqrt() / GtG.det().sqrt();
                        e += h * 0.5 * ((mu * lambda) / (2 * mu + lambda)) * ((J - 1).powN(2) + (J.inv() - 1).powN(2));
                    }

                    // Terms for the full shell model with curved rest configurations
                    if (term == StrainEnergyTerm::CurvedMembrane1)
                        e += (h - K * (h.powN(3) / 12.0)) * W_m_quadratic(E);
                    if (term == StrainEnergyTerm::CurvedMembrane2)
                        e += ((h.powN(3) / 12.0) - K * (h.powN(5) / 80.0)) * W_m_quadratic(E * b + c * Gamma);
                    if (term == StrainEnergyTerm::CurvedMembrane3)
                        e += (h.powN(3) / 6.0) * W_m_bilinear(E, c * Gamma * b - 2.0 * H * c * Gamma);
                    if (term == StrainEnergyTerm::CurvedMembrane4)
                        e += (h.powN(5) / 80.0) * W_mp(E * b * b + c * Gamma * b);
                    if (term == StrainEnergyTerm::CurvedBending1)
                        e += (h - K * (h.powN(3) / 12.0)) * W_curv(Gamma);
                    if (term == StrainEnergyTerm::CurvedBending2)
                        e += ((h.powN(3) / 12.0) - K * (h.powN(5) / 80.0)) * W_curv(Gamma * b);
                    if (term == StrainEnergyTerm::CurvedBending3)
                        e += (h.powN(5) / 80.0) * W_curv(Gamma * b * b);
                }

                symx::Scalar E_tot = e * q_weight * form(G) * DispEl::ref_area;
                energy.set(E_tot);
            }
        };
    };

    const auto add_energy_per_discretization = [&](const std::string& name, auto&& energy) {
        stark.global_energy.add_energy(fmt::format("EnergyMpShellTri3Tri3{}", name), this->conn_disp_tri3_rot_tri3.connectivity, energy(this->conn_disp_tri3_rot_tri3));
        stark.global_energy.add_energy(fmt::format("EnergyMpShellTri6Tri3{}", name), this->conn_disp_tri6_rot_tri3.connectivity, energy(this->conn_disp_tri6_rot_tri3));
        stark.global_energy.add_energy(fmt::format("EnergyMpShellTri6Tri6{}", name), this->conn_disp_tri6_rot_tri6.connectivity, energy(this->conn_disp_tri6_rot_tri6));
        stark.global_energy.add_energy(fmt::format("EnergyMpShellTri10Tri3{}", name), this->conn_disp_tri10_rot_tri3.connectivity, energy(this->conn_disp_tri10_rot_tri3));
        stark.global_energy.add_energy(fmt::format("EnergyMpShellTri10Tri6{}", name), this->conn_disp_tri10_rot_tri6.connectivity, energy(this->conn_disp_tri10_rot_tri6));
        stark.global_energy.add_energy(fmt::format("EnergyMpShellQuad4Quad4{}", name), this->conn_disp_quad4_rot_quad4.connectivity, energy(this->conn_disp_quad4_rot_quad4));
        stark.global_energy.add_energy(fmt::format("EnergyMpShellQuad9Quad4{}", name), this->conn_disp_quad9_rot_quad4.connectivity, energy(this->conn_disp_quad9_rot_quad4));
        stark.global_energy.add_energy(fmt::format("EnergyMpShellQuad9Quad9{}", name), this->conn_disp_quad9_rot_quad9.connectivity, energy(this->conn_disp_quad9_rot_quad9));
    };

    const std::string gamma_label = stark.settings.models.enable_model_mp_shell_use_quaternion_gamma ? "QuatGamma" : "MatGamma";
    const std::string corot_label = stark.settings.models.enable_model_mp_shell_use_corot_rotation_interpolation ? "Corot" : "";

    if (!stark.settings.models.enable_model_mp_shell_full_rest_curvature_terms) {
        if (!stark.settings.models.enable_model_mp_shell_use_nonlinear_volume_terms) {
            add_energy_per_discretization(fmt::format("{}{}StrainFlatMembrane", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::FlatMembrane, conn); });
            add_energy_per_discretization(fmt::format("{}{}StrainFlatBending", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::FlatBending, conn); });
            add_energy_per_discretization(fmt::format("{}{}StrainFlatCurvature", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::FlatCurvature, conn); });
        } else {
            add_energy_per_discretization(fmt::format("{}{}StrainFlatMembraneNoVol", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::FlatMembraneNoVol, conn); });
            add_energy_per_discretization(fmt::format("{}{}StrainFlatNonlinVol", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::FlatNonlinVol, conn); });
            add_energy_per_discretization(fmt::format("{}{}StrainFlatBending", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::FlatBending, conn); });
            add_energy_per_discretization(fmt::format("{}{}StrainFlatCurvature", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::FlatCurvature, conn); });
        }
    } else {
        add_energy_per_discretization(fmt::format("{}{}StrainCurvedMembrane1", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::CurvedMembrane1, conn); });
        add_energy_per_discretization(fmt::format("{}{}StrainCurvedMembrane2", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::CurvedMembrane2, conn); });
        add_energy_per_discretization(fmt::format("{}{}StrainCurvedMembrane3", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::CurvedMembrane3, conn); });
        add_energy_per_discretization(fmt::format("{}{}StrainCurvedMembrane4", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::CurvedMembrane4, conn); });
        add_energy_per_discretization(fmt::format("{}{}StrainCurvedBending1", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::CurvedMembrane1, conn); });
        add_energy_per_discretization(fmt::format("{}{}StrainCurvedBending2", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::CurvedMembrane2, conn); });
        add_energy_per_discretization(fmt::format("{}{}StrainCurvedBending3", gamma_label, corot_label), [&](auto& conn) { return make_strain_energy(StrainEnergyTerm::CurvedMembrane3, conn); });
    }
}

void stark::EnergyMicropolarShells::_initialize_bc_energies(stark::core::Stark& stark)
{
    const auto prescribed_ang_vel_energy = [&](symx::Energy& energy, symx::Element& VERTEX_BC)
    {
        // Create symbols
        symx::Vector w1 = energy.make_dof_vector(this->micropolar_dof, this->w1, VERTEX_BC["dof"]);
        symx::Vector w1_prescribed = energy.make_vector(this->bc_target_ang_vel, VERTEX_BC["group"]);
        symx::Vector prescribed_active = energy.make_vector(this->bc_target_ang_vel_active, VERTEX_BC["group"]);
        symx::Scalar k = energy.make_scalar(this->bc_target_ang_vel_stiffness, VERTEX_BC["group"]);

        // Energy
        symx::Scalar E = energy.make_zero();
        E += prescribed_active[0] * (w1[0] - w1_prescribed[0]).powN(2);
        E += prescribed_active[1] * (w1[1] - w1_prescribed[1]).powN(2);
        E += prescribed_active[2] * (w1[2] - w1_prescribed[2]).powN(2);
        E *= 0.5 * k;
        energy.set(E);
    };

    const auto prescribed_rot_energy = [&](symx::Energy& energy, symx::Element& VERTEX_BC) {
        // Create symbols
        const symx::Vector w1 = energy.make_dof_vector(this->micropolar_dof, this->w1, VERTEX_BC["dof"]);
        const symx::Matrix R1_prescribed = energy.make_matrix(this->bc_target_rot, {3, 3}, VERTEX_BC["group"]);
        const symx::Vector prescribed_active = energy.make_vector(this->bc_target_rot_active, VERTEX_BC["group"]);
        const symx::Scalar k = energy.make_scalar(this->bc_target_rot_stiffness, VERTEX_BC["group"]);
        const symx::Scalar dt = energy.make_scalar(stark.dt);

        const symx::Vector q0_vec = energy.make_vector(this->q0, VERTEX_BC["dof"]);
        const Quaternion q0 = Quaternion::from_vec4(q0_vec);
        const Quaternion q1 = q0.time_integration_with_global_w(w1, dt).normalized();
        const symx::Matrix R1 = q1.to_rotation(energy);

        // Energy
        symx::Scalar E = energy.make_zero();
        E += prescribed_active[0] * (R1.col(0) - R1_prescribed.col(0)).squared_norm();
        E += prescribed_active[1] * (R1.col(1) - R1_prescribed.col(1)).squared_norm();
        E += prescribed_active[2] * (R1.col(2) - R1_prescribed.col(2)).squared_norm();
        E *= 0.5 * k;
        energy.set(E);
    };

    stark.global_energy.add_energy("EnergyMpShellPrescribedAngVel", this->conn_prescribed_ang_vel, prescribed_ang_vel_energy);
    stark.global_energy.add_energy("EnergyMpShellPrescribedRotation", this->conn_prescribed_rot, prescribed_rot_energy);
}

void stark::EnergyMicropolarShells::_before_simulation(stark::core::Stark &stark)
{
    // Initialize rest rotations and rest tensors at quadrature points
    for (int group = 0; group < this->labels.size(); ++group) {
        this->_initialize_rest_state_quantities(group);
    }

    // Print information about the micropolar energy
    auto print = [&](const std::string& string) {
        stark.console.print(string, stark::core::ConsoleVerbosity::TimeSteps);
    };

    print(fmt::format("\n---\nMicropolar Shell model\n"));

    {
        print(fmt::format("General settings:\n"));
        print(fmt::format("\tenable_model_mp_shell: {}\n", stark.settings.models.enable_model_mp_shell));
        print(fmt::format("\tenable_model_mp_shell_full_rest_curvature_terms: {}\n", stark.settings.models.enable_model_mp_shell_full_rest_curvature_terms));
        print(fmt::format("\tenable_model_mp_shell_use_corot_rotation_interpolation: {}\n", stark.settings.models.enable_model_mp_shell_use_corot_rotation_interpolation));
        print(fmt::format("\tenable_model_mp_shell_use_quaternion_gamma: {}\n", stark.settings.models.enable_model_mp_shell_use_quaternion_gamma));
        print(fmt::format("\tenable_model_mp_shell_use_nonlinear_volume_terms: {}\n", stark.settings.models.enable_model_mp_shell_use_nonlinear_volume_terms));
        print(fmt::format("\tnever_project_mp_shell: {}\n", stark.settings.models.never_project_mp_shell));
        print(fmt::format("\tWRITE_QUADRATURE_MESHES: {}\n", this->WRITE_QUADRATURE_MESHES));
        print(fmt::format("\tWRITE_SUBDIVIDED_MESHES: {}\n", this->WRITE_SUBDIVIDED_MESHES));
        print(fmt::format("\tMESH_SUBDIVISION_LEVEL: {}\n", this->MESH_SUBDIVISION_LEVEL));
    }

    // Print quadrature information
    {
        print(fmt::format("Discretization and quadrature:\n"));
        auto print_quadrature_name = [&]<typename D, typename R>(const Fem::MixedFemEnergyConnectivity<D, R>& connectivity) {
            print(fmt::format("\tdisplacement: {}, rotation: {} -> quadrature: {} ({} points)\n", Fem::element_type_to_string(D::element_type), Fem::element_type_to_string(R::element_type), connectivity.quadrature.name, connectivity.quadrature_points.size()));
        };

        this->for_each_connectivity(print_quadrature_name);
    }

    // Print information about all shell objects
    {
        print(fmt::format("Meshes:\n"));
        for (int group = 0; group < this->labels.size(); ++group) {
            const auto& set = this->point_sets.at(group);
            int n_micropolar_dofs = this->micropolar_dof_offsets.at(group+1) - this->micropolar_dof_offsets.at(group+0);
            const std::string type_displacement = Fem::element_type_to_string(this->displacement_meshes.element_type.at(group));
            const std::string type_rotation = Fem::element_type_to_string(this->rotation_meshes.element_type.at(group));
            int n_elements_displacement = this->displacement_meshes.num_elements(group);
            int n_elements_rotation = this->rotation_meshes.num_elements(group);
            print(fmt::format("\tMesh #{}: {} displacement elements ({}, {} total nodes), {} rotation elements ({}, {} total nodes)\n", group, n_elements_displacement, type_displacement, set.size(), n_elements_rotation, type_rotation, n_micropolar_dofs));

            print(fmt::format("\t\tthickness: {:.6e}\n", this->thickness.at(group)));
            print(fmt::format("\t\tdensity: {:.6e}, angular inertia density: {:.6e}, inertia damping: {:.6e}\n", this->density.at(group), this->angular_inertia_density.at(group), this->inertia_damping.at(group)));
            print(fmt::format("\t\tyoungs modulus: {:.6e}, poissons ratio: {:.6e}\n", this->youngs_modulus.at(group), this->poissons_ratio.at(group)));
            print(fmt::format("\t\tmu_c_factor: {:.6e}, L_c: {:.6e}\n", this->mu_c_factor.at(group), this->length_scale.at(group)));
            print(fmt::format("\t\tshear_correction: {:.6e}, bending_correction: {:.6e}\n", this->shear_correction.at(group), this->bending_correction.at(group)));
        }
        print(fmt::format("---\n\n"));
    }
}

void stark::EnergyMicropolarShells::_before_time_step(stark::core::Stark& stark)
{
    // Reset initial guess for angular velocities
    std::fill(this->w1.begin(), this->w1.end(), Eigen::Vector3d::Zero());
}

void stark::EnergyMicropolarShells::_after_time_step(stark::core::Stark& stark)
{
    // Reset angular velocities for constrained nodes
    for (int i = 0; i < this->conn_prescribed_ang_vel.size(); ++i) {
        const int bc_group_i = this->conn_prescribed_ang_vel[i][1];
        const int dof_i = this->conn_prescribed_ang_vel[i][2];
        
        for (int j = 0; j < 3; ++j) {
            this->w1[dof_i][j] = this->bc_target_ang_vel_active[bc_group_i][j] ? this->bc_target_ang_vel[bc_group_i][j] : this->w1[dof_i][j];
        }
    }

    // Important to use the same integration rule as in energy to avoid non-equilibrium state for next time step
    const auto quaternion_time_integration_global_w = [&](const Eigen::Quaterniond& q_start, const Eigen::Vector3d& w_glob, const double dt, const bool normalize) -> Eigen::Quaterniond {
        // https://stackoverflow.com/questions/46908345/integrate-angular-velocity-as-quaternion-rotation?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
        Eigen::Quaterniond w_ = Eigen::Quaterniond(0, w_glob.x(), w_glob.y(), w_glob.z());
        Eigen::Quaterniond q_end = q_start;
        q_end.coeffs() += 0.5 * dt * (w_ * q_start).coeffs();
        if (normalize) {
            q_end.normalize();
        }
        return q_end;
    };

    // Time integration of angular velocity DOFs
    const double dt = stark.dt;
    for (int i = 0; i < this->q1.size(); i++) {
        const Eigen::Quaterniond q_start = this->q1[i];
        const Eigen::Quaterniond q_end = quaternion_time_integration_global_w(q_start, this->w1[i], dt, true);
        this->q1[i] = q_end;
    }

    // Time integration of angular velocities at quadrature points
    auto update_quadrature_point_rotations = [&]<typename D, typename R>(const Fem::MixedFemEnergyConnectivity<D, R>& connectivity)
    {
        using DispEl = Fem::Element<D>;
        using RotEl = Fem::Element<R>;

        for (int group = 0; group < this->labels.size(); ++group) {
            // Skip meshes of different element types
            if (!this->has_element_types<DispEl, RotEl>(group)) {
                continue;
            }

            const int micropolar_dof_offset = this->micropolar_dof_offsets.at(group);
            const int quadrature_data_offset = this->quadrature_data.offsets.at(group);

            const auto& rot_mesh = this->rotation_meshes.get_mesh<RotEl>(group);
            for (int ele_i = 0; ele_i < rot_mesh.size(); ++ele_i) {
                const std::array<int, RotEl::num_nodes> ele = rot_mesh[ele_i];
                // Collect nodal angular velocities
                std::array<Eigen::Vector3d, RotEl::num_nodes> w1;
                for (int i = 0; i < RotEl::num_nodes; ++i) {
                    w1[i] = this->w1[micropolar_dof_offset + ele[i]];
                }
                // Interpolate to quadrature point and perform time integration
                const int num_quadrature_points = connectivity.quadrature_points.size();
                for (int i = 0; i < num_quadrature_points; ++i) {
                    const Eigen::Vector3d w1_quad = RotEl::interpolate(w1, connectivity.quadrature_points[i]);
                    const Eigen::Quaterniond q_start = this->quadrature_data.q0_quad.at(quadrature_data_offset + ele_i * num_quadrature_points + i);
                    const Eigen::Quaterniond q_end = quaternion_time_integration_global_w(q_start, w1_quad, dt, true);
                    this->quadrature_data.q0_quad.at(quadrature_data_offset + ele_i * num_quadrature_points + i) = q_end;
                }
            }
        }
    };

    this->for_each_connectivity(update_quadrature_point_rotations);

    // Update previous state
    this->q0 = this->q1;
    this->w0 = this->w1;
}

void stark::EnergyMicropolarShells::_write_frame(stark::core::Stark& stark)
{
    // Rotation DOF mesh: interpolate positions of nodes of rotation mesh and export micropolar DOF
    for (int group = 0; group < this->labels.size(); ++group) {
        vtkio::VTKFile vtk_file;

        const auto& label = this->labels.at(group);
        const std::string path = stark.get_frame_path(label + "_rotations") + ".vtk";

        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3d> angular_velocities;
        std::vector<Eigen::Vector3d> dir_x, dir_y, dir_z;

        const int mp_dof_begin = this->micropolar_dof_offsets.at(group);
        const int mp_dof_end = this->micropolar_dof_offsets.at(group+1);
        const int n_micropolar_dofs = mp_dof_end - mp_dof_begin;

        vertices.resize(n_micropolar_dofs, Eigen::Vector3d::Zero());
        angular_velocities.insert(angular_velocities.end(), this->w1.begin() + mp_dof_begin, this->w1.begin() + mp_dof_end);
        
        dir_x.resize(n_micropolar_dofs, Eigen::Vector3d::UnitX());
        dir_y.resize(n_micropolar_dofs, Eigen::Vector3d::UnitY());
        dir_z.resize(n_micropolar_dofs, Eigen::Vector3d::UnitZ());

        for (int dof_i = mp_dof_begin; dof_i < mp_dof_end; ++dof_i) {
            const Eigen::Matrix3d R = this->q1[dof_i].to_rotation_matrix();
            dir_x[dof_i - mp_dof_begin] = R.col(0);
            dir_y[dof_i - mp_dof_begin] = R.col(1);
            dir_z[dof_i - mp_dof_begin] = R.col(2);
        }

        const auto interpolate_displacements = [&]<typename D, typename R>(const Fem::MixedFemEnergyConnectivity<D, R>& _connectivity, const int group)
        {
            using DispEl = Fem::Element<D>;
            using RotEl = Fem::Element<R>;

            const auto& point_set = this->point_sets.at(group);

            // Skip meshes of other element types
            if (!this->has_element_types<DispEl, RotEl>(group)) {
                return;
            }

            std::vector<std::array<int, DispEl::num_nodes>>& displacement_mesh = this->displacement_meshes.get_mesh<DispEl>(group);
            std::vector<std::array<int, RotEl::num_nodes>>& rotation_mesh = this->rotation_meshes.get_mesh<RotEl>(group);
            // Get the node coordinates of the rotation mesh reference element
            const Eigen::Matrix<double, 2, RotEl::num_nodes> rotation_ref_nodes = RotEl::nodes();

            for (int ele_i = 0; ele_i < displacement_mesh.size(); ele_i++) {
                const std::array<int, DispEl::num_nodes>& ele_disp = displacement_mesh[ele_i];
                const std::array<int, DispEl::num_nodes> ele_disp_glob = point_set.get_global_indices(ele_disp);

                Eigen::Matrix<double, 3, DispEl::num_nodes> x1_mat;
                for (int i = 0; i < DispEl::num_nodes; ++i) {
                    x1_mat.col(i) = this->dyn->x1.get(ele_disp_glob[i]);
                }
                // Interpolate position for each node of the rotation mesh element
                const std::array<int, RotEl::num_nodes>& ele_rot = rotation_mesh[ele_i];
                for (int node_i = 0; node_i < RotEl::num_nodes; ++node_i) {
                    const Eigen::Vector2d node = rotation_ref_nodes.col(node_i);
                    const Eigen::Vector<double, DispEl::num_nodes> phi = DispEl::basis(node);
                    const Eigen::Vector3d node_deformed = x1_mat * phi;
                    // Overwrite any previously interpolated positions, should be identical for conforming meshes
                    vertices.at(ele_rot[node_i]) = node_deformed;
                }
            }
        };

        this->for_each_connectivity([&](const auto& connectivity) { interpolate_displacements(connectivity, group); });

        if (vertices.empty()) {
            vtk_file.write_empty(path);
        } else {
            vtk_file.set_points_from_twice_indexable(vertices);
            vtk_file.set_cells_as_particles(vertices.size());
            vtk_file.set_point_data_from_twice_indexable("angular_velocity", angular_velocities, vtkio::AttributeType::Vectors);
            vtk_file.set_point_data_from_twice_indexable("x_axis", dir_x, vtkio::AttributeType::Vectors);
            vtk_file.set_point_data_from_twice_indexable("y_axis", dir_y, vtkio::AttributeType::Vectors);
            vtk_file.set_point_data_from_twice_indexable("z_axis", dir_z, vtkio::AttributeType::Vectors);
            vtk_file.write(path);
        }
    }

    // Nonlinear mesh subdivision
    if (this->WRITE_SUBDIVIDED_MESHES) {
        for (int group = 0; group < this->labels.size(); ++group) {
            vtkio::VTKFile vtk_file;
            vtkio::VTKFile vtk_file_edges;

            const auto& label = this->labels.at(group);
            const std::string path = stark.get_frame_path(label + "_surface_subdiv") + ".vtk";
            const std::string path_edges = stark.get_frame_path(label + "_surface_edges") + ".vtk";

            std::vector<Eigen::Vector3d> vertices_subdiv;
            std::vector<Eigen::Vector3d> normals_subdiv;
            std::vector<std::array<int, 3>> triangles_subdiv;

            std::vector<std::vector<double>> scalar_fields;
            scalar_fields.resize(this->scalar_field_callbacks.size());

            std::vector<Eigen::Vector3d> vertices_edges;
            std::vector<std::array<int, 2>> lines_edges;

            const auto extract_subdivided_mesh = [&]<typename D, typename R>(
                    const Fem::MixedFemEnergyConnectivity<D, R>& connectivity,
                    int group)
            {
                using DispEl = Fem::Element<D>;
                using RotEl = Fem::Element<R>;
                
                if (!this->has_element_types<DispEl, RotEl>(group)) {
                    return;
                }

                const auto& set = this->point_sets.at(group);
                const auto& mesh = this->displacement_meshes.get_mesh<DispEl>(group);

                const auto subdivided = this->subdivided_meshes.at(group);
                const auto& triangle_soup = subdivided.first;
                const auto& vertex_data = subdivided.second;

                if (triangle_soup.empty()) {
                    return;
                }

                vertices_subdiv.reserve(vertices_subdiv.size() + triangle_soup.size() * 3);
                normals_subdiv.reserve(normals_subdiv.size() + triangle_soup.size() * 3);
                triangles_subdiv.reserve(triangles_subdiv.size() + triangle_soup.size());
                for (auto& field : scalar_fields) {
                    field.reserve(field.size() + triangle_soup.size() * 3);
                }
                // Upper-bound for edges contains all edges/vertices of subdivided triangles
                vertices_edges.reserve(vertices_edges.size() + triangle_soup.size() * 6);
                lines_edges.reserve(lines_edges.size() + triangle_soup.size() * 3);

                for (int ele_i = 0; ele_i < triangle_soup.size(); ++ele_i) {
                    const auto& local_tri = triangle_soup[ele_i];
                    const int origin_element = vertex_data[local_tri[0]].origin_tri;

                    const std::array<int, DispEl::num_nodes>& element_conn = mesh.at(origin_element);
                    const std::array<int, DispEl::num_nodes> element_conn_glob = set.get_global_indices(element_conn);

                    Eigen::Matrix<double, 3, DispEl::num_nodes> X_mat;
                    for (int i = 0; i < DispEl::num_nodes; ++i) {
                        X_mat.col(i) = this->dyn->X.get(element_conn_glob[i]);
                    }

                    Eigen::Matrix<double, 3, DispEl::num_nodes> x1_mat;
                    for (int i = 0; i < DispEl::num_nodes; ++i) {
                        x1_mat.col(i) = this->dyn->x1.get(element_conn_glob[i]);
                    }

                    const int offset = (int) vertices_subdiv.size();
                    triangles_subdiv.push_back({offset + 0, offset + 1, offset + 2});

                    // Add triangle vertices
                    for (int vert_i = 0; vert_i < 3; ++vert_i) {
                        const auto& xi = vertex_data[local_tri[vert_i]].local_coords;
                        const Eigen::Matrix<double, DispEl::num_nodes, 1> phi = DispEl::basis(xi);
                        const Eigen::Matrix<double, DispEl::num_nodes, 2> dphidxi = DispEl::basis_grad(xi).transpose();
                        const Eigen::Vector3d x1_local = x1_mat * phi;
                        vertices_subdiv.push_back(x1_local);

                        const Eigen::Matrix<double, 3, 2> g = x1_mat * dphidxi;
                        const Eigen::Vector3d n = g.col(0).cross(g.col(1)).normalized();
                        normals_subdiv.push_back(n);

                        for (int field_i = 0; field_i < scalar_fields.size(); ++field_i) {
                            const Eigen::Vector3d X_local = X_mat * phi;
                            scalar_fields[field_i].push_back(this->scalar_field_callbacks.at(field_i).callback(X_local));
                        }
                    }

                    // Add exterior edges
                    for (int i = 0; i < 3; ++i) {
                        const int i_next = (i + 1) % 3;

                        const auto& xi = vertex_data[local_tri[i]].local_coords;
                        const auto& xi_next = vertex_data[local_tri[(i + 1) % 3]].local_coords;

                        if constexpr (DispEl::ref_shape == Fem::ElementShape::Tri2d) {
                            if ((xi.x() == 0.0 && xi_next.x() == 0.0) || (xi.y() == 0.0 && xi_next.y() == 0.0) ||
                                ((xi.x() + xi.y()) == 1.0 && (xi_next.x() + xi_next.y()) == 1.0)) {
                                const int offset_edge = (int) vertices_edges.size();
                                lines_edges.push_back({offset_edge + 0, offset_edge + 1});
                                vertices_edges.push_back(vertices_subdiv[offset + i]);
                                vertices_edges.push_back(vertices_subdiv[offset + i_next]);
                            }
                        } else if constexpr (DispEl::ref_shape == Fem::ElementShape::Quad2d) {
                            if ((xi.x() == 0.0 && xi_next.x() == 0.0) || (xi.y() == 0.0 && xi_next.y() == 0.0) ||
                                    (xi.x() == 1.0 && xi_next.x() == 1.0) || (xi.y() == 1.0 && xi_next.y() == 1.0)) {
                                const int offset_edge = (int) vertices_edges.size();
                                lines_edges.push_back({offset_edge + 0, offset_edge + 1});
                                vertices_edges.push_back(vertices_subdiv[offset + i]);
                                vertices_edges.push_back(vertices_subdiv[offset + i_next]);
                            }
                        } else {
                            static_assert(!std::is_same_v<D, D>, "unsupported shape");
                        }
                    }
                }
            };

            this->for_each_connectivity([&](const auto& conn) { extract_subdivided_mesh(conn, group); });

            if (vertices_subdiv.empty()) {
                vtk_file.write_empty(path);
                vtk_file_edges.write_empty(path_edges);
            } else {
                vtk_file.set_points_from_twice_indexable(vertices_subdiv);
                vtk_file.set_cells_from_twice_indexable(triangles_subdiv, vtkio::CellType::Triangle);
                vtk_file.set_point_data_from_twice_indexable("normals", normals_subdiv, vtkio::AttributeType::Vectors);
                for (int i = 0; i < scalar_fields.size(); ++i) {
                    vtk_file.set_point_data_from_indexable(this->scalar_field_callbacks.at(i).name, scalar_fields[i], vtkio::AttributeType::Scalars);
                }
                vtk_file.write(path);

                vtk_file_edges.set_points_from_twice_indexable(vertices_edges);
                vtk_file_edges.set_cells_from_twice_indexable(lines_edges, vtkio::CellType::Line);
                vtk_file_edges.write(path_edges);
            }
        }
    }

    // Quadrature point mesh
    if (this->WRITE_QUADRATURE_MESHES) {
        for (int group = 0; group < this->labels.size(); ++group) {
            vtkio::VTKFile vtk_file;

            const auto& label = this->labels.at(group);
            const std::string path = stark.get_frame_path(label + "_quadrature") + ".vtk";

            struct QuadratureOutputData {
                std::vector<Eigen::Vector3d> quadrature_points;
                std::vector<int> face_id;
                std::vector<Eigen::Vector3d> a1_cov, a2_cov, n0;
                std::vector<Eigen::Vector3d> a1_conv, a2_conv;
                std::vector<Eigen::Vector3d> dir_x, dir_y, dir_z;
                std::vector<double> Gamma11, Gamma12, Gamma13, Gamma21, Gamma22, Gamma23, Gamma31, Gamma32, Gamma33;
                std::vector<double> H, K;
            };

            QuadratureOutputData data;

            const int num_total_quadrature_points = this->quadrature_data.q0_quad.size();
            data.quadrature_points.resize(num_total_quadrature_points);
            data.face_id.resize(num_total_quadrature_points);
            data.a1_cov.resize(num_total_quadrature_points);
            data.a2_cov.resize(num_total_quadrature_points);
            data.n0.resize(num_total_quadrature_points);
            data.a1_conv.resize(num_total_quadrature_points);
            data.a2_conv.resize(num_total_quadrature_points);
            data.dir_x.resize(num_total_quadrature_points);
            data.dir_y.resize(num_total_quadrature_points);
            data.dir_z.resize(num_total_quadrature_points);
            data.H.resize(num_total_quadrature_points);
            data.K.resize(num_total_quadrature_points);
            data.Gamma11.resize(num_total_quadrature_points);
            data.Gamma12.resize(num_total_quadrature_points);
            data.Gamma13.resize(num_total_quadrature_points);
            data.Gamma21.resize(num_total_quadrature_points);
            data.Gamma22.resize(num_total_quadrature_points);
            data.Gamma23.resize(num_total_quadrature_points);
            data.Gamma31.resize(num_total_quadrature_points);
            data.Gamma32.resize(num_total_quadrature_points);
            data.Gamma33.resize(num_total_quadrature_points);

            const auto extract_quadrature_points = [&]<typename D, typename R>(
                    const Fem::MixedFemEnergyConnectivity<D, R>& connectivity,
                    int group)
            {
                using DispEl = Fem::Element<D>;
                using RotEl = Fem::Element<R>;

                if (!this->has_element_types<DispEl, RotEl>(group)) {
                    return;
                }

                const auto& point_set = this->point_sets.at(group);
                const auto& mesh = this->displacement_meshes.get_mesh<DispEl>(group);
                const auto& rot_mesh = this->rotation_meshes.get_mesh<RotEl>(group);

                const int num_quadrature_points = connectivity.quadrature_points.size();
                const int quadrature_offset = this->quadrature_data.offsets.at(group);

                for (int ele_i = 0; ele_i < mesh.size(); ++ele_i) {
                    const std::array<int, DispEl::num_nodes>& ele = mesh[ele_i];
                    const std::array<int, RotEl::num_nodes>& rot_ele = rot_mesh[ele_i];
                    const std::array<int, RotEl::num_nodes>& rot_ele_glob = this->micropolar_get_global_indices(group, rot_ele);

                    Eigen::Matrix<double, 3, DispEl::num_nodes> X_mat;
                    for (int i = 0; i < DispEl::num_nodes; ++i) {
                        X_mat.col(i) = point_set.get_rest_position(ele[i]);
                    }

                    Eigen::Matrix<double, 3, DispEl::num_nodes> x1_mat;
                    for (int i = 0; i < DispEl::num_nodes; ++i) {
                        x1_mat.col(i) = point_set.get_position(ele[i]);
                    }

                    std::array<Eigen::Matrix3d, RotEl::num_nodes> Ri;
                    for (int i = 0; i < RotEl::num_nodes; ++i) {
                        const Eigen::Quaterniond q = this->q1[rot_ele_glob[i]];
                        Ri[i] = q.toRotationMatrix();
                    }

                    for (int quad_i = 0; quad_i < num_quadrature_points; ++quad_i) {
                        const int global_quad_i = quadrature_offset + ele_i * num_quadrature_points + quad_i;

                        const Eigen::Vector2d q_point = connectivity.quadrature.points[quad_i];
                        const Eigen::Vector<double, DispEl::num_nodes> phi = DispEl::basis(q_point);
                        const Eigen::Vector3d x_quad = x1_mat * phi;

                        data.face_id[global_quad_i] = ele_i;
                        data.quadrature_points[global_quad_i] = x_quad;

                        const Eigen::Matrix<double, DispEl::num_nodes, 2> dphidxi = DispEl::basis_grad(q_point).transpose();		// Nx2
                        const Eigen::Matrix<double, 3, 2> G = X_mat * dphidxi;										// 3x2
                        const Eigen::Matrix<double, 3, 2> g = x1_mat * dphidxi;										// 3x2
                        const Eigen::Matrix<double, 2, 3> Ginv = (G.transpose() * G).inverse() * G.transpose();     // 2x3

                        const Eigen::Vector3d a1_cov = G.col(0);
                        const Eigen::Vector3d a2_cov = G.col(1);
                        const Eigen::Vector3d a1_conv = Ginv.transpose().col(0);
                        const Eigen::Vector3d a2_conv = Ginv.transpose().col(1);
                        const Eigen::Vector3d n0 = a1_cov.cross(a2_cov).normalized();

                        data.a1_cov[global_quad_i] = a1_cov;
                        data.a2_cov[global_quad_i] = a2_cov;
                        data.n0[global_quad_i] = n0;
                        data.a1_conv[global_quad_i] = a1_conv;
                        data.a2_conv[global_quad_i] = a2_conv;

                        const Eigen::Matrix3d Q = this->quadrature_data.q0_quad[global_quad_i].to_rotation_matrix();
                        data.dir_x[global_quad_i] = Q.col(0);
                        data.dir_y[global_quad_i] = Q.col(1);
                        data.dir_z[global_quad_i] = Q.col(2);

                        const Eigen::Matrix2d S = this->quadrature_data.S_rest_quad[global_quad_i].to_matrix();
                        const Eigen::Matrix3d b =  this->quadrature_data.b_rest_quad[global_quad_i].to_matrix();
                        const double H = 0.5 * b.trace();
                        const double K = 0.5 * (b.trace()*b.trace() - (b*b).trace());

                        data.H[global_quad_i] = H;
                        data.K[global_quad_i] = K;

                        const auto axl = [](const Eigen::Matrix3d& m) { return Eigen::Vector3d(m(2,1), m(0,2), m(1,0)); };
                        const Eigen::Matrix<double, RotEl::num_nodes, 2> dphidxi_rot = RotEl::basis_grad(q_point).transpose();		// Nx2
                        Eigen::Matrix3d grad_R_alpha = Eigen::Matrix3d::Zero();
                        Eigen::Matrix3d grad_R_beta = Eigen::Matrix3d::Zero();
                        for (int i = 0; i < RotEl::num_nodes; ++i) {
                            grad_R_alpha += dphidxi_rot(i, 0) * Ri[i];
                            grad_R_beta += dphidxi_rot(i, 1) * Ri[i];
                        }
                        const Eigen::Vector3d Gamma_alpha = axl(Q.transpose() * grad_R_alpha);
                        const Eigen::Vector3d Gamma_beta = axl(Q.transpose() * grad_R_beta);

                        Eigen::Matrix<double, 3, 2> Gamma_omega;
                        Gamma_omega.col(0) = Gamma_alpha;
                        Gamma_omega.col(1) = Gamma_beta;
                        Eigen::Matrix3d Gamma = Gamma_omega * Ginv;

                        data.Gamma11[global_quad_i] = Gamma(0, 0);
                        data.Gamma12[global_quad_i] = Gamma(0, 1);
                        data.Gamma13[global_quad_i] = Gamma(0, 2);
                        data.Gamma21[global_quad_i] = Gamma(1, 0);
                        data.Gamma22[global_quad_i] = Gamma(1, 1);
                        data.Gamma23[global_quad_i] = Gamma(1, 2);
                        data.Gamma31[global_quad_i] = Gamma(2, 0);
                        data.Gamma32[global_quad_i] = Gamma(2, 1);
                        data.Gamma33[global_quad_i] = Gamma(2, 2);
                    }
                }
            };

            this->for_each_connectivity([&](const auto& conn) { extract_quadrature_points(conn, group); });

            if (false) {
                vtk_file.write_empty(path);
            } else {
                vtk_file.set_points_from_twice_indexable(data.quadrature_points);
                vtk_file.set_cells_as_particles(data.quadrature_points.size());
                vtk_file.set_point_data_from_indexable("face_id", data.face_id, vtkio::AttributeType::Scalars);

                vtk_file.set_point_data_from_twice_indexable("a1_cov", data.a1_cov, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_twice_indexable("a2_cov", data.a2_cov, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_twice_indexable("n0", data.n0, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_twice_indexable("a1_conv", data.a1_conv, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_twice_indexable("a2_conv", data.a2_conv, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_indexable("H0", data.H, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("K0", data.K, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_twice_indexable("x_axis", data.dir_x, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_twice_indexable("y_axis", data.dir_y, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_twice_indexable("z_axis", data.dir_z, vtkio::AttributeType::Vectors);
                vtk_file.set_point_data_from_indexable("Gamma11", data.Gamma11, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma12", data.Gamma12, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma13", data.Gamma13, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma21", data.Gamma21, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma22", data.Gamma22, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma23", data.Gamma23, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma31", data.Gamma31, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma32", data.Gamma32, vtkio::AttributeType::Scalars);
                vtk_file.set_point_data_from_indexable("Gamma33", data.Gamma33, vtkio::AttributeType::Scalars);

                vtk_file.write(path);
            }
        }
    }
}

void stark::EnergyMicropolarShells::_initialize_rest_state_quantities(int group)
{
    const auto compute = [&]<typename D, typename R>(
            const Fem::MixedFemEnergyConnectivity<D, R>& connectivity,
            const int group)
    {
        using DispEl = Fem::Element<D>;
        using RotEl = Fem::Element<R>;

        if (!this->has_element_types<DispEl, RotEl>(group)) {
            return;
        }

        const int num_quadrature_points = connectivity.quadrature_points.size();
        const int quadrature_begin = this->quadrature_data.offsets.at(group);

        const auto& set = this->point_sets.at(group);
        const auto& displacement_mesh = this->displacement_meshes.get_mesh<DispEl>(group);

        const Eigen::Matrix<double, 2, DispEl::num_nodes> ref_nodes = DispEl::nodes();

        for (int ele_i = 0; ele_i < displacement_mesh.size(); ++ele_i) {
            const std::array<int, DispEl::num_nodes>& conn = displacement_mesh[ele_i];
            const std::array<int, DispEl::num_nodes> conn_tri_glob = set.get_global_indices(conn);

            Eigen::Matrix<double, 3, DispEl::num_nodes> X_mat;
            for (int i = 0; i < DispEl::num_nodes; ++i) {
                X_mat.col(i) = this->dyn->X.get(conn_tri_glob[i]);
            }

            Eigen::Matrix<double, 3, DispEl::num_nodes> n_mat;
            for (int i = 0; i < DispEl::num_nodes; ++i) {
                const Eigen::Vector2d& x_ref = ref_nodes.col(i);
                const Eigen::Matrix<double, DispEl::num_nodes, 2> dphidxi = DispEl::basis_grad(x_ref).transpose();		    // Nx2
                const Eigen::Matrix<double, 3, 2> G = X_mat * dphidxi;										// 3x2

                const Eigen::Vector3d a1_cov = G.col(0);
                const Eigen::Vector3d a2_cov = G.col(1);
                const Eigen::Vector3d n = a1_cov.cross(a2_cov).normalized();

                n_mat.col(i) = n;
            }

            for (int quad_i = 0; quad_i < num_quadrature_points; ++quad_i) {
                const int global_quad_i = quadrature_begin + ele_i * num_quadrature_points + quad_i;
                const Eigen::Vector2d q_point = connectivity.quadrature.points[quad_i];

                const Eigen::Matrix<double, DispEl::num_nodes, 2> dphidxi = DispEl::basis_grad(q_point).transpose();		// Nx2
                const Eigen::Matrix<double, 3, 2> G = X_mat * dphidxi;										// 3x2
                const Eigen::Matrix<double, 2, 3> Ginv = (G.transpose() * G).inverse() * G.transpose();

                const double g_inv = std::sqrt((G.transpose() * G).inverse().determinant());
                if (!std::isfinite(g_inv)) {
                    fmt::print("Error: Non finite inverse rest metric encountered for element {}\n", ele_i);
                }

                const Eigen::Vector3d a1_cov = G.col(0);
                const Eigen::Vector3d a2_cov = G.col(1);
                const Eigen::Vector3d n = a1_cov.cross(a2_cov).normalized();

                Eigen::Matrix3d F;
                F.col(0) = a1_cov;
                F.col(1) = a2_cov;
                F.col(2) = n;

                Eigen::Matrix3d R0 = polar_decomposition_3d(F);
                this->quadrature_data.R_rest_quad[global_quad_i] = R0;

                const Eigen::Matrix<double, 3, 2> grad_n = n_mat * dphidxi;
                const Eigen::Matrix2d II_dn = -1.0 * grad_n.transpose() * G;
                const Eigen::Matrix2d II = II_dn;
                // Shape operator 2x2
                const Eigen::Matrix2d S = (G.transpose() * G).inverse() * II;
                // "Second fundamental tensor" 3x3
                const Eigen::Matrix3d b = -1.0 * grad_n * Ginv;

                const Eigen::Matrix3d a = G * Ginv;

                this->quadrature_data.a_rest_quad[global_quad_i] = a;
                this->quadrature_data.b_rest_quad[global_quad_i] = b;
                this->quadrature_data.S_rest_quad[global_quad_i] = S;
            }
        }
    };

    this->for_each_connectivity([&](const auto& conn) { compute(conn, group); });
}

template<typename D, typename R>
stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl(
        const std::string& label,
        const PointSetHandler& set,
        const std::vector<std::array<int, D::num_nodes>>& displacement_elements,
        const Mesh<R::num_nodes>& rotation_mesh,
        const Params& params)
{
    using DispEl = Fem::Element<D>;
    using RotEl = Fem::Element<R>;

    const auto& rotation_elements = rotation_mesh.conn;

    if (displacement_elements.size() != rotation_elements.size()) {
        throw std::runtime_error(fmt::format("EnergyMicropolarShells::add_mesh_impl: Number of displacement elements ({}) must match number of rotation elements ({})", displacement_elements.size(), rotation_elements.size()));
    }

    const int group = (int)this->labels.size();
    this->labels.push_back(label);
    this->point_sets.push_back(set);

    // Store parameters
    this->thickness.push_back(params.thickness);
    this->youngs_modulus.push_back(params.youngs_modulus);
    this->poissons_ratio.push_back(params.poissons_ratio);
    this->density.push_back(params.density);
    this->angular_inertia_density.push_back(params.angular_inertia_density);
    this->inertia_damping.push_back(params.inertia_damping);
    this->mu_c_factor.push_back(params.mu_c_factor);
    this->length_scale.push_back(params.length_scale);
    this->shear_correction.push_back(params.shear_correction);
    this->bending_correction.push_back(params.bending_correction);

    // Store meshes
    this->displacement_meshes.add_mesh<DispEl>(displacement_elements);
    this->rotation_meshes.add_mesh<RotEl>(rotation_elements);
    this->rotation_mesh_vertices.push_back(rotation_mesh.vertices);

    // Initialize micropolar DOF
    const int num_micropolar_dof = rotation_mesh.vertices.size();
    const int micropolar_dof_offset = this->micropolar_dof_offsets.back();
    this->micropolar_dof_offsets.push_back(micropolar_dof_offset + num_micropolar_dof);
    this->q0.resize(this->q0.size() + num_micropolar_dof, QuaternionWrapperd::identity());
    this->q1.resize(this->q1.size() + num_micropolar_dof, QuaternionWrapperd::identity());
    this->w0.resize(this->w0.size() + num_micropolar_dof, Eigen::Vector3d::Zero());
    this->w1.resize(this->w1.size() + num_micropolar_dof, Eigen::Vector3d::Zero());
    this->aa.resize(this->aa.size() + num_micropolar_dof, Eigen::Vector3d::Zero());
    this->t.resize(this->t.size() + num_micropolar_dof, Eigen::Vector3d::Zero());

    // Compute nodal lumped area and create lumped inertia energy connectivity
    {
        std::vector<double> lumped_area_disp(set.size(), 0.0);
        std::vector<double> lumped_area_rot(num_micropolar_dof, 0.0);

        for (int ele_i = 0; ele_i < displacement_elements.size(); ele_i++) {
            std::array<Eigen::Vector3d, DispEl::num_nodes> nodes;

            const std::array<int, DispEl::num_nodes>& conn_disp = displacement_elements[ele_i];
            const std::array<int, RotEl::num_nodes>& conn_rot = rotation_elements[ele_i];

            for (int i = 0; i < DispEl::num_nodes; i++) {
                nodes[i] = set.get_rest_position(conn_disp[i]);
            }

            const double total_area = DispEl::area(nodes);
            for (int i = 0; i < conn_disp.size(); i++) {
                lumped_area_disp.at(conn_disp[i]) += total_area / DispEl::num_nodes;
            }
            for (int i = 0; i < conn_rot.size(); i++) {
                lumped_area_rot.at(conn_rot[i]) += total_area / RotEl::num_nodes;
            }
        }
        this->nodal_lumped_area_disp.insert(this->nodal_lumped_area_disp.end(), lumped_area_disp.begin(), lumped_area_disp.end());
        this->nodal_lumped_area_rot.insert(this->nodal_lumped_area_rot.end(), lumped_area_rot.begin(), lumped_area_rot.end());

        for (int i = 0; i < lumped_area_disp.size(); i++) {
            this->conn_lumped_disp_inertia.numbered_push_back({group, set.get_global_index(i)});
        }

        for (int i = 0; i < lumped_area_rot.size(); i++) {
            this->conn_lumped_rot_inertia.numbered_push_back({group, micropolar_dof_offset + i});
        }
    }

    auto add_micropolar_connectivity = [this, &set, group]<typename D2, typename R2>(
            Fem::MixedFemEnergyConnectivity<D2, R2>& connectivity,
            const std::vector<std::array<int, D2::num_nodes>>& disp_elements,
            const std::vector<std::array<int, R2::num_nodes>>& rot_elements)
    {
        using DispEl = Fem::Element<D2>;
        using RotEl = Fem::Element<R2>;

        auto& quadrature = connectivity.quadrature;

        const int num_elements = disp_elements.size();
        const int num_quad_points = quadrature.length();

        const int quadrature_data_offset = this->quadrature_data.offsets.back();
        const int num_total_quad_points = num_quad_points * num_elements;

        // Initialize connectivity
        for (int ele_i = 0; ele_i < num_elements; ele_i++) {
            const std::array<int, DispEl::num_nodes>& conn_disp = disp_elements[ele_i];
            const std::array<int, DispEl::num_nodes> conn_disp_glob = set.get_global_indices(conn_disp);

            const std::array<int, RotEl::num_nodes>& conn_rot = rot_elements[ele_i];
            const std::array<int, RotEl::num_nodes> conn_rot_glob = this->micropolar_get_global_indices(group, conn_rot);

            for (int quad_i = 0; quad_i < num_quad_points; ++quad_i) {
                const int q_idx_global = quadrature_data_offset + ele_i * num_quad_points + quad_i;
                connectivity.push_element(ele_i, group, conn_disp_glob, conn_rot_glob, quad_i, q_idx_global);
            }
        }

        // Initialize per-quadrature-point data
        this->quadrature_data.offsets.push_back(quadrature_data_offset + num_total_quad_points);
        this->quadrature_data.q0_quad.resize(this->quadrature_data.q0_quad.size() + num_total_quad_points, QuaternionWrapperd::identity());
        this->quadrature_data.R_rest_quad.resize(this->quadrature_data.a_rest_quad.size() + num_total_quad_points, MatrixWrapper3d::identity());
        this->quadrature_data.S_rest_quad.resize(this->quadrature_data.a_rest_quad.size() + num_total_quad_points, MatrixWrapper2d::identity());
        this->quadrature_data.a_rest_quad.resize(this->quadrature_data.a_rest_quad.size() + num_total_quad_points, MatrixWrapper3d::identity());
        this->quadrature_data.b_rest_quad.resize(this->quadrature_data.a_rest_quad.size() + num_total_quad_points, MatrixWrapper3d::identity());
        this->quadrature_data.Gamma_rest_quad.resize(this->quadrature_data.Gamma_rest_quad.size() + num_total_quad_points, MatrixWrapper3d::zero());
    };

    if constexpr (DispEl::element_type == Fem::ElementType::Tri3 && RotEl::element_type == Fem::ElementType::Tri3) {
        add_micropolar_connectivity(this->conn_disp_tri3_rot_tri3, displacement_elements, rotation_elements);
    } else if constexpr (DispEl::element_type == Fem::ElementType::Tri6 && RotEl::element_type == Fem::ElementType::Tri3) {
        add_micropolar_connectivity(this->conn_disp_tri6_rot_tri3, displacement_elements, rotation_elements);
    } else if constexpr (DispEl::element_type == Fem::ElementType::Tri6 && RotEl::element_type == Fem::ElementType::Tri6) {
        add_micropolar_connectivity(this->conn_disp_tri6_rot_tri6, displacement_elements, rotation_elements);
    } else if constexpr (DispEl::element_type == Fem::ElementType::Tri10 && RotEl::element_type == Fem::ElementType::Tri3) {
        add_micropolar_connectivity(this->conn_disp_tri10_rot_tri3, displacement_elements, rotation_elements);
    } else if constexpr (DispEl::element_type == Fem::ElementType::Tri10 && RotEl::element_type == Fem::ElementType::Tri6) {
        add_micropolar_connectivity(this->conn_disp_tri10_rot_tri6, displacement_elements, rotation_elements);
    } else if constexpr (DispEl::element_type == Fem::ElementType::Quad4 && RotEl::element_type == Fem::ElementType::Quad4) {
        add_micropolar_connectivity(this->conn_disp_quad4_rot_quad4, displacement_elements, rotation_elements);
    } else if constexpr (DispEl::element_type == Fem::ElementType::Quad9 && RotEl::element_type == Fem::ElementType::Quad4) {
        add_micropolar_connectivity(this->conn_disp_quad9_rot_quad4, displacement_elements, rotation_elements);
    } else if constexpr (DispEl::element_type == Fem::ElementType::Quad9 && RotEl::element_type == Fem::ElementType::Quad9) {
        add_micropolar_connectivity(this->conn_disp_quad9_rot_quad9, displacement_elements, rotation_elements);
    } else {
        static_assert(!std::is_same_v<D, D>, "Unsupported element combination");
    }

    // Create subdivided mesh for visualization
    {
        if constexpr (DispEl::element_type == Fem::ElementType::Tri3) {
            // No need to subdivide Tri3 meshes, just output the original mesh with local coordinates for consistency
            std::vector<std::array<int, 3>> triangle_soup;
            std::vector<VertexLocalCoords> local_coords;
            triangle_soup.resize(displacement_elements.size());
            local_coords.resize(3 * displacement_elements.size());
            for (int i = 0; i < displacement_elements.size(); ++i) {
                triangle_soup[i] = {3*i+0, 3*i+1, 3*i+2};
                local_coords[3 * i + 0] = {i, Eigen::Vector2d(0, 0)};
                local_coords[3 * i + 1] = {i, Eigen::Vector2d(0, 1)};
                local_coords[3 * i + 2] = {i, Eigen::Vector2d(1, 0)};
            }
            this->subdivided_meshes.emplace_back(triangle_soup, local_coords);
        } else if constexpr (DispEl::element_type == Fem::ElementType::Tri6) {
            this->subdivided_meshes.push_back(tri6_to_tri3_subdivide_with_local_coords(displacement_elements, this->MESH_SUBDIVISION_LEVEL));
        } else if constexpr (DispEl::element_type == Fem::ElementType::Tri10) {
            this->subdivided_meshes.push_back(tri10_to_tri3_subdivide_with_local_coords(displacement_elements, this->MESH_SUBDIVISION_LEVEL));
        } else if constexpr (DispEl::element_type == Fem::ElementType::Quad4) {
            this->subdivided_meshes.push_back(quad4_to_tri3_subdivide_with_local_coords(displacement_elements, this->MESH_SUBDIVISION_LEVEL));
        } else if constexpr (DispEl::element_type == Fem::ElementType::Quad9) {
            this->subdivided_meshes.push_back(quad9_to_tri3_subdivide_with_local_coords(displacement_elements, this->MESH_SUBDIVISION_LEVEL));
        } else {
            this->subdivided_meshes.emplace_back();
        }
    }

    return EnergyMicropolarShells::Handler(this, group);
}

template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisTri3, stark::Fem::BasisTri3>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisTri3::num_nodes>>&, const Mesh<Fem::BasisTri3::num_nodes>&, const stark::EnergyMicropolarShells::Params&);
template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisTri6, stark::Fem::BasisTri3>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisTri6::num_nodes>>&, const Mesh<Fem::BasisTri3::num_nodes>&, const stark::EnergyMicropolarShells::Params&);
template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisTri6, stark::Fem::BasisTri6>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisTri6::num_nodes>>&, const Mesh<Fem::BasisTri6::num_nodes>&, const stark::EnergyMicropolarShells::Params&);
template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisTri10, stark::Fem::BasisTri3>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisTri10::num_nodes>>&, const Mesh<Fem::BasisTri3::num_nodes>&, const stark::EnergyMicropolarShells::Params&);
template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisTri10, stark::Fem::BasisTri6>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisTri10::num_nodes>>&, const Mesh<Fem::BasisTri6::num_nodes>&, const stark::EnergyMicropolarShells::Params&);
template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisQuad4, stark::Fem::BasisQuad4>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisQuad4::num_nodes>>&, const Mesh<Fem::BasisQuad4::num_nodes>&, const stark::EnergyMicropolarShells::Params&);
template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisQuad9, stark::Fem::BasisQuad4>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisQuad9::num_nodes>>&, const Mesh<Fem::BasisQuad4::num_nodes>&, const stark::EnergyMicropolarShells::Params&);
template stark::EnergyMicropolarShells::Handler stark::EnergyMicropolarShells::add_mesh_impl<stark::Fem::BasisQuad9, stark::Fem::BasisQuad9>(const std::string&, const PointSetHandler&, const std::vector<std::array<int, Fem::BasisQuad9::num_nodes>>&, const Mesh<Fem::BasisQuad9::num_nodes>&, const stark::EnergyMicropolarShells::Params&);

stark::EnergyMicropolarShells::Params stark::EnergyMicropolarShells::get_params(const Handler& handler) const
{
    handler.exit_if_not_valid("EnergyMicropolarShells::get_params");

    const int group = handler.get_idx();

    Params params;
    params.thickness = this->thickness.at(group);
    params.youngs_modulus = this->youngs_modulus.at(group);
    params.poissons_ratio = this->poissons_ratio.at(group);
    params.density = this->density.at(group);
    params.angular_inertia_density = this->angular_inertia_density.at(group);
    params.inertia_damping = this->inertia_damping.at(group);
    params.mu_c_factor = this->mu_c_factor.at(group);
    params.length_scale = this->length_scale.at(group);
    params.shear_correction = this->shear_correction.at(group);
    params.bending_correction = this->bending_correction.at(group);
    return params;
}

void stark::EnergyMicropolarShells::set_params(const Handler& handler, const Params& params)
{
    handler.exit_if_not_valid("EnergyMicropolarShells::set_params");

    const int group = handler.get_idx();
    this->thickness.at(group) = params.thickness;
    this->youngs_modulus.at(group) = params.youngs_modulus;
    this->poissons_ratio.at(group) = params.poissons_ratio;
    this->density.at(group) = params.density;
    this->angular_inertia_density.at(group) = params.angular_inertia_density;
    this->inertia_damping.at(group) = params.inertia_damping;
    this->mu_c_factor.at(group) = params.mu_c_factor;
    this->length_scale.at(group) = params.length_scale;
    this->shear_correction.at(group) = params.shear_correction;
    this->bending_correction.at(group) = params.bending_correction;
}

void stark::EnergyMicropolarShells::set_all_rotations(int group, const Eigen::Quaterniond& q)
{
    const int micropolar_begin = this->micropolar_dof_offsets.at(group);
    const int micropolar_end = this->micropolar_dof_offsets.at(group+1);

    std::fill(this->q0.begin() + micropolar_begin, this->q0.begin() + micropolar_end, QuaternionWrapperd::from_quaternion(q));
    std::fill(this->q1.begin() + micropolar_begin, this->q1.begin() + micropolar_end, QuaternionWrapperd::from_quaternion(q));

    const int quadrature_begin = this->quadrature_data.offsets.at(group);
    const int quadrature_end = this->quadrature_data.offsets.at(group+1);
    std::fill(this->quadrature_data.q0_quad.begin() + quadrature_begin, this->quadrature_data.q0_quad.begin() + quadrature_end, QuaternionWrapperd::from_quaternion(q));
}

std::vector<Eigen::Vector3d> stark::EnergyMicropolarShells::interpolate_quadrature_points_mesh(int group)
{
    std::vector<Eigen::Vector3d> x1_quad;

    const auto& point_set = this->point_sets.at(group);

    const auto interpolate_quadrature_points = [&]<typename D, typename R>(const Fem::MixedFemEnergyConnectivity<D, R>& connectivity, const int group)
    {
        using DispEl = Fem::Element<D>;
        using RotEl = Fem::Element<R>;

        if (!this->has_element_types<DispEl, RotEl>(group)) {
            return;
        }

        const auto& mesh = this->displacement_meshes.get_mesh<DispEl>(group);

        const int num_quadrature_points = connectivity.quadrature_points.size();
        x1_quad.resize(mesh.size() * num_quadrature_points, Eigen::Vector3d::Zero());

        for (int ele_i = 0; ele_i < mesh.size(); ++ele_i) {
            const std::array<int, DispEl::num_nodes>& ele = mesh[ele_i];

            std::array<Eigen::Vector3d, DispEl::num_nodes> nodal_x1;
            for (int node_i = 0; node_i < DispEl::num_nodes; ++node_i) {
                nodal_x1[node_i] = point_set.get_position(ele[node_i]);
            }

            for (int quad_i = 0; quad_i < num_quadrature_points; ++quad_i) {
                const Eigen::Vector2d q_point = connectivity.quadrature.points[quad_i];
                const Eigen::Vector3d x1 = DispEl::interpolate(nodal_x1, q_point);
                x1_quad[ele_i * num_quadrature_points + quad_i] = x1;
            }
        }
    };

    this->for_each_connectivity([&](const auto& conn) { interpolate_quadrature_points(conn, group); });

    return x1_quad;
}

Eigen::Matrix3d stark::EnergyMicropolarShells::get_rest_deformation_gradient(int group, int quad_i) const
{
    const int quadrature_begin = this->quadrature_data.offsets.at(group);
    const int quadrature_end = this->quadrature_data.offsets.at(group + 1);

    if (quadrature_begin + quad_i >= quadrature_end) {
        throw std::runtime_error("stark::EnergyMicropolarShells::Handler::set_rest_deformation_gradient: quadrature index out of bounds");
    }

    return this->quadrature_data.a_rest_quad[quadrature_begin + quad_i].to_matrix();
}

void stark::EnergyMicropolarShells::set_rest_deformation_gradient(int group, int quad_i, const Eigen::Matrix3d& F0)
{
    const int quadrature_begin = this->quadrature_data.offsets.at(group);
    const int quadrature_end = this->quadrature_data.offsets.at(group + 1);

    if (quadrature_begin + quad_i >= quadrature_end) {
        throw std::runtime_error("stark::EnergyMicropolarShells::Handler::set_rest_deformation_gradient: quadrature index out of bounds");
    }

    this->quadrature_data.a_rest_quad[quadrature_begin + quad_i] = MatrixWrapper3d(F0);
}

void stark::EnergyMicropolarShells::set_rest_curvature(int group, int quad_i, const Eigen::Matrix3d& curvature)
{
    const int quadrature_begin = this->quadrature_data.offsets.at(group);
    const int quadrature_end = this->quadrature_data.offsets.at(group + 1);

    if (quadrature_begin + quad_i >= quadrature_end) {
        throw std::runtime_error("stark::EnergyMicropolarShells::Handler::set_rest_curvature: quadrature index out of bounds");
    }

    this->quadrature_data.Gamma_rest_quad[quadrature_begin + quad_i] = MatrixWrapper3d(curvature);
}
void stark::EnergyMicropolarShells::prescribe_angular_velocity_inside_aabb(
        int dof_group,
        const Eigen::Vector3d& center,
        const Eigen::Vector3d& size,
        const Eigen::Vector3d& omega_bc,
        const std::array<bool, 3> active,
        const double stiffness)
{
    const auto& point_set = this->point_sets.at(dof_group);
    const Eigen::AlignedBox3d aabb = Eigen::AlignedBox3d(center - 0.5 * size, center + 0.5 * size);

    const int bc_group = this->bc_target_ang_vel.size();
    this->bc_target_ang_vel.push_back(omega_bc);
    this->bc_target_ang_vel_active.emplace_back(static_cast<double>(active[0]), static_cast<double>(active[1]), static_cast<double>(active[2]));
    this->bc_target_ang_vel_stiffness.push_back(stiffness);

    this->add_bc_connectivity(this->conn_prescribed_ang_vel, dof_group, bc_group, aabb);
}
int stark::EnergyMicropolarShells::prescribe_rotation_inside_aabb(
        int dof_group,
        const Eigen::Vector3d& center,
        const Eigen::Vector3d& size,
        const Eigen::Matrix3d& rotation_bc,
        const std::array<bool, 3> active_axis,
        const double stiffness)
{
    const Eigen::AlignedBox3d aabb = Eigen::AlignedBox3d(center - 0.5 * size, center + 0.5 * size);

    const int bc_group = this->bc_target_rot.size();
    this->bc_target_rot.push_back(rotation_bc);
    this->bc_target_rot_active.emplace_back(static_cast<double>(active_axis[0]), static_cast<double>(active_axis[1]), static_cast<double>(active_axis[2]));
    this->bc_target_rot_stiffness.push_back(stiffness);

    this->add_bc_connectivity(this->conn_prescribed_rot, dof_group, bc_group, aabb);

    return bc_group;
}
void stark::EnergyMicropolarShells::add_bc_connectivity(
        symx::LabelledConnectivity<3>& bc_conn,
        int dof_group,
        int bc_group,
        const Eigen::AlignedBox3d& aabb)
{
    const int micropolar_dof_offset = this->micropolar_dof_offsets.at(dof_group);
    const auto& vertices = this->rotation_mesh_vertices.at(dof_group);

    for (int i = 0; i < vertices.size(); ++i) {
        if (aabb.contains(vertices[i])) {
            bc_conn.numbered_push_back({bc_group, micropolar_dof_offset + i});
        }
    }
}
void stark::EnergyMicropolarShells::set_prescribed_rotation(int bc_idx, const Eigen::Matrix3d& rotation_bc)
{
    this->bc_target_rot.at(bc_idx) = rotation_bc;
}

void stark::EnergyMicropolarShells::Handler::update_rest_configuration(const bool reset_rotations)
{
    this->model->_initialize_rest_state_quantities(this->get_idx());
    if (reset_rotations) {
        this->model->set_all_rotations(this->get_idx(), Eigen::Quaterniond::Identity());
    }
}
std::vector<Eigen::Vector3d> stark::EnergyMicropolarShells::Handler::interpolate_quadrature_points_mesh() const
{
    return this->model->interpolate_quadrature_points_mesh(this->get_idx());
}
Eigen::Matrix3d stark::EnergyMicropolarShells::Handler::get_rest_deformation_gradient(int quad_i) const
{
    return this->model->get_rest_deformation_gradient(this->get_idx(), quad_i);
}
void stark::EnergyMicropolarShells::Handler::set_rest_deformation_gradient(int quad_i, const Eigen::Matrix3d& F0)
{
    this->model->set_rest_deformation_gradient(this->get_idx(), quad_i, F0);
}
void stark::EnergyMicropolarShells::Handler::set_rest_curvature(int quad_i, const Eigen::Matrix3d& curvature)
{
    this->model->set_rest_curvature(this->get_idx(), quad_i, curvature);
}
void stark::EnergyMicropolarShells::Handler::prescribe_angular_velocity_inside_aabb(const Eigen::Vector3d& center, const Eigen::Vector3d& size, const Eigen::Vector3d& omega_bc, std::array<bool, 3> active, double stiffness)
{
    this->model->prescribe_angular_velocity_inside_aabb(this->get_idx(), center, size, omega_bc, active, stiffness);
}
int stark::EnergyMicropolarShells::Handler::prescribe_rotation_inside_aabb(const Eigen::Vector3d& center, const Eigen::Vector3d& size, const Eigen::Matrix3d& rotation_bc, std::array<bool, 3> axis_active, double stiffness)
{
    return this->model->prescribe_rotation_inside_aabb(this->get_idx(), center, size, rotation_bc, axis_active, stiffness);
}
void stark::EnergyMicropolarShells::Handler::set_prescribed_rotation(int bc_idx, const Eigen::Matrix3d& rotation_bc)
{
    this->model->set_prescribed_rotation(bc_idx, rotation_bc);
}
void stark::EnergyMicropolarShells::Handler::set_prescribed_rotation(int bc_idx, const double angle_deg, const Eigen::Vector3d& axis)
{
    this->model->set_prescribed_rotation(bc_idx, Eigen::AngleAxisd(deg2rad(angle_deg), axis.normalized()).toRotationMatrix());
}

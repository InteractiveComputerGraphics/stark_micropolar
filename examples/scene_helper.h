#pragma once

#include <stark>

#include "utils.h"

/// Helper type to dynamically construct shells with different element types from parameters
struct SurfaceMaterialBuilder
{
    enum class StrainModel {
        None,
        Kim20,
        Wen23,
        Micropolar,
    };

    enum class BendingModel {
        None,
        DiscreteShells
    };

    struct GeneralSurfaceHandler {
        Option<stark::PointSetHandler> point_set;
        Option<stark::EnergyLumpedInertia::Handler> inertia;
        Option<stark::EnergyFrictionalContact::Handler> contact;
        // Strain energies
        Option<stark::EnergyTriangleStrainKim20::Handler> strain_kim20;
        Option<stark::EnergyTriangleStrainWen23::Handler> strain_wen23;
        Option<stark::EnergyMicropolarShells::Handler> strain_micropolar;
        // Bending energies
        Option<stark::EnergyDiscreteShells::Handler> bending_discrete_shells;
    };

    bool elasticity_only = true;

    Option<bool> enable_contact;
    Option<double> contact_distance;

    Option<double> inertia_area_density;
    Option<double> inertia_damping;

    Option<double> thickness;
    Option<double> strain_young_modulus;
    Option<double> strain_poissons_ratio;

    Option<double> strain_damping;
    Option<double> strain_limit;
    Option<double> strain_limit_stiffness;

    Option<double> bending_stiffness;
    Option<double> bending_damping;
    Option<double> edge_bending_stiffness;

    Option<stark::Fem::ElementType> micropolar_element_type_displacement;
    Option<stark::Fem::ElementType> micropolar_element_type_rotation;
    Option<double> micropolar_mu_c_factor;
    Option<double> micropolar_length_scale;
    Option<double> micropolar_shear_correction;
    Option<double> micropolar_bending_correction;
    Option<double> micropolar_angular_inertia;

    StrainModel strain_model = StrainModel::None;
    BendingModel bending_model = BendingModel::None;

    GeneralSurfaceHandler add_surface_grid(stark::Simulation& sim, const std::string& label, const Eigen::Vector2d& center, const Eigen::Vector2d& dimensions, const std::array<int, 2>& n_quads_per_dim, const double z = 0.0) const {
        const auto element = this->micropolar_element_type_displacement.value_or(stark::Fem::ElementType::Tri3);
        const auto shape = stark::Fem::element_type_to_shape(element);

        if (shape == stark::Fem::ElementShape::Tri2d) {
            auto [vertices, triangles] = stark::generate_triangle_grid(center, dimensions, n_quads_per_dim, z);
            return this->add_surface_from_tri3(sim, label, vertices, triangles);
        } else if (shape == stark::Fem::ElementShape::Quad2d) {
            auto [vertices, quads] = stark::generate_quad_grid(center, dimensions, n_quads_per_dim, z);
            return this->add_surface_from_quad4(sim, label, vertices, quads);
        } else {
            throw std::runtime_error("SurfaceMaterialBuilder::add_surface_grid: unsupported element shape");
        }
    }

    /// Adds a surface to the simulator from a triangle mesh, automatically performs subdivision for higher order elements
    GeneralSurfaceHandler add_surface_from_tri3(stark::Simulation& sim, const std::string& label, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles) const {
        GeneralSurfaceHandler handler;

        std::vector<std::array<int, 3>> contact_mesh = triangles;
        std::vector<std::array<int, 3>> output_mesh = triangles;

        const auto default_models = [&]() {
            handler.point_set = sim.deformables->point_sets->add(vertices);
            handler.inertia = this->add_lumped_inertia_tri3(sim, handler.point_set, triangles);
        };

        switch (this->strain_model) {
            case StrainModel::Kim20: {
                default_models();
                handler.strain_kim20 = this->add_kim20(sim, handler.point_set, triangles);
            } break;
            case StrainModel::Wen23: {
                default_models();
                handler.strain_wen23 = this->add_wen23(sim, handler.point_set, triangles);
            } break;
            case StrainModel::Micropolar: {
                auto& disp_el = this->micropolar_element_type_displacement;
                auto& rot_el = this->micropolar_element_type_rotation;

                auto invalid_combination = [&]() { throw std::runtime_error(fmt::format("SurfaceMaterialBuilder::add_surface_from_tri3: unsupported element type combination ({} + {})", stark::Fem::element_type_to_string(disp_el), stark::Fem::element_type_to_string(rot_el))); };
                auto unimplemented_combination = [&]() { throw std::runtime_error(fmt::format("SurfaceMaterialBuilder::add_surface_from_tri3: unimplemented element type combination ({} + {})", stark::Fem::element_type_to_string(disp_el), stark::Fem::element_type_to_string(rot_el))); };

                if (!disp_el.has_value() || !rot_el.has_value()) {
                    throw std::runtime_error("SurfaceMaterialBuilder::add_surface_from_tri3: micropolar element types not set");
                }

                switch (disp_el) {
                    case stark::Fem::ElementType::Tri3: {
                        handler.point_set = sim.deformables->point_sets->add(vertices);
                        if (rot_el == stark::Fem::ElementType::Tri3) {
                            auto rot_mesh = stark::Mesh<3> { vertices, triangles };
                            handler.strain_micropolar = this->add_mp<stark::Fem::BasisTri3, stark::Fem::BasisTri3>(sim, label, handler.point_set, triangles, rot_mesh);
                        } else {
                            invalid_combination();
                        }
                        break;
                    }
                    case stark::Fem::ElementType::Tri6: {
                        auto [vertices_quadratic, triangles_quadratic] = stark::tri3_to_tri6(vertices, triangles);
                        auto [vertices_quadratic_subdiv, triangle_quadratic_subdiv] = stark::tri6_to_tri3_subdivide(vertices_quadratic, triangles_quadratic);
                        handler.point_set = sim.deformables->point_sets->add(vertices_quadratic);

                        contact_mesh = triangle_quadratic_subdiv;
                        output_mesh = triangle_quadratic_subdiv;

                        if (rot_el == stark::Fem::ElementType::Tri3) {
                            auto rot_mesh = stark::Mesh<3> { vertices, triangles };
                            handler.strain_micropolar = this->add_mp<stark::Fem::BasisTri6, stark::Fem::BasisTri3>(sim, label, handler.point_set, triangles_quadratic, rot_mesh);
                        } else if (rot_el == stark::Fem::ElementType::Tri6) {
                            auto rot_mesh = stark::Mesh<6> { vertices_quadratic, triangles_quadratic};
                            handler.strain_micropolar = this->add_mp<stark::Fem::BasisTri6, stark::Fem::BasisTri6>(sim, label, handler.point_set, triangles_quadratic, rot_mesh);
                        } else {
                            invalid_combination();
                        }
                        break;
                    }
                        // TODO: Remaining element combinations
                    case stark::Fem::ElementType::Tri10: {
                        auto [vertices_quadratic, triangles_quadratic] = stark::tri3_to_tri6(vertices, triangles);
                        auto [vertices_cubic, triangles_cubic] = stark::tri3_to_tri10(vertices, triangles);
                        auto [vertices_cubic_subdiv, triangle_cubic_subdiv] = stark::tri10_to_tri3_subdivide(vertices_cubic, triangles_cubic);
                        handler.point_set = sim.deformables->point_sets->add(vertices_cubic);

                        contact_mesh = triangle_cubic_subdiv;
                        output_mesh = triangle_cubic_subdiv;

                        if (rot_el == stark::Fem::ElementType::Tri3) {
                            auto rot_mesh = stark::Mesh<3> { vertices, triangles };
                            handler.strain_micropolar = this->add_mp<stark::Fem::BasisTri10, stark::Fem::BasisTri3>(sim, label, handler.point_set, triangles_cubic, rot_mesh);
                        } else if (rot_el == stark::Fem::ElementType::Tri6) {
                            auto rot_mesh = stark::Mesh<6> { vertices_quadratic, triangles_quadratic};
                            handler.strain_micropolar = this->add_mp<stark::Fem::BasisTri10, stark::Fem::BasisTri6>(sim, label, handler.point_set, triangles_cubic, rot_mesh);
                        } else {
                            invalid_combination();
                        }
                        break;
                    }
                    default: {
                        invalid_combination();
                    }
                }
            } break;
            case StrainModel::None:
                break;
        };

        // Add output mesh
        if (this->strain_model != StrainModel::None && !label.empty()) {
            sim.deformables->output->add_triangle_mesh(label, handler.point_set, output_mesh);
        }

        // Add contact mesh
        if (this->strain_model != StrainModel::None && enable_contact.value_or(false)) {
            handler.contact = this->add_contact_tri3(sim, handler.point_set, contact_mesh);
            sim.deformables->output->add_triangle_mesh(label + "_collision_mesh", handler.point_set, contact_mesh);
        }

        // Add bending model
        switch (bending_model) {
            case BendingModel::DiscreteShells: {
                handler.bending_discrete_shells = this->add_discrete_shells(sim, *handler.point_set, triangles);
            } break;
            case BendingModel::None:
                break;
        };

        return handler;
    }

    GeneralSurfaceHandler add_surface_from_tri6(stark::Simulation& sim, const std::string& label, const std::vector<Eigen::Vector3d>& quadratic_vertices, const std::vector<std::array<int, 6>>& quadratic_triangles) const
    {
        GeneralSurfaceHandler handler;

        auto [subdiv_vertices, subdiv_triangles] = stark::tri6_to_tri3_subdivide(quadratic_vertices, quadratic_triangles);

        if (this->strain_model == StrainModel::Micropolar) {
            auto& disp_el = this->micropolar_element_type_displacement;
            auto& rot_el = this->micropolar_element_type_rotation;

            auto invalid_combination = [&]() { throw std::runtime_error(fmt::format("SurfaceMaterialBuilder::add_surface_from_tri6: unsupported element type combination ({} + {})", stark::Fem::element_type_to_string(disp_el), stark::Fem::element_type_to_string(rot_el))); };
            auto unimplemented_combination = [&]() { throw std::runtime_error(fmt::format("SurfaceMaterialBuilder::add_surface_from_tri6: unimplemented element type combination ({} + {})", stark::Fem::element_type_to_string(disp_el), stark::Fem::element_type_to_string(rot_el))); };

            if (!disp_el.has_value() || !rot_el.has_value()) {
                throw std::runtime_error("SurfaceMaterialBuilder::add_surface_from_tri6: micropolar element types not set");
            }

            if (disp_el == stark::Fem::ElementType::Tri3) {
                return this->add_surface_from_tri3(sim, label, subdiv_vertices, subdiv_triangles);
            } else if (disp_el == stark::Fem::ElementType::Tri6) {
                if (rot_el == stark::Fem::ElementType::Tri3) {
                    auto [linear_vertices, linear_triangles] = stark::tri6_to_tri3_coarsen(quadratic_vertices, quadratic_triangles);
                    auto rot_mesh = stark::Mesh<3> { linear_vertices, linear_triangles };
                    handler.point_set = sim.deformables->point_sets->add(quadratic_vertices);
                    handler.strain_micropolar = this->add_mp<stark::Fem::BasisTri6, stark::Fem::BasisTri3>(sim, label, handler.point_set, quadratic_triangles, rot_mesh);
                } else if (rot_el == stark::Fem::ElementType::Tri6) {
                    auto rot_mesh = stark::Mesh<6> { quadratic_vertices, quadratic_triangles };
                    handler.point_set = sim.deformables->point_sets->add(quadratic_vertices);
                    handler.strain_micropolar = this->add_mp<stark::Fem::BasisTri6, stark::Fem::BasisTri6>(sim, label, handler.point_set, quadratic_triangles, rot_mesh);
                } else {
                    invalid_combination();
                }
            } else {
                unimplemented_combination();
            }
        } else {
            return this->add_surface_from_tri3(sim, label, subdiv_vertices, subdiv_triangles);
        }

        // Add output mesh
        if (this->strain_model != StrainModel::None && !label.empty()) {
            sim.deformables->output->add_triangle_mesh(label, handler.point_set, subdiv_triangles);
        }

        // Add contact mesh
        if (this->strain_model != StrainModel::None && enable_contact.value_or(false)) {
            handler.contact = this->add_contact_tri3(sim, handler.point_set, subdiv_triangles);
            sim.deformables->output->add_triangle_mesh(label + "_collision_mesh", handler.point_set, subdiv_triangles);
        }

        return handler;
    }

    GeneralSurfaceHandler add_surface_from_quad4(stark::Simulation& sim, const std::string& label, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 4>>& quads) const
    {
        GeneralSurfaceHandler handler;

        // TODO: Implement quad4 to tri3 subdivision
        std::vector<std::array<int, 3>> output_mesh;
        std::vector<std::array<int, 3>> contact_mesh;

        if (this->strain_model == StrainModel::Micropolar) {
            auto& disp_el = this->micropolar_element_type_displacement;
            auto& rot_el = this->micropolar_element_type_rotation;

            auto invalid_combination = [&]() { throw std::runtime_error(fmt::format("SurfaceMaterialBuilder::add_surface_from_quad4: unsupported element type combination ({} + {})", stark::Fem::element_type_to_string(disp_el), stark::Fem::element_type_to_string(rot_el))); };
            auto unimplemented_combination = [&]() { throw std::runtime_error(fmt::format("SurfaceMaterialBuilder::add_surface_from_quad4: unimplemented element type combination ({} + {})", stark::Fem::element_type_to_string(disp_el), stark::Fem::element_type_to_string(rot_el))); };

            if (!disp_el.has_value() || !rot_el.has_value()) {
                throw std::runtime_error("SurfaceMaterialBuilder::add_surface_from_quad4: micropolar element types not set");
            }

            if (disp_el == stark::Fem::ElementType::Quad4) {
                handler.point_set = sim.deformables->point_sets->add(vertices);
                if (rot_el == stark::Fem::ElementType::Quad4) {
                    auto tri3_mesh = stark::quad4_to_tri3(vertices, quads);
                    output_mesh = tri3_mesh.conn;
                    contact_mesh = tri3_mesh.conn;

                    auto rot_mesh = stark::Mesh<4>{vertices, quads};
                    handler.strain_micropolar = this->add_mp<stark::Fem::BasisQuad4, stark::Fem::BasisQuad4>(sim, label, handler.point_set, quads, rot_mesh);
                } else {
                    invalid_combination();
                }
            } else if (disp_el == stark::Fem::ElementType::Quad9) {
                auto [vertices_quadratic, quads_quadratic] = stark::quad4_to_quad9(vertices, quads);
                auto [vertices_quadratic_subdiv, quads_quadratic_subdiv] = stark::quad9_to_quad4_subdivide(vertices_quadratic, quads_quadratic);
                auto [vertices_quadratic_subdiv_tri, tris_quadratic_subdiv] = stark::quad4_to_tri3(vertices_quadratic_subdiv, quads_quadratic_subdiv);
                handler.point_set = sim.deformables->point_sets->add(vertices_quadratic);

                contact_mesh = tris_quadratic_subdiv;
                output_mesh = tris_quadratic_subdiv;

                if (rot_el == stark::Fem::ElementType::Quad4) {
                    auto rot_mesh = stark::Mesh<4>{vertices, quads};
                    handler.strain_micropolar = this->add_mp<stark::Fem::BasisQuad9, stark::Fem::BasisQuad4>(sim, label,
                                                                                                             handler.point_set,
                                                                                                             quads_quadratic,
                                                                                                             rot_mesh);
                } else if (rot_el == stark::Fem::ElementType::Quad9) {
                    auto rot_mesh = stark::Mesh<9>{vertices_quadratic, quads_quadratic};
                    handler.strain_micropolar = this->add_mp<stark::Fem::BasisQuad9, stark::Fem::BasisQuad9>(sim, label,
                                                                                                             handler.point_set,
                                                                                                             quads_quadratic,
                                                                                                             rot_mesh);
                } else {
                    invalid_combination();
                }
            } else {
                unimplemented_combination();
            }
        } else {
            throw std::runtime_error("SurfaceMaterialBuilder::add_surface_from_quad4: only micropolar strain model supported");
        }

        // Add output mesh
        if (this->strain_model != StrainModel::None && !label.empty()) {
            sim.deformables->output->add_triangle_mesh(label, handler.point_set, output_mesh);
        }

        // Add contact mesh
        if (this->strain_model != StrainModel::None && enable_contact.value_or(false)) {
            handler.contact = this->add_contact_tri3(sim, handler.point_set, contact_mesh);
            sim.deformables->output->add_triangle_mesh(label + "_collision_mesh", handler.point_set, contact_mesh);
        }

        return handler;
    }

    stark::EnergyLumpedInertia::Handler add_lumped_inertia_tri3(stark::Simulation& sim, const stark::PointSetHandler& point_set, const std::vector<std::array<int, 3>>& triangles) const {
        return sim.deformables->lumped_inertia->add(point_set, triangles,
                                                    stark::EnergyLumpedInertia::Params()
                                                            .set_density(inertia_area_density)
                                                            .set_damping(inertia_damping)
        );
    }

    stark::EnergyFrictionalContact::Handler add_contact_tri3(stark::Simulation& sim, const stark::PointSetHandler& point_set, const std::vector<std::array<int, 3>>& triangles) const {
        return sim.interactions->contact->add_triangles(point_set, triangles,
                                                        stark::EnergyFrictionalContact::Params()
                                                                .set_contact_thickness(this->contact_distance)
        );
    }

    stark::EnergyDiscreteShells::Handler add_discrete_shells(stark::Simulation& sim, const stark::PointSetHandler& point_set, const std::vector<std::array<int, 3>>& triangles) const {
        return sim.deformables->discrete_shells->add(point_set, triangles,
                                                     stark::EnergyDiscreteShells::Params()
                                                             .set_elasticity_only(true)
                                                             .set_stiffness(bending_stiffness)
        );
    }

    stark::EnergyTriangleStrainKim20::Handler add_kim20(stark::Simulation& sim, const stark::PointSetHandler& point_set, const std::vector<std::array<int, 3>>& triangles) const {
        return sim.deformables->strain_kim_20->add(point_set, triangles,
                                                   stark::EnergyTriangleStrainKim20::Params()
                                                           .set_thickness(thickness)
                                                           .set_stiffness(strain_young_modulus)
                                                           .set_strain_limit(strain_limit)
                                                           .set_strain_limit_stiffness(strain_limit_stiffness)
        );
    }

    stark::EnergyTriangleStrainWen23::Handler add_wen23(stark::Simulation& sim, const stark::PointSetHandler& point_set, const std::vector<std::array<int, 3>>& triangles) const {
        return sim.deformables->strain_wen_23->add(point_set, triangles,
                                                   stark::EnergyTriangleStrainWen23::Params()
                                                           .set_thickness(thickness)
                                                           .set_youngs_modulus(strain_young_modulus)
                                                           .set_poissons_ratio(strain_poissons_ratio)
        );
    }

    template<typename DispEl, typename RotEl>
    stark::EnergyMicropolarShells::Handler add_mp(stark::Simulation& sim, const std::string& label, const stark::PointSetHandler& point_set, const std::vector<std::array<int, DispEl::num_nodes>>& elements, stark::Mesh<RotEl::num_nodes>& rotation_mesh) const {
        return sim.deformables->strain_micropolar_shells->add_mesh<DispEl, RotEl>(label, point_set, elements, rotation_mesh,
                                                                                  stark::EnergyMicropolarShells::Params()
                                                                                          .set_thickness(thickness)
                                                                                          .set_youngs_modulus(strain_young_modulus)
                                                                                          .set_poissons_ratio(strain_poissons_ratio)
                                                                                          .set_mu_c_factor(this->micropolar_mu_c_factor)
                                                                                          .set_length_scale(this->micropolar_length_scale)
                                                                                          .set_shear_correction(this->micropolar_shear_correction)
                                                                                          .set_bending_correction(this->micropolar_bending_correction)
                                                                                          .set_angular_inertia_density(this->micropolar_angular_inertia)
                                                                                          .set_density(this->inertia_area_density / this->thickness)
        );
    }
};

/// Helper type to initialize the simulator with default settings and override scene specific settings
struct ShellSceneSettings {
    ShellSceneSettings() = default;
    explicit ShellSceneSettings(const std::string& scene_name) : scene_name(scene_name) {}

    std::string scene_name = "micropolar_scene";
    std::string file_prefix = "";

    bool subfolder_per_simulation = false;

    Option<int> n_threads;
    Option<int> output_fps;
    Option<bool> adaptive_time_step;
    Option<double> max_timestep;
    Option<bool> project_to_pd;
    Option<double> residual_tol;
    Option<int> max_newton_iter;
    Option<stark::LinearSystemSolver> solver;
    Option<SurfaceMaterialBuilder::StrainModel> shell_model;

    Option<bool> mp_never_project_to_pd;
    Option<bool> mp_use_full_model;
    Option<bool> mp_use_quaternion_gamma;
    Option<bool> mp_use_nonlinear_volume;

    Option<double> thickness;
    Option<double> youngs_modulus;

    Option<double> mp_mu_c_factor;
    Option<double> mp_length_scale;
    Option<double> mp_shear_correction;
    Option<double> mp_bending_correction;

    Option<stark::Fem::ElementType> mp_element_displacement;
    Option<stark::Fem::ElementType> mp_element_rotation;

    bool write_quadrature_mesh = false;
    bool write_subdivided_mesh = false;
    int mesh_subdivision_level = 3;

    /// Returns the output directory for the scene
    std::string get_output_directory(const std::string& custom_scene_name = "") const {
        if (this->subfolder_per_simulation) {
            return fmt::format("{}/{}/{}", OUTPUT_PATH, scene_name, this->get_simulation_name(custom_scene_name));
        } else {
            return fmt::format("{}/{}", OUTPUT_PATH, scene_name);
        }
    }

    /// Returns the simulation name for the scene
    std::string get_simulation_name(const std::string& custom_scene_name = "") const {
        return fmt::format(
                "{}{}{}",
                file_prefix,
                file_prefix.empty() || custom_scene_name.empty() ? "" : "_",
                custom_scene_name
        );
    }

    /// Returns a string representing of the element types specified by the settings
    std::string get_element_types_string() const {
        if (!this->mp_element_displacement.has_value() || !this->mp_element_rotation.has_value()) {
            throw std::runtime_error("ShellSceneSettings::get_element_types_string: called without overriding element types");
        }
        return stark::Fem::element_type_to_string(this->mp_element_displacement) + stark::Fem::element_type_to_string(this->mp_element_rotation);
    }

    /// Returns simulator settings initialized with default values
    stark::Settings simulator_settings(const std::string& custom_scene_name = "") const {
        stark::Settings settings = stark::Settings();
        // Default settings
        settings.output.simulation_name = this->get_simulation_name(custom_scene_name);
        settings.output.output_directory = this->get_output_directory(custom_scene_name);
        settings.output.codegen_directory = COMPILE_PATH;
        settings.output.console_verbosity = stark::ConsoleVerbosity::TimeSteps;
        settings.output.fps = 30;
        settings.debug.symx_check_for_NaNs = true;
        settings.debug.symx_finite_difference_check = false;
        settings.simulation.init_frictional_contact = false;
        settings.newton.project_to_PD = false;
        settings.newton.linear_system_solver = stark::LinearSystemSolver::DirectIndefinite;
        settings.newton.residual = { stark::ResidualType::Force, 1e-7 };
        settings.newton.max_line_search_iterations = 20;
        settings.newton.enable_noise_resistant_line_search = true;
        settings.newton.enable_flipping_on_non_descent = true;

        settings.models.never_project_mp_shell = true;
        settings.models.never_project_tri_wen23 = true;

        if (this->shell_model.has_value()) {
            switch (this->shell_model) {
                case SurfaceMaterialBuilder::StrainModel::Micropolar: {
                    settings.models.enable_model_mp_shell = true;
                    settings.models.enable_model_wen23 = false;
                    break;
                }
                case SurfaceMaterialBuilder::StrainModel::Wen23: {
                    settings.models.enable_model_mp_shell = false;
                    settings.models.enable_model_wen23 = true;
                    break;
                }
            }
        }


        return settings;
    }

    /// Applies the settings overrides to the simulator settings
    void apply_overrides_to(stark::Settings& settings) const {
        this->n_threads.apply_to(settings.execution.n_threads);
        this->output_fps.apply_to(settings.output.fps);
        this->adaptive_time_step.apply_to(settings.simulation.use_adaptive_time_step);
        this->max_timestep.apply_to(settings.simulation.max_time_step_size);
        this->project_to_pd.apply_to(settings.newton.project_to_PD);
        this->solver.apply_to(settings.newton.linear_system_solver);
        this->residual_tol.apply_to(settings.newton.residual.tolerance);
        this->max_newton_iter.apply_to(settings.newton.max_newton_iterations);
        this->mp_never_project_to_pd.apply_to(settings.models.never_project_mp_shell);

        this->mp_use_quaternion_gamma.apply_to(settings.models.enable_model_mp_shell_use_quaternion_gamma);
        this->mp_use_nonlinear_volume.apply_to(settings.models.enable_model_mp_shell_use_nonlinear_volume_terms);
        this->mp_use_full_model.apply_to(settings.models.enable_model_mp_shell_full_rest_curvature_terms);
    }

    /// Applies the settings overrides to the micropolar shell model
    void apply_overrides_to(stark::EnergyMicropolarShells* mp) const {
        if (mp != nullptr) {
            mp->WRITE_QUADRATURE_MESHES = this->write_quadrature_mesh;
            mp->WRITE_SUBDIVIDED_MESHES = this->write_subdivided_mesh;
            mp->MESH_SUBDIVISION_LEVEL = this->mesh_subdivision_level;
        }
    }

    void apply_overrides_to(SurfaceMaterialBuilder& material) const {
        this->thickness.apply_to(material.thickness);
        this->youngs_modulus.apply_to(material.strain_young_modulus);
        this->mp_mu_c_factor.apply_to(material.micropolar_mu_c_factor);
        this->mp_length_scale.apply_to(material.micropolar_length_scale);
        this->mp_shear_correction.apply_to(material.micropolar_shear_correction);
        this->mp_bending_correction.apply_to(material.micropolar_bending_correction);

        this->mp_element_displacement.apply_to(material.micropolar_element_type_displacement);
        this->mp_element_rotation.apply_to(material.micropolar_element_type_rotation);
    }
};

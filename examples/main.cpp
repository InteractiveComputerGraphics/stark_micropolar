#include <fmt/format.h>

#include <stark>

#include "paths.h"
#include "scene_helper.h"

/// Default boundary condition parameters
constexpr auto bc_params = [](const double stiffness = 1e6){
    return stark::EnergyPrescribedPositions::Params().set_stiffness(stiffness);
};

/// Main scenes from the paper "Curved Three-Director Cosserat Shells with Strong Coupling", 2024
/// by F. Löschner, J. A. Fernández-Fernández, S. R. Jeske and J. Bender
namespace paper_scenes
{
void circle_growth_mp()
{
    ShellSceneSettings scene_settings("circle_growth_no_contacts_mp");

    enum class GrowthOrigin {
        Inner,
        Outer,
    };

    struct Parametrization {
        std::string file_prefix;
        GrowthOrigin growth_origin;
        int mesh_index;
        int sim_subdivision;
        bool quadratic_elements;
        int render_subdivision;
    };

    std::vector<Parametrization> parametrizations = {
            { "", GrowthOrigin::Inner, 3, 0, true, 3 },
            { "", GrowthOrigin::Outer, 3, 0, true, 3 },
    };

    //scene_settings.write_quadrature_mesh = true;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.mp_use_quaternion_gamma = true;
    scene_settings.mp_use_nonlinear_volume = false;
    scene_settings.mp_use_full_model = false;

    //scene_settings.mp_never_project_to_pd = false;
    //scene_settings.project_to_pd = true;

    const auto scene = [&](const ShellSceneSettings& scene_settings, GrowthOrigin growth_origin, const stark::Mesh<3>& mesh, int subdivision){
        stark::Settings settings = scene_settings.simulator_settings(fmt::format("{}_{}tris_n{}", growth_origin == GrowthOrigin::Inner ? "inner" : "outer", mesh.conn.size(), subdivision));

        settings.models.enable_model_mp_shell = true;
        settings.execution.end_simulation_time = 8.0;
        //settings.simulation.init_frictional_contact = true;
        //settings.simulation.gravity = {0.0, 0.0, -9.81};
        settings.simulation.gravity = {0.0, 0.0, 0.0};

        //settings.simulation.max_time_step_size = 1e-3;
        //settings.simulation.adaptive_time_step.set(0.0, 0.001, 0.020);

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);

        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Global Contact Settings
        /*
        const double contact_distance  = 0.0001;
        const double contact_mu = 1.0;

        simulation.interactions->contact->set_global_params(
        stark::EnergyFrictionalContact::GlobalParams()
            .set_default_contact_thickness(contact_distance)
            .set_friction_stick_slide_threshold(0.01)
            .set_min_contact_stiffness(1e7)
        );
        */

        // Mesh
        const double r = 0.025;
        auto [vertices, triangles] = mesh;

        if (subdivision > 0) {
            stark::Mesh<3> subdivided_mesh;
            subdivided_mesh.vertices = vertices;
            subdivided_mesh.conn = triangles;

            for (int i = 0; i < subdivision; ++i) {
                const auto tri6_mesh = stark::tri3_to_tri6(subdivided_mesh.vertices, subdivided_mesh.conn);
                subdivided_mesh = stark::tri6_to_tri3_subdivide(tri6_mesh.vertices, tri6_mesh.conn);
            }

            vertices = subdivided_mesh.vertices;
            triangles = subdivided_mesh.conn;
        }

        const auto edges = stark::find_edges_from_simplices(triangles, vertices.size());
        double dx = std::numeric_limits<double>::max();
        for (const auto& e: edges) {
            const auto& v1 = vertices[e[0]];
            const auto& v2 = vertices[e[1]];
            dx = std::min(dx, (v2-v1).norm());
        }

        fmt::print("Shortest edge (dx): {:.e}\n", dx);

        auto material = SurfaceMaterialBuilder();
        material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;
        //material.enable_contact = true;
        //material.contact_distance = contact_distance;

        material.thickness = 3e-4;
        material.inertia_area_density = 250.0 * material.thickness;
        material.inertia_damping = 2.0;

        material.strain_young_modulus = 1e6;
        material.strain_poissons_ratio = 0.4;

        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_length_scale = 1e-3 * material.thickness;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        material.micropolar_element_type_displacement = stark::Fem::ElementType::Tri6;
        material.micropolar_element_type_rotation = stark::Fem::ElementType::Tri3;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);

        // Project to sphere and plane
        for (int i = 0; i < shell.point_set->size(); ++i) {
            auto v = shell.point_set->get_rest_position(i);
            Eigen::Vector3d v_new = v.normalized() * r;
            v_new.z() = 0.0;
            shell.point_set->set_rest_position(i, v_new);
        }

        // Project boundary to circle
        std::vector<int> boundary;

        {
            auto [quadratic_vertices, quadratic_triangles] = stark::tri3_to_tri6(vertices, triangles);
            auto [subdiv_vertices, subdiv_triangles] = stark::tri6_to_tri3_subdivide(quadratic_vertices, quadratic_triangles);

            auto [boundary_edges, map] = stark::find_perimeter_edges(subdiv_triangles, subdiv_vertices.size());
            for (const auto& e : boundary_edges) {
                for (const int i : e) {
                    boundary.push_back(map[i]);
                    auto v = shell.point_set->get_rest_position(map[i]);
                    Eigen::Vector2d v_new_2d = v.segment(0, 2).normalized() * r;
                    Eigen::Vector3d v_new = {v_new_2d.x(), v_new_2d.y(), 0.0};
                    shell.point_set->set_rest_position(map[i], v_new);
                }
            }
        }

        // Initial transformation
        for (int i = 0; i < shell.point_set->size(); ++i) {
            auto x_rest = shell.point_set->get_rest_position(i);
            shell.point_set->set_position(i, x_rest);
        }
        shell.point_set->add_rotation(1.5, Eigen::Vector3d::UnitX(), Eigen::Vector3d::Zero(), false);
        shell.point_set->add_displacement(0.01 * Eigen::Vector3d::UnitZ(), false);

        std::vector<Eigen::Vector3d> X_quad = shell.strain_micropolar->interpolate_quadrature_points_mesh();
        std::vector<Eigen::Matrix3d> F0_orig;

        /*
        // Floor
        const double floor_size = 0.5;

        auto [V, T] = stark::generate_triangle_grid({ 0.0, 0.0 }, { floor_size, floor_size }, { 10, 10 });
        auto HF = simulation.deformables->point_sets->add(V);
        auto BC = simulation.deformables->prescribed_positions->add(HF, HF.all(),
            stark::EnergyPrescribedPositions::Params()
                .set_stiffness(1e6)
                //.set_tolerance(0.002)
        );
        simulation.deformables->output->add_triangle_mesh("floor", HF, T);
        auto contact_handle_floor = simulation.interactions->contact->add_triangles(HF, T, { contact_distance });

        simulation.interactions->contact->set_friction(shell.contact, contact_handle_floor, contact_mu);
        */

        const double t_end = 5;

        const auto compute_growth = [&simulation,r,growth_origin,t_end](const Eigen::Vector3d& X) -> double {
            const double t = simulation.get_time();
            const double factor = stark::blend(0.0, 1.0, 1.0, t_end, t, stark::BlendType::Linear);
            const double max_growth = 0.15;
            const double current_growth = factor * max_growth;

            const bool from_middle = growth_origin == GrowthOrigin::Inner;

            const double stop_factor = from_middle ? 0.9 : 0.4;
            const double r_stop = from_middle ? r * stop_factor : r * (1.0 - stop_factor);

            const double r_i = X.segment(0, 2).norm();
            const double r_i_rel = from_middle ? std::min(r_i, r_stop) / r_stop : std::max(r_i - r_stop, 0.0) / (r - r_stop);
            const double alpha = from_middle ? std::max(1.0 - r_i_rel, 0.0) : r_i_rel;

            const double s = current_growth * alpha;
            const double growth = std::pow(2.71828, 2.0 * s);

            return growth;
        };

        mp->scalar_field_callbacks.push_back({
                "growth",
                compute_growth
        });

        simulation.add_time_event(1.0, t_end, [&](double t) {
            for (int i = 0; i < X_quad.size(); ++i) {
                const double growth = compute_growth(X_quad[i]);
                //fmt::print("{} r_i {} r_i_rel {} alpha {} s {}\n", i, r_i, r_i_rel, alpha, s);
                const Eigen::Matrix3d F0 = F0_orig[i];
                shell.strain_micropolar->set_rest_deformation_gradient(i, growth * F0);
            }
        });

        // Run
        simulation.run([&] {
            if (simulation.get_time() > 0) {
                //std::abort();
            }

            if (F0_orig.empty()) {
                for (int i = 0; i < X_quad.size(); ++i) {
                    F0_orig.push_back(shell.strain_micropolar->get_rest_deformation_gradient(i));
                }
            }
        });
    };

    for (const auto& params : parametrizations) {
        scene_settings.file_prefix = params.file_prefix;
        scene_settings.mesh_subdivision_level = params.render_subdivision;
        const auto mesh = stark::load_obj(fmt::format("{}/icosphere_cap_{}.obj", MODELS_PATH, params.mesh_index)).back();
        scene(scene_settings, params.growth_origin, mesh, params.sim_subdivision);
    }
}

void curvature_modes_mp()
{
    ShellSceneSettings scene_settings("curvature_modes_mp");
    scene_settings.write_quadrature_mesh = false;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.mesh_subdivision_level = 4;
    scene_settings.output_fps = 60.0;

    std::vector<double> curvature_factors = {
            1.0, 1.0, 1.0,
            3.0, 2.0, 2.0,
    };

    std::vector<int> axes_to_run = {
            0,
            1,
            2,
            3,
            4,
            5,
    };

    const auto scene = [](const ShellSceneSettings& scene_settings, int n, int axis, double curvature_factor) { ;
        stark::Settings settings = scene_settings.simulator_settings();
        settings.output.simulation_name = fmt::format("{}_n_{}_axis_{}", scene_settings.get_simulation_name(), n, axis);

        settings.models.enable_model_mp_shell = true;
        settings.execution.end_simulation_time = 5.0;
        settings.simulation.gravity = {0.0, 0.0, 0.0};

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);

        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Mesh
        const double w = 0.125;
        const int length_factor = 8;
        const double l = length_factor * w;
        auto [vertices, triangles] = stark::generate_triangle_grid(
                { 0.0, 0.0},
                { l, w },
                { length_factor*n, n }
        );
        const double dx = w / n;

        auto material = SurfaceMaterialBuilder();
        material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

        material.thickness = 1e-3;
        material.inertia_area_density = 1000.0 * material.thickness;
        material.inertia_damping = 2.0;
        material.strain_young_modulus = 1e4;
        material.strain_poissons_ratio = 0.3;
        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_length_scale = 1;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        material.micropolar_element_type_displacement = stark::Fem::ElementType::Tri6;
        material.micropolar_element_type_rotation = stark::Fem::ElementType::Tri3;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);

        // BC: Fix central line
        const Eigen::Vector3d aabb_pos = {0.0, 0.0, 0.0};
        const Eigen::Vector3d aabb_size = {0.1*dx, 1.1*w, 0.1*dx};
        simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, aabb_pos, aabb_size, bc_params());

        std::vector<Eigen::Vector3d> X_quad = shell.strain_micropolar->interpolate_quadrature_points_mesh();

        // Run
        simulation.run([&]() {
            const double t_end = 2.5;
            const double t = simulation.get_time();

            const double factor = std::min(t, t_end) / t_end;
            const double max_curvature = 2.0*M_PI;
            const double current_curvature = factor * max_curvature;

            for (int i = 0; i < X_quad.size(); ++i) {
                const auto& X_q = X_quad[i];
                const double x = X_q.x();
                const double x_rel = x/(0.5 * l);
                //double x_rel = 0.0;

                // Set only dRy/dx
                //const double c = std::min(1.0-std::abs(x_rel), 1.0) * current_curvature;
                const double c = current_curvature * (x_rel > 0 ? 1.0 : -1.0);

                std::array<double, 9> array = {
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                };
                array[axis] = curvature_factor * c;

                stark::MatrixWrapper3d curvature = stark::MatrixWrapper3d::from_matrix(stark::MatrixWrapper3d::from_array(array).to_matrix().transpose());
                shell.strain_micropolar->set_rest_curvature(i, curvature.to_matrix());
            }

            //fmt::print("current_curvature: {}\n", current_curvature);
        });
    };

    for (int axis : axes_to_run) {
        scene(scene_settings, 2, axis, curvature_factors[axis]);
    }
    //scene(scene_settings, 20);
}

void armadillo_twisting_mp()
{
    ShellSceneSettings scene_settings("armadillo_twisting_mp");

    struct Parametrization {
        std::string file_prefix;
        std::string mesh;
        int curvature_index;
        int sim_subdivision;
        bool full_model;
        bool use_quadratic_elements;
        int render_subdivision;
    };

    std::vector<Parametrization> parametrizations = {
            { "", "armadillo_lvl2_verts7743", 2, 0, false, true, 1 },
            { "", "armadillo_lvl2_verts7743", 5, 0, false, true, 1 },
            { "", "armadillo_lvl2_verts7743", 8, 0, false, true, 1 },
    };

    scene_settings.file_prefix = "";
    //scene_settings.write_quadrature_mesh = true;
    scene_settings.output_fps = 60.0;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.mp_use_quaternion_gamma = true;
    scene_settings.mp_use_nonlinear_volume = false;
    scene_settings.subfolder_per_simulation = true;
    scene_settings.residual_tol = 1e-7;

    const bool with_contact = false;

    const auto scene = [&](const ShellSceneSettings& scene_settings, const std::string& mesh, int subdivision, int curvature_index){
        stark::Settings settings = scene_settings.simulator_settings(fmt::format("c{}_{}_n{}", curvature_index, mesh, subdivision));

        settings.execution.end_simulation_time = 4.0;
        settings.models.enable_model_mp_shell = true;
        settings.simulation.gravity = {0.0, 0.0, 0.0};

        settings.simulation.init_frictional_contact = with_contact;
        settings.newton.project_to_PD = with_contact;

        settings.simulation.max_time_step_size = 5e-3;
        //settings.simulation.adaptive_time_step.set(0.0, 0.001, 0.020);

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);
        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Global Contact Settings
        const double contact_distance  = 1e-4;
        const double contact_mu = 1.0;

        simulation.interactions->contact->set_global_params(
                stark::EnergyFrictionalContact::GlobalParams()
                        .set_default_contact_thickness(contact_distance)
                        .set_friction_stick_slide_threshold(0.01)
                        .set_min_contact_stiffness(1e6)
        );

        // Mesh
        const double height = 0.2;
        auto meshes = stark::load_obj(fmt::format("{}/{}.obj", MODELS_PATH, mesh));
        auto [vertices, triangles] = meshes.back();

        const auto rot = Eigen::AngleAxisd(stark::deg2rad(90), Eigen::Vector3d::UnitX());
        Eigen::AlignedBox3d aabb;
        for (auto& v : vertices) {
            v = rot * v;
            aabb.extend(v);
        }

        const Eigen::Vector3d center = aabb.center();
        const double scaling = height / aabb.sizes().z();
        aabb = Eigen::AlignedBox3d();

        for (auto& v : vertices) {
            v -= center;
            v = scaling * v;
            aabb.extend(v);
        }

        if (subdivision > 0) {
            stark::Mesh<3> subdivided_mesh;
            subdivided_mesh.vertices = vertices;
            subdivided_mesh.conn = triangles;

            for (int i = 0; i < subdivision; ++i) {
                const auto tri6_mesh = stark::tri3_to_tri6(subdivided_mesh.vertices, subdivided_mesh.conn);
                subdivided_mesh = stark::tri6_to_tri3_subdivide(tri6_mesh.vertices, tri6_mesh.conn);
            }

            vertices = subdivided_mesh.vertices;
            triangles = subdivided_mesh.conn;
        }

        const auto edges = stark::find_edges_from_simplices(triangles, vertices.size());
        double dx = std::numeric_limits<double>::max();
        for (const auto& e: edges) {
            const auto& v1 = vertices[e[0]];
            const auto& v2 = vertices[e[1]];
            dx = std::min(dx, (v2-v1).norm());
        }

        fmt::print("Shortest edge (dx): {:.e}\n", dx);

        auto material = SurfaceMaterialBuilder();
        material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;
        if (with_contact) {
            material.enable_contact = true;
            material.contact_distance = contact_distance;
        }

        material.thickness = 1e-3;
        material.inertia_area_density = 1000.0 * material.thickness;
        material.inertia_damping = 0.0;

        material.strain_young_modulus = 1e5;
        material.strain_poissons_ratio = 0.4;

        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_length_scale = 0.1;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);
        //simulation.interactions->contact->set_friction(shell.contact, shell.contact, contact_mu);

        const Eigen::Vector3d bc_center = {0.0, 0.0, aabb.min().z()};
        const Eigen::Vector3d bc_size = {1.1*aabb.sizes().x(), 1.1*aabb.sizes().y(), 0.01};

        Eigen::AlignedBox3d bc_aabb(bc_center - 0.5*bc_size, bc_center + 0.5*bc_size);
        std::vector<int> bc_points;
        for (int i = 0; i < shell.point_set->size(); i++) {
            if (bc_aabb.contains(shell.point_set->get_position(i))) {
                bc_points.push_back(i);
            }
        }
        //simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, bc_center, bc_size, bc_params(1e7));

        // Floor
        const double size = 5.0 * height;
        const double floor_height = 0.01;
        auto [floor_vertices, floor_triangles, floor] = simulation.presets->rigidbodies->add_box("floor", 1.0, { size, size, 0.01 });
        floor.rigidbody.set_translation({ 0.0, 0.0, aabb.min().z() - 0.5 * floor_height - 2.0*contact_distance });
        simulation.rigidbodies->add_constraint_fix(floor.rigidbody);

        //shell.point_set->add_rotation(45.0, Eigen::Vector3d::UnitY(), Eigen::Vector3d::Zero(), false);
        simulation.deformables->prescribed_positions->add(shell.point_set, bc_points, bc_params(1e7));

        //simulation.interactions->contact->set_friction(shell.contact, contact_handle_floor, contact_mu);

        std::vector<Eigen::Vector3d> X_quad = shell.strain_micropolar->interpolate_quadrature_points_mesh();

        const double t_end = 2.0;
        const double c_max = 1 * M_PI / 0.2;

        simulation.add_time_event(0.0, t_end, [&](double t) {
            const double c = stark::blend(0.0, c_max, 0.0, t_end, t, stark::BlendType::Linear);

            for (int i = 0; i < X_quad.size(); ++i) {
                //const double x = X_q.x();
                //const double x_rel = (x + l) / (2.0 * l);

                // Set only dRy/dx
                std::array<double, 9> curvature_arr = {
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                };
                curvature_arr[curvature_index] = -c;
                stark::MatrixWrapper3d curvature = stark::MatrixWrapper3d::from_array(curvature_arr);
                shell.strain_micropolar->set_rest_curvature(i, curvature.to_matrix());
            }
        });

        // Run
        simulation.run([&] {
            if (simulation.get_time() > 0) {
                //std::abort();
            }
        });
    };

    for (const auto& params : parametrizations) {
        scene_settings.file_prefix = params.file_prefix;
        scene_settings.mp_use_full_model = params.full_model;
        scene_settings.mesh_subdivision_level = params.render_subdivision;
        scene_settings.mp_element_displacement = params.use_quadratic_elements ? stark::Fem::ElementType::Tri6 : stark::Fem::ElementType::Tri3;
        scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
        scene(scene_settings, params.mesh, params.sim_subdivision, params.curvature_index);
    }
}

void plate_twisting_mp_refinement(const ShellSceneSettings& scene_settings, int ny, bool symmetric)
{
    stark::Settings settings = scene_settings.simulator_settings(fmt::format("{}_ny_{}", scene_settings.get_element_types_string(), ny));

    settings.models.enable_model_mp_shell = true;
    settings.execution.end_simulation_time = 4.0;
    // Turn projection off by default. Causes problems for this scene.
    settings.newton.project_to_PD = false;
    settings.simulation.max_time_step_size = 0.01;
    settings.simulation.use_adaptive_time_step = false;
    settings.simulation.gravity = {0.0, 0.0, 0.0};

    scene_settings.apply_overrides_to(settings);
    stark::Simulation simulation(settings);
    auto mp = simulation.deformables->strain_micropolar_shells.get();
    scene_settings.apply_overrides_to(mp);

    // Mesh
    const double l = 100 * 1e-3;
    const double w = 10 * 1e-3;
    const double dx = w/ny;
    auto [vertices, triangles] = stark::generate_triangle_grid(
            { 0.0, 0.0 },
            { l, w },
            { 10 * ny, ny });

    auto material = SurfaceMaterialBuilder();
    material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

    // Original parameters from SNB16 twisting
    material.thickness = 2.0 * 1e-3;
    material.inertia_area_density = 1000.0 * material.thickness;
    material.inertia_damping = 5.0;

    material.strain_young_modulus = 12.86 * 1e9;
    material.strain_poissons_ratio = 0.1393;
    material.micropolar_mu_c_factor = 1.0;
    material.micropolar_shear_correction = 1.0;
    material.micropolar_length_scale = material.thickness * 1e-3;
    material.micropolar_bending_correction = 1.0;
    material.micropolar_angular_inertia = 0.0;

    // Try out softer parameters
    if (false) {
        material.thickness = 2.0 * 1e-3;
        material.inertia_area_density = 1000.0 * material.thickness;

        material.strain_young_modulus = 1e7;
        material.strain_poissons_ratio = 0.3;
        material.micropolar_length_scale = material.thickness * 1e-3;
    }

    // Try out parameters of aluminium. Source: Wolfram
    if (true) {
        material.inertia_area_density = material.thickness * 2700;
        material.strain_young_modulus = 70 * 1e9;
        material.strain_poissons_ratio = 0.35;
    }

    scene_settings.apply_overrides_to(material);
    auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);

    const Eigen::Vector3d aabb_center_left = { -l * 0.5, 0.0, 0.0 };
    const Eigen::Vector3d aabb_center_right = { l * 0.5, 0.0, 0.0 };
    const Eigen::Vector3d aabb_size = {dx, 1.1 * w, dx};

    const double omega_deg = (symmetric) ? 180.0 : 360.0;
    const double omega_rad = omega_deg / 180.0 * M_PI;

    // Needed because of high stiffness of material
    const double bc_stiffness = 1.5e10;
    // Lower stiffness appears to be fine for rotations
    const double bc_stiffness_rot = 1e4;

    auto bc_left = simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set,
                                                                                 aabb_center_left,
                                                                                 aabb_size,
                                                                                 stark::EnergyPrescribedPositions::Params().set_stiffness(bc_stiffness));

    auto bc_right = simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set,
                                                                                  aabb_center_right,
                                                                                  aabb_size,
                                                                                  stark::EnergyPrescribedPositions::Params().set_stiffness(bc_stiffness));


    auto bc_left_rot = shell.strain_micropolar->prescribe_rotation_inside_aabb(aabb_center_left, aabb_size, Eigen::Matrix3d::Identity(), {false, false, true}, bc_stiffness_rot);
    auto bc_right_rot = shell.strain_micropolar->prescribe_rotation_inside_aabb(aabb_center_right, aabb_size, Eigen::Matrix3d::Identity(), {false, false, true}, bc_stiffness_rot);

    // Run
    simulation.run([&]() {
        // Make sure to apply the transformation for the *end* of the time step
        const double t = simulation.get_time() + simulation.get_time_step_size();
        const double angle = t*omega_rad;

        Eigen::AngleAxisd rot_aa = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX());
        Eigen::Matrix3d rot = rot_aa.toRotationMatrix();

        bc_right.set_transformation(Eigen::Vector3d::Zero(), rot);
        shell.strain_micropolar->set_prescribed_rotation(bc_right_rot, rot);

        if (symmetric) {
            bc_left.set_transformation(Eigen::Vector3d::Zero(), rot.transpose());
            shell.strain_micropolar->set_prescribed_rotation(bc_left_rot, rot.transpose());
        }
    });
}
void plate_twisting_mp_refinement_study() {
    ShellSceneSettings scene_settings("plate_twisting_mp_refinement");
    scene_settings.file_prefix = "twisting";

    scene_settings.output_fps = 60;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.write_quadrature_mesh = false;
    scene_settings.subfolder_per_simulation = true;

    struct Parametrization {
        bool use_quadratic_elements;
        int ny;
        int render_subdivisions;
    };

    std::vector<Parametrization> parametrizations = {
            //{true, 1, 4},
            //{true, 2, 4},
            {true, 4, 2},
            // {true, 8, 3},
            // {true, 16, 2},
            // {true, 32, 1},
            //{false, 1, 0},
            // {false, 2, 0},
            // {false, 4, 0},
            // {false, 8, 0},
            // {false, 16, 0},
            // {false, 32, 0},
            // {false, 64, 0},
            // {false, 128, 0},
    };

    for (const auto& params : parametrizations ) {
        scene_settings.mp_element_displacement = params.use_quadratic_elements ? stark::Fem::ElementType::Tri6 : stark::Fem::ElementType::Tri3;
        scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
        scene_settings.mesh_subdivision_level = params.render_subdivisions;
        plate_twisting_mp_refinement(scene_settings, params.ny, false);
    }
}

void plate_loading_mp_refinement(const ShellSceneSettings& scene_settings, int nx, int ny)
{
    stark::Settings settings = scene_settings.simulator_settings(fmt::format("{}_nx_{}_ny_{}", scene_settings.get_element_types_string(), nx, ny));

    settings.models.enable_model_mp_shell = true;
    settings.execution.end_simulation_time = 2.0;
    // Adaptive timestep and projection not needed for this scene
    settings.newton.project_to_PD = false;
    settings.simulation.use_adaptive_time_step = false;
    settings.simulation.max_time_step_size = 0.01;
    settings.simulation.gravity = {0.0, 0.0, 0.0};

    scene_settings.apply_overrides_to(settings);
    stark::Simulation simulation(settings);
    auto mp = simulation.deformables->strain_micropolar_shells.get();
    scene_settings.apply_overrides_to(mp);

    // Mesh
    const double l = 0.1;
    const double w = 0.01;
    auto [vertices, triangles] = stark::generate_triangle_grid(
            { 0.0, 0.0},
            { l, w },
            { nx, ny }
    );
    const double dx = l / nx;

    auto material = SurfaceMaterialBuilder();
    material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

    material.inertia_damping = 1.0;

    // Parameters of L-shaped beam in SNB16
    {
        material.thickness = 0.6 * 1e-3;						// 0.6 mm
        material.inertia_area_density = 1000.0 * material.thickness;
        material.strain_young_modulus = 71240 * 1e6;			// 71240 N/mm^2 = 71GPa ~ Young's modulus of aluminium
        material.strain_poissons_ratio = 0.31;
        //material.micropolar_mu_c_factor = 0.0;
        material.micropolar_length_scale = 0.6 * 1e-3 * 1e-3;	// 0.6 * 1e-3 mm
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        material.micropolar_mu_c_factor = 1.0;
    }

    if (true) {
        // Aluminium plate
        material.thickness = 2.0 * 1e-3;
        material.inertia_area_density = 2700 * material.thickness;
        material.strain_young_modulus = 70 * 1e9;
        material.strain_poissons_ratio = 0.35;
        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_length_scale = material.thickness * 1e-3;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;
    }

    scene_settings.apply_overrides_to(material);
    auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);

    std::vector<int> force_vertices;
    for (int i = 0; i < shell.point_set->size(); ++i) {
        if (shell.point_set->get_rest_position(i).x() > 0.5*l - 0.1*dx) {
            force_vertices.push_back(i);
        }
    }
    int n_force_verts = force_vertices.size();
    fmt::print("Number of vertices where force will be applied: {}\n", n_force_verts);

    // BC: Clamp on left side
    const Eigen::Vector3d aabb_center = { -0.5*l, 0.0, 0.0 };
    const Eigen::Vector3d aabb_dim = {0.1*dx, 1.1*w, 0.1*dx};
    simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, aabb_center, aabb_dim, bc_params());
    shell.strain_micropolar->prescribe_rotation_inside_aabb(aabb_center, aabb_dim, Eigen::Matrix3d::Identity(), {false, false, true}, 1e4);

    // Run
    fmt::print("\n");
    simulation.run([&]() {
        const double target_force = 50.0;
        //const double target_force = 1.80;
        const double t_end = 0.5 * settings.execution.end_simulation_time;
        //const double t_end = 5.0;
        const double t = simulation.get_time();

        const double factor = std::min(t, t_end) / t_end;
        const double total_force = factor * target_force;
        const double nodal_force = total_force / n_force_verts;
        //fmt::print("nodal_force={:.6e}", nodal_force);

        for (const int i : force_vertices) {
            shell.point_set->set_force(i, {0.0, 0.0, nodal_force});
        }
    });
}
void plate_loading_mp_refinement_study() {
    ShellSceneSettings scene_settings("plate_loading_mp_refinement");
    scene_settings.file_prefix = "loading";

    scene_settings.output_fps = 60;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.write_quadrature_mesh = false;
    scene_settings.subfolder_per_simulation = true;

    struct Parametrization {
        int ny = 1;
        bool use_quadratic_elements = true;
        int mesh_subdivision_level;
    };

    std::vector<Parametrization> parametrizations = {
            {  1,  true, 3},
            {  2,  true, 3},
            {  4,  true, 2},
            {  8,  true, 1},
            { 16,  true, 1},
            { 32,  true, 0},
            { 64,  true, 0},
            // {128,  true, 0},
            {  1, false, 0},
            {  2, false, 0},
            {  4, false, 0},
            {  8, false, 0},
            { 16, false, 0},
            { 32, false, 0},
            { 64, false, 0},
            // {128, false, 0},
    };

    for (const auto& params : parametrizations) {
        scene_settings.mp_element_displacement = params.use_quadratic_elements ? stark::Fem::ElementType::Tri6 : stark::Fem::ElementType::Tri3;
        scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
        scene_settings.mesh_subdivision_level = params.mesh_subdivision_level;
        plate_loading_mp_refinement(scene_settings, 10 * params.ny, params.ny);
    }
}

void plate_roll_mp_refinement(const ShellSceneSettings& scene_settings, int nx, int ny)
{
    stark::Settings settings = scene_settings.simulator_settings(fmt::format("{}_nx_{}_ny_{}", scene_settings.get_element_types_string(), nx, ny));

    settings.models.enable_model_mp_shell = true;
    settings.execution.end_simulation_time = 12.0;
    // Adaptive timestep and projection not needed for this scene
    settings.newton.project_to_PD = false;
    settings.simulation.use_adaptive_time_step = false;
    settings.simulation.max_time_step_size = 0.01;
    settings.simulation.gravity = {0.0, 0.0, 0.0};

    scene_settings.apply_overrides_to(settings);
    stark::Simulation simulation(settings);

    auto mp = simulation.deformables->strain_micropolar_shells.get();
    scene_settings.apply_overrides_to(mp);

    // Mesh
    const double l = 0.5;
    const double w = 0.125;
    const double dx = l/nx;
    auto [vertices, triangles] = stark::generate_triangle_grid(
            { 0.0, 0.0 },
            {l, w},
            {nx, ny});

    auto material = SurfaceMaterialBuilder();
    material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

    material.thickness = 1.0e-2;
    material.inertia_area_density = 1000.0 * material.thickness;
    material.inertia_damping = 5.0;

    material.strain_young_modulus = 1.0e6;
    material.strain_poissons_ratio = 0.3;
    material.micropolar_mu_c_factor = 1.0;
    material.micropolar_shear_correction = 1.0;
    material.micropolar_length_scale = 1.0e-2;
    material.micropolar_bending_correction = 1.0;
    material.micropolar_angular_inertia = 0.0;

    scene_settings.apply_overrides_to(material);
    auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);

    const Eigen::Vector3d aabb_center = { -0.5*l, 0.0, 0.0 };
    const Eigen::Vector3d aabb_dim = {0.1 * dx, 1.1 * w, 0.1 * dx};

    simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, aabb_center, aabb_dim, bc_params());
    shell.strain_micropolar->prescribe_rotation_inside_aabb(aabb_center, aabb_dim, Eigen::Matrix3d::Identity(), {false, false, true}, 1e4);

    std::vector<Eigen::Vector3d> X_quad = shell.strain_micropolar->interpolate_quadrature_points_mesh();

    // Run
    simulation.run([&]() {
        const double t_end = 2.5;
        const double t = simulation.get_time();

        const double factor = std::min(t, t_end) / t_end;
        const double max_curvature = 4.0 * M_PI;
        const double current_curvature = factor * max_curvature;

        for (int i = 0; i < X_quad.size(); ++i) {
            //const double x = X_q.x();
            //const double x_rel = (x + l) / (2.0 * l);

            // Set only dRy/dx
            const double c = current_curvature;
            stark::MatrixWrapper3d curvature = stark::MatrixWrapper3d::from_array({
                                                                                          0.0, 0.0, 0.0,
                                                                                          c, 0.0, 0.0,
                                                                                          0.0, 0.0, 0.0,
                                                                                  });

            shell.strain_micropolar->set_rest_curvature(i, curvature.to_matrix());
        }

        //fmt::print("current_curvature: {}\n", current_curvature);
    });
}
void plate_roll_mp_refinement_study() {
    ShellSceneSettings scene_settings("plate_roll_mp_refinement_cylinder");
    scene_settings.file_prefix = "rolling";

    scene_settings.output_fps = 60;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.write_quadrature_mesh = false;
    scene_settings.subfolder_per_simulation = true;

    struct Parametrization {
        int ny = 1;
        bool use_quadratic_elements = true;
        int mesh_subdivision_level;
    };

    std::vector<Parametrization> parametrizations = {
            {  1,false, 0},
            {  2,false, 0},
            {  4,false, 0},
            {  8,false, 0},
            { 16,false, 0},
            { 32,false, 0},
            { 64,false, 0},
            // {128,false, 0},
            {  1, true, 3},
            {  2, true, 3},
            {  4, true, 2},
            {  8, true, 1},
            { 16, true, 1},
            { 32, true, 0},
            { 64, true, 0},
            // {128, true, 0},
    };

    for (const auto& params : parametrizations ) {
        scene_settings.mp_element_displacement = params.use_quadratic_elements ? stark::Fem::ElementType::Tri6 : stark::Fem::ElementType::Tri3;
        scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
        scene_settings.mesh_subdivision_level = params.mesh_subdivision_level;
        plate_roll_mp_refinement(scene_settings, 4 * params.ny, params.ny);
    }
}

void lotus_mp_refinement_study()
{
    ShellSceneSettings scene_settings("lotus_mp_refinement_study");
    scene_settings.mp_element_displacement = stark::Fem::ElementType::Tri6;
    scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
    scene_settings.write_quadrature_mesh = false;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.subfolder_per_simulation = true;
    scene_settings.output_fps = 60;

    struct Parametrization {
        std::string file_prefix;
        int mesh_subdivision;
        int n;
        bool with_rot;
    };

    std::vector<Parametrization> parametrizations = {
            { "render", 4, 4, false },
            { "render", 3, 8, false },
            { "render", 2, 16, false },
            { "render", 1, 32, false },
            { "render", 0, 64, false },
            // { "comparison", 1, 4, false },
            // { "comparison", 1, 8, false },
            // { "comparison", 1, 16, false },
            // { "comparison", 1, 32, false },
            // { "comparison", 1, 64, false },
    };

    const auto scene = [](const ShellSceneSettings& scene_settings, int n, double bending, double rot_z) {
        std::string bending_value = fmt::format("{:.2f}", bending);
        std::replace( bending_value.begin(), bending_value.end(), '.', '_');
        stark::Settings settings = scene_settings.simulator_settings(fmt::format("n_{}_bending_{}", n, bending_value));
        //stark::Settings settings = scene_settings.simulator_settings(fmt::format("n_{}_subdiv_{}_bending_{}", n, scene_settings.mesh_subdivision_level, bending_value));

        settings.models.enable_model_mp_shell = true;
        settings.models.enable_model_mp_shell_use_quaternion_gamma = true;

        settings.execution.end_simulation_time = 10.0;
        settings.simulation.gravity = {0.0, 0.0, 0.0};

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);

        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Cloth
        const double l = 0.5;
        auto [vertices, triangles] = stark::generate_triangle_grid(
                { 0.0, 0.0},
                { 2.0*l, 2.0*l },
                { n, n }
        );
        const double dx = (l * 2) / n;

        auto material = SurfaceMaterialBuilder();
        material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

        material.thickness = 1e-2;
        material.inertia_area_density = 3000.0 * material.thickness;
        material.inertia_damping = 0.0;
        material.strain_young_modulus = 5e5;
        material.strain_poissons_ratio = 0.3;
        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_length_scale = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);

        const auto rot = Eigen::AngleAxis(stark::deg2rad(rot_z), Eigen::Vector3d::UnitZ()).toRotationMatrix();
        const Eigen::Vector3d aabb_pos = {0.0, 0.0, 0.0};
        const Eigen::Vector3d aabb_size = 0.1 * dx * Eigen::Vector3d::Ones();

        // BC: Fix central vertex
        simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, aabb_pos, aabb_size, bc_params());
        //shell.strain_micropolar->prescribe_rotation_inside_aabb(aabb_pos, aabb_size, rot, {true, true, true}, 1e4);

        std::vector<Eigen::Vector3d> X_quad = shell.strain_micropolar->interpolate_quadrature_points_mesh();

        shell.point_set->add_rotation(rot_z, Eigen::Vector3d::UnitZ(), {0,0,0}, false);


        // Run
        simulation.run([&]() {
            const double t_end = 6.0;
            const double t = simulation.get_time();

            const double alpha_max = 10.0;
            const double alpha = std::min(t / t_end, 1.0) * alpha_max;

            auto f = [&](double x, double t) {
                double y = std::abs(x*std::min(t + 1.0, 1.8));
                return bending*(-0.2314286 + 10.89881*y - 14.28571*y*y + 4.427083*y*y*y);
            };

            for (int i = 0; i < X_quad.size(); ++i) {
                const auto quad_i = X_quad[i];

                const double x_rel = quad_i.x() / l;
                const double y_rel = quad_i.y() / l;
                const double sx = f(x_rel, t / t_end);
                const double sy = f(y_rel, t / t_end);

                stark::MatrixWrapper3d curvature = stark::MatrixWrapper3d(stark::MatrixWrapper3d::from_array({
                                                                                                                     0.0, alpha*sy, 0.0,
                                                                                                                     -alpha*sx, 0.0, 0.0,
                                                                                                                     0.0, 0.0, 0.0
                                                                                                             }).to_matrix());

                shell.strain_micropolar->set_rest_curvature(i, curvature.to_matrix());
            }
        });
    };

    for (const auto& params : parametrizations) {
        scene_settings.file_prefix = params.file_prefix;
        scene_settings.mesh_subdivision_level = params.mesh_subdivision;

        scene(scene_settings, params.n, 0.1, 0.0);
        scene(scene_settings, params.n, 0.2, params.with_rot ? 30.0 : 0.0);
        scene(scene_settings, params.n, 0.35, params.with_rot ? 72.0 : 0.0);
    }
}

void bunny_comparison_mp()
{
    using StrainModel = SurfaceMaterialBuilder::StrainModel;
    ShellSceneSettings scene_settings("bunny_comparison_mp");

    struct Parametrization {
        std::string file_prefix;
        std::string mesh;
        StrainModel model;
        double thickness;
        double max_dt;
        int sim_subdivision;
        bool full_model;
        bool use_quadratic_elements;
        int render_subdivision;
    };

    std::vector<Parametrization> parametrizations = {
            { "", "bunny_lvl3_verts3308", StrainModel::Wen23, 1.00e-2, 1e-2, 0, false, true, 2 },
            { "lc5e-1h", "bunny_lvl3_verts3308", StrainModel::Micropolar, 1.0e-2, 1e-2, 0, false, true, 2 },
            { "", "bunny_lvl3_verts3308", StrainModel::Wen23, 1.5e-3, 1e-3, 0, false, true, 2 },
            { "lc5e-1h", "bunny_lvl3_verts3308", StrainModel::Micropolar, 1.5e-3, 1e-3, 0, false, true, 2 },
            { "", "bunny_lvl3_verts3308", StrainModel::Wen23, 3.75e-3, 1e-3, 0, false, true, 2 },
            { "lc5e-1h", "bunny_lvl3_verts3308", StrainModel::Micropolar, 3.75e-3, 1e-3, 0, false, true, 2 },
    };

    scene_settings.file_prefix = "";
    //scene_settings.write_quadrature_mesh = true;
    scene_settings.output_fps = 60.0;
    //scene_settings.write_subdivided_mesh = true;
    scene_settings.mp_use_quaternion_gamma = true;
    scene_settings.mp_use_nonlinear_volume = true;
    scene_settings.subfolder_per_simulation = true;
    scene_settings.residual_tol = 1e-7;

    const auto scene = [&](const ShellSceneSettings& scene_settings, const std::string& mesh, double thickness, int subdivision){
        std::string thickness_value = fmt::format("{:.1e}", thickness);
        std::replace( thickness_value.begin(), thickness_value.end(), '.', '_');

        stark::Settings settings = scene_settings.simulator_settings(fmt::format("{}_h{}_{}_n{}", scene_settings.shell_model == StrainModel::Wen23 ? "Wen23" : "Micropolar", thickness_value, mesh, subdivision));

        settings.execution.end_simulation_time = 5.0;
        //settings.simulation.gravity = {0.0, 0.0, 0.0};

        settings.simulation.init_frictional_contact = true;
        settings.newton.project_to_PD = true;

        //settings.simulation.max_time_step_size = 1e-3;
        //settings.simulation.adaptive_time_step.set(0.0, 0.001, 0.020);

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);
        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Global Contact Settings
        const double contact_distance  = 5e-5;
        const double contact_mu = 1.0;

        simulation.interactions->contact->set_global_params(
                stark::EnergyFrictionalContact::GlobalParams()
                        .set_default_contact_thickness(contact_distance)
                        .set_friction_stick_slide_threshold(0.01)
                        .set_min_contact_stiffness(1e6)
        );

        // Mesh
        const double height = 0.2;
        auto meshes = stark::load_obj(fmt::format("{}/{}.obj", MODELS_PATH, mesh));
        auto [vertices, triangles] = meshes.back();

        const auto rot = Eigen::AngleAxisd(stark::deg2rad(90), Eigen::Vector3d::UnitX());
        Eigen::AlignedBox3d aabb;
        for (auto& v : vertices) {
            v = rot * v;
            aabb.extend(v);
        }

        const Eigen::Vector3d center = aabb.center();
        const double scaling = height / aabb.sizes().z();
        aabb = Eigen::AlignedBox3d();

        for (auto& v : vertices) {
            v -= center;
            v = scaling * v;
            aabb.extend(v);
        }

        if (subdivision > 0) {
            stark::Mesh<3> subdivided_mesh;
            subdivided_mesh.vertices = vertices;
            subdivided_mesh.conn = triangles;

            for (int i = 0; i < subdivision; ++i) {
                const auto tri6_mesh = stark::tri3_to_tri6(subdivided_mesh.vertices, subdivided_mesh.conn);
                subdivided_mesh = stark::tri6_to_tri3_subdivide(tri6_mesh.vertices, tri6_mesh.conn);
            }

            vertices = subdivided_mesh.vertices;
            triangles = subdivided_mesh.conn;
        }

        const auto edges = stark::find_edges_from_simplices(triangles, vertices.size());
        double dx = std::numeric_limits<double>::max();
        for (const auto& e: edges) {
            const auto& v1 = vertices[e[0]];
            const auto& v2 = vertices[e[1]];
            dx = std::min(dx, (v2-v1).norm());
        }

        fmt::print("Shortest edge (dx): {:.e}\n", dx);

        auto material = SurfaceMaterialBuilder();
        material.strain_model = scene_settings.shell_model;

        material.enable_contact = true;
        material.contact_distance = contact_distance;

        material.thickness = thickness;
        material.inertia_area_density = 500.0 * material.thickness;
        material.inertia_damping = 1.0;

        material.strain_young_modulus = 1e5;
        material.strain_poissons_ratio = 0.4;

        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_length_scale = thickness * 5e-1;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "", vertices, triangles);
        simulation.interactions->contact->set_friction(shell.contact, shell.contact, contact_mu);

        const Eigen::Vector3d bc_center = {0.0, 0.0, aabb.min().z()};
        const Eigen::Vector3d bc_size_equator = {1.1*aabb.sizes().x(), 1.1*aabb.sizes().y(), 0.01};
        simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, bc_center, bc_size_equator, bc_params(1e7));

        // Floor
        const double size = 5.0 * height;
        const double floor_height = 0.01;
        auto [floor_vertices, floor_triangles, floor] = simulation.presets->rigidbodies->add_box("floor", 1.0, { size, size, 0.01 });
        floor.rigidbody.set_translation({ 0.0, 0.0, aabb.min().z() - 0.5 * floor_height - 2.0*contact_distance });
        simulation.rigidbodies->add_constraint_fix(floor.rigidbody);

        //simulation.interactions->contact->set_friction(shell.contact, contact_handle_floor, contact_mu);

        // Run
        simulation.run([&] {
            if (simulation.get_time() > 0) {
                //std::abort();
            }
        });
    };

    for (const auto& params : parametrizations) {
        scene_settings.shell_model = params.model;
        scene_settings.file_prefix = params.file_prefix;
        scene_settings.mp_use_full_model = params.full_model;
        scene_settings.mesh_subdivision_level = params.render_subdivision;
        scene_settings.max_timestep = params.max_dt;
        scene_settings.mp_element_displacement = params.use_quadratic_elements ? stark::Fem::ElementType::Tri6 : stark::Fem::ElementType::Tri3;
        scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
        scene(scene_settings, params.mesh, params.thickness, params.sim_subdivision);
    }
}

void curved_shell_mp_halfsphere_twisting()
{
    ShellSceneSettings scene_settings("curved_shell_mp_halfsphere_twisting");

    struct Parametrization {
        std::string file_prefix;
        int sim_subdivision;
        bool full_model;
        bool use_quadratic_elements;
        int render_subdivision;
    };

    std::vector<Parametrization> parametrizations = {
            //{ "reduced", 0, false, true, 4 },
            //{ "reduced", 1, false, true, 4 },
            //{ "reduced", 2, false, true, 3 },
            { "reduced", 3, false, true, 2 },
            //{ "reduced", 4, false, true, 1 },
            //{ "full", 0, true, true, 4 },
            //{ "full", 1, true, true, 4 },
            //{ "full", 2, true, true, 3 },
            { "full", 3, true, true, 2 },
            //{ "full", 4, true, true, 1 },
    };

    std::vector<double> thicknesses = {1e-3};

    scene_settings.file_prefix = "";
    //scene_settings.write_quadrature_mesh = true;
    scene_settings.output_fps = 60.0;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.mp_use_quaternion_gamma = true;
    scene_settings.mp_use_nonlinear_volume = false;
    scene_settings.subfolder_per_simulation = true;
    scene_settings.residual_tol = 1e-7;
    //scene_settings.quadrature = stark::TriQuadrature::tri_p4();
    //scene_settings.set_tri3_p2p2();

    const auto scene = [&](const ShellSceneSettings& scene_settings, double thickness, int subdivision){
        std::string thickness_value = fmt::format("{:.1e}", thickness);
        std::replace( thickness_value.begin(), thickness_value.end(), '.', '_');

        stark::Settings settings = scene_settings.simulator_settings(fmt::format("h{}_n{}", thickness_value, subdivision));

        settings.models.enable_model_mp_shell = true;
        settings.execution.end_simulation_time = 6.0;
        //settings.simulation.gravity = {0.0, 0.0, 10};
        settings.simulation.gravity = {0.0, 0.0, 0.0};

        settings.simulation.init_frictional_contact = true;

        //settings.simulation.max_time_step_size = 1e-3;
        //settings.simulation.adaptive_time_step.set(0.0, 0.001, 0.020);

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);
        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Global Contact Settings
        const double contact_distance  = 0.005;
        const double contact_mu = 1.0;

        simulation.interactions->contact->set_global_params(
                stark::EnergyFrictionalContact::GlobalParams()
                        .set_default_contact_thickness(contact_distance)
                        .set_friction_stick_slide_threshold(0.01)
                        .set_min_contact_stiffness(1e6)
        );

        // Mesh
        const double r = 1.0;
        auto meshes = stark::load_obj(fmt::format("{}/uvsphere_half_8x4.obj", MODELS_PATH));
        auto [vertices, triangles] = meshes.back();

        if (subdivision > 0) {
            stark::Mesh<3> subdivided_mesh;
            subdivided_mesh.vertices = vertices;
            subdivided_mesh.conn = triangles;

            for (int i = 0; i < subdivision; ++i) {
                const auto tri6_mesh = stark::tri3_to_tri6(subdivided_mesh.vertices, subdivided_mesh.conn);
                subdivided_mesh = stark::tri6_to_tri3_subdivide(tri6_mesh.vertices, tri6_mesh.conn);
            }

            vertices = subdivided_mesh.vertices;
            triangles = subdivided_mesh.conn;
        }

        const auto edges = stark::find_edges_from_simplices(triangles, vertices.size());
        double dx = std::numeric_limits<double>::max();
        for (const auto& e: edges) {
            const auto& v1 = vertices[e[0]];
            const auto& v2 = vertices[e[1]];
            dx = std::min(dx, (v2-v1).norm());
        }

        fmt::print("Shortest edge (dx): {:.e}\n", dx);

        auto material = SurfaceMaterialBuilder();
        material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

        material.enable_contact = true;
        material.contact_distance = contact_distance;

        material.thickness = thickness;
        material.inertia_area_density = 1000.0 * material.thickness;
        material.inertia_damping = 4.0;

        material.strain_young_modulus = 10000;
        material.strain_poissons_ratio = 0.4;

        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_length_scale = 1e-5;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "", vertices, triangles);
        simulation.interactions->contact->set_friction(shell.contact, shell.contact, contact_mu);

        const Eigen::Vector3d bc_center = {0.0, 0.0, 0.0};
        const Eigen::Vector3d bc_size_equator = {2.1*r, 2.1*r, 0.1*dx};
        const Eigen::Vector3d bc_size_northpole = {0.1*dx, 0.1*dx, 2.1*r};
        //const Eigen::AlignedBox3d aabb = Eigen::AlignedBox3d(bc_center - 0.5 * bc_size, bc_center + 0.5 * bc_size);

        for (int i = 0; i < shell.point_set->size(); ++i) {
            auto v = shell.point_set->get_rest_position(i);
            const Eigen::Vector3d v_new = v.normalized() * r;
            shell.point_set->set_rest_position(i, v_new);
        }

        // Curved initial mesh
        for (int i = 0; i < shell.point_set->size(); ++i) {
            shell.point_set->set_position(i, shell.point_set->get_rest_position(i));
        }

        // BC: Fix deformation at equator
        simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, bc_center, bc_size_equator, bc_params(1e7));

        // Twisting
        const double omega_deg = 90.0;
        const double omega_rad = omega_deg / 180.0 * M_PI;
        //shell.strain_micropolar->prescribe_omega_inside_aabb(bc_center, bc_size_northpole, {0.0, 0.0, omega_rad}, {true, true, true}, 1e3);
        auto bc_rot = shell.strain_micropolar->prescribe_rotation_inside_aabb(bc_center, bc_size_northpole, Eigen::Matrix3d::Identity(), {true, true, true}, 1e3);

        // Run
        simulation.run([&] {
            if (simulation.get_time() > 0) {
                //std::abort();
            }

            // Make sure to apply the transformation for the *end* of the time step
            const double t_stop = 4.0;
            const double t = std::min(simulation.get_time() + simulation.get_time_step_size(), t_stop);
            const double angle = t*omega_rad;

            Eigen::AngleAxisd rot_aa = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
            Eigen::Matrix3d rot = rot_aa.toRotationMatrix();

            shell.strain_micropolar->set_prescribed_rotation(bc_rot, rot);
        });
    };

    for (double h : thicknesses) {
        for (const auto& params : parametrizations) {
            scene_settings.file_prefix = params.file_prefix;
            scene_settings.mp_use_full_model = params.full_model;
            scene_settings.mesh_subdivision_level = params.render_subdivision;
            scene_settings.mp_element_displacement = params.use_quadratic_elements ? stark::Fem::ElementType::Tri6 : stark::Fem::ElementType::Tri3;
            scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
            scene(scene_settings, h, params.sim_subdivision);
        }
    }
}

void moebius_strip_mp(bool generate_meshes = true, bool run_simulation = true)
{
    ShellSceneSettings scene_settings("moebius_strip_mp_make_mesh");
    scene_settings.file_prefix = "";
    scene_settings.mp_element_displacement = stark::Fem::ElementType::Tri6;
    scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;

    scene_settings.mp_use_quaternion_gamma = true;

    const auto scene_mesh = [](const ShellSceneSettings& scene_settings, int ny) {
        stark::Settings settings = scene_settings.simulator_settings(fmt::format("nx_{}_ny_{}", 10 * ny, ny));

        settings.models.enable_model_mp_shell = true;
        settings.execution.end_simulation_time = 6.0;
        settings.newton.project_to_PD = false;
        settings.simulation.use_adaptive_time_step = false;
        settings.simulation.max_time_step_size = 0.01;
        settings.simulation.gravity = {0.0, 0.0, 0.0};

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);
        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Mesh
        const double l = 100 * 1e-3;
        const double w = 10 * 1e-3;
        const double dx = w/ny;
        auto [vertices, triangles] = stark::generate_triangle_grid(
                { 0.0, 0.0 },
                { l, w },
                { 10 * ny, ny });

        auto material = SurfaceMaterialBuilder();
        material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

        material.thickness = 0.5 * 1e-3;
        material.inertia_area_density = 1000.0 * material.thickness;
        material.inertia_damping = 5.0;

        material.strain_young_modulus = 1e6;
        material.strain_poissons_ratio = 0.3;
        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_length_scale = material.thickness * 1e-3;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "", vertices, triangles);

        const Eigen::Vector3d aabb_center_left = { -l * 0.5, 0.0, 0.0 };
        const Eigen::Vector3d aabb_center_right = { l * 0.5, 0.0, 0.0 };
        const Eigen::Vector3d aabb_size = {0.1*dx, 1.1 * w, 0.1*dx};

        // Needed because of high stiffness of material
        const double bc_stiffness = 1e8;
        // Lower stiffness appears to be fine for rotations
        const double bc_stiffness_rot = 1e4;

        std::vector<int> verts_left;
        std::vector<int> verts_right;
        const auto aabb_left = Eigen::AlignedBox3d(aabb_center_left - 0.5*aabb_size, aabb_center_left + 0.5*aabb_size);
        const auto aabb_right = Eigen::AlignedBox3d(aabb_center_right - 0.5*aabb_size, aabb_center_right + 0.5*aabb_size);
        for (int i = 0; i < shell.point_set->size(); ++i) {
            const auto v = shell.point_set->get_position(i);
            if (aabb_left.contains(v)) {
                verts_left.push_back(i);
            } else if (aabb_right.contains(v)) {
                verts_right.push_back(i);
            }
        }

        const auto cmp = [&](const auto& v1, const auto& v2) { return shell.point_set->get_position(v1).y() < shell.point_set->get_position(v2).y(); };
        std::sort(verts_left.begin(), verts_left.end(), cmp);
        std::sort(verts_right.begin(), verts_right.end(), cmp);
        std::reverse(verts_right.begin(), verts_right.end());

        auto bc_left = simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set,
                                                                                     aabb_center_left,
                                                                                     aabb_size,
                                                                                     stark::EnergyPrescribedPositions::Params().set_stiffness(bc_stiffness));

        auto bc_right = simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set,
                                                                                      aabb_center_right,
                                                                                      aabb_size,
                                                                                      stark::EnergyPrescribedPositions::Params().set_stiffness(bc_stiffness));

        const auto dt = [&]() { return simulation.get_time_step_size(); };

        const double twist_end_time = 1.0;
        const double reset_time = 1.5;
        const double curl_start = 2.0;
        const double curl_end = 4.0;
        const double attach_time = 4.5;
        const double final_time = 5.0;
        const double save_time = 5.5;

        const double target_twist = 180.0;
        const Eigen::Vector3d axis = Eigen::Vector3d::UnitX();

        std::vector<Eigen::Vector3d> X_quad = shell.strain_micropolar->interpolate_quadrature_points_mesh();

        simulation.add_time_event(0.0, twist_end_time, [&, dt](double t) {
            const double angle = stark::blend(0.0, target_twist/2.0, 0.0, twist_end_time, t + dt(), stark::BlendType::Linear);
            bc_right.set_transformation(Eigen::Vector3d::Zero(), angle, axis);
            bc_left.set_transformation(Eigen::Vector3d::Zero(), -angle, axis);
        });

        simulation.add_time_event(curl_start, curl_end, [&, dt](double t) {
            const double max_curvature = 2.0 * M_PI / l;
            const double curvature = stark::blend(0.0, max_curvature, curl_start, curl_end, t + dt(), stark::BlendType::Linear);

            for (int i = 0; i < X_quad.size(); ++i) {
                stark::MatrixWrapper3d gamma = stark::MatrixWrapper3d::from_array({
                                                                                          0.0, 0.0, 0.0,
                                                                                          0.0, 0.0, 0.0,
                                                                                          curvature, 0.0, 0.0
                                                                                  });
                shell.strain_micropolar->set_rest_curvature(i, gamma.to_matrix());
            }
        });


        // Run
        bool twist_done = false;
        bool curl_done = false;
        bool final_reset_done = false;
        bool saved = false;
        simulation.run([&]() {
            const double t = simulation.get_time();
            if (!twist_done && t > reset_time) {
                twist_done = true;
                for (int i = 0; i < shell.point_set->size(); ++i) {
                    shell.point_set->set_rest_position(i, shell.point_set->get_position(i));
                }
                shell.strain_micropolar->update_rest_configuration(false);

                auto params = shell.strain_micropolar->get_params();
                params.youngs_modulus = 1e4;
                params.length_scale = 1.0;
                shell.strain_micropolar->set_params(params);

                bc_right.set_params(stark::EnergyPrescribedPositions::Params().set_stiffness(0.0));
                bc_left.set_params(stark::EnergyPrescribedPositions::Params().set_stiffness(0.0));
            }

            if (!curl_done && t > attach_time) {
                curl_done = true;
                auto& set = shell.point_set;
                simulation.interactions->attachments->add(set, set, verts_left, verts_right, stark::EnergyAttachments::Params().set_stiffness(1e6));
            }

            if (!final_reset_done && t > final_time) {
                final_reset_done = true;

                for (int i = 0; i < X_quad.size(); ++i) {
                    shell.strain_micropolar->set_rest_curvature(i, Eigen::Matrix3d::Zero());
                }

                for (int i = 0; i < shell.point_set->size(); ++i) {
                    shell.point_set->set_rest_position(i, shell.point_set->get_position(i));
                }
                shell.strain_micropolar->update_rest_configuration(true);
            }

            if (!saved && t > save_time) {
                saved = true;

                std::vector<std::array<int, 6>> quadratic_triangles = shell.strain_micropolar->get_mesh<stark::Fem::BasisTri6>();
                std::vector<Eigen::Vector3d> quadratic_vertices;
                for (int i = 0; i < shell.point_set->size(); ++i) {
                    quadratic_vertices.push_back(shell.point_set->get_position(i));
                }

                // Replace left side vertices by right side in elements
                for (auto& tri : quadratic_triangles) {
                    for (int i = 0; i < verts_left.size(); ++i) {
                        for (int j = 0; j < tri.size(); ++j) {
                            if (tri[j] == verts_left[i]) {
                                tri[j] = verts_right[i];
                                break;
                            }
                        }
                    }
                }

                std::vector<int> vertex_map(quadratic_vertices.size(), -1);
                std::vector<Eigen::Vector3d> filtered_vertices;
                for (auto& tri : quadratic_triangles) {
                    for (int i = 0; i < tri.size(); ++i) {
                        if (vertex_map[tri[i]] == -1) {
                            vertex_map[tri[i]] = filtered_vertices.size();
                            filtered_vertices.push_back(quadratic_vertices[tri[i]]);
                        }
                        tri[i] = vertex_map[tri[i]];
                    }
                }

                std::vector<std::array<int, 6>> vtk_quadratic_triangles;
                for (const std::array<int, 6>& tri : quadratic_triangles) {
                    vtk_quadratic_triangles.push_back({
                        // VTK node ordering
                        tri[0],
                        tri[2],
                        tri[4],
                        tri[1],
                        tri[3],
                        tri[5]
                    });
                }

                const std::string path = fmt::format("{}/moebius_strip_tri6_n{}.vtk", OUTPUT_PATH, ny);
                vtkio::VTKFile vtk_file;

                if (filtered_vertices.empty()) {
                    vtk_file.write_empty(path);
                } else {
                    vtk_file.set_points_from_twice_indexable(filtered_vertices);
                    vtk_file.set_cells_from_twice_indexable(vtk_quadratic_triangles, vtkio::CellType::Triangle6);
                    vtk_file.write(path);
                }
            }
        });
    };

    if (generate_meshes) {
        scene_mesh(scene_settings, 1);
        scene_mesh(scene_settings, 2);
        scene_mesh(scene_settings, 4);
        scene_mesh(scene_settings, 8);
        scene_mesh(scene_settings, 16);
        scene_mesh(scene_settings, 32);
    }

    auto load_moebius_mesh = [](int n) {
        auto vtkfile = vtkio::VTKFile();
        vtkfile.read(fmt::format("{}/moebius_strip_tri6_n{}.vtk", OUTPUT_PATH, n));

        std::vector<Eigen::Vector3d> vertices;
        vtkfile.get_points_to_twice_indexable(vertices);

        std::vector<std::array<int, 6>> triangles;
        auto cell_type = vtkfile.get_cells_to_twice_indexable(triangles);
        if (cell_type != vtkio::CellType::Triangle6) {
            throw std::runtime_error("Cell type is not Triangle6");
        }

        for (auto& tri : triangles) {
            // Unscramble VTK ordering
            const auto tri_copy = tri;
            tri[0] = tri_copy[0];
            tri[1] = tri_copy[3];
            tri[2] = tri_copy[1];
            tri[3] = tri_copy[4];
            tri[4] = tri_copy[2];
            tri[5] = tri_copy[5];
        }

        return std::make_pair(vertices, triangles);
    };

    scene_settings.scene_name = "moebius_strip_mp_sim";
    scene_settings.write_subdivided_mesh = true;
    scene_settings.mesh_subdivision_level = 1;
    //scene_settings.write_quadrature_mesh = true;
    scene_settings.output_fps = 30;
    scene_settings.subfolder_per_simulation = true;
    scene_settings.mp_use_nonlinear_volume = false;
    scene_settings.mp_use_quaternion_gamma = false;

    enum class MaterialModel {
        Micropolar,
        Wen23,
    };

    const auto model_name = [](MaterialModel model) -> std::string {
        switch (model) {
            case MaterialModel::Micropolar:
                return "mp";
            case MaterialModel::Wen23:
                return "wen23";
            default:
                throw std::runtime_error("not implemented");
        }
    };

    struct Parametrization {
        std::string file_prefix;
        MaterialModel model;
        int mesh_n;
        bool use_full_model;
    };

    std::vector<Parametrization> parametrizations = {
            {"plate", MaterialModel::Micropolar, 2, false },
            {"shell", MaterialModel::Micropolar, 2, true },
            // {"plate", MaterialModel::Micropolar, 4, false },
            // {"shell", MaterialModel::Micropolar, 4, true },
            // {"plate", MaterialModel::Micropolar, 8, false },
            // {"shell", MaterialModel::Micropolar, 8, true },
            // {"plate", MaterialModel::Micropolar, 16, false },
            // {"shell", MaterialModel::Micropolar, 16, true },
            // {"plate", MaterialModel::Micropolar, 32, false },
            // {"shell", MaterialModel::Micropolar, 32, true },
    };

    const auto scene_sim = [&](const ShellSceneSettings& scene_settings, MaterialModel model, int ny) {
        stark::Settings settings = scene_settings.simulator_settings(fmt::format("{}_nx_{}_ny_{}", model_name(model), 10 * ny, ny));

        settings.execution.end_simulation_time = 4.0;
        settings.models.enable_model_mp_shell = true;

        settings.newton.project_to_PD = false;
        settings.simulation.max_time_step_size = 0.001;
        //settings.simulation.gravity = {0.0, 0.0, 0.0};

        settings.simulation.use_adaptive_time_step = false;

        switch (model) {
            case MaterialModel::Micropolar: {
                settings.models.enable_model_mp_shell = true;
                settings.models.enable_model_wen23 = false;
                break;
            }
            case MaterialModel::Wen23: {
                settings.models.enable_model_mp_shell = false;
                settings.models.enable_model_wen23 = true;
                break;
            }
        }

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);
        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Load mesh
        auto [vertices_tri6, triangles_tri6] = load_moebius_mesh(ny);

        auto material = SurfaceMaterialBuilder();
        switch (model) {
            case MaterialModel::Micropolar: {
                material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;
                break;
            }
            case MaterialModel::Wen23: {
                material.strain_model = SurfaceMaterialBuilder::StrainModel::Wen23;
                break;
            }
        }

        material.thickness = 0.5 * 1e-3;
        material.inertia_area_density = 2000.0 * material.thickness;
        material.inertia_damping = 2.0;

        material.strain_young_modulus = 1e6;
        material.strain_poissons_ratio = 0.4;
        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_length_scale = material.thickness * 1e-3;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri6(simulation, "", vertices_tri6, triangles_tri6);

        const Eigen::Vector3d aabb_center = { 0.0, -0.015, 0.0 };
        const Eigen::Vector3d aabb_size = {0.0001, 0.015, 0.0001};
        simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, aabb_center, aabb_size, bc_params());
        if (model == MaterialModel::Micropolar) {
            shell.strain_micropolar->prescribe_rotation_inside_aabb(aabb_center, aabb_size, Eigen::Matrix3d::Identity(), {false, false, true}, 1e4);
        }

        simulation.run([&]() {
            if (simulation.get_time() > 0) {
                //std::abort();
            }
        });
    };

    if (run_simulation) {
        for (const auto& params: parametrizations) {
            scene_settings.file_prefix = params.file_prefix;
            if (params.use_full_model) {
                scene_settings.mp_use_full_model = true;
            }
            scene_sim(scene_settings, params.model, params.mesh_n);
        }
    }
}

void scroll()
{
    bool simulate_paper = true;

    // Set up
    //// Micropolar
    ShellSceneSettings scene_settings("scroll");
    scene_settings.file_prefix = "scroll";

    scene_settings.output_fps = 30;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.write_quadrature_mesh = false;
    scene_settings.subfolder_per_simulation = false;
    scene_settings.mp_use_quaternion_gamma = true;

    //// Stark global
    stark::Settings settings = scene_settings.simulator_settings();
    //settings.output.simulation_name = "scroll";
    settings.models.enable_model_mp_shell = true;
    settings.newton.residual.tolerance = 1e-3;
    settings.newton.enable_flipping_on_non_descent = false;
    settings.newton.project_to_PD = true;
    settings.simulation.max_time_step_size = 0.01;
    settings.simulation.gravity = { 0.0, 0.0, -9.81 };
    settings.simulation.init_frictional_contact = true;
    scene_settings.apply_overrides_to(settings);
    stark::Simulation simulation(settings);

    auto mp = simulation.deformables->strain_micropolar_shells.get();
    scene_settings.apply_overrides_to(mp);

    //// Contact
    double contact_thickness = 0.0005;
    simulation.interactions->contact->set_global_params(
            stark::EnergyFrictionalContact::GlobalParams()
                    .set_default_contact_thickness(contact_thickness)
                    .set_friction_stick_slide_threshold(0.01)
                    .set_min_contact_stiffness(1e9)
    );

    // Floor
    auto [floor_vertices, floor_triangles, floor] = simulation.presets->rigidbodies->add_box("floor", 1.0, { 2.0, 2.0, 0.01 });
    floor.rigidbody.set_translation({ 0.5, 0.5, -0.005 -2.0*contact_thickness -0.005 });
    simulation.rigidbodies->add_constraint_fix(floor.rigidbody);

    // Paper
    //// Geometry
    double lx = 0.5;
    double ly = 0.25;
    int ny = 10;
    int nx = std::ceil(lx/ly * (double)ny);
    auto [vertices, triangles] = stark::generate_triangle_grid({ 0.0, 0.0 }, { lx, ly }, { nx, ny });

    //// Material
    auto material = SurfaceMaterialBuilder();
    material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;
    material.micropolar_element_type_displacement = stark::Fem::ElementType::Tri6;
    material.micropolar_element_type_rotation = stark::Fem::ElementType::Tri3;

    material.thickness = 0.0002;
    material.inertia_area_density = 0.08;
    material.inertia_damping = 1.0;

    material.strain_young_modulus = 1.0e6;
    material.strain_poissons_ratio = 0.3;
    material.micropolar_mu_c_factor = 1.0;
    material.micropolar_shear_correction = 1.0;
    material.micropolar_length_scale = 0.002;  // Stiffness of the turn that prescribes the curvature
    material.micropolar_bending_correction = 1.0;
    material.micropolar_angular_inertia = 0.0;

    material.enable_contact = true;
    material.contact_distance = contact_thickness;

    // Paper
    //// Add
    scene_settings.apply_overrides_to(material);
    auto shell = material.add_surface_from_tri3(simulation, "paper", vertices, triangles);

    //// Waves
    {
        auto& points = shell.point_set;
        const int n_verts = points->size();

        for (int i = 0; i < n_verts; ++i) {
            Eigen::Vector3d v = points->get_rest_position(i);
            v.z() = 0.1;
            points->set_position(i, v);
        }

        struct Wave {
            double frequency;
            double amplitude;
        };

        std::vector<Wave> waves_x{
                { 2, 0.005 },
                { 8.0, 0.001 }
        };

        std::vector<Wave> waves_y{
                { 3.5, 0.001 },
        };

        for (const auto& wave : waves_x) {
            for (int i = 0; i < n_verts; ++i) {
                Eigen::Vector3d v = points->get_rest_position(i);
                const double x_rel = v.x() / (0.5 * lx);
                v.z() += wave.amplitude * std::cos(wave.frequency * x_rel * M_PI);
                points->set_rest_position(i, v);
                points->set_position(i, v);
            }
        }

        for (const auto& wave : waves_y) {
            for (int i = 0; i < n_verts; ++i) {
                Eigen::Vector3d v = points->get_rest_position(i);
                const double y_rel = v.y() / (0.5 * ly);
                v.z() += wave.amplitude * std::sin(wave.frequency * y_rel * M_PI);
                points->set_rest_position(i, v);
                points->set_position(i, v);
            }
        }
    }

    //// BC
    simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, { 0.0, 0.5*ly, 0.0 }, { 0.0001, 0.0001, 1.0 }, bc_params());
    simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, { 0.0, -0.5*ly, 0.0 }, { 0.0001, 0.0001, 1.0 }, bc_params());
    //shell.strain_micropolar->prescribe_rotation_inside_aabb({ 0.0, 0.0, 0.0 }, { 0.0001, 2.0*ly, 2.0*ly }, Eigen::Matrix3d::Identity(), { false, false, true }, 1e4);

    //// Friction
    shell.contact->set_friction(floor.contact, 1.0);

    // Obstacles
    double obs_size = 0.03;
    double obs_init_height = 0.5;
    double obs_floor_height = 0.5 * obs_size + 2.0 * contact_thickness;
    std::vector<Eigen::Vector3d> obstacle_init_positions = { { 0.003688, 0.0666, obs_init_height }, { 0.0, 0.0, obs_init_height }, { -0.003688, -0.0666, obs_init_height } };
    std::vector<Eigen::Vector3d> obstacle_mid_positions = { { 0.003688, 0.0666, obs_floor_height }, { 0.0, 0.0, obs_floor_height }, { -0.003688, -0.0666, obs_floor_height } };
    std::vector<Eigen::Vector3d> obstacle_end_positions = { { -0.15362, 0.0666, obs_floor_height }, { 0.12337, 0.0, obs_floor_height }, { -0.2163, -0.0666, obs_floor_height } };
    std::vector<stark::RBCFixHandler> fixes;
    for (auto& p : obstacle_init_positions) {
        auto mesh = stark::load_obj(fmt::format("{}/rounded_box.obj", MODELS_PATH))[0];
        stark::scale(mesh.vertices, 0.03);

        auto box = simulation.presets->rigidbodies->add("boxes", 0.5, stark::inertia_tensor_box(0.5, { 0.03, 0.03, 0.03 }), mesh.vertices, mesh.conn);
        box.rigidbody.set_translation(p);
        fixes.push_back(simulation.rigidbodies->add_constraint_fix(box.rigidbody));
    }

    // Script
    std::vector<Eigen::Vector3d> X_quad = shell.strain_micropolar->interpolate_quadrature_points_mesh();
    double closing_duration = 5.0;
    if (simulate_paper) {
        // Curvature
        simulation.add_time_event(0.0, closing_duration, [&](double t)
        {
            double magnitude_blend = stark::blend(0.0, 2e3 * M_PI, 0.0, closing_duration, t, stark::BlendType::Linear);
            double closing_in_blend = stark::blend(0.0, 0.05, 0.0, closing_duration, t, stark::BlendType::Linear);
            double magnitude_cap = 100.0;
            double max = 0.0;
            for (int i = 0; i < X_quad.size(); ++i) {
                double x = (closing_in_blend + std::abs(X_quad[i].x())) / lx;
                double p3_bend = std::pow(x, 3) * magnitude_blend;
                double cx = std::min(p3_bend, magnitude_cap);

                double y = (closing_in_blend + std::abs(X_quad[i].y())) / ly;
                double p3_bend_y = std::pow(y, 2) * magnitude_blend;
                double cy = 0.5 * std::min(p3_bend_y, 10.0);

                double cz = 0.0;
                if (X_quad[i].x() < 0.0) {
                    cz *= -1.0;
                }

                stark::MatrixWrapper3d curvature = stark::MatrixWrapper3d::from_array({
                    0.0, -cy, 0.0,
                    -cx, 0.0, 0.0,
                    cz, 0.0, 0.0,
                });
                Eigen::Matrix3d K = curvature.to_matrix();
                Eigen::Matrix3d R = Eigen::AngleAxisd(stark::deg2rad(3.0), Eigen::Vector3d::UnitZ()).toRotationMatrix();

                shell.strain_micropolar->set_rest_curvature(i, R.transpose() * K * R);
                max = std::max(max, std::abs(cx));
            }
            std::cout << " cx = " << max << " ";
        });
    }

    // Boxes descend
    double descend_begin = 6.0;
    double descend_end = 8.0;
    simulation.add_time_event(descend_begin, descend_end, [&](double t)
    {
        for (int box = 0; box < 3; box++) {
            Eigen::Vector3d target;
            for (int dim = 0; dim < 3; dim++) {
                target[dim] = stark::blend(obstacle_init_positions[box][dim], obstacle_mid_positions[box][dim], descend_begin, descend_end, t, stark::BlendType::Linear);
            }
            fixes[box].set_transformation(target, Eigen::Matrix3d::Identity());
        }
    });

    // Paper opens
    double open_begin = 9.0;
    double open_end = 12.0;
    simulation.add_time_event(open_begin, open_end, [&](double t)
    {
        for (int box = 0; box < 3; box++) {
            Eigen::Vector3d target;
            for (int dim = 0; dim < 3; dim++) {
                target[dim] = stark::blend(obstacle_mid_positions[box][dim], obstacle_end_positions[box][dim], open_begin, open_end, t, stark::BlendType::Linear);
            }
            fixes[box].set_transformation(target, Eigen::Matrix3d::Identity());
        }
    });

    // Run
    double duration = 15.0;
    simulation.run(duration, [&]() {
        std::cout << "  fc = " << fixes[1].get_anchor_point().get_violation_in_m_and_force()[1];
    });
}
}

namespace extra_scenes {
void plate_angular_inertia_mp() {
    ShellSceneSettings scene_settings("plate_angular_inertia_mp");
    scene_settings.output_fps = 60;

    const auto scene = [&](const ShellSceneSettings& scene_settings, const bool enable_angular_inertia) {
        stark::Settings settings = scene_settings.simulator_settings();

        settings.models.enable_model_mp_shell = true;
        settings.execution.end_simulation_time = 8.0;

        //settings.simulation.max_time_step_size = 1.0/60.0;   // 60 Hz
        settings.simulation.max_time_step_size = 0.005;

        scene_settings.apply_overrides_to(settings);
        stark::Simulation simulation(settings);

        auto mp = simulation.deformables->strain_micropolar_shells.get();
        scene_settings.apply_overrides_to(mp);

        // Mesh
        const int ny = 4;
        const double l = 20 * 1e-2;
        const double w = 2 * 1e-2;
        const double dx = w/ny;
        auto [vertices, triangles] = stark::generate_triangle_grid(
                { 0.0, 0.0 },
                { l, w },
                { 10 * ny, ny });

        auto material = SurfaceMaterialBuilder();
        material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

        const double density = 1000.0;

        material.thickness = 1.0 * 1e-3;
        material.inertia_area_density = density * material.thickness;
        //material.inertia_damping = 0.0;

        material.strain_young_modulus = 1.0 * 1e6;
        material.strain_poissons_ratio = 0.4;
        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_length_scale = 0.0;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;

        if (enable_angular_inertia) {
            material.micropolar_angular_inertia = 1e-1 * density;
        } else {
            material.micropolar_angular_inertia = 0.0;
        }

        material.micropolar_element_type_displacement = stark::Fem::ElementType::Tri6;
        material.micropolar_element_type_rotation = stark::Fem::ElementType::Tri3;

        scene_settings.apply_overrides_to(material);
        auto shell = material.add_surface_from_tri3(simulation, "shell", vertices, triangles);

        const Eigen::Vector3d aabb_center = { -0.5*l, 0.0, 0.0 };
        const Eigen::Vector3d aabb_dim = {0.1 * dx, 1.1 * w, 1.1 * w};

        simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, aabb_center, aabb_dim, bc_params());
        //shell.prescribe_rotation_inside_aabb(aabb_center, aabb_dim, Eigen::Matrix3d::Identity(), {false, false, true}, 1e4);

        simulation.run([&]() {});
    };

    scene_settings.file_prefix = "ang_inertia_on";
    scene(scene_settings, false);
    scene_settings.file_prefix = "ang_inertia_off";
    scene(scene_settings, true);
}

void plate_loading_mp_refinement_with_quads(const ShellSceneSettings& scene_settings, int nx, int ny)
{
    stark::Settings settings = scene_settings.simulator_settings(fmt::format("{}_nx_{}_ny_{}", scene_settings.get_element_types_string(), nx, ny));

    settings.models.enable_model_mp_shell = true;
    settings.execution.end_simulation_time = 2.0;
    // Adaptive timestep and projection not needed for this scene
    settings.newton.project_to_PD = false;
    settings.simulation.use_adaptive_time_step = false;
    settings.simulation.max_time_step_size = 0.01;
    settings.simulation.gravity = {0.0, 0.0, 0.0};

    scene_settings.apply_overrides_to(settings);
    stark::Simulation simulation(settings);
    auto mp = simulation.deformables->strain_micropolar_shells.get();
    scene_settings.apply_overrides_to(mp);

    // Mesh
    const double l = 0.1;
    const double w = 0.01;
    const double dx = l / nx;

    auto material = SurfaceMaterialBuilder();
    material.strain_model = SurfaceMaterialBuilder::StrainModel::Micropolar;

    material.inertia_damping = 1.0;

    // Parameters of L-shaped beam in SNB16
    {
        material.thickness = 0.6 * 1e-3;						// 0.6 mm
        material.inertia_area_density = 1000.0 * material.thickness;
        material.strain_young_modulus = 71240 * 1e6;			// 71240 N/mm^2 = 71GPa ~ Young's modulus of aluminium
        material.strain_poissons_ratio = 0.31;
        //material.micropolar_mu_c_factor = 0.0;
        material.micropolar_length_scale = 0.6 * 1e-3 * 1e-3;	// 0.6 * 1e-3 mm
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;

        material.micropolar_mu_c_factor = 1.0;
    }

    if (true) {
        // Aluminium plate
        material.thickness = 2.0 * 1e-3;
        material.inertia_area_density = 2700 * material.thickness;
        material.strain_young_modulus = 70 * 1e9;
        material.strain_poissons_ratio = 0.35;
        material.micropolar_mu_c_factor = 1.0;
        material.micropolar_length_scale = material.thickness * 1e-3;
        material.micropolar_shear_correction = 1.0;
        material.micropolar_bending_correction = 1.0;
        material.micropolar_angular_inertia = 0.0;
    }

    scene_settings.apply_overrides_to(material);
    auto shell = material.add_surface_grid(simulation, "shell", { 0.0, 0.0}, { l, w }, { nx, ny });

    std::vector<int> force_vertices;
    for (int i = 0; i < shell.point_set->size(); ++i) {
        if (shell.point_set->get_rest_position(i).x() > 0.5*l - 0.1*dx) {
            force_vertices.push_back(i);
        }
    }
    int n_force_verts = force_vertices.size();
    fmt::print("Number of vertices where force will be applied: {}\n", n_force_verts);

    // BC: Clamp on left side
    const Eigen::Vector3d aabb_center = { -0.5*l, 0.0, 0.0 };
    const Eigen::Vector3d aabb_dim = {0.1*dx, 1.1*w, 0.1*dx};
    simulation.deformables->prescribed_positions->add_inside_aabb(shell.point_set, aabb_center, aabb_dim, bc_params());
    shell.strain_micropolar->prescribe_rotation_inside_aabb(aabb_center, aabb_dim, Eigen::Matrix3d::Identity(), {false, false, true}, 1e4);

    // Run
    fmt::print("\n");
    simulation.run([&]() {
        const double target_force = 50.0;
        //const double target_force = 1.80;
        const double t_end = 0.5 * settings.execution.end_simulation_time;
        //const double t_end = 5.0;
        const double t = simulation.get_time();

        const double factor = std::min(t, t_end) / t_end;
        const double total_force = factor * target_force;
        const double nodal_force = total_force / n_force_verts;
        //fmt::print("nodal_force={:.6e}", nodal_force);

        for (const int i : force_vertices) {
            shell.point_set->set_force(i, {0.0, 0.0, nodal_force});
        }
    });
}
void plate_loading_mp_refinement_study_with_quads() {
    ShellSceneSettings scene_settings("plate_loading_mp_refinement_with_quads");
    scene_settings.file_prefix = "loading";

    scene_settings.output_fps = 60;
    scene_settings.write_subdivided_mesh = true;
    scene_settings.write_quadrature_mesh = false;
    scene_settings.subfolder_per_simulation = true;

    struct Parametrization {
        int ny = 1;
        stark::Fem::ElementType element_type;
        int mesh_subdivision_level;
        bool quadratic_rotations = false;
    };

    std::vector<Parametrization> parametrizations = {
            {  1, stark::Fem::ElementType::Quad4, 3},
            {  2, stark::Fem::ElementType::Quad4, 3},
            {  4, stark::Fem::ElementType::Quad4, 3},
            {  8, stark::Fem::ElementType::Quad4, 3},
            {  16, stark::Fem::ElementType::Quad4, 3},
            {  32, stark::Fem::ElementType::Quad4, 3},
            {  64, stark::Fem::ElementType::Quad4, 3},
            {  1, stark::Fem::ElementType::Tri6, 3},
            {  2, stark::Fem::ElementType::Tri6, 3},
            {  4, stark::Fem::ElementType::Tri6, 2},
            {  8, stark::Fem::ElementType::Tri6, 1},
            { 16, stark::Fem::ElementType::Tri6, 1},
            { 32, stark::Fem::ElementType::Tri6, 0},
            { 64, stark::Fem::ElementType::Tri6, 0},
            {  1, stark::Fem::ElementType::Tri3, 0},
            {  2, stark::Fem::ElementType::Tri3, 0},
            {  4, stark::Fem::ElementType::Tri3, 0},
            {  8, stark::Fem::ElementType::Tri3, 0},
            { 16, stark::Fem::ElementType::Tri3, 0},
            { 32, stark::Fem::ElementType::Tri3, 0},
            { 64, stark::Fem::ElementType::Tri3, 0},
    };

    for (const auto& params : parametrizations) {
        scene_settings.mp_element_displacement = params.element_type;
        if (stark::Fem::element_type_to_shape(params.element_type) == stark::Fem::ElementShape::Tri2d) {
            if (params.quadratic_rotations) {
                scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri6;
            } else {
                scene_settings.mp_element_rotation = stark::Fem::ElementType::Tri3;
            }
        } else {
            if (params.quadratic_rotations) {
                scene_settings.mp_element_rotation = stark::Fem::ElementType::Quad9;
            } else {
                scene_settings.mp_element_rotation = stark::Fem::ElementType::Quad4;
            }
        }
        scene_settings.mesh_subdivision_level = params.mesh_subdivision_level;
        plate_loading_mp_refinement_with_quads(scene_settings, 10 * params.ny, params.ny);
    }
}
}

int main() {
    //paper_scenes::circle_growth_mp();
    //paper_scenes::curvature_modes_mp();
    //paper_scenes::armadillo_twisting_mp();
    paper_scenes::plate_twisting_mp_refinement_study();
    //paper_scenes::plate_loading_mp_refinement_study();
    //paper_scenes::plate_roll_mp_refinement_study();
    //paper_scenes::lotus_mp_refinement_study();
    //paper_scenes::bunny_comparison_mp();
    //paper_scenes::curved_shell_mp_halfsphere_twisting();
    //paper_scenes::moebius_strip_mp(true, false);
    //paper_scenes::moebius_strip_mp(false, true);
    //paper_scenes::scroll();

    //extra_scenes::plate_angular_inertia_mp();
    //extra_scenes::plate_loading_mp_refinement_study_with_quads();

    fmt::println("Exiting.");
    return 0;
}

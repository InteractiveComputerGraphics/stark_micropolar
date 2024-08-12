#pragma once
#include <string>
#include <limits>

#include <Eigen/Dense>
#include <symx>

#include "Console.h"

namespace stark
{
	// User facing enums
	enum class ResidualType { Force, Acceleration };
	enum class LinearSystemSolver { CG, EigenCG, DirectLU };
	struct Residual { ResidualType type; double tolerance; };

	// Wrapper to set SymX compiler
	inline void set_compiler_command(const std::string& str)
	{
		symx::compiler_command = str;
	}
}

namespace stark::core
{
	struct Settings
	{
		struct Output
		{
			std::string simulation_name = "";
			std::string output_directory = "";
			std::string codegen_directory = "";
			std::string time_stamp = "SET_TO_CREATION_TIME_BY_DEFAULT";
			int fps = 30;
			ConsoleVerbosity console_verbosity = ConsoleVerbosity::TimeSteps;
			ConsoleOutputTo console_output_to = ConsoleOutputTo::FileAndConsole;
			bool enable_output = true;
		};
		struct Simulation
		{
			Eigen::Vector3d gravity = { 0.0, 0.0, -9.81 };
			bool init_frictional_contact = true;
			double max_time_step_size = 0.01;
			bool use_adaptive_time_step = true;
			double time_step_size_success_muliplier = 1.05;
			double time_step_size_lower_bound = 1e-6;
		};
		struct NewtonsMethod
		{
			Residual residual = { ResidualType::Acceleration, 1.0 };
			LinearSystemSolver linear_system_solver = LinearSystemSolver::CG;
			bool project_to_PD = true;

			int max_newton_iterations = 100;
			int max_line_search_iterations = 10;
            bool enable_noise_resistant_line_search = true;
            bool enable_flipping_on_non_descent = true;
			double line_search_multiplier = 0.5;
			double cg_max_iterations_multiplier = 1.0;
			double epsilon_residual = 1e-12;  // Does not apply the correction if the residual is below this value. Avoids numerical instability in CG.
		};
		struct Execution
		{
			double allowed_execution_time = std::numeric_limits<double>::max();
			double end_simulation_time = std::numeric_limits<double>::max();
			int end_frame = std::numeric_limits<int>::max();
			int n_threads = -1;
		};
		struct Debug
		{
			bool symx_check_for_NaNs = false;
			bool symx_suppress_compiler_output = true;
			bool symx_force_compilation = false;
			bool symx_force_load = false;
			bool symx_finite_difference_check = false;
			bool line_search_output = false;
		};
        struct Models
        {
            /// Flag to enable/disable compilation of the default triangle strain model
            bool enable_default_tri_strain = false;
            /// Flag to enable/disable compilation of the default tet strain model
            bool enable_default_tet_strain = false;

            /// Flag to enable/disable compilation of the Wen23 Kirchhoff-Love shell model
            bool enable_model_wen23 = false;
            /// Flag to always disable PD projection of the Wen23 model, even when PD projection is enabled
        	bool never_project_tri_wen23 = false;

            /// Flag to enable/disable compilation of the Micropolar shell model
            bool enable_model_mp_shell = false;
            /// Flag to enable/disable use of corotated rotation interpolation instead of angular velocity based interpolation
        	bool enable_model_mp_shell_use_corot_rotation_interpolation = false;
            /// Flag to enable/disable use of quaternion vs rotation matrix based curvature tensor Gamma
            bool enable_model_mp_shell_use_quaternion_gamma = false;
            /// Flag to enable/disable use of full rest curvature terms in the micropolar shell model
            bool enable_model_mp_shell_full_rest_curvature_terms = false;
            /// Flag to enable/disable use of nonlinear volume term in the micropolar shell model
            bool enable_model_mp_shell_use_nonlinear_volume_terms = false;
            // Flag to always disable PD projection of the micropolar shell model, even when PD projection is enabled
            bool never_project_mp_shell = false;
        };

		Output output;
		Simulation simulation;
		NewtonsMethod newton;
		Execution execution;
		Debug debug;
        Models models;

		/* Methods */
		Settings();
		std::string as_string() const;
	};
}

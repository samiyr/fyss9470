#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Defaults {
	static bool pythia_printing 			= true;
	static bool include_decayed 			= true;
	static bool mpi 						= true;
	static bool use_biasing					= false;
	static double bias_power 				= 4.0;
	static double bias_reference 			= 10.0;
	static int random_seed					= -1;
	static bool variable_seed				= false;
	static int pythia_next 					= 10000;
	static bool parallelize					= false;
	static int status_threshold 			= 100;
	static double sigma_eff					= 1.0;
	static double m 						= 0.0;
	static string histogram_file_extension 	= ".csv";
	static string run_data_file_extension 	= ".txt";
}

#include "Analyzer.cc"

namespace Constants {
	/*
		Misc
	 */
	static string cmnd_input = "pp.cmnd";
	/*
		Particle groups
	 */
	static std::vector<int> pions = {111};

	namespace RunParameters {
		static std::vector<Analyzer::Parameters> STAR7 = {
			Analyzer::Parameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_1420_mpi_unity.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_2024_mpi_unity.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_2428_mpi_unity.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_2850_mpi_unity.csv"),

			Analyzer::Parameters(
				1.4, 2.0,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1420_2024_mpi_unity.csv"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1420_2428_mpi_unity.csv"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1420_2850_mpi_unity.csv"),

			Analyzer::Parameters(
				2.0, 2.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_2024_2428_mpi_unity.csv"),
			Analyzer::Parameters(
				2.0, 2.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_2024_2850_mpi_unity.csv"),

			Analyzer::Parameters(
				2.4, 2.8,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_2428_2850_mpi_unity.csv"),
		};
		static std::vector<Analyzer::Parameters> STAR8 = {
			Analyzer::Parameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_1420_mpi__.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_2024_mpi__.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_2428_mpi__.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_2850_mpi__.csv"),

			Analyzer::Parameters(
				1.4, 2.0,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1420_2024_mpi__.csv"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1420_2428_mpi__.csv"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1420_2850_mpi__.csv"),

			Analyzer::Parameters(
				2.0, 2.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_2024_2428_mpi__.csv"),
			Analyzer::Parameters(
				2.0, 2.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_2024_2850_mpi__.csv"),

			Analyzer::Parameters(
				2.4, 2.8,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_2428_2850_mpi__.csv"),
		};

		static std::vector<Analyzer::Parameters> STAR9 = {
			Analyzer::Parameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_1420_mpi__.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_2024_mpi__.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_2428_mpi__.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_2850_mpi__.csv"),

			Analyzer::Parameters(
				1.4, 2.0,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1420_2024_mpi__.csv"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1420_2428_mpi__.csv"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1420_2850_mpi__.csv"),

			Analyzer::Parameters(
				2.0, 2.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_2024_2428_mpi__.csv"),
			Analyzer::Parameters(
				2.0, 2.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_2024_2850_mpi__.csv"),

			Analyzer::Parameters(
				2.4, 2.8,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_2428_2850_mpi__.csv"),
		};
		static std::vector<Analyzer::Parameters> STAR = {
			Analyzer::Parameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_1014_1420"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_1014_2024"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_1014_2428"),
			Analyzer::Parameters(
				1.0, 1.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_1014_2850"),

			Analyzer::Parameters(
				1.4, 2.0,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_1420_2024"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_1420_2428"),
			Analyzer::Parameters(
				1.4, 2.0,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_1420_2850"),

			Analyzer::Parameters(
				2.0, 2.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_2024_2428"),
			Analyzer::Parameters(
				2.0, 2.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_2024_2850"),

			Analyzer::Parameters(
				2.4, 2.8,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_2641_2641_2428_2850"),
		};
		static std::vector<Analyzer::Parameters> birapidity_window_test = {
			Analyzer::Parameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				std::nullopt, std::nullopt,
				"delta_phi_1e7_2641__1014_1420_mpi__.csv"),
			Analyzer::Parameters(
				OptionalRange<double>(1.0, 1.4),
				OptionalRange<double>(1.4, 2.0),
				OptionalRange<double>(2.6, 4.1),
				OptionalRange<double>(2.6, 4.1, true),
				"delta_phi_1e7_2641_2641c_1014_1420_mpi__.csv"),
			Analyzer::Parameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_1e7_2641_2641_1014_1420_mpi__.csv"),
		};
	}
}

#endif // CONSTANTS_H
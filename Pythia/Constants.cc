#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Defaults {
	static bool pythia_printing 	= true;
	static bool include_decayed 	= true;
	static bool mpi 				= true;
	static bool use_biasing			= false;
	static double bias_power 		= 4.0;
	static double bias_reference 	= 10.0;
	static int random_seed			= -1;
	static bool variable_seed		= false;
	static int pythia_next 			= 10990;
	static bool parallelize			= false;
	static int status_threshold 	= 100;
}

#include "AnalysisParameters.cc"

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
		static std::vector<AnalysisParameters> STAR7 = {
			AnalysisParameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_1420_mpi_unity.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_2024_mpi_unity.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_2428_mpi_unity.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1014_2850_mpi_unity.csv"),

			AnalysisParameters(
				1.4, 2.0,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1420_2024_mpi_unity.csv"),
			AnalysisParameters(
				1.4, 2.0,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1420_2428_mpi_unity.csv"),
			AnalysisParameters(
				1.4, 2.0,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_1420_2850_mpi_unity.csv"),

			AnalysisParameters(
				2.0, 2.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_2024_2428_mpi_unity.csv"),
			AnalysisParameters(
				2.0, 2.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_2024_2850_mpi_unity.csv"),

			AnalysisParameters(
				2.4, 2.8,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR7/delta_phi_1e7_2641_2641_2428_2850_mpi_unity.csv"),
		};
		static std::vector<AnalysisParameters> STAR8 = {
			AnalysisParameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_1420_mpi__.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_2024_mpi__.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_2428_mpi__.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1014_2850_mpi__.csv"),

			AnalysisParameters(
				1.4, 2.0,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1420_2024_mpi__.csv"),
			AnalysisParameters(
				1.4, 2.0,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1420_2428_mpi__.csv"),
			AnalysisParameters(
				1.4, 2.0,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_1420_2850_mpi__.csv"),

			AnalysisParameters(
				2.0, 2.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_2024_2428_mpi__.csv"),
			AnalysisParameters(
				2.0, 2.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_2024_2850_mpi__.csv"),

			AnalysisParameters(
				2.4, 2.8,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR8/delta_phi_1e8_2641_2641_2428_2850_mpi__.csv"),
		};

		static std::vector<AnalysisParameters> STAR9 = {
			AnalysisParameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_1420_mpi__.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_2024_mpi__.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_2428_mpi__.csv"),
			AnalysisParameters(
				1.0, 1.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1014_2850_mpi__.csv"),

			AnalysisParameters(
				1.4, 2.0,
				2.0, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1420_2024_mpi__.csv"),
			AnalysisParameters(
				1.4, 2.0,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1420_2428_mpi__.csv"),
			AnalysisParameters(
				1.4, 2.0,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_1420_2850_mpi__.csv"),

			AnalysisParameters(
				2.0, 2.4,
				2.4, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_2024_2428_mpi__.csv"),
			AnalysisParameters(
				2.0, 2.4,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_2024_2850_mpi__.csv"),

			AnalysisParameters(
				2.4, 2.8,
				2.8, 5.0,
				2.6, 4.1,
				2.6, 4.1,
				"STAR9/delta_phi_1e9_2641_2641_2428_2850_mpi__.csv"),
		};
		static std::vector<AnalysisParameters> birapidity_window_test = {
			AnalysisParameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				std::nullopt, std::nullopt,
				"delta_phi_1e7_2641__1014_1420_mpi__.csv"),
			AnalysisParameters(
				OptionalRange<double>(1.0, 1.4),
				OptionalRange<double>(1.4, 2.0),
				OptionalRange<double>(2.6, 4.1),
				OptionalRange<double>(2.6, 4.1, true),
				"delta_phi_1e7_2641_2641c_1014_1420_mpi__.csv"),
			AnalysisParameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"delta_phi_1e7_2641_2641_1014_1420_mpi__.csv"),
		};
	}
}

#endif // CONSTANTS_H
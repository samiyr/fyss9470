#ifndef CONSTANTS_H
#define CONSTANTS_H

#define EVENT_COUNT_TYPE unsigned long long

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
		static std::vector<Analyzer::Parameters> moving_pT_window_25 = {
			Analyzer::Parameters(
				1.0, 1.4,
				1.4, 2.0,
				2.6, 4.1,
				2.6, 4.1,
				"data1",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.1, 1.5,
				1.5, 2.1,
				2.6, 4.1,
				2.6, 4.1,
				"data2",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.2, 1.6,
				1.6, 2.2,
				2.6, 4.1,
				2.6, 4.1,
				"data3",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.3, 1.7,
				1.7, 2.3,
				2.6, 4.1,
				2.6, 4.1,
				"data4",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.4, 1.8,
				1.8, 2.4,
				2.6, 4.1,
				2.6, 4.1,
				"data5",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.5, 1.9,
				1.9, 2.5,
				2.6, 4.1,
				2.6, 4.1,
				"data6",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.6, 2.0,
				2.0, 2.6,
				2.6, 4.1,
				2.6, 4.1,
				"data7",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.7, 2.1,
				2.1, 2.7,
				2.6, 4.1,
				2.6, 4.1,
				"data8",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.8, 2.2,
				2.2, 2.8,
				2.6, 4.1,
				2.6, 4.1,
				"data9",
				1.0,
				25.0),
			Analyzer::Parameters(
				1.9, 2.3,
				2.3, 2.9,
				2.6, 4.1,
				2.6, 4.1,
				"data10",
				1.0,
				25.0),
			Analyzer::Parameters(
				2.0, 2.4,
				2.4, 3.0,
				2.6, 4.1,
				2.6, 4.1,
				"data11",
				1.0,
				25.0),
		};
	}
}

#endif // CONSTANTS_H
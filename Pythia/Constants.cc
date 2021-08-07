#ifndef CONSTANTS_H
#define CONSTANTS_H

/// Defines the type used for indicating the number of events to be generated.
/// Needed because int might overflow for large enough event counts.
#define EVENT_COUNT_TYPE unsigned long long

/// Provides a list of default parameter values.
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
}

/// Provides a list of constant values.
namespace Constants {
	/// The file name for the Pythia cmnd input file.
	static string cmnd_input = "pp.cmnd";
	/// The ID for neutral pions.
	static std::vector<int> pions = {111};
}

#endif // CONSTANTS_H
#ifndef GENERATOR_PARAMETERS_H
#define GENERATOR_PARAMETERS_H

#include "Helpers.cc"

struct GeneratorParameters {
	double cm_energy;
	int event_count;
	bool pythia_printing;

	bool include_decayed;

	bool mpi;

	OptionalRange<double> pT_range;
	OptionalRange<double> y_range;

	bool use_biasing;
	double bias_power;
	double bias_reference;

	int random_seed;

	std::vector<int> particle_ids;
};

#endif // GENERATOR_PARAMETER_H
#ifndef GENERATOR_PARAMETERS_H
#define GENERATOR_PARAMETERS_H

#include "Helpers.cc"
#include "Beam.cc"

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

	Beam beam_A;
	Beam beam_B;
};

#endif // GENERATOR_PARAMETER_H
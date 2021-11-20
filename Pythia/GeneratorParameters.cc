#ifndef GENERATOR_PARAMETERS_H
#define GENERATOR_PARAMETERS_H

#include "Helpers.cc"
#include "Beam.cc"

/// An encapsulation of generator parameters.
struct GeneratorParameters {
	/// The center-of-mass energy, in GeV.
	double cm_energy;
	/// The number of events to be generated.
	EVENT_COUNT_TYPE event_count;
	/// A flag determining whether Pythia should print information prior to event generation.
	bool pythia_printing;
	/// Whether to include decayed particles.
	bool include_decayed;
	/// Whether MPI is enabled in Pythia.
	bool mpi;
	/// The acceptance range of transverse momentum.
	OptionalRange<double> pT_range;
	/// The acceptance range of rapidity.
	OptionalRange<double> y_range;
	/// Whether to use phase space biasing.
	bool use_biasing;
	/// The exponent used in phase space biasing. If use_biasing is set to false, this value has no effect.
	double bias_power;
	/// The reference bias used in phase space biasing. If use_biasing is set to false, this value has no effect.
	double bias_reference;
	/// The random seed for the PRNG used in Pythia.
	int random_seed;
	/// The list of accepted particle identifiers.
	std::vector<int> particle_ids;
	/// Beam A parameters.
	Beam beam_A;
	/// Beam B parameters.
	Beam beam_B;
	/// Type of process used in Pythia.
	Process process;
	/// Enable ncoll mode.
	bool use_ncoll;
	/// Number of retries per event in ncoll mode.
	EVENT_COUNT_TYPE ncoll_retries;
	/// ncoll index offset and multiplier, i.e. the requested ncoll number is given by 
	/// i_event + ncoll_multiplier * ncoll_offset. Therefore ncoll_multiplier = 0 disabled offsetting.
	EVENT_COUNT_TYPE ncoll_offset;
	EVENT_COUNT_TYPE ncoll_multiplier;
};

#endif // GENERATOR_PARAMETER_H
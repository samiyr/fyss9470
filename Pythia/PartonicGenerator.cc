#ifndef PARTONIC_GENERATOR_H
#define PARTONIC_GENERATOR_H

#include "Pythia8/Pythia.h"
#include "ParticleGenerator.cc"

using namespace Pythia8;

class PartonicGenerator {
public:
	double cm_energy;
	int event_count;
	bool pythia_printing = true;

	bool include_decayed = true;

	bool mpi = true;

	OptionalRange<double> pT_range;
	OptionalRange<double> y_range;
	OptionalRange<double> pT_hat_range;

	bool use_biasing = false;
	double bias_power = 4.0;
	double bias_reference = 10.0;

	bool variable_seed = false;
	int random_seed = -1;

	std::vector<int> particle_ids = {111};

	bool parallelize = false;
	std::vector<OptionalRange<double>> pT_hat_bins;

	PartonicGenerator(double energy, int count, std::vector<OptionalRange<double>> bins) {
		cm_energy = energy;
		event_count = count;
		pT_hat_bins = bins;
	}

	template <typename F, typename G>
	void generate(F &&lambda, G &&completion) {
		#pragma omp parallel for if(parallelize)
		for (std::vector<OptionalRange<double>>::size_type i = 0; i < pT_hat_bins.size(); i++) {
			const auto range = pT_hat_bins[i];
			ParticleGenerator generator(cm_energy, event_count);

			generator.particle_ids = particle_ids;
			generator.include_decayed = include_decayed;
			generator.y_range = y_range;
			generator.pT_range = pT_range;
			generator.pT_hat_range = range;
			generator.use_biasing = use_biasing;
			generator.bias_power = bias_power;
			generator.pythia_printing = pythia_printing;
			generator.mpi = mpi;

			if (variable_seed) {
				generator.random_seed = i + random_seed;
			} else {
				generator.random_seed = random_seed;
			}

			generator.initialize();

			const std::vector<std::vector<ParticleContainer>> particles = generator.generate();
			lambda(particles, &generator);
		}
		completion();
	}
};

#endif // PARTONIC_GENERATOR_H
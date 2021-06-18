#ifndef PARTONIC_GENERATOR_H
#define PARTONIC_GENERATOR_H

#include "Pythia8/Pythia.h"
#include "ParticleGenerator.cc"

using namespace Pythia8;

class PartonicGenerator {
public:
	double cm_energy;
	int event_count;
	bool pythia_printing = Defaults::pythia_printing;

	bool include_decayed = Defaults::include_decayed;

	bool mpi = Defaults::mpi;

	OptionalRange<double> pT_range;
	OptionalRange<double> y_range;
	OptionalRange<double> pT_hat_range;

	bool use_biasing = Defaults::use_biasing;
	double bias_power = Defaults::bias_power;
	double bias_reference = Defaults::bias_reference;

	bool variable_seed = Defaults::variable_seed;
	int random_seed = Defaults::random_seed;

	std::vector<int> particle_ids = Constants::pions;

	bool parallelize = Defaults::parallelize;
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

	template <typename F, typename G>
	void generate_at_loop(F &&lambda, G &&completion) {
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

			generator.generate([&lambda, &generator](std::vector<ParticleContainer> particles) {
				lambda(particles, &generator);
			});
		}
		completion();
	}
};

#endif // PARTONIC_GENERATOR_H
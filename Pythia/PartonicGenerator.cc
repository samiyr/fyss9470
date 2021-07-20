#ifndef PARTONIC_GENERATOR_H
#define PARTONIC_GENERATOR_H

#include "Pythia8/Pythia.h"
#include "ParticleGenerator.cc"

using namespace Pythia8;

class PartonicGenerator {
public:
	// double cm_energy;
	// int event_count;
	// bool pythia_printing = Defaults::pythia_printing;

	// bool include_decayed = Defaults::include_decayed;

	// bool mpi = Defaults::mpi;

	// OptionalRange<double> pT_range;
	// OptionalRange<double> y_range;
	// OptionalRange<double> pT_hat_range;

	// bool use_biasing = Defaults::use_biasing;
	// double bias_power = Defaults::bias_power;
	// double bias_reference = Defaults::bias_reference;

	// bool variable_seed = Defaults::variable_seed;
	// int random_seed = Defaults::random_seed;

	// std::vector<int> particle_ids = Constants::pions;
	GeneratorParameters params;

	bool parallelize;
	bool variable_seed;
	std::vector<OptionalRange<double>> pT_hat_bins;

	PartonicGenerator(GeneratorParameters p, std::vector<OptionalRange<double>> bins) {
		params = p;
		pT_hat_bins = bins;
	}

	template <typename F, typename G>
	void generate(F &&lambda, G &&completion) {
		#pragma omp parallel for if(parallelize)
		for (std::vector<OptionalRange<double>>::size_type i = 0; i < pT_hat_bins.size(); i++) {
			const auto range = pT_hat_bins[i];

			ParticleGenerator generator(params);

			generator.pT_hat_range = range;

			if (variable_seed) {
				generator.params.random_seed = i + params.random_seed;
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

			ParticleGenerator generator(params);

			generator.pT_hat_range = range;

			if (variable_seed) {
				generator.params.random_seed = i + params.random_seed;
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
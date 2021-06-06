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

	OptionalRange<double> pT_range = OptionalRange<double>();
	OptionalRange<double> y_range = OptionalRange<double>();
	OptionalRange<double> pT_hat_range = OptionalRange<double>();

	bool use_biasing = false;
	double bias_power = 4.0;
	double bias_reference = 10.0;

	std::vector<int> particle_ids = {111};

	bool parallelize = false;
	std::vector<std::optional<double>> pT_hat_bins;

	PartonicGenerator(double energy, int count, std::vector<std::optional<double>> bins) {
		cm_energy = energy;
		event_count = count;
		pT_hat_bins = bins;
	}

	template <typename F, typename G>
	void start(F &&lambda, G &&completion) {
		#pragma omp parallel for if(parallelize)
		for (std::vector<double>::size_type i = 0; i < pT_hat_bins.size() - 1; i++) {
			const double pT_hat_min = *pT_hat_bins[i];
			const double pT_hat_max = *pT_hat_bins[i + 1];

			ParticleGenerator generator(cm_energy, event_count);

			generator.particle_ids = particle_ids;
			generator.include_decayed = include_decayed;
			generator.y_range = y_range;
			generator.pT_hat_range = OptionalRange<double>(pT_hat_min, pT_hat_max);
			generator.use_biasing = use_biasing;
			generator.bias_power = bias_power;

			generator.initialize();

			const std::vector<ParticleContainer> particles = generator.generate();
			lambda(particles, &generator);
		}
		completion();
	}
};

#endif // PARTONIC_GENERATOR_H
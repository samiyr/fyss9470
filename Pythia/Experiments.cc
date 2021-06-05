#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Pythia8/Pythia.h"
#include <string>
#include "ParticleGenerator.cc"
#include "Histogram.cc"
#include "Helpers.cc"
#include <math.h>
#include <execution>

using namespace Pythia8;

class CrossSectionExperiment {
public:
	double energy;
	int count;
	std::vector<double> bins;
	std::vector<std::optional<double>> pT_hat_bins;
	OptionalRange<double> y_range;
	bool include_decayed;
	bool use_biasing;
	double bias_power = 4.0;
	bool parallelize;
	std::string filename = "pT_histogram.csv";

	void run() {
		std::vector<ValueHistogram<double>> containers;
		std::vector<double> weights;

		#pragma omp parallel for if(parallelize)
		for (std::vector<double>::size_type i = 0; i < pT_hat_bins.size() - 1; i++) {
			const double pT_hat_min = *pT_hat_bins[i];
			const double pT_hat_max = *pT_hat_bins[i + 1];

			ParticleGenerator generator(energy, count);

			generator.include_decayed = include_decayed;
			generator.y_range = y_range;
			generator.pT_hat_range = OptionalRange<double>(pT_hat_min, pT_hat_max);
			generator.use_biasing = use_biasing;
			generator.bias_power = bias_power;

			generator.initialize();

			const std::vector<ParticleContainer> pions = generator.generate();
			const double sigma = generator.sigma();

			const std::vector<double> pTs = find_pTs(pions);
			const std::vector<double> event_weights = find_event_weights(pions);

			Histogram<double> partial = Histogram<double>(bins);
			partial.fill(pTs, event_weights);

			const auto partial_container = partial.normalize(generator.total_weight(), sigma, y_range, use_biasing);
			containers.push_back(partial_container);

			weights.push_back(1.0);
		}

		const auto combined = combine(containers, weights);

		cout << "Normalized pT histogram" << "\n";
		combined.print();
		cout << "\n";
		combined.export_histogram(filename);
	}
};

class AzimuthCorrelationExperiment {
public:
	double energy;
	int count;
	OptionalRange<double> y_range;
	bool include_decayed;
	bool parallelize;

	void run() {
		cout << "Starting experiment with E = " << energy << ", N = " << count << "\n";
		cout << "Generating pions" << "\n";

		ParticleGenerator generator(energy, count);

		generator.include_decayed = include_decayed;
		generator.y_range = y_range;

		generator.initialize();

		const std::vector<ParticleContainer> pions = generator.generate();
		cout << "Generated " << pions.size() << " pions" << "\n";

		Histogram<double> hist(10, 0, M_PI);
		#pragma omp parallel if(parallelize)
		{
			std::vector<double> _deltas;
			// TODO: reserve space for _deltas
			#pragma omp for collapse(2)
			for (std::vector<ParticleContainer>::size_type i = 0; i < pions.size(); i++) {
				const Particle p1 = pions[i].particle;
				const double phi1 = p1.phi();
				for (std::vector<Particle>::size_type j = i + 1; j < pions.size(); j++) {
					const Particle p2 = pions[j].particle;
					const double phi2 = p2.phi();
					const double delta_phi = abs(phi1 - phi2);
					_deltas.push_back(min(delta_phi, 2 * M_PI - delta_phi));
				}
			}
			#pragma omp critical
			hist.fill(_deltas);
		}

		ValueHistogram<double> exported = hist.export_to_values();
		exported.export_histogram("delta_phi.csv");
	}
};


int main() {
	CrossSectionExperiment cs;

	cs.energy = 200;
	cs.count = 5'000'000;
	cs.bins = {
		1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0
	};
	cs.pT_hat_bins = {
		2.0, 5.0, 10.0, 40.0, std::nullopt
	};
	cs.y_range = OptionalRange<double>(-0.35, 0.35);
	cs.include_decayed = true;
	cs.use_biasing = true;
	cs.bias_power = 4.0;
	cs.parallelize = true;
	cs.filename = "pT_histogram.csv";

	cs.run();
	
	AzimuthCorrelationExperiment ac;

	ac.energy = 200;
	ac.count = 10000;
	ac.y_range = OptionalRange<double>(-0.35, 0.35);
	ac.include_decayed = true;
	ac.parallelize = true;

	// ac.run();

	return 0;
}

#endif // EXPERIMENTS_H
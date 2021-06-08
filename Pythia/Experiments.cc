#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Pythia8/Pythia.h"
#include "PartonicGenerator.cc"
#include "Histogram.cc"

using namespace Pythia8;

/**
 * Represents a generic Pythia experiment.
 * Must be subclassed, calling `run()` directly
 * will kill the execution.
 */
class Experiment {
public:
	/// Center-of-mass energy in GeV.
	double energy;
	/// Number of events per partonic bin.
	int count;
	/// Histogram bins.
	std::vector<double> bins;
	/// Partonic bins for generating high-pT particles.
	std::vector<std::optional<double>> pT_hat_bins;
	/// Allowed rapidity range.
	OptionalRange<double> y_range;
	/// Whether to include decayed particles.
	bool include_decayed;
	/// Whether to use partonic pT biasing.
	bool use_biasing;
	/// The power to which partonic pT will be raised when biasing.
	/// Ignored if `use_biasing` is set to `false`.
	double bias_power;
	/// Use OpenMP multithreading if enabled.
	bool parallelize;
	/// Filename of the exported file. 
	/// If set to `nullopt`, no data is exported.
	std::optional<std::string> filename;
	/// Controls the 'Print:quiet' Pythia flag.
	bool pythia_printing;
	/// Runs the experiment. Subclasses must implement this method.
	void run() {
		abort();
	}
};

class CrossSectionExperiment : public Experiment {
public:
	void run() {
		std::vector<ValueHistogram<double>> containers;

		PartonicGenerator generator(energy, count, pT_hat_bins);

		generator.include_decayed = include_decayed;
		generator.y_range = y_range;
		generator.use_biasing = use_biasing;
		generator.bias_power = bias_power;
		generator.parallelize = parallelize;
		generator.pythia_printing = pythia_printing;

		generator.generate([&containers, this](std::vector<ParticleContainer> pions, ParticleGenerator *particle_generator) {
			const std::vector<double> pTs = find_pTs(pions);
			const std::vector<double> event_weights = find_event_weights(pions);

			Histogram<double> partial(bins);
			partial.fill(pTs, event_weights);

			const double sigma = particle_generator->sigma();
			const double total_weight = particle_generator->total_weight();
			const auto partial_container = partial.normalize(total_weight, sigma, y_range, use_biasing);
			containers.push_back(partial_container);
		}, [&containers, this]{
			const auto combined = combine(containers);

			cout << "Normalized pT histogram" << "\n";
			combined.print();
			cout << "\n";
			if (filename) {
				combined.export_histogram(*filename);
			}
		});
	}
};

class AzimuthCorrelationExperiment : public Experiment {
public:
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
		if (filename) {
			exported.export_histogram(*filename);
		}
	}
};


int main() {
	CrossSectionExperiment cs;

	cs.energy = 200;
	cs.count = 10'000;
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
	cs.pythia_printing = false;

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
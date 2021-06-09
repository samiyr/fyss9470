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
	std::vector<OptionalRange<double>> pT_hat_bins;
	/// Allowed rapidity range.
	OptionalRange<double> y_range;
	/// Allowed transverse momentum range.
	OptionalRange<double> pT_range;
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
	/// Initializes an instance of `PartonicGenerator`
	/// with the appropriate parameters set.
	PartonicGenerator create_generator() {
		PartonicGenerator generator(energy, count, pT_hat_bins);

		generator.include_decayed = include_decayed;
		generator.y_range = y_range;
		generator.pT_range = pT_range;
		generator.use_biasing = use_biasing;
		generator.bias_power = bias_power;
		generator.parallelize = parallelize;
		generator.pythia_printing = pythia_printing;

		return generator;
	}
};

class CrossSectionExperiment : public Experiment {
public:
	void run() {
		std::vector<ValueHistogram<double>> containers;
		PartonicGenerator generator = create_generator();

		generator.generate([&containers, this](std::vector<std::vector<ParticleContainer>> pions_by_event, ParticleGenerator *particle_generator) {
			const auto pions = flatten(pions_by_event);
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
	Range<double> pT_1;
	Range<double> pT_2;
	void run() {
		std::vector<std::vector<ParticleContainer>> pions;
		PartonicGenerator generator = create_generator();

		generator.generate([&pions](std::vector<std::vector<ParticleContainer>> particles, [[maybe_unused]] ParticleGenerator *particle_generator) {
			#pragma omp critical
			pions.insert(pions.end(), particles.begin(), particles.end());
		}, [&pions, this]() {
			std::vector<std::vector<ValueHistogram<unsigned long long int>>> histograms;
			#pragma omp parallel if(parallelize)
			{
				std::vector<ValueHistogram<unsigned long long int>> _histograms;
				#pragma omp for nowait
				for (auto &list : pions) {
					ValueHistogram<unsigned long long int> _hist(bins);
					const auto N = list.size();
					for (std::vector<ParticleContainer>::size_type i = 0; i < N; i++) {
						const Particle particle1 = list[i].particle;
						const bool check11 = pT_1.in_range(particle1.pT());
						const bool check12 = pT_2.in_range(particle1.pT());
						if (!(check11 || check12)) {
							continue;
						}
						const double phi1 = particle1.phi();
						for (std::vector<ParticleContainer>::size_type j = i + 1; j < N; j++) {
							const Particle particle2 = list[j].particle;
							const bool check21 = pT_1.in_range(particle2.pT());
							const bool check22 = pT_2.in_range(particle2.pT());
							if (!((check21 && !check11) || (check22 && !check12))) {
								continue;
							}
							const double phi2 = particle2.phi();

							const double delta_phi = abs(phi1 - phi2);
							const double value = min(delta_phi, 2 * M_PI - delta_phi);

							_hist.fill(value);	
						}
					}
					#pragma omp critical
					_histograms.push_back(_hist);
				}
				#pragma omp critical
				histograms.push_back(_histograms);
			}
			const auto combined = combine(flatten(histograms));
			const auto normalized = combined.normalize_to_unity();
			cout << "Azimuth histogram" << "\n";
			normalized.print_with_bars();
			cout << "\n";
			if (filename) {
				normalized.export_histogram(*filename);
			}			
		});
	}
};


int main() {
	CrossSectionExperiment cs;

	cs.energy = 200;
	cs.count = 10000;
	cs.bins = {
		1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0
	};
	cs.pT_hat_bins = {
		OptionalRange<double>(2.0, 5.0),
		OptionalRange<double>(5.0, 10.0),
		OptionalRange<double>(10.0, 40.0),
		OptionalRange<double>(40.0, std::nullopt)
	};
	cs.y_range = OptionalRange<double>(-0.35, 0.35);
	cs.pT_range = OptionalRange<double>();
	cs.include_decayed = true;
	cs.use_biasing = true;
	cs.bias_power = 4.0;
	cs.parallelize = true;
	cs.filename = "pT_histogram.csv";
	cs.pythia_printing = false;

	//cs.run();

	
	AzimuthCorrelationExperiment ac;

	ac.energy = 200;
	ac.count = 100'000 / 16;
	ac.bins = fixed_range(0.0, M_PI, 20);
	ac.pT_hat_bins = std::vector<OptionalRange<double>>(16, OptionalRange<double>(1.0, std::nullopt));
	ac.y_range = OptionalRange<double>();
	ac.pT_range = OptionalRange<double>(1.0, 2.0);
	ac.pT_1 = Range<double>(1.0, 1.4);
	ac.pT_2 = Range<double>(1.4, 2.0);
	ac.include_decayed = true;
	ac.use_biasing = true;
	ac.bias_power = 4.0;
	ac.parallelize = true;
	ac.filename = "delta_phi.csv";
	ac.pythia_printing = false;

	ac.run();

	return 0;
}

#endif // EXPERIMENTS_H
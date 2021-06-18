#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Pythia8/Pythia.h"
#include "PartonicGenerator.cc"
#include "Histogram.cc"
#include "CorrelationAnalyzer.cc"
#include "Constants.cc"

using namespace Pythia8;

/**
 * Represents a generic Pythia experiment.
 * Must be subclassed, calling `run()` directly
 * will kill the execution.
 */
class Experiment {
public:
	enum class Normalization {
		None, Unity, Count
	};
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
	/// Controls the 'Print:quiet' Pythia flag.
	bool pythia_printing;
	/// The base random seed used in Pythia. See also `variable_seed`.
	int random_seed = -1;
	/// If set, the random seed forwarded to Pythia will be
	/// `random_seed` + i, where i is the iterator index of
	/// the partonic pT bins `pT_hat_bins`.
	bool variable_seed = false;
	/// Whether to use multiparton interactions.
	bool mpi = true;
	/// Normalization type to use.
	Normalization normalization;
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
		generator.variable_seed = variable_seed;
		generator.random_seed = random_seed;

		return generator;
	}

	template<typename T>
	ValueHistogram<T> normalize(ValueHistogram<T> hist) {
		switch (normalization) {
			case Normalization::None:
				return hist;
				break;
			case Normalization::Unity:
				return hist.normalize_to_unity();
				break;
			case Normalization::Count:
				return hist.normalize_by(hist.total());
				break;
		}
	}
};

class CrossSectionExperiment : public Experiment {
public:
	std::optional<string> filename;

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
			const ValueHistogram<double> combined = normalize(combine(containers));

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
	std::vector<CorrelationAnalyzerParameters> runs;

	void run() {
		PartonicGenerator generator = create_generator();

		std::vector<CorrelationAnalyzer> analyzers;

		for (auto &run : runs) {
			analyzers.emplace_back(run, bins);
		}

		generator.generate_at_loop([&analyzers](std::vector<ParticleContainer> particles, [[maybe_unused]] ParticleGenerator *particle_generator) {
			for (auto &analyzer : analyzers) {
				analyzer.book(&particles);
			}
		}, [&analyzers, this] {
			for (auto &analyzer : analyzers) {
				const auto normalized = normalize(analyzer.histogram);
				cout << "\nAzimuth histogram\n";
				normalized.print_with_bars();
				cout << "\n";
				if (analyzer.parameters.filename) {
					normalized.export_histogram(*analyzer.parameters.filename);
				}
			}
		});

		// std::vector<std::vector<ParticleContainer>> pions;
		// PartonicGenerator generator = create_generator();

		// generator.generate([&pions](std::vector<std::vector<ParticleContainer>> particles, [[maybe_unused]] ParticleGenerator *particle_generator) {
		// 	#pragma omp critical
		// 	pions.insert(pions.end(), particles.begin(), particles.end());
		// }, [&pions, this]() {
		// 	const auto run_count = runs.size();
		// 	for (std::vector<CorrelationAnalyzerParameters>::size_type i = 0; i < run_count; i++) {
		// 		const auto params = runs[i];

		// 		CorrelationAnalyzer analyzer(params, &pions);
		// 		analyzer.bins = bins;
		// 		analyzer.run_index = i;
		// 		analyzer.run_count = run_count;

		// 		const auto hist = analyzer.analyze();
		// 		const auto normalized = normalize(hist);
		// 		cout << "\nAzimuth histogram\n";
		// 		normalized.print_with_bars();
		// 		cout << "\n";
		// 		if (params.filename) {
		// 			normalized.export_histogram(*params.filename);
		// 		}
		// 	}
		// });
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
	ac.pT_range = OptionalRange<double>(1.0, 5.0);

	ac.include_decayed = true;
	ac.mpi = true;

	ac.use_biasing = true;
	ac.bias_power = 4.0;

	ac.parallelize = true;
	ac.pythia_printing = false;
	
	ac.variable_seed = true;
	ac.random_seed = 1;

	ac.normalization = Experiment::Normalization::None;

	// ac.runs = Constants::CorrelationRuns::birapidity_window_test;
	ac.runs = Constants::CorrelationRuns::STAR7;

	ac.run();

	return 0;
}

#endif // EXPERIMENTS_H
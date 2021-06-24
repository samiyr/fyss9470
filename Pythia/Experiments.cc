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
		None, Unity, Count, Trigger
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
	bool include_decayed = Defaults::include_decayed;
	/// Whether to use partonic pT biasing.
	bool use_biasing = Defaults::use_biasing;
	/// The power to which partonic pT will be raised when biasing.
	/// Ignored if `use_biasing` is set to `false`.
	double bias_power = Defaults::bias_power;

	double bias_reference = Defaults::bias_reference;
	/// Use OpenMP multithreading if enabled.
	bool parallelize = Defaults::parallelize;
	/// Controls the 'Print:quiet' Pythia flag.
	bool pythia_printing = Defaults::pythia_printing;
	/// The base random seed used in Pythia. See also `variable_seed`.
	int random_seed = Defaults::random_seed;
	/// If set, the random seed forwarded to Pythia will be
	/// `random_seed` + i, where i is the iterator index of
	/// the partonic pT bins `pT_hat_bins`.
	bool variable_seed = Defaults::variable_seed;
	/// Whether to use multiparton interactions.
	bool mpi = Defaults::mpi;
	/// Normalization type to use.
	Normalization normalization;
	/// Runs the experiment. Subclasses must implement this method.
	void run() {
		abort();
	}
	/// Initializes an instance of `PartonicGenerator`
	/// with the appropriate parameters set.
	PartonicGenerator create_generator() {
		GeneratorParameters params;

		params.cm_energy = energy;
		params.event_count = count;
		params.include_decayed = include_decayed;
		params.y_range = y_range;
		params.pT_range = pT_range;
		params.use_biasing = use_biasing;
		params.bias_power = bias_power;
		params.bias_reference = bias_reference;
		params.pythia_printing = pythia_printing;
		params.random_seed = random_seed;
		params.mpi = mpi;
		params.particle_ids = Constants::pions;

		PartonicGenerator generator(params, pT_hat_bins);
		generator.parallelize = parallelize;
		generator.variable_seed = variable_seed;

		return generator;
	}

	template<typename T>
	ValueHistogram<T> normalize(ValueHistogram<T> hist) {
		switch (normalization) {
			case Normalization::Unity:
				return hist.normalize_to_unity();
				break;
			case Normalization::Count:
				return hist.normalize_by(hist.total());
				break;
			default:
				return hist;
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
			const ValueHistogram<double> combined = normalize(ValueHistogram<double>::combine(containers));

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

		generator.generate_at_loop([&analyzers](std::vector<ParticleContainer> particles, [[maybe_unused]] ParticleGenerator *particle_generator, [[maybe_unused]] bool last_event) {
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
	}
};

class DPSExperiment : public Experiment {
public:
	OptionalRange<double> pT_1;
	OptionalRange<double> pT_2;

	OptionalRange<double> y_1;
	OptionalRange<double> y_2;

	double m;
	int A;
	int B;
	double sigma_eff;

	void run() {
		std::vector<ValueHistogram<double>> containers;
		PartonicGenerator generator = create_generator();

		double sigma_sps_1 = 0.0;
		double sigma_sps_2 = 0.0;

		ParticleFilter filter_1;

		filter_1.allowed_particle_ids = Constants::pions;
		filter_1.include_decayed = include_decayed;
		filter_1.pT_range = pT_1;
		filter_1.y_range = y_1;

		ParticleFilter filter_2;

		filter_2.allowed_particle_ids = Constants::pions;
		filter_2.include_decayed = include_decayed;
		filter_2.pT_range = pT_2;
		filter_2.y_range = y_2;

		generator.generate_at_loop([&filter_1, &filter_2, &sigma_sps_1, &sigma_sps_2]
			(std::vector<ParticleContainer> particles, [[maybe_unused]] ParticleGenerator *particle_generator, bool last_event) {
			for (auto &particle : particles) {
				// pT and/or y ranges are assumed to be non-overlapping, i.e. a particle passes at most one filter
				if (filter_1.is_allowed(particle)) {
					sigma_sps_1 += particle.event_weight;
				} else if (filter_2.is_allowed(particle)) {
					sigma_sps_2 += particle.event_weight;
				}
			}
			if (last_event) {
				sigma_sps_1 *= particle_generator->sigma() / particle_generator->total_weight();
				sigma_sps_2 *= particle_generator->sigma() / particle_generator->total_weight();
			}
		}, [&sigma_sps_1, &sigma_sps_2] {
			cout << "sps 1: " << sigma_sps_1 << "\n";
			cout << "sps 2: " << sigma_sps_2 << "\n";
		});
	}
};

int main() {
	DPSExperiment dps;

	dps.energy = 200;
	dps.count = 1000000;
	dps.mpi = false;
	dps.pT_hat_bins = {
		OptionalRange<double>(1.0, std::nullopt)
	};

	dps.y_1 = OptionalRange<double>(2.6, 4.1);
	dps.y_2 = OptionalRange<double>(2.6, 4.1);
	dps.y_range = OptionalRange<double>(2.6, 4.1);

	dps.pT_1 = OptionalRange<double>(1.0, 1.4);
	dps.pT_2 = OptionalRange<double>(1.4, 2.0);
	dps.pT_range = OptionalRange<double>(1.0, 2.0);

	dps.include_decayed = true;
	dps.use_biasing = false;
	dps.parallelize = false;

	dps.run();


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
	ac.count = 10'000'000 / 16;
	ac.bins = fixed_range(0.0, M_PI, 20);

	ac.pT_hat_bins = std::vector<OptionalRange<double>>(16, OptionalRange<double>(1.0, std::nullopt));
	ac.y_range = OptionalRange<double>();
	ac.pT_range = OptionalRange<double>(1.0, 2.0);

	ac.include_decayed = true;
	ac.mpi = false;

	ac.use_biasing = true;
	ac.bias_power = 4.0;

	ac.parallelize = true;
	ac.pythia_printing = false;
	
	ac.variable_seed = true;
	ac.random_seed = 1;

	ac.normalization = Experiment::Normalization::Unity;

	ac.runs = {
		CorrelationAnalyzerParameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			std::nullopt, std::nullopt,
			"Tests/delta_phi_1e7_2641__1014_1420__unity.csv"),
		CorrelationAnalyzerParameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			-5.0, 5.0,
			"Tests/delta_phi_1e7_2641_5050_1014_1420__unity.csv"),
		CorrelationAnalyzerParameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			0.0, 5.0,
			"Tests/delta_phi_1e7_2641_0050_1014_1420__unity.csv"),
		CorrelationAnalyzerParameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			1.0, 4.5,
			"Tests/delta_phi_1e7_1045_2641_1014_1420__unity.csv"),
		CorrelationAnalyzerParameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			2.6, 4.1,
			"Tests/delta_phi_1e7_2641_2641_1014_1420__unity.csv"),
	};

	// ac.run();

	return 0;
}

#endif // EXPERIMENTS_H
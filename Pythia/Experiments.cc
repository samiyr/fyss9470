#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Pythia8/Pythia.h"
#include "PartonicGenerator.cc"
#include "Histogram.cc"
#include "Analyzer.cc"
#include "Constants.cc"
#include "EventGenerator.cc"

using namespace Pythia8;

/**
 * Represents a generic Pythia experiment.
 * Must be subclassed.
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
	virtual void run() = 0;

	GeneratorParameters create_parameters() {
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

		return params;
	}

	/// Initializes an instance of `PartonicGenerator`
	/// with the appropriate parameters set.
	PartonicGenerator create_generator() {
		GeneratorParameters params = create_parameters();
		// GeneratorParameters params;

		// params.cm_energy = energy;
		// params.event_count = count;
		// params.include_decayed = include_decayed;
		// params.y_range = y_range;
		// params.pT_range = pT_range;
		// params.use_biasing = use_biasing;
		// params.bias_power = bias_power;
		// params.bias_reference = bias_reference;
		// params.pythia_printing = pythia_printing;
		// params.random_seed = random_seed;
		// params.mpi = mpi;
		// params.particle_ids = Constants::pions;

		PartonicGenerator generator(params, pT_hat_bins);
		generator.parallelize = parallelize;
		generator.variable_seed = variable_seed;

		return generator;
	}

	template<typename T>
	ValueHistogram<T> normalize(Normalization _normalization, ValueHistogram<T> hist) {
		switch (_normalization) {
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

	template<typename T>
	ValueHistogram<T> normalize(ValueHistogram<T> hist) {
		return normalize(normalization, hist);
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
	std::vector<AnalysisParameters> runs;

	void run() {
		PartonicGenerator generator = create_generator();

		std::vector<Analyzer> analyzers;

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
	}
};

class DPSExperiment : public Experiment {
public:
	enum class MPIStrategy {
		Disabled, PythiaMPI, DPS
	};
	std::vector<AnalysisParameters> runs;
	MPIStrategy mpi_strategy;

	double m;
	int A;
	int B;
	double sigma_eff;

	void run() {
		switch(mpi_strategy) {
			case MPIStrategy::Disabled:
				mpi = false;
				break;
			case MPIStrategy::PythiaMPI:
				mpi = true;
				break;
			case MPIStrategy::DPS:
				mpi = false;
				break;
		}

		std::vector<std::vector<EventGenerator::Result>> per_thread_results;
		#pragma omp parallel if(parallelize)
		{
			std::vector<std::vector<EventGenerator::Result>> _results;
			#pragma omp for nowait
			for (std::vector<OptionalRange<double>>::size_type i = 0; i < pT_hat_bins.size(); i++) {
				auto params = create_parameters();
				if (variable_seed) {
					params.random_seed += i;
				}
				auto const range = pT_hat_bins[i];
				EventGenerator gen(params, bins, range, runs);

				std::vector<EventGenerator::Result> result = gen.run();
				_results.push_back(result);
			}
			#pragma omp critical
			per_thread_results.insert(per_thread_results.end(), _results.begin(), _results.end());
		}

		std::vector<EventGenerator::Result> results = EventGenerator::Result::combine(per_thread_results);

		for (auto &result : results) {
			auto normalized = normalize(Normalization::Unity, result.histogram);
			normalized.print_with_bars();

			if (mpi_strategy == MPIStrategy::DPS) {
				const double sps1 = result.sigma_sps_1 / pT_hat_bins.size();
				const double sps2 = result.sigma_sps_2 / pT_hat_bins.size();
				const double ssps = result.sigma_sps / pT_hat_bins.size();

				const double dps = (m * A * B / sigma_eff) * sps1 * sps2;
				const double sps = A * B * ssps;

				const double den = sps + dps;
				const double alpha = sps / den;
				const double beta = dps / den;

				cout << "sps = " << sps << "\n";
				cout << "dps = " << dps << "\n";
				cout << "alpha = ";
				print_with_precision(alpha, 12);
				cout << "beta = ";
				print_with_precision(beta, 12);

				normalized *= alpha;
				normalized += beta / M_PI;

				normalized.print_with_bars();
			}
		}
	}
};

int main() {
	DPSExperiment dps;

	dps.energy = 200;
	dps.count = 1'000'000 / 16;
	dps.mpi_strategy = DPSExperiment::MPIStrategy::DPS;
	dps.bins = fixed_range(0.0, M_PI, 20);
	dps.pT_hat_bins = std::vector<OptionalRange<double>>(16, OptionalRange<double>(1.0, std::nullopt));

	dps.runs = {
		AnalysisParameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			2.6, 4.1,
			std::nullopt)
	};

	dps.pT_range = OptionalRange<double>(1.0, 2.0);
	dps.y_range = OptionalRange<double>(2.6, 4.1);

	dps.include_decayed = true;
	dps.use_biasing = true;
	dps.parallelize = true;
	dps.pythia_printing = false;

	dps.variable_seed = true;
	dps.random_seed = 1;

	dps.m = 1;
	dps.A = 1;
	dps.B = 1;
	dps.sigma_eff = 10;

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
	ac.count = 100'000 / 16;
	ac.bins = fixed_range(0.0, M_PI, 20);

	ac.pT_hat_bins = std::vector<OptionalRange<double>>(16, OptionalRange<double>(1.0, std::nullopt));
	ac.y_range = OptionalRange<double>();
	ac.pT_range = OptionalRange<double>(1.0, 2.0);

	ac.include_decayed = true;
	ac.mpi = true;

	ac.use_biasing = true;
	ac.bias_power = 4.0;

	ac.parallelize = true;
	ac.pythia_printing = false;
	
	ac.variable_seed = true;
	ac.random_seed = 1;

	ac.normalization = Experiment::Normalization::Unity;

	ac.runs = {
		AnalysisParameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			2.6, 4.1,
			std::nullopt)
	};

	// ac.run();

	return 0;
}

#endif // EXPERIMENTS_H
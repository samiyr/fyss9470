#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Pythia8/Pythia.h"
#include "PartonicGenerator.cc"
#include "Histogram.cc"
#include "Analyzer.cc"
#include "Constants.cc"
#include "EventGenerator.cc"
#include "Beam.cc"
#include "Around.cc"

using namespace Pythia8;

/**
 * Represents a generic Pythia experiment.
 * Must be subclassed.
 */
class Experiment {
public:
	/**
	 * Types of normalization to apply in data analysis.
	 */
	enum class Normalization {
		/// No normalization
		None, 
		/// Normalize by integral
		Unity, 
		/// Normalize by event count
		Count
	};
	/// Center-of-mass energy in GeV. Required.
	double energy;
	/// Number of events per partonic bin. Required.
	int count;
	/// Histogram bins. Required, must be non-empty.
	std::vector<double> bins;
	/// Partonic bins for generating high-pT particles. Required, must be non-empty.
	std::vector<OptionalRange<double>> pT_hat_bins;
	/// Allowed rapidity range. Defaults to (-inf, inf).
	OptionalRange<double> y_range = OptionalRange<double>();
	/// Allowed transverse momentum range. Defaults to (-inf, inf).
	OptionalRange<double> pT_range = OptionalRange<double>();
	/// Whether to include decayed particles. Default specified in Constants.cc.
	bool include_decayed = Defaults::include_decayed;
	/// Whether to use partonic pT biasing. Default specified in Constants.cc.
	bool use_biasing = Defaults::use_biasing;
	/// The power to which partonic pT will be raised when biasing.
	/// Default specified in Constants.cc.
	/// Ignored if `use_biasing` is set to `false`.
	double bias_power = Defaults::bias_power;
	/// Reference bias used in biasing.
	/// Default specified in Constants.cc.
	/// See Pythia's documentation for more.
	/// Ignored if `use_biasing` is set to `false`.
	double bias_reference = Defaults::bias_reference;
	/// Use OpenMP multithreading if enabled. Default specified in Constants.cc.
	bool parallelize = Defaults::parallelize;
	/// Controls the 'Print:quiet' Pythia flag. Default specified in Constants.cc.
	bool pythia_printing = Defaults::pythia_printing;
	/// The base random seed used in Pythia. See also `variable_seed`.
	/// Default specified in Constants.cc.
	int random_seed = Defaults::random_seed;
	/// If set, the random seed forwarded to Pythia will be
	/// `random_seed` + i, where i is the iterator index of
	/// the partonic pT bins `pT_hat_bins`.
	/// Default specified in Constants.cc.
	bool variable_seed = Defaults::variable_seed;
	/// Whether to use Pythia's multiparton interactions. 
	/// Default specified in Constants.cc.
	/// Provided value superseded by `mpi_strategy` in `DPSExperiment`. 
	bool mpi = Defaults::mpi;
	/// Normalization type to use.
	Normalization normalization;
	/// Beam A parameters.
	Beam beam_A;
	/// Beam B parameters.
	Beam beam_B;

	/// Runs the experiment. Subclasses must implement this method.
	virtual void run() = 0;

	/// Packages the given parameters to a `GeneratorParameters` instance.
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
		params.beam_A = beam_A;
		params.beam_B = beam_B;

		return params;
	}

	/**
	 * Initializes an instance of `PartonicGenerator`
	 * with the appropriate parameters set.
	 */
	PartonicGenerator create_generator() {
		GeneratorParameters params = create_parameters();
		PartonicGenerator generator(params, pT_hat_bins);
		generator.parallelize = parallelize;
		generator.variable_seed = variable_seed;

		return generator;
	}

	/**
	 * Normalizes a histogram according to the given normalization type.
	 * \param _normalization The requested normalization, of type `Normalization`.
	 * \param hist Input histogram.
	 */
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
	/**
	 * Normalizes a histogram according to the normalization type specified in member `normalization`.
	 * \param hist Input histogram.
	 */
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
	std::vector<Analyzer::Parameters> runs;

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

/**
 * DPS experiment.
 */
class DPSExperiment : public Experiment {
public:
	/**
	 * Enum encapsulating the strategy with multiparton interactions.
	 */
	enum class MPIStrategy {
		/// No MPI
		Disabled, 
		/// Use Pythia's MPI model
		PythiaMPI, 
		/// Use an analytic DPS model
		DPS
	};
	/**
	 * A list of analysis parameters,
	 * each of which will produce its own
	 * analysis result using the same events.
	 */
	std::vector<Analyzer::Parameters> runs;
	/**
	 * Strategy for dealing with multiparton interactions.
	 * This setting will automatically change the `mpi` field.
	 * Manually changing the `mpi` field has no effect.
	 */
	MPIStrategy mpi_strategy;

	void run() {
		// Set the `mpi` field to the appropriate field.
		switch(mpi_strategy) {
			// MPI disabled ==> mpi = false;
			case MPIStrategy::Disabled:
				mpi = false;
				break;
			// Use Pythia's model ==> mpi = true;
			case MPIStrategy::PythiaMPI:
				mpi = true;
				break;
			// Use DPS model ==> mpi = false to avoid double counting.
			case MPIStrategy::DPS:
				mpi = false;
				break;
		}

		// Initalize a shared vector collecting the results from each thread.
		std::vector<std::vector<EventGenerator::Result>> per_thread_results;
		#pragma omp parallel if(parallelize)
		{
			// A local variable for collecting an individual thread's results.
			std::vector<std::vector<EventGenerator::Result>> _results;
			// Loop over partonic pT bins in parallel (if enabled)
			#pragma omp for nowait
			for (std::vector<OptionalRange<double>>::size_type i = 0; i < pT_hat_bins.size(); i++) {
				auto params = create_parameters();
				// Change seeds if needed
				if (variable_seed) {
					params.random_seed += i;
				}
				auto const range = pT_hat_bins[i];
				// Create an event generator with the given parameters
				EventGenerator gen(params, bins, range, runs);
				// Generate events and collect results...
				std::vector<EventGenerator::Result> result = gen.run();
				// ...and append to the local variable
				_results.push_back(result);
			}
			// Collect per-thread results to the shared variable, one thread at a time to avoid data races
			#pragma omp critical
			per_thread_results.insert(per_thread_results.end(), _results.begin(), _results.end());
		}
		// Combine the per-thread results to a single list of results, each corresponding to a single run
		std::vector<EventGenerator::Result> results = EventGenerator::Result::combine(per_thread_results);
		// Analyze the results, one run at a time
		for (std::vector<EventGenerator::Result>::size_type run_index = 0; run_index < results.size(); run_index++) {
			auto result = results[run_index];
			auto normalized = normalize(Normalization::Unity, result.histogram);

			cout << "\n";
			cout << "*------------------------------------------------------------------------------------*" << "\n";
			cout << "|                                                                                    |" << "\n";
			cout << "|                                       RUN #" << run_index + 1 << "                                       |\n";
			cout << "|                                                                                    |" << "\n";
			cout << "*------------------------------------------------------------------------------------*" << "\n";
			cout << "\n";

			if (mpi_strategy == MPIStrategy::DPS) {
				const int A = beam_A.nucleus.mass_number;
				const int B = beam_B.nucleus.mass_number;

				const Around<double> sps1 = Around(result.sigma_sps_1);
				const Around<double> sps2 = Around(result.sigma_sps_2);
				const Around<double> ssps = Around(result.sigma_sps);

				const double m = result.parameters.m;
				const double sigma_pp = result.parameters.sigma_eff;
				const double sigma_eff = calculate_sigma_eff(beam_B, sigma_pp);

				const Around<double> dps = (m * A * B / sigma_eff) * sps1 * sps2;
				const Around<double> sps = (double)A * (double)B * ssps;

				const Around<double> den = sps + dps;
				const Around<double> alpha = sps / den;
				const Around<double> beta = dps / den;

				cout << "pT_1\t\t= " << result.parameters.pT_small.extent() << "\n";
				cout << "pT_2\t\t= " << result.parameters.pT_large.extent() << "\n";
				cout << "y_1\t\t= " << result.parameters.y_small.extent() << "\n";
				cout << "y_2\t\t= " << result.parameters.y_large.extent() << "\n";
				cout << "m\t\t= " << m << "\n";
				cout << "sigma_eff\t= " << sigma_eff << "\n";
				cout << "sps1\t\t= "; print_with_precision(sps1, 6);
				cout << "sps2\t\t= "; print_with_precision(sps2, 6);
				cout << "ssps\t\t= "; print_with_precision(ssps, 6);
				cout << "sps\t\t= "; print_with_precision(sps, 6);
				cout << "dps\t\t= "; print_with_precision(dps, 6);
				cout << "alpha\t\t= "; print_with_precision(alpha, 6);
				cout << "beta\t\t= "; print_with_precision(beta, 6);

				normalized *= alpha;
				normalized += beta / M_PI;
			}
			normalized.print_with_bars();

			if (result.parameters.filename) {
				normalized.export_histogram(*result.parameters.filename);
			}
		}
	}
};

int main() {
	/* --- DPS --- */

	DPSExperiment dps;

	dps.energy = 200;
	dps.count = 10'000'000 / 16;
	dps.mpi_strategy = DPSExperiment::MPIStrategy::DPS;
	dps.bins = fixed_range(0.0, M_PI, 20);
	dps.pT_hat_bins = std::vector<OptionalRange<double>>(16, OptionalRange<double>(10.0, std::nullopt));

	dps.runs = {
		Analyzer::Parameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			2.6, 4.1,
			"STAR7/delta_phi_1e7_1014_1420_2641_2641_250_dps_pp_100.csv",
			1.0,
			25.0),
	};

	dps.pT_range = OptionalRange<double>(1.0, 2.0);
	dps.y_range = OptionalRange<double>(2.6, 4.1);

	dps.include_decayed = true;
	dps.use_biasing = true;
	dps.parallelize = true;
	dps.pythia_printing = false;

	dps.variable_seed = true;
	dps.random_seed = 1;

	dps.beam_A = Beam();
	dps.beam_B = Beam();
	// dps.beam_B = Beam(13, 27, Beam::NuclearPDF::EPPS16NLO, false);

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
		Analyzer::Parameters(
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
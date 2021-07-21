#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "PartonicGenerator.cc"
#include "Histogram.cc"
#include "Analyzer.cc"
#include "Constants.cc"
#include "EventGenerator.cc"
#include "Beam.cc"
#include "Around.cc"
#include <chrono>

/**
 * Represents a generic Pythia experiment.
 * Must be subclassed.
 */
class Experiment {
public:
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
	/// Should be enabled when parallelizing Pythia by specifying
	/// multiple pT_hat bins.
	bool variable_seed = Defaults::variable_seed;
	/// Whether to use Pythia's multiparton interactions. 
	/// Default specified in Constants.cc.
	/// Provided value superseded by `mpi_strategy` in `DPSExperiment`. 
	bool mpi = Defaults::mpi;
	/// Normalization type to use.
	Normalization normalization;
	/// Beam A parameters. Should be a proton since nuclear structure is not
	/// taken into account when calculating the effective cross section (sigma_eff).
	Beam beam_A = Beam();
	/// Beam B parameters. Required.
	Beam beam_B;
	/// Type of process, e.g. Process::HardQCD or Process::SoftQCDNonDiffractive. Required.
	Process process;
	/// Path component appended to the filename during export.
	/// Either the filename or the working directory has to include
	/// a path component separator '/' since one isn't automatically
	/// appended to the export file. Defaults to nullopt.
	std::optional<string> working_directory = std::nullopt;
	/// The file extension used in the exported histogram. 
	/// Default specified in Constants.cc. Must include '.', e.g. '.csv'.
	string histogram_file_extension = Defaults::histogram_file_extension;
	/// The file extension used in the exported parameter data file. 
	/// Default specified in Constants.cc. Must include '.', e.g. '.txt'.
	string run_data_file_extension = Defaults::run_data_file_extension;

	bool cross_section_error;
	bool histogram_fluctuation_error;

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
		params.process = process;

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
	ValueHistogram<T> normalize(Normalization _normalization, ValueHistogram<T> hist, double param = 0.0) {
		switch (_normalization) {
			case Normalization::Unity:
				return hist.normalize_to_unity();
				break;
			case Normalization::Count:
				return hist.normalize_by(hist.total());
				break;
			case Normalization::STARC:
				return hist.normalize_to_star_C(param);
				break;
			default:
				return hist;
				break;
		}
	}

	ValueHistogram<double> normalize(Normalization _normalization, EventGenerator::Result result) {
		const auto histogram = ValueHistogram<double>::combine(result.histograms, histogram_fluctuation_error);
		switch (_normalization) {
			case Normalization::Unity:
				return histogram.normalize_to_unity();
				break;
			case Normalization::Count:
				return histogram.normalize_by(histogram.total());
				break;
			case Normalization::STARC:
				return histogram.normalize_to_star_C(result.N_trigger);
				break;
			default:
				return histogram;
				break;
		}
	}
	/**
	 * Normalizes a histogram according to the normalization type specified in member `normalization`.
	 * \param hist Input histogram.
	 */
	template<typename T>
	ValueHistogram<T> normalize(ValueHistogram<T> hist, double param = 0.0) {
		return normalize(normalization, hist, param);
	}

	ValueHistogram<double> normalize(EventGenerator::Result result) {
		return normalize(normalization, result);
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
		// Start the clock.
		const auto start_time = std::chrono::high_resolution_clock::now();
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
			auto normalized = normalize(result);
			// Stop the clock
			const auto end_time = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double> elapsed_duration = end_time - start_time;
			const double elapsed_time = elapsed_duration.count();
			// Prefix for the export file path.
			const string prefix = working_directory ? *working_directory : "";
			// Initialize data parameter output file.
			ofstream file;
			if (result.parameters.filename) {
				// Construct data output file path.
				const string filename = prefix + *result.parameters.filename + run_data_file_extension;
				file.open(filename);
				file << std::setprecision(6);

				file << "elapsed_time\t= " << elapsed_time << " s\n\n";
			}

			if (mpi_strategy == MPIStrategy::DPS) {
				const int A = beam_A.nucleus.mass_number;
				const int B = beam_B.nucleus.mass_number;

				Around<double> sps1 = Around(result.sigma_sps_1);
				Around<double> sps2 = Around(result.sigma_sps_2);
				Around<double> ssps = Around(result.sigma_sps);

				if (!cross_section_error) {
					sps1.error = std::nullopt;
					sps2.error = std::nullopt;
					ssps.error = std::nullopt;
				}

				const double m = result.parameters.m;
				const double sigma_pp = result.parameters.sigma_eff;
				const double sigma_eff = calculate_sigma_eff(beam_B, sigma_pp);

				const Around<double> dps = (m * A * B / sigma_eff) * sps1 * sps2;
				const Around<double> sps = (double)A * (double)B * ssps;

				const Around<double> den = sps + dps;
				const Around<double> alpha = sps / den;
				const Around<double> beta = dps / den;

				if (result.parameters.filename) {
					file << "eCM\t\t= " << energy << "\n";
					file << "count\t\t= " << count * pT_hat_bins.size() << "\n\n";

					file << "process\t\t= " << to_string(process) << "\n";
					file << "include_decayed\t= " << bool_to_string(include_decayed) << "\n";
					file << "mpi\t\t= " << to_string(mpi_strategy) << "\n\n";
					file << "normalization\t= " << to_string(normalization) << "\n\n";

					file << "beam_A\t\t= " << beam_A << "\n";
					file << "beam_B\t\t= " << beam_B << "\n\n";

					file << "use_biasing\t= " << bool_to_string(use_biasing) << "\n";
					file << "bias_power\t= " << bias_power << "\n";
					file << "bias_reference\t= " << bias_reference << "\n\n";

					file << "pT_1\t\t= " << result.parameters.pT_small.extent() << "\n";
					file << "pT_2\t\t= " << result.parameters.pT_large.extent() << "\n";
					file << "y_1\t\t= " << result.parameters.y_small.extent() << "\n";
					file << "y_2\t\t= " << result.parameters.y_large.extent() << "\n\n";

					file << "m\t\t= " << m << "\n";
					file << "sigma_pp\t= " << sigma_pp << "\n";
					file << "sigma_eff\t= " << sigma_eff << "\n\n";

					file << "sps1\t\t= " << sps1 << "\n";
					file << "sps2\t\t= " << sps2 << "\n";
					file << "ssps\t\t= " << ssps << "\n\n";

					file << "sps\t\t= " << sps << "\n";
					file << "dps\t\t= " << dps << "\n\n";

					file << "alpha\t\t= " << alpha << "\n";
					file << "beta\t\t= " << beta << "\n\n";
				}

				normalized *= alpha;
				normalized += beta / M_PI;
			}
			cout << normalized;

			if (result.parameters.filename) {
				file << "azimuthal histogram:";
				file << normalized;
				file << "\n";
				file << "pT_hat_bins:\n";
				for (auto bin : pT_hat_bins) {
					file << bin.extent() << "\n";
				}
				file.close();
				normalized.export_histogram(prefix + *result.parameters.filename + histogram_file_extension);
			}
		}
	}
};

int main() {
	ValueHistogram<double> h1({0.0, 1.0, 2.0});
	h1.fill(0.5, 3.0);
	h1.fill(1.5, 1.0);

	ValueHistogram<double> h2({0.0, 1.0, 2.0});
	h2.fill(0.7, 3.0);
	h2.fill(1.3, 1.0);

	ValueHistogram<double> h3({0.0, 1.0, 2.0});
	h3.fill(0.5, 2.0);
	h3.fill(1.5, 2.0);

	cout << ValueHistogram<double>::combine({h1, h2});
	cout << ValueHistogram<double>::combine({h2, h3});

	/* --- DPS --- */

	DPSExperiment dps;

	dps.process = Process::HardQCD;
	dps.energy = 200;
	dps.count = 100'000'000 / 16;
	dps.mpi_strategy = MPIStrategy::DPS;
	dps.normalization = Normalization::Unity;
	dps.bins = fixed_range(0.0, M_PI, 20);
	dps.pT_hat_bins = std::vector<OptionalRange<double>>(16, OptionalRange<double>(1.5, std::nullopt));
	dps.cross_section_error = true;
	dps.histogram_fluctuation_error = true;

	dps.runs = {
		Analyzer::Parameters(
			1.0, 1.4,
			1.4, 2.0,
			2.6, 4.1,
			2.6, 4.1,
			"data1",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.1, 1.5,
			1.5, 2.1,
			2.6, 4.1,
			2.6, 4.1,
			"data2",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.2, 1.6,
			1.6, 2.2,
			2.6, 4.1,
			2.6, 4.1,
			"data3",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.3, 1.7,
			1.7, 2.3,
			2.6, 4.1,
			2.6, 4.1,
			"data4",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.4, 1.8,
			1.8, 2.4,
			2.6, 4.1,
			2.6, 4.1,
			"data5",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.5, 1.9,
			1.9, 2.5,
			2.6, 4.1,
			2.6, 4.1,
			"data6",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.6, 2.0,
			2.0, 2.6,
			2.6, 4.1,
			2.6, 4.1,
			"data7",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.7, 2.1,
			2.1, 2.7,
			2.6, 4.1,
			2.6, 4.1,
			"data8",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.8, 2.2,
			2.2, 2.8,
			2.6, 4.1,
			2.6, 4.1,
			"data9",
			1.0,
			25.0),
		Analyzer::Parameters(
			1.9, 2.3,
			2.3, 2.9,
			2.6, 4.1,
			2.6, 4.1,
			"data10",
			1.0,
			25.0),
		Analyzer::Parameters(
			2.0, 2.4,
			2.4, 3.0,
			2.6, 4.1,
			2.6, 4.1,
			"data11",
			1.0,
			25.0),
	};

	dps.pT_range = OptionalRange<double>(1.0, 2.0);
	dps.y_range = OptionalRange<double>(2.6, 4.1);

	dps.include_decayed = true;
	dps.use_biasing = false;
	dps.parallelize = true;
	dps.pythia_printing = false;

	dps.variable_seed = true;
	dps.random_seed = 1;

	dps.beam_A = Beam();
	dps.beam_B = Beam();
	// dps.beam_B = Beam(13, 27, Beam::NuclearPDF::EPPS16NLO, true);
	// dps.beam_B = Beam(97, 197, Beam::NuclearPDF::EPPS16NLO, true);

	dps.working_directory = "Tests/pT/HardQCD/";
	dps.histogram_file_extension = ".csv";
	dps.run_data_file_extension = ".txt";

	// dps.run();


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

	ac.normalization = Normalization::Unity;

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
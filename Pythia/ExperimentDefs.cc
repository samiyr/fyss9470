#ifndef EXPERIMENT_DEFS_H
#define EXPERIMENT_DEFS_H

#include "PartonicGenerator.cc"
#include "Histogram.cc"
#include "Analyzer.cc"
#include "Constants.cc"
#include "EventGenerator.cc"
#include "Beam.cc"
#include "Around.cc"
#include <chrono>
#include <optional>
#include <filesystem>

extern int THREAD_COUNT;

/**
 * Represents a generic Pythia experiment.
 * Must be subclassed.
 */
class Experiment {
public:
	/// Center-of-mass energy in GeV. Required.
	double energy;
	/// Number of events per partonic bin. Required.
	EVENT_COUNT_TYPE count;
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
	std::optional<std::filesystem::path> working_directory = std::nullopt;
	/// The file extension used in the exported histogram. 
	/// Default specified in Constants.cc. Must include '.', e.g. '.csv'.
	string histogram_file_extension = Defaults::histogram_file_extension;
	/// The file extension used in the exported parameter data file. 
	/// Default specified in Constants.cc. Must include '.', e.g. '.txt'.
	string run_data_file_extension = Defaults::run_data_file_extension;

	bool cross_section_error;
	bool histogram_fluctuation_error;
	bool experimental_histogram_error;

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
		// Start the clock.
		const auto start_time = std::chrono::high_resolution_clock::now();

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
		}, [&containers, this, start_time]{
			ValueHistogram<double> combined = normalize(ValueHistogram<double>::combine(containers, histogram_fluctuation_error));
			if (experimental_histogram_error) {
				combined = ValueHistogram<double>::calculate_statistical_error(combined);
			}
			cout << "Normalized pT histogram" << "\n";
			combined.print();
			cout << "\n";

			// Stop the clock
			const auto end_time = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double> elapsed_duration = end_time - start_time;
			const double elapsed_time = elapsed_duration.count();

			if (filename) {
				const auto histogram_path = construct_path(working_directory, *filename, histogram_file_extension);
				const auto data_path = construct_path(working_directory, *filename, run_data_file_extension);
				ofstream file;
				file.open(data_path);

				file << std::setprecision(6);

				file << "elapsed_time\t= " << elapsed_time << " s\n\n";

				file << "eCM\t\t= " << energy << "\n";
				file << "count\t\t= " << count * (EVENT_COUNT_TYPE)pT_hat_bins.size() << "\n\n";

				file << "process\t\t= " << to_string(process) << "\n";
				file << "include_decayed\t= " << bool_to_string(include_decayed) << "\n";
				file << "mpi\t\t= " << bool_to_string(mpi) << "\n\n";
				file << "normalization\t= " << "cross section" << "\n\n";

				file << "beam_A\t\t= " << beam_A << "\n";
				file << "beam_B\t\t= " << beam_B << "\n\n";

				file << "use_biasing\t= " << bool_to_string(use_biasing) << "\n";
				file << "bias_power\t= " << bias_power << "\n";
				file << "bias_reference\t= " << bias_reference << "\n\n";

				file << "pT\t\t= " << pT_range.extent() << "\n";
				file << "y\t\t= " << y_range.extent() << "\n";

				file << "pT histogram:";
				file << combined;
				file << "\n";
				file << "pT_hat_bins:\n";
				for (auto bin : pT_hat_bins) {
					file << bin.extent() << "\n";
				}
				file.close();

				combined.export_histogram(histogram_path);
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

			if (experimental_histogram_error) {
				normalized = ValueHistogram<double>::calculate_statistical_error(normalized);
			}
			// Stop the clock
			const auto end_time = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double> elapsed_duration = end_time - start_time;
			const double elapsed_time = elapsed_duration.count();
			// Initialize data parameter output file.
			ofstream file;
			if (result.parameters.filename) {
				// Construct data output file path.
				const auto data_path = construct_path(working_directory, *result.parameters.filename, run_data_file_extension);
				file.open(data_path);
				file << std::setprecision(6);

				file << "elapsed_time\t= " << elapsed_time << " s\n\n";

				file << "eCM\t\t= " << energy << "\n";
				file << "count\t\t= " << count * (EVENT_COUNT_TYPE)pT_hat_bins.size() << "\n\n";

				file << "process\t\t= " << to_string(process) << "\n";
				file << "include_decayed\t= " << bool_to_string(include_decayed) << "\n";
				file << "mpi\t\t= " << to_string(mpi_strategy) << "\n\n";
				file << "normalization\t= " << to_string(normalization) << "\n\n";

				file << "beam_A\t\t= " << beam_A << "\n";
				file << "beam_B\t\t= " << beam_B << "\n\n";

				file << "use_biasing\t= " << bool_to_string(use_biasing) << "\n";
				file << "bias_power\t= " << bias_power << "\n";
				file << "bias_reference\t= " << bias_reference << "\n\n";

				file << "pT\t\t= " << pT_range.extent() << "\n";
				file << "y\t\t= " << y_range.extent() << "\n";
				file << "pT_1\t\t= " << result.parameters.pT_small.extent() << "\n";
				file << "pT_2\t\t= " << result.parameters.pT_large.extent() << "\n";
				file << "y_1\t\t= " << result.parameters.y_small.extent() << "\n";
				file << "y_2\t\t= " << result.parameters.y_large.extent() << "\n\n";
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
				const auto histogram_path = construct_path(working_directory, *result.parameters.filename, histogram_file_extension);
				normalized.export_histogram(histogram_path);
			}
		}
	}
};

// --- pT cross section ---

CrossSectionExperiment pT_template(EVENT_COUNT_TYPE count, bool mpi = Defaults::mpi) {
	CrossSectionExperiment cs;

	cs.process = Process::HardQCD;
	cs.mpi = mpi;
	cs.energy = 200;
	cs.count = count;
	cs.bins = {
		1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0
	};

	cs.y_range = OptionalRange<double>(-0.35, 0.35);
	cs.pT_range = OptionalRange<double>();
	cs.include_decayed = true;
	cs.bias_power = 4.0;
	cs.parallelize = true;
	cs.pythia_printing = false;
	cs.normalization = Normalization::None;

	return cs;
}

void pT_cross_section(EVENT_COUNT_TYPE count, bool biasing, bool subdivision, string fn, bool mpi = Defaults::mpi, int seed = -1, bool experimental_histogram_error = false) {
	CrossSectionExperiment cs = pT_template(count, mpi);

	if (subdivision) {
		cs.pT_hat_bins = {
			OptionalRange<double>(2.0, 5.0),
			OptionalRange<double>(5.0, 10.0),
			OptionalRange<double>(10.0, 20.0),
			OptionalRange<double>(20.0, 40.0),
			OptionalRange<double>(40.0, std::nullopt)
		};		
	} else {
		cs.pT_hat_bins = {OptionalRange<double>(2.0, std::nullopt)};
	}

	cs.use_biasing = biasing;
	cs.working_directory = "Data/pp/pT";
	cs.filename = fn;
	cs.random_seed = seed;
	cs.cross_section_error = false;
	cs.histogram_fluctuation_error = false;
	cs.experimental_histogram_error = experimental_histogram_error;

	cs.run();
}

void run_pT_experiment() {
	const auto t1 = std::chrono::high_resolution_clock::now();
	// pT_cross_section(10'000'000, true, true, "pT.csv");
	const auto t2 = std::chrono::high_resolution_clock::now();
	pT_cross_section(100'000, false, false, "MPI/No bias no subdivision/1e5/data.csv");	
	const auto t3 = std::chrono::high_resolution_clock::now();
	pT_cross_section(100'000, true, false, "MPI/Bias no subdivision/1e5/data.csv");
	const auto t4 = std::chrono::high_resolution_clock::now();
	pT_cross_section(100'000 / 5, false, true, "MPI/No bias subdivision/1e5/data.csv");
	const auto t5 = std::chrono::high_resolution_clock::now();
	pT_cross_section(100'000 / 5, true, true, "MPI/Bias subdivision/1e5/data.csv");
	const auto t6 = std::chrono::high_resolution_clock::now();

	const std::chrono::duration<double> duration1 = t2 - t1;
	const std::chrono::duration<double> duration2 = t3 - t2;
	const std::chrono::duration<double> duration3 = t4 - t3;
	const std::chrono::duration<double> duration4 = t5 - t4;
	const std::chrono::duration<double> duration5 = t6 - t5;
	
	cout << "pT1 = " << duration1.count() << "\n";	
	cout << "pT2 = " << duration2.count() << "\n";	
	cout << "pT3 = " << duration3.count() << "\n";	
	cout << "pT4 = " << duration4.count() << "\n";	
	cout << "pT5 = " << duration5.count() << "\n";
}

// --- Azimuth correlation ---

DPSExperiment dps_template(EVENT_COUNT_TYPE count, Process process, MPIStrategy mpi, double pT_hat_min, Beam b, string wd) {
	DPSExperiment dps;

	dps.energy = 200;
	dps.count = count / THREAD_COUNT;
	dps.mpi_strategy = mpi;
	dps.process = process;
	dps.normalization = Normalization::Unity;
	dps.bins = fixed_range(0.0, M_PI, 20);
	dps.pT_hat_bins = std::vector<OptionalRange<double>>(THREAD_COUNT, OptionalRange<double>(pT_hat_min, std::nullopt));
	dps.cross_section_error = false;
	dps.histogram_fluctuation_error = false;
	dps.experimental_histogram_error = true;

	dps.include_decayed = true;
	dps.use_biasing = false;
	dps.parallelize = true;
	dps.pythia_printing = false;

	dps.variable_seed = true;
	dps.random_seed = 1;

	dps.working_directory = wd;

	dps.beam_A = Beam();
	dps.beam_B = b;

	dps.runs = {
		Analyzer::Parameters(1.0, 1.4, 1.4, 2.0, 2.6, 4.1, 2.6, 4.1, "1014_1420_dps10", 0.5, 10.0),
		Analyzer::Parameters(1.0, 1.4, 1.4, 2.0, 2.6, 4.1, 2.6, 4.1, "1014_1420_dps25", 0.5, 25.0),

		Analyzer::Parameters(1.0, 1.4, 2.0, 2.4, 2.6, 4.1, 2.6, 4.1, "1014_2024_dps10", 0.5, 10.0),
		Analyzer::Parameters(1.0, 1.4, 2.0, 2.4, 2.6, 4.1, 2.6, 4.1, "1014_2024_dps25", 0.5, 25.0),
		
		Analyzer::Parameters(1.0, 1.4, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "1014_2428_dps10", 0.5, 10.0),
		Analyzer::Parameters(1.0, 1.4, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "1014_2428_dps25", 0.5, 25.0),
		
		Analyzer::Parameters(1.0, 1.4, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "1014_2850_dps10", 0.5, 10.0),
		Analyzer::Parameters(1.0, 1.4, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "1014_2850_dps25", 0.5, 25.0),		


		Analyzer::Parameters(1.4, 2.0, 2.0, 2.4, 2.6, 4.1, 2.6, 4.1, "1420_2024_dps10", 0.5, 10.0),
		Analyzer::Parameters(1.4, 2.0, 2.0, 2.4, 2.6, 4.1, 2.6, 4.1, "1420_2024_dps25", 0.5, 25.0),
		
		Analyzer::Parameters(1.4, 2.0, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "1420_2428_dps10", 0.5, 10.0),
		Analyzer::Parameters(1.4, 2.0, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "1420_2428_dps25", 0.5, 25.0),
		
		Analyzer::Parameters(1.4, 2.0, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "1420_2850_dps10", 0.5, 10.0),
		Analyzer::Parameters(1.4, 2.0, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "1420_2850_dps25", 0.5, 25.0),	


		Analyzer::Parameters(2.0, 2.4, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "2024_2428_dps10", 0.5, 10.0),
		Analyzer::Parameters(2.0, 2.4, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "2024_2428_dps25", 0.5, 25.0),
		
		Analyzer::Parameters(2.0, 2.4, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "2024_2850_dps10", 0.5, 10.0),
		Analyzer::Parameters(2.0, 2.4, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "2024_2850_dps25", 0.5, 25.0),	


		Analyzer::Parameters(2.4, 2.8, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "2428_2850_dps10", 0.5, 10.0),
		Analyzer::Parameters(2.4, 2.8, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "2424_2850_dps25", 0.5, 25.0),	
	};

	return dps;
}

DPSExperiment mpi_template(EVENT_COUNT_TYPE count, Process process, MPIStrategy mpi, double pT_hat_min, Beam b, string wd) {
	DPSExperiment dps;

	dps.energy = 200;
	dps.count = count / THREAD_COUNT;
	dps.mpi_strategy = mpi;
	dps.process = process;
	dps.normalization = Normalization::Unity;
	dps.bins = fixed_range(0.0, M_PI, 20);
	dps.pT_hat_bins = std::vector<OptionalRange<double>>(THREAD_COUNT, OptionalRange<double>(pT_hat_min, std::nullopt));
	dps.cross_section_error = false;
	dps.histogram_fluctuation_error = false;
	dps.experimental_histogram_error = true;

	dps.include_decayed = true;
	dps.use_biasing = false;
	dps.parallelize = true;
	dps.pythia_printing = false;

	dps.variable_seed = true;
	dps.random_seed = 1;

	dps.working_directory = wd;

	dps.beam_A = Beam();
	dps.beam_B = b;

	dps.runs = {
		Analyzer::Parameters(1.0, 1.4, 1.4, 2.0, 2.6, 4.1, 2.6, 4.1, "1014_1420"),
		Analyzer::Parameters(1.0, 1.4, 2.0, 2.4, 2.6, 4.1, 2.6, 4.1, "1014_2024"),
		Analyzer::Parameters(1.0, 1.4, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "1014_2428"),
		Analyzer::Parameters(1.0, 1.4, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "1014_2850"),

		Analyzer::Parameters(1.4, 2.0, 2.0, 2.4, 2.6, 4.1, 2.6, 4.1, "1420_2024"),
		Analyzer::Parameters(1.4, 2.0, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "1420_2428"),
		Analyzer::Parameters(1.4, 2.0, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "1420_2850"),

		Analyzer::Parameters(2.0, 2.4, 2.4, 2.8, 2.6, 4.1, 2.6, 4.1, "2024_2428"),
		Analyzer::Parameters(2.0, 2.4, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "2024_2850"),

		Analyzer::Parameters(2.4, 2.8, 2.8, 5.0, 2.6, 4.1, 2.6, 4.1, "2428_2850"),
	};

	return dps;
}

void pp_dps_run(EVENT_COUNT_TYPE count, Process process, MPIStrategy mpi, double pT_hat_min, string wd) {
	DPSExperiment dps = dps_template(count, process, mpi, pT_hat_min, Beam(), wd);
	dps.run();
}

void pp_mpi_run(EVENT_COUNT_TYPE count, Process process, MPIStrategy mpi, double pT_hat_min, string wd) {
	DPSExperiment dps = mpi_template(count, process, mpi, pT_hat_min, Beam(), wd);
	dps.run();
}

void Al_run(EVENT_COUNT_TYPE count, Process process, MPIStrategy mpi, double pT_hat_min, string wd, bool nPDF = false) {
	DPSExperiment dps = dps_template(count, process, mpi, pT_hat_min, Beam(13, 27, Beam::NuclearPDF::EPPS16NLO, nPDF), wd);
	dps.run();
}

void Au_run(EVENT_COUNT_TYPE count, Process process, MPIStrategy mpi, double pT_hat_min, string wd, bool nPDF = false) {
	DPSExperiment dps = dps_template(count, process, mpi, pT_hat_min, Beam(97, 197, Beam::NuclearPDF::EPPS16NLO, nPDF), wd);
	dps.run();
}

#endif // EXPERIMENT_DEFS_H
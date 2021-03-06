#ifndef EVENT_GENERATOR_H
#define EVENT_GENERATOR_H

/// An abstraction for generating events for azimuthal correlations.
class EventGenerator {
public:
	/// A container for the event generation results.
	struct Result {
		/// A list of histograms, each corresponding to a pT_hat range.
		std::vector<ValueHistogram<double>> histograms;
		/// A list of nuclear histograms, each correspnding to a pT_hat range.
		std::vector<ValueHistogram<double>> nuclear_histograms;

		std::vector<ValueHistogram<double>> partial_histograms;
		/// A list of momentum fraction histograms.
		std::vector<ValueHistogram<double>> x1_pre_histograms;
		std::vector<ValueHistogram<double>> x2_pre_histograms;
		std::vector<ValueHistogram<double>> x1_post_histograms;
		std::vector<ValueHistogram<double>> x2_post_histograms;
		/// The total weight of trigger particles.
		double N_trigger;
		/// The total weight of associated particles.
		double N_assoc;
		/// THe total weight of paired particles.
		double N_pair;

		double partial_N_trigger;
		double partial_N_assoc;
		double partial_N_pair;

		/// A list of SPS integrals, each corresponding to a pT_hat range.
		std::vector<double> sigma_sps_1;
		std::vector<double> sigma_sps_2;
		std::vector<double> sigma_sps;
		/// The total cross section and weight, given by Pythia.
		double sigma_gen;
		double total_weight;
		/// The parameter set associated with this run.
		Analyzer::Parameters parameters;
		/// The index of the current run.
		std::vector<Analyzer>::size_type run_index;

		EVENT_COUNT_TYPE next_calls;
		EVENT_COUNT_TYPE next_skips;
		EVENT_COUNT_TYPE partial_events;
		ValueHistogram<double> partial_ncolls;
		/// Overloads the += operator for combining multiple Result objects.
		Result& operator+=(Result rhs) {
			histograms.insert(histograms.end(), rhs.histograms.begin(), rhs.histograms.end());
			nuclear_histograms.insert(nuclear_histograms.end(), rhs.nuclear_histograms.begin(), rhs.nuclear_histograms.end());
			partial_histograms.insert(partial_histograms.end(), rhs.partial_histograms.begin(), rhs.partial_histograms.end());

			x1_pre_histograms.insert(x1_pre_histograms.end(), rhs.x1_pre_histograms.begin(), rhs.x1_pre_histograms.end());
			x2_pre_histograms.insert(x2_pre_histograms.end(), rhs.x2_pre_histograms.begin(), rhs.x2_pre_histograms.end());
			x1_post_histograms.insert(x1_post_histograms.end(), rhs.x1_post_histograms.begin(), rhs.x1_post_histograms.end());
			x2_post_histograms.insert(x2_post_histograms.end(), rhs.x2_post_histograms.begin(), rhs.x2_post_histograms.end());

			N_trigger += rhs.N_trigger;
			N_assoc += rhs.N_assoc;
			N_pair += rhs.N_pair;

			partial_N_trigger += rhs.partial_N_trigger;
			partial_N_assoc += rhs.partial_N_assoc;
			partial_N_pair += rhs.partial_N_pair;

			sigma_sps_1.insert(sigma_sps_1.end(), rhs.sigma_sps_1.begin(), rhs.sigma_sps_1.end());
			sigma_sps_2.insert(sigma_sps_2.end(), rhs.sigma_sps_2.begin(), rhs.sigma_sps_2.end());
			sigma_sps.insert(sigma_sps.end(), rhs.sigma_sps.begin(), rhs.sigma_sps.end());
			sigma_gen += rhs.sigma_gen;
			total_weight += rhs.total_weight;
			next_calls += rhs.next_calls;
			next_skips += rhs.next_skips;
			partial_events += rhs.partial_events;
			partial_ncolls += rhs.partial_ncolls;

			return *this;
		}
		/// Combines a list of Result objects, each corresponding to a pT_hat range.
		static std::vector<Result> combine(std::vector<std::vector<Result>> input) {
			std::vector<Result> results;
			for (std::vector<Result>::size_type i = 0; i < input[0].size(); i++) {
				Result result = input[0][i];
				for (std::vector<std::vector<Result>>::size_type j = 1; j < input.size(); j++) {
					result += input[j][i];
				}
				results.push_back(result);
			}
			return results;
		}
	};
	/// Parameters for the particle generator.
	GeneratorParameters params;
	/// The bins of the delta_phi histogram.
	std::vector<double> bins;
	/// The pT_hat range, defining the phase space cuts.
	OptionalRange<double> pT_hat_range;
	/// A list of run parameters.
	std::vector<Analyzer::Parameters> runs;
	/// A list of Analyzer instances.
	std::vector<Analyzer> analyzers;
	/// Initalizes an event generator with the given parameters, histogram bins,
	/// phase space cuts and run parameters.
	EventGenerator(
		GeneratorParameters _params, 
		std::vector<double> _bins, 
		OptionalRange<double> _pT_hat_range, 
		std::vector<Analyzer::Parameters> _runs)
		:params(_params), bins(_bins), pT_hat_range(_pT_hat_range), runs(_runs) {
		for (auto &run : runs) {
			analyzers.emplace_back(run, bins);
		}
	}
	/// Generate and analyze events and return a list of Results.
	std::vector<Result> run() {
		/// Initialize the particle generator.
		ParticleGenerator generator(params, pT_hat_range);
		generator.initialize();
		/// Run the particle generator.
		generator.generate([this](std::vector<ParticleContainer> particles, bool is_partial) {
			/// For each event, send the generated particles to all analyzers.
			for (auto &analyzer : analyzers) {
				analyzer.book(&particles/*, &info*/, is_partial);
			}
		}, [this](ParticleGeneratorInfo info) {
			for (auto &analyzer : analyzers) {
				analyzer.next_calls = info.next_calls;
				analyzer.next_skips = info.next_skips;
				analyzer.partial_events = info.partial_events;
				analyzer.partial_ncolls = info.partial_ncolls;
			}
		});
		/// Get the total cross section and weight.
		const double sigma_gen = generator.sigma();
		const double total_weight = generator.total_weight();
		/// Initalize the result list.
		std::vector<Result> result_vector;
		/// Loop over each analyzer.
		for (std::vector<Analyzer>::size_type i = 0; i < analyzers.size(); i++) {
			/// For each analyzer, construct a Result.
			auto analyzer = analyzers[i];
			const auto factor = sigma_gen / total_weight;

			Result result;
			
			result.histograms = {analyzer.histogram};
			result.nuclear_histograms = {analyzer.nuclear_histogram};
			result.partial_histograms = {analyzer.partial_histogram};

			result.x1_pre_histograms = {analyzer.x1_pre_histogram};
			result.x2_pre_histograms = {analyzer.x2_pre_histogram};
			result.x1_post_histograms = {analyzer.x1_post_histogram};
			result.x2_post_histograms = {analyzer.x2_post_histogram};

			result.N_trigger = analyzer.N_trigger;
			result.N_assoc = analyzer.N_assoc;
			result.N_pair = analyzer.N_pair;

			result.partial_N_trigger = analyzer.partial_N_trigger;
			result.partial_N_assoc = analyzer.partial_N_assoc;
			result.partial_N_pair = analyzer.partial_N_pair;
			
			result.sigma_sps_1 = {analyzer.N_trigger * factor};
			result.sigma_sps_2 = {analyzer.N_assoc * factor};
			result.sigma_sps = {analyzer.N_pair * factor};
			
			result.sigma_gen = sigma_gen;
			result.total_weight = total_weight;
			
			result.parameters = analyzer.parameters;
			result.run_index = i;

			result.next_calls = analyzer.next_calls;
			result.next_skips = analyzer.next_skips;
			result.partial_events = analyzer.partial_events;
			result.partial_ncolls = analyzer.partial_ncolls;
			
			result_vector.push_back(result);
		}

		return result_vector;
	}	
};

#endif // EVENT_GENERATOR_H
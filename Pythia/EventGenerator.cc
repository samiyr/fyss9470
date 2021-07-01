#ifndef EVENT_GENERATOR_H
#define EVENT_GENERATOR_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

class EventGenerator {
public:
	struct Result {
		ValueHistogram<double> histogram;
		double sigma_sps_1;
		double sigma_sps_2;
		double sigma_sps;

		double sigma_gen;
		double total_weight;

		AnalysisParameters parameters;
		std::vector<Analyzer>::size_type run_index;

		Result& operator+=(Result rhs) {
			histogram += rhs.histogram;
			sigma_sps_1 += rhs.sigma_sps_1;
			sigma_sps_2 += rhs.sigma_sps_2;
			sigma_sps += rhs.sigma_sps;
			sigma_gen += rhs.sigma_gen;
			total_weight += rhs.total_weight;
			return *this;
		}

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
	GeneratorParameters params;
	std::vector<double> bins;
	OptionalRange<double> pT_hat_range;
	std::vector<AnalysisParameters> runs;
	std::vector<Analyzer> analyzers;
	Pythia pythia;

	EventGenerator(GeneratorParameters _params, std::vector<double> _bins, OptionalRange<double> _pT_hat_range, std::vector<AnalysisParameters> _runs) {
		params = _params;
		bins = _bins;
		pT_hat_range = _pT_hat_range;
		runs = _runs;
		for (auto &run : runs) {
			analyzers.emplace_back(run, bins);
		}
		initialize();
	}
	void initialize() {
		Settings &settings = pythia.settings;

		pythia.readFile(Constants::cmnd_input);
		settings.parm("Beams:eCM", params.cm_energy);
		settings.parm("PhaseSpace:pTHatMin", pT_hat_range.start.has_value() ? *pT_hat_range.start : 0);
		settings.parm("PhaseSpace:pTHatMax", pT_hat_range.end.has_value() ? *pT_hat_range.end : -1);
        settings.flag("PhaseSpace:bias2Selection", params.use_biasing);
        settings.parm("PhaseSpace:bias2SelectionPow", params.bias_power);
        settings.parm("PhaseSpace:bias2SelectionRef", params.bias_reference);
        settings.flag("Print:quiet", !params.pythia_printing);
        settings.mode("Next:numberCount", Defaults::pythia_next);
        settings.mode("Random:seed", params.random_seed);
        settings.flag("PartonLevel:MPI", params.mpi);

        params.beam_A.apply_to(settings, "A");
        params.beam_B.apply_to(settings, "B");

		pythia.init();
	}

	std::vector<Result> run() {
		ParticleGenerator generator(params);
		generator.initialize();
		generator.generate([this](std::vector<ParticleContainer> particles) {			
			for (auto &analyzer : analyzers) {
				analyzer.book(&particles);
			}
		});

		const double sigma_gen = generator.sigma();
		const double total_weight = generator.total_weight();

		std::vector<Result> result_vector;

		for (std::vector<Analyzer>::size_type i = 0; i < analyzers.size(); i++) {
			auto analyzer = analyzers[i];
			const auto factor = sigma_gen / total_weight;
			Result result;
			result.histogram = analyzer.histogram;
			result.sigma_sps_1 = analyzer.sigma_sps_1 * factor;
			result.sigma_sps_2 = analyzer.sigma_sps_2 * factor;
			result.sigma_sps = analyzer.sigma_sps * factor;
			result.sigma_gen = sigma_gen;
			result.total_weight = total_weight;
			result.parameters = analyzer.parameters;
			result.run_index = i;
			result_vector.push_back(result);
		}

		return result_vector;
	}	
};

#endif // EVENT_GENERATOR_H
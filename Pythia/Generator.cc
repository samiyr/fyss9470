#ifndef GENERATOR_H
#define GENERATOR_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

struct Result {
	ValueHistogram<double> histogram;
	double sigma_sps_1;
	double sigma_sps_2;
	double sigma_sps;
};

class Generator {
public:
	GeneratorParameters params;
	std::vector<double> bins;
	OptionalRange<double> pT_hat_range;
	std::vector<AnalysisParameters> runs;
	std::vector<Analyzer> analyzers;
	Pythia pythia;

	Generator(GeneratorParameters _params, std::vector<double> _bins, OptionalRange<double> _pT_hat_range, std::vector<AnalysisParameters> _runs) {
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

		pythia.init();
	}

	std::vector<Result> run() {
		ParticleGenerator generator(params);
		generator.initialize();
		generator.generate([this](std::vector<ParticleContainer> particles, [[maybe_unused]] bool last_event) {			
			for (auto &analyzer : analyzers) {
				analyzer.book(&particles);
			}
		});

		std::vector<Result> result_vector;

		for (auto &analyzer : analyzers) {
			Result result;
			result.histogram = analyzer.histogram;
			result.sigma_sps_1 = analyzer.sigma_sps_1;
			result.sigma_sps_2 = analyzer.sigma_sps_2;
			result.sigma_sps = analyzer.sigma_sps;
			result_vector.push_back(result);
		}

		return result_vector;
	}	
};

#endif // GENERATOR_H
#ifndef PARTICLE_GENERATOR_H
#define PARTICLE_GENERATOR_H

#include "Pythia8/Pythia.h"
#include "ParticleFilter.cc"
#include "Constants.cc"
#include "GeneratorParameters.cc"

using namespace Pythia8;

class ParticleGenerator {
public:
	GeneratorParameters params;
	OptionalRange<double> pT_hat_range;
	Pythia *pythia;

	ParticleGenerator(GeneratorParameters p) : params(p) {}
	
	void initialize() {
		pythia = new Pythia("../share/Pythia8/xmldoc", params.pythia_printing);
		Settings &settings = pythia->settings;

		pythia->readFile(Constants::cmnd_input);
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

        switch(params.process) {
        	case Process::HardQCD:
        		settings.flag("HardQCD:all", true);
        		break;
        	case Process::SoftQCDNonDiffractive:
        		settings.flag("SoftQCD:nonDiffractive", true);
        		break;
        }

        params.beam_A.apply_to(settings, "A");
        params.beam_B.apply_to(settings, "B");

		pythia->init();
	}

	template <typename F>
	void generate(F lambda) {
		ParticleFilter filter;
		filter.allowed_particle_ids = params.particle_ids;
		filter.include_decayed = params.include_decayed;
		filter.pT_range = params.pT_range;
		filter.y_range = params.y_range;

		for (int i = 0; i < params.event_count; ++i) {
			if (!pythia->next()) {
				continue;
			}
			const Event &event = pythia->event;
			const Info &info = pythia->info;

			const int particle_count = event.size();
			std::vector<ParticleContainer> particles;
			particles.reserve(particle_count);

			for (int j = 0; j < particle_count; j++) {
				Particle particle = event[j];
				if (filter.is_allowed(particle)) {
					particles.emplace_back(particle, info.weight());
				}
			}
			lambda(particles);
		}
	}

	std::vector<std::vector<ParticleContainer>> generate() {
		std::vector<std::vector<ParticleContainer>> particles(params.event_count, std::vector<ParticleContainer>());

		generate([&particles](std::vector<ParticleContainer> generated) {
			particles.push_back(generated);
		});

		return particles;
	}
	double sigma() const {
		return pythia->info.sigmaGen();
	}
	double total_weight() const {
		return pythia->info.weightSum();
	}
};

#endif // PARTICLE_GENERATOR_H
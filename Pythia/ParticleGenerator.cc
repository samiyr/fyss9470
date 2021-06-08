#ifndef PARTICLE_GENERATOR_H
#define PARTICLE_GENERATOR_H

#include "Pythia8/Pythia.h"
#include "ParticleFilter.cc"

using namespace Pythia8;

// Constants
const std::string cmnd_input = "pp.cmnd";

class ParticleGenerator {
public:
	double cm_energy;
	int event_count;
	bool pythia_printing = true;

	bool include_decayed = true;

	OptionalRange<double> pT_range;
	OptionalRange<double> y_range;
	OptionalRange<double> pT_hat_range;

	bool use_biasing = false;
	double bias_power = 4.0;
	double bias_reference = 10.0;

	std::vector<int> particle_ids = {111};

	Pythia pythia;

	ParticleGenerator(double energy, int count) {
		cm_energy = energy;
		event_count = count;
	}
	void initialize() {
		Settings &settings = pythia.settings;

		pythia.readFile(cmnd_input);
		settings.parm("Beams:eCM", cm_energy);
		settings.parm("PhaseSpace:pTHatMin", pT_hat_range.start.has_value() ? *pT_hat_range.start : 0);
		settings.parm("PhaseSpace:pTHatMax", pT_hat_range.end.has_value() ? *pT_hat_range.end : -1);
        settings.flag("PhaseSpace:bias2Selection", use_biasing);
        settings.parm("PhaseSpace:bias2SelectionPow", bias_power);
        settings.parm("PhaseSpace:bias2SelectionRef", bias_reference);
        settings.flag("Print:quiet", !pythia_printing);

		pythia.init();
	}
	std::vector<ParticleContainer> generate() {
		std::vector<ParticleContainer> particles;

		ParticleFilter filter;
		filter.allowed_particle_ids = particle_ids;
		filter.include_decayed = include_decayed;
		filter.pT_range = pT_range;
		filter.y_range = y_range;

		for (int i = 0; i < event_count; ++i) {
			if (!pythia.next()) {
				continue;
			}
			const Event &event = pythia.event;
			const Info &info = pythia.info;
			const double pT_hat = info.pTHat();

			const int particle_count = event.size();

			for (int j = 0; j < particle_count; j++) {
				Particle particle = event[j];
				if (filter.is_allowed(particle)) {
					particles.emplace_back(particle, pT_hat, info.weight());
				}
			}
		}	
		return particles;
	}
	double sigma() const {
		return pythia.info.sigmaGen();
	}
	double total_weight() const {
		return pythia.info.weightSum();
	}
};

#endif // PARTICLE_GENERATOR_H
#ifndef PARTICLE_GENERATOR_H
#define PARTICLE_GENERATOR_H

#include "Pythia8/Pythia.h"
#include <string>
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

	double pT_min = NAN;
	double pT_max = NAN;
	
	double y_min = NAN;
	double y_max = NAN;

	double pT_hat_min = 0;
	double pT_hat_max = -1;

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
		settings.parm("PhaseSpace:pTHatMin", pT_hat_min);
		settings.parm("PhaseSpace:pTHatMax", pT_hat_max);
		if (!pythia_printing) {
			pythia.readString("Print:quiet = on");
		}
		pythia.init();
	}
	std::vector<Particle> generate() {
		std::vector<Event> events;

		for (int i = 0; i < event_count; ++i) {
			if (!pythia.next()) {
				continue;
			}

			if (i != 0) {
				events.push_back(pythia.event);
			}
		}	

		const std::vector<Particle> particles = all_particles(events);

		ParticleFilter filter;
		filter.allowed_particle_ids = particle_ids;
		filter.include_decayed = include_decayed;
		filter.pT_min = pT_min;
		filter.pT_max = pT_max;
		filter.y_min = y_min;
		filter.y_max = y_max;

		const std::vector<Particle> filtered = filter.filter(particles);

		return filtered;
	}
	double sigma() {
		return pythia.info.sigmaGen();
	}
};

#endif // PARTICLE_GENERATOR_H
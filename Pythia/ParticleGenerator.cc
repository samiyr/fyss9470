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

	Range<double> pT_range = Range<double>(NAN, NAN);
	Range<double> y_range = Range<double>(NAN, NAN);
	Range<double> pT_hat_range = Range<double>(0, -1);

	std::vector<int> particle_ids = {111};

	double pT_hat_total = 0;

	Pythia pythia;

	ParticleGenerator(double energy, int count) {
		cm_energy = energy;
		event_count = count;
	}
	void initialize() {
		Settings &settings = pythia.settings;

		pythia.readFile(cmnd_input);
		settings.parm("Beams:eCM", cm_energy);
		settings.parm("PhaseSpace:pTHatMin", pT_hat_range.start);
		settings.parm("PhaseSpace:pTHatMax", pT_hat_range.end);
		if (!pythia_printing) {
			pythia.readString("Print:quiet = on");
		}
		pythia.init();
	}
	std::vector<ParticleContainer> generate() {
		std::vector<Event> events;
		events.reserve(event_count);
		std::vector<ParticleContainer> particles;

		for (int i = 0; i < event_count; ++i) {
			if (!pythia.next()) {
				continue;
			}

			if (i != 0) {
				events.push_back(pythia.event);
				const double pT_hat = pythia.info.pTHat();
				pT_hat_total += pT_hat;

				const int particle_count = pythia.event.size();

				for (int j = 0; j < particle_count; j++) {
					const ParticleContainer container = ParticleContainer(pythia.event[j], pT_hat);
					particles.push_back(container);
				}
			}
		}	

		ParticleFilter filter;
		filter.allowed_particle_ids = particle_ids;
		filter.include_decayed = include_decayed;
		filter.pT_range = pT_range;
		filter.y_range = y_range;

		const std::vector<ParticleContainer> filtered = filter.filter(particles);

		return filtered;
	}
	double sigma() {
		return pythia.info.sigmaGen();
	}
};

#endif // PARTICLE_GENERATOR_H
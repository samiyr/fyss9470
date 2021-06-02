#include "Pythia8/Pythia.h"
#include <string>
#include "ParticleFilter.cc"

using namespace Pythia8;

// Constants
const std::string cmnd_input = "pp.cmnd";

class PionGenerator {
public:
	double cm_energy;
	int event_count;
	bool pythia_printing = true;

	bool pi_pos = false;
	bool pi_neg = false;
	bool pi_0 = true;

	bool include_decayed = true;

	double pT_min = NAN;
	double pT_max = NAN;
	
	double y_min = NAN;
	double y_max = NAN;

	std::vector<int> pion_ids;

	Pythia pythia;

	PionGenerator(double energy, int count) {
		cm_energy = energy;
		event_count = count;
	}
	void initialize() {
		if (pi_pos) {
			pion_ids.push_back(211);
		}
		if (pi_neg) {
			pion_ids.push_back(-211);
		}
		if (pi_0) {
			pion_ids.push_back(111);
		}

		pythia.readFile(cmnd_input);
		pythia.readString("Beams:eCM = " + std::to_string(cm_energy));
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
		filter.allowed_particle_ids = pion_ids;
		filter.include_decayed = include_decayed;
		filter.pT_min = pT_min;
		filter.pT_max = pT_max;
		filter.y_min = y_min;
		filter.y_max = y_max;

		const std::vector<Particle> pions = filter.filter(particles);

		return pions;
	}
	double sigma() {
		return pythia.info.sigmaGen();
	}
};


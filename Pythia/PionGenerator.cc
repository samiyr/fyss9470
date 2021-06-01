#include "Pythia8/Pythia.h"
#include <string>
#include "Helpers.cc"

using namespace Pythia8;

// Constants
const std::string cmnd_input = "pp.cmnd";
const std::vector<int> pion_ids {-211, 211, 111};

class PionGenerator {
public:
	double cm_energy;
	int event_count;
	bool pythia_printing;
	Pythia pythia;

	PionGenerator(double energy, int count, bool printing = true) {
		cm_energy = energy;
		event_count = count;
		pythia_printing = printing;
	}
	void initialize() {
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
		const std::vector<Particle> pions = filter_particles(particles, pion_ids);

		return pions;
	}
};


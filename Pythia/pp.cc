#include "Pythia8/Pythia.h"
#include <string>
#include "PionGenerator.cc"

using namespace Pythia8;

void create_histogram(int energy, int count) {
	Pythia pythia;
	pythia.readFile(cmnd_input);
	pythia.readString("Beams::eCM = " + std::to_string(energy));
	pythia.init();

	std::vector<Event> events;
	Hist azimuth_hist("azimuth angles", 50, -4, 4);

	for (int i = 0; i < count; ++i) {
		if (!pythia.next()) {
			continue;
		}

		events.push_back(pythia.event);
	}
	const std::vector<Particle> particles = all_particles(events);
	const std::vector<Particle> pions = filter_particles(particles, pion_ids);
	const std::vector<double> azimuths = find_azimuths(pions);

	for (double azimuth : azimuths) {
		azimuth_hist.fill(azimuth);
	}

	cout << azimuth_hist << "\n";

	const std::string python_plot_file = "azimuth_angles";
	HistPlot hpl(python_plot_file);
	hpl.frame("azimuth_angles", "Azimuth angles", "$\\phi$", "N");
	hpl.add(azimuth_hist);
	hpl.plot(); 

	pythia.stat();
}

int pion_production(double energy, int count) {
	PionGenerator generator(energy, count);
	const std::vector<Particle> pions = generator.generate();

	return pions.size();
}

void test_pion_production() {
	std::vector<std::pair<int, int>> pion_counts;
	for (int energy = 1; energy <= 10000; energy *= 10) {
		const int pions = pion_production(energy, 100);
		pion_counts.push_back(std::make_pair(pions, energy));
	}
	for (std::pair<int, int> e : pion_counts) {
		cout << "pions at " << e.second << " GeV: " << e.first << "\n";
	}
}

int main() {
	// create_histogram(8000, 100);

	test_pion_production();

	return 0;
}


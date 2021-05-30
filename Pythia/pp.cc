#include "Pythia8/Pythia.h"
#include <algorithm>
#include <string>
#include <utility>

using namespace Pythia8;

// Constants
const std::string input = "pp.cmnd";
const std::vector<int> pion_ids {-211, 211, 111};

template <typename T>
bool contains(std::vector<T> vec, T element) {
	return std::find(vec.begin(), vec.end(), element) != vec.end();
}

std::vector<Particle> all_particles(std::vector<Event> events) {
	std::vector<Particle> particles;
	for (Event event : events) {
		const int particle_count = event.size();
		for (int i = 0; i < particle_count; ++i) {
			const Particle particle = event[i];
			particles.push_back(particle);
		}
	}
	return particles;
}

std::vector<Particle> filter_particles(std::vector<Particle> particles, std::vector<int> types, bool keep_decayed = false) {
	particles.erase(std::remove_if(particles.begin(), particles.end(), [types, keep_decayed](Particle particle) {
		if (!contains(types, particle.id())) {
			return true;
		}
		if (!particle.isFinal() && !keep_decayed) {
			return true;
		}
		return false;
	}), particles.end());
	return particles;
}

std::vector<double> find_azimuths(std::vector<Particle> particles) {
	std::vector<double> phis;
	for (Particle particle : particles) {
		phis.push_back(particle.phi());
	}
	return phis;
}

void create_histogram(int energy, int count) {
	Pythia pythia;
	pythia.readFile(input);
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
	Pythia pythia;
	pythia.readFile(input);
	pythia.readString("Beams::eCM = " + std::to_string(energy));
	pythia.init();

	std::vector<Event> events;

	for (int i = 0; i < count; ++i) {
		if (!pythia.next()) {
			continue;
		}

		if (i != 0) {
			events.push_back(pythia.event);
		}
	}

	const std::vector<Particle> particles = all_particles(events);
	const std::vector<Particle> pions = filter_particles(particles, pion_ids);

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


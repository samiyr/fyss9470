#include "Pythia8/Pythia.h"
#include <algorithm>
#include <string>
#include <utility>

using namespace Pythia8;

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

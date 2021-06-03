#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include "Pythia8/Pythia.h"
#include <string>
#include "Helpers.cc"
#include <functional>
#include <math.h>

using namespace Pythia8;

class ParticleFilter {
private:
	std::vector<Particle> apply_filter(std::vector<Particle> particles, bool (ParticleFilter::*filter)(Particle)) {
		particles.erase(std::remove_if(particles.begin(), particles.end(), [this, filter](Particle particle) {
			const bool value = (this->*filter)(particle);
			return !value;
		}), particles.end());
		return particles;
	}
	std::vector<Particle> apply_filter_chain(std::vector<Particle> particles, std::vector<bool (ParticleFilter::*)(Particle)> chain) {
		std::vector<Particle> current = particles;
		for (bool (ParticleFilter::*filter)(Particle) : chain) {
			current = apply_filter(current, filter);
		}
		return current;
	}

	bool id_filter(Particle p) {
		return contains(allowed_particle_ids, p.id());
	}
	bool decay_filter(Particle p) {
		if (!include_decayed) {
			return p.isFinal();
		}
		return true;
	}
	bool pT_filter(Particle p) {
		return in_range(p.pT(), pT_min, pT_max);
	}
	bool rapidity_filter(Particle p) {
		return in_range(p.y(), y_min, y_max);
	}
	bool in_range(double value, double min, double max) {
		if (!std::isnan(min) && value < min) {
			return false;
		}
		if (!std::isnan(max) && value > max) {
			return false;
		}
		return true;
	}

public:
	std::vector<int> allowed_particle_ids;
	bool include_decayed;

	double pT_min;
	double pT_max;
	
	double y_min;
	double y_max;

	std::vector<bool (ParticleFilter::*)(Particle)> filters = {
		&ParticleFilter::id_filter,
		&ParticleFilter::decay_filter,
		&ParticleFilter::pT_filter,
		&ParticleFilter::rapidity_filter,
	};
			
	std::vector<Particle> filter(const std::vector<Particle> particles) {
		filters.push_back(&ParticleFilter::id_filter);

		return apply_filter_chain(particles, filters);
	}
};

#endif // PARTICLE_FILTER_H
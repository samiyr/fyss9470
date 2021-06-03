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
	std::vector<ParticleContainer> apply_filter(std::vector<ParticleContainer> particles, bool (ParticleFilter::*filter)(ParticleContainer)) {
		particles.erase(std::remove_if(particles.begin(), particles.end(), [this, filter](ParticleContainer particle) {
			const bool value = (this->*filter)(particle);
			return !value;
		}), particles.end());
		return particles;
	}
	std::vector<ParticleContainer> apply_filter_chain(std::vector<ParticleContainer> particles, std::vector<bool (ParticleFilter::*)(ParticleContainer)> chain) {
		std::vector<ParticleContainer> current = particles;
		for (bool (ParticleFilter::*filter)(ParticleContainer) : chain) {
			current = apply_filter(current, filter);
		}
		return current;
	}

	bool id_filter(ParticleContainer p) {
		return contains(allowed_particle_ids, p.particle.id());
	}
	bool decay_filter(ParticleContainer p) {
		if (!include_decayed) {
			return p.particle.isFinal();
		}
		return true;
	}
	bool pT_filter(ParticleContainer p) {
		return in_range(p.particle.pT(), pT_min, pT_max);
	}
	bool rapidity_filter(ParticleContainer p) {
		return in_range(p.particle.y(), y_min, y_max);
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

	std::vector<bool (ParticleFilter::*)(ParticleContainer)> filters = {
		&ParticleFilter::id_filter,
		&ParticleFilter::decay_filter,
		&ParticleFilter::pT_filter,
		&ParticleFilter::rapidity_filter,
	};
			
	std::vector<ParticleContainer> filter(const std::vector<ParticleContainer> particles) {
		filters.push_back(&ParticleFilter::id_filter);

		return apply_filter_chain(particles, filters);
	}
};

#endif // PARTICLE_FILTER_H
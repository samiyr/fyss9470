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
	bool check_filter(ParticleContainer particle, bool (ParticleFilter::*filter)(ParticleContainer)) {
		return (this->*filter)(particle);
	}
	std::vector<ParticleContainer> apply_filter(std::vector<ParticleContainer> particles, bool (ParticleFilter::*filter)(ParticleContainer)) {
		particles.erase(std::remove_if(particles.begin(), particles.end(), [this, filter](ParticleContainer particle) {
			return !check_filter(particle, filter);
		}), particles.end());
		return particles;
	}
	bool check_filter_chain(ParticleContainer particle, std::vector<bool (ParticleFilter::*)(ParticleContainer)> chain) {
		for (bool (ParticleFilter::*filter)(ParticleContainer) : chain) {
			if(!check_filter(particle, filter)) {
				return false;
			}
		}
		return true;
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
		return pT_range.in_range(p.particle.pT());
	}
	bool rapidity_filter(ParticleContainer p) {
		return y_range.in_range(p.particle.y());
	}

public:
	std::vector<int> allowed_particle_ids;
	bool include_decayed;

	OptionalRange<double> pT_range = OptionalRange<double>();
	OptionalRange<double> y_range = OptionalRange<double>();

	std::vector<bool (ParticleFilter::*)(ParticleContainer)> filters = {
		&ParticleFilter::id_filter,
		&ParticleFilter::decay_filter,
		&ParticleFilter::pT_filter,
		&ParticleFilter::rapidity_filter,
	};
			
	std::vector<ParticleContainer> filter(const std::vector<ParticleContainer> particles) {
		return apply_filter_chain(particles, filters);
	}
	bool is_allowed(const ParticleContainer particle) {
		return check_filter_chain(particle, filters);
	}
};

#endif // PARTICLE_FILTER_H
#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include "Pythia8/Pythia.h"
#include "Helpers.cc"

using namespace Pythia8;

class ParticleFilter {
private:
	bool check_filter(ParticleContainer particle, bool (ParticleFilter::*f)(ParticleContainer)) {
		return (this->*f)(particle);
	}
	// std::vector<Particle> apply_filter(std::vector<Particle> particles, bool (ParticleFilter::*f)(Particle)) {
	// 	particles.erase(std::remove_if(particles.begin(), particles.end(), [this, f](Particle particle) {
	// 		return !check_filter(particle, f);
	// 	}), particles.end());
	// 	return particles;
	// }
	bool check_filter_chain(ParticleContainer particle, std::vector<bool (ParticleFilter::*)(ParticleContainer)> chain) {
		for (bool (ParticleFilter::*f)(ParticleContainer) : chain) {
			if (!check_filter(particle, f)) {
				return false;
			}
		}
		return true;
	}

	// std::vector<Particle> apply_filter_chain(std::vector<Particle> particles, std::vector<bool (ParticleFilter::*)(Particle)> chain) {
	// 	std::vector<Particle> current = particles;
	// 	for (bool (ParticleFilter::*f)(Particle) : chain) {
	// 		current = apply_filter(current, f);
	// 	}
	// 	return current;
	// }
	bool id_filter(ParticleContainer p) {
		return contains(&allowed_particle_ids, p.id);
	}
	bool decay_filter(ParticleContainer p) {
		if (!include_decayed) {
			return p.is_final;
		}
		return true;
	}
	bool pT_filter(ParticleContainer p) {
		return pT_range.in_range(p.pT);
	}
	bool rapidity_filter(ParticleContainer p) {
		return y_range.in_range(p.y);
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
			
	bool is_allowed(const ParticleContainer particle) {
		return check_filter_chain(particle, filters);
	}
	bool is_allowed(const Particle particle) {
		return check_filter_chain(ParticleContainer(particle, 0.0), filters);
	}
};

#endif // PARTICLE_FILTER_H
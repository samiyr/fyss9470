#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include "Pythia8/Pythia.h"
#include "Helpers.cc"

using namespace Pythia8;

class ParticleFilter {
private:
	bool check_filter(Particle particle, bool (ParticleFilter::*f)(Particle)) {
		return (this->*f)(particle);
	}
	bool check_filter(ParticleContainer container, bool (ParticleFilter::*f)(Particle)) {
		return check_filter(container.particle, f);
	}
	std::vector<Particle> apply_filter(std::vector<Particle> particles, bool (ParticleFilter::*f)(Particle)) {
		particles.erase(std::remove_if(particles.begin(), particles.end(), [this, f](Particle particle) {
			return !check_filter(particle, f);
		}), particles.end());
		return particles;
	}
	std::vector<ParticleContainer> apply_filter(std::vector<ParticleContainer> particles, bool (ParticleFilter::*f)(Particle)) {
		particles.erase(std::remove_if(particles.begin(), particles.end(), [this, f](ParticleContainer particle) {
			return !check_filter(particle, f);
		}), particles.end());
		return particles;
	}
	bool check_filter_chain(Particle particle, std::vector<bool (ParticleFilter::*)(Particle)> chain) {
		for (bool (ParticleFilter::*f)(Particle) : chain) {
			if (!check_filter(particle, f)) {
				return false;
			}
		}
		return true;
	}
	bool check_filter_chain(ParticleContainer container, std::vector<bool (ParticleFilter::*)(Particle)> chain) {
		return check_filter_chain(container.particle, chain);
	}

	std::vector<Particle> apply_filter_chain(std::vector<Particle> particles, std::vector<bool (ParticleFilter::*)(Particle)> chain) {
		std::vector<Particle> current = particles;
		for (bool (ParticleFilter::*f)(Particle) : chain) {
			current = apply_filter(current, f);
		}
		return current;
	}
	std::vector<ParticleContainer> apply_filter_chain(std::vector<ParticleContainer> particles, std::vector<bool (ParticleFilter::*)(Particle)> chain) {
		std::vector<ParticleContainer> current = particles;
		for (bool (ParticleFilter::*f)(Particle) : chain) {
			current = apply_filter(current, f);
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
		return pT_range.in_range(p.pT());
	}
	bool rapidity_filter(Particle p) {
		return y_range.in_range(p.y());
	}

public:
	std::vector<int> allowed_particle_ids;
	bool include_decayed;

	OptionalRange<double> pT_range = OptionalRange<double>();
	OptionalRange<double> y_range = OptionalRange<double>();

	std::vector<bool (ParticleFilter::*)(Particle)> filters = {
		&ParticleFilter::id_filter,
		&ParticleFilter::decay_filter,
		&ParticleFilter::pT_filter,
		&ParticleFilter::rapidity_filter,
	};
			
	std::vector<ParticleContainer> filter(const std::vector<ParticleContainer> particles) {
		return apply_filter_chain(particles, filters);
	}
	bool is_allowed(const Particle particle) {
		return check_filter_chain(particle, filters);
	}
};

#endif // PARTICLE_FILTER_H
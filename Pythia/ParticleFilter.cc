#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include "Pythia8/Pythia.h"
#include "Helpers.cc"

using namespace Pythia8;

/// A general particle filtration system.
class ParticleFilter {
private:
	/// Applies the given filter function to a particle and returns the result of the filtration.
	bool check_filter(ParticleContainer particle, bool (ParticleFilter::*f)(ParticleContainer)) {
		return (this->*f)(particle);
	}
	/// Applies all the given filter functions to a particle and returns the result of the filtration. 
	/// Returns immediately if a filter is not passed.
	bool check_filter_chain(ParticleContainer particle, std::vector<bool (ParticleFilter::*)(ParticleContainer)> chain) {
		for (bool (ParticleFilter::*f)(ParticleContainer) : chain) {
			if (!check_filter(particle, f)) {
				return false;
			}
		}
		return true;
	}
	/// Checks if the particle container has the allowed particle ID.
	bool id_filter(ParticleContainer p) {
		return contains(&allowed_particle_ids, p.id);
	}
	/// Checks if the particle container has decayed and whether that's allowed.
	bool decay_filter(ParticleContainer p) {
		if (!include_decayed) {
			return p.is_final;
		}
		return true;
	}
	/// Checks if the particle's transverse momentum is within range.
	bool pT_filter(ParticleContainer p) {
		return pT_range.in_range(p.pT);
	}
	/// Checks if the particle's rapidity is within range.
	bool rapidity_filter(ParticleContainer p) {
		return y_range.in_range(p.y);
	}

public:
	/// A list of allowed particle IDs, as given by Pythia.
	std::vector<int> allowed_particle_ids;
	/// Whether to allow decayed particles.
	bool include_decayed;
	/// The transverse momentum acceptance range. Defaults to (-inf, inf).
	OptionalRange<double> pT_range = OptionalRange<double>();
	/// The rapidity acceptance range. Defaults to (-inf, inf),
	OptionalRange<double> y_range = OptionalRange<double>();
	/// A list of built-in filter functions.
	std::vector<bool (ParticleFilter::*)(ParticleContainer)> filters = {
		&ParticleFilter::id_filter,
		&ParticleFilter::decay_filter,
		&ParticleFilter::pT_filter,
		&ParticleFilter::rapidity_filter,
	};
	/// Checks whether a particle container passes all the filters.
	bool is_allowed(const ParticleContainer particle) {
		return check_filter_chain(particle, filters);
	}
	/// Checks whether a particle passes all the filters.
	bool is_allowed(const Particle particle) {
		return check_filter_chain(ParticleContainer(particle, 0.0), filters);
	}
};

#endif // PARTICLE_FILTER_H
#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

/**
 * Simple container for particles and associated properties not stored directly in Pythia8::Particle.
 */
struct ParticleContainer {
	Particle particle;
	double pT_hat;
	double event_weight;
	int event_id;

	ParticleContainer(Particle p, double h, double w, int e) {
		particle = p;
		pT_hat = h;
		event_weight = w;
		event_id = e;
	}
};

#endif // PARTICLE_CONTAINER_H
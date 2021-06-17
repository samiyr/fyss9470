#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

/**
 * Simple container for particles and associated properties not stored directly in Pythia8::Particle.
 */
struct ParticleContainer {
	double pT;
	double y;
	double phi;
	double event_weight;

	ParticleContainer(Particle p, double w) {
		pT = p.pT();
		y = p.y();
		phi = p.phi();
		event_weight = w;
	}
};

#endif // PARTICLE_CONTAINER_H
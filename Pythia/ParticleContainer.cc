#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

struct ParticleContainer {
	Particle particle;
	double pT_hat;

	ParticleContainer(Particle p, double h) {
		particle = p;
		pT_hat = h;
	}
};

#endif // PARTICLE_CONTAINER_H
#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

/**
 * Simple container for particles and associated properties not stored directly in Pythia8::Particle.
 */
struct ParticleContainer {
	/// Transverse momentum.
	double pT;
	/// Rapidity.
	double y;
	/// Azimuthal angle.
	double phi;
	/// Associated event weight. If biasing is disabled, equals 1.
	double event_weight;
	/// Particle ID.
	int id;
	/// True if the particle hasn't decayed.
	bool is_final;
	/// Constructs a particle container from the Pythia Particle class and with the event weight.
	ParticleContainer(Particle p, double w) 
	: pT(p.pT()), y(p.y()), phi(p.phi()), event_weight(w), id(p.id()), is_final(p.isFinal()) {}
};

#endif // PARTICLE_CONTAINER_H
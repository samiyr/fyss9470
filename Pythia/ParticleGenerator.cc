#ifndef PARTICLE_GENERATOR_H
#define PARTICLE_GENERATOR_H

#include "Pythia8/Pythia.h"
#include "ParticleFilter.cc"
#include "Constants.cc"

using namespace Pythia8;

class ParticleGenerator {
public:
	double cm_energy;
	int event_count;
	bool pythia_printing;

	bool include_decayed;

	bool mpi;

	OptionalRange<double> pT_range;
	OptionalRange<double> y_range;
	OptionalRange<double> pT_hat_range;

	bool use_biasing;
	double bias_power;
	double bias_reference;

	int random_seed;

	std::vector<int> particle_ids;

	Pythia pythia;

	ParticleGenerator(double energy, int count) {
		cm_energy = energy;
		event_count = count;
	}
	void initialize() {
		Settings &settings = pythia.settings;

		pythia.readFile(Constants::cmnd_input);
		settings.parm("Beams:eCM", cm_energy);
		settings.parm("PhaseSpace:pTHatMin", pT_hat_range.start.has_value() ? *pT_hat_range.start : 0);
		settings.parm("PhaseSpace:pTHatMax", pT_hat_range.end.has_value() ? *pT_hat_range.end : -1);
        settings.flag("PhaseSpace:bias2Selection", use_biasing);
        settings.parm("PhaseSpace:bias2SelectionPow", bias_power);
        settings.parm("PhaseSpace:bias2SelectionRef", bias_reference);
        settings.flag("Print:quiet", !pythia_printing);
        settings.mode("Next:numberCount", Defaults::pythia_next);
        settings.mode("Random:seed", random_seed);
        settings.flag("PartonLevel:MPI", mpi);

		pythia.init();
	}

	template <typename F>
	void generate(F lambda) {
		ParticleFilter filter;
		filter.allowed_particle_ids = particle_ids;
		filter.include_decayed = include_decayed;
		filter.pT_range = pT_range;
		filter.y_range = y_range;

		for (int i = 0; i < event_count; ++i) {
			if (!pythia.next()) {
				continue;
			}
			const Event &event = pythia.event;
			const Info &info = pythia.info;

			const int particle_count = event.size();
			std::vector<ParticleContainer> particles;
			particles.reserve(particle_count);

			for (int j = 0; j < particle_count; j++) {
				Particle particle = event[j];
				if (filter.is_allowed(particle)) {
					particles.emplace_back(particle, info.weight());
				}
			}

			lambda(particles);
		}
	}

	std::vector<std::vector<ParticleContainer>> generate() {
		std::vector<std::vector<ParticleContainer>> particles(event_count, std::vector<ParticleContainer>());

		generate([&particles](std::vector<ParticleContainer> generated) {
			particles.push_back(generated);
		});

		return particles;

		// ParticleFilter filter;
		// filter.allowed_particle_ids = particle_ids;
		// filter.include_decayed = include_decayed;
		// filter.pT_range = pT_range;
		// filter.y_range = y_range;

		// for (int i = 0; i < event_count; ++i) {
		// 	if (!pythia.next()) {
		// 		continue;
		// 	}
		// 	const Event &event = pythia.event;
		// 	const Info &info = pythia.info;

		// 	const int particle_count = event.size();

		// 	for (int j = 0; j < particle_count; j++) {
		// 		Particle particle = event[j];
		// 		if (filter.is_allowed(particle)) {
		// 			particles[i].emplace_back(particle, info.weight());
		// 		}
		// 	}
		// }
		// return particles;
	}
	double sigma() const {
		return pythia.info.sigmaGen();
	}
	double total_weight() const {
		return pythia.info.weightSum();
	}
};

#endif // PARTICLE_GENERATOR_H
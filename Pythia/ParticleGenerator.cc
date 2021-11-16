#ifndef PARTICLE_GENERATOR_H
#define PARTICLE_GENERATOR_H

#include "Pythia8/Pythia.h"
#include "ParticleFilter.cc"
#include "Constants.cc"
#include "GeneratorParameters.cc"

using namespace Pythia8;

struct ParticleGeneratorInfo {
	Pythia *pythia;
	EVENT_COUNT_TYPE next_calls;
	EVENT_COUNT_TYPE next_skips;
	EVENT_COUNT_TYPE rejections;

	ParticleGeneratorInfo(Pythia *_pythia, EVENT_COUNT_TYPE _next_calls, EVENT_COUNT_TYPE _next_skips, EVENT_COUNT_TYPE _rejections)
	:pythia(_pythia), next_calls(_next_calls), next_skips(_next_skips), rejections(_rejections) {}
};

/// The most basic generator, responsible for generating particles and interfacing with Pythia.
class ParticleGenerator {
public:
	/// Generator parameters.
	GeneratorParameters params;
	/// A list of phase space cuts.
	OptionalRange<double> pT_hat_range;
	/// A pointer to the Pythia instance associated with the generator.
	Pythia *pythia;
	/// Constructs a particle generator from the parameters and phase space cuts.
	ParticleGenerator(GeneratorParameters p, OptionalRange<double> _pT_hat_range) : params(p), pT_hat_range(_pT_hat_range) {}
	/// Initializes the generator and Pythia, and applies the settings.
	void initialize() {
		pythia = new Pythia("../share/Pythia8/xmldoc", params.pythia_printing);
		Settings &settings = pythia->settings;

		pythia->readFile(Constants::cmnd_input);
		settings.parm("Beams:eCM", params.cm_energy);
		settings.parm("PhaseSpace:pTHatMin", pT_hat_range.start.has_value() ? *pT_hat_range.start : 0);
		settings.parm("PhaseSpace:pTHatMax", pT_hat_range.end.has_value() ? *pT_hat_range.end : -1);
        settings.flag("PhaseSpace:bias2Selection", params.use_biasing);
        settings.parm("PhaseSpace:bias2SelectionPow", params.bias_power);
        settings.parm("PhaseSpace:bias2SelectionRef", params.bias_reference);
        settings.flag("Print:quiet", !params.pythia_printing);
        settings.mode("Next:numberCount", Defaults::pythia_next);
        settings.mode("Random:seed", params.random_seed);
        settings.flag("PartonLevel:MPI", params.mpi);

        switch(params.process) {
        	case Process::HardQCD:
        		settings.flag("HardQCD:all", true);
        		break;
        	case Process::SoftQCDNonDiffractive:
        		settings.flag("SoftQCD:nonDiffractive", true);
        		break;
        }

        params.beam_A.apply_to(settings, "A");
        params.beam_B.apply_to(settings, "B");

		pythia->init();
	}

	template <typename F, typename G>
	/// Start generating events in the event loop. The lambda
	/// expression is called with the particles generated from each event.
	void generate(F lambda, G completion) {
		ParticleFilter filter;
		filter.allowed_particle_ids = params.particle_ids;
		filter.include_decayed = params.include_decayed;
		filter.pT_range = params.pT_range;
		filter.y_range = params.y_range;

		EVENT_COUNT_TYPE next_calls = 0;
		EVENT_COUNT_TYPE next_skips = 0;
		EVENT_COUNT_TYPE rejections = 0;

		if (params.use_ncoll) {
			// Non-standard event generation loop to include ncoll effects for nuclear beams.
			EVENT_COUNT_TYPE retry_counter = 0;
			const auto ncoll_list = params.beam_B.ncoll_list();
			for (EVENT_COUNT_TYPE i = 0; i < params.event_count; ++i) {
				double x1_tot = 0.0;
				const int n_coll = ncoll_list.has_value() ? (*ncoll_list)(i) : 1;
				std::vector<ParticleContainer> particles;
				// Collision loop.
				for (int i_coll = 0; i_coll < n_coll; i_coll++) {
					next_calls++;
					if (!pythia->next()) {
						next_skips++;
						continue;
					}
					const Event &event = pythia->event;
					const Info &info = pythia->info;

					const int particle_count = event.size();

					for (int j = 0; j < particle_count; j++) {
						Particle particle = event[j];
						if (filter.is_allowed(particle)) {
							particles.emplace_back(particle, info.weight());
						}
					}
					const int n_MPI = info.nMPI();
					double x1_sum = 0.0;

					for (int i_MPI = 0; i_MPI < n_MPI; i_MPI++) {
						const double x1 = (*info.beamAPtr)[i_MPI].x();
						x1_sum += x1;
					}
					x1_tot += x1_sum;
				} // End of collision loop

				// Check conservation of momentum.
				if (x1_tot > 1.0) {
					retry_counter++;
					if (retry_counter > params.ncoll_retries) {
						rejections++;
					} else {
						i--;
					}
					continue;
				}
				retry_counter = 0;
				lambda(particles);
			}
			completion(ParticleGeneratorInfo(pythia, next_calls, next_skips, rejections));
		} else {
			// Standard event generation loop.
			for (EVENT_COUNT_TYPE i = 0; i < params.event_count; ++i) {
				next_calls++;
				if (!pythia->next()) {
					next_skips++;
					continue;
				}
				const Event &event = pythia->event;
				const Info &info = pythia->info;

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
			completion(ParticleGeneratorInfo(pythia, next_calls, next_skips, rejections));
		}
	}
	/// Generates, stores and returns a 2D matrix of all the particles generated in the event loop, event-by-event.
	std::vector<std::vector<ParticleContainer>> generate() {
		std::vector<std::vector<ParticleContainer>> particles(params.event_count, std::vector<ParticleContainer>());

		generate([&particles](std::vector<ParticleContainer> generated) {
			particles.push_back(generated);
		}, []([[maybe_unused]] ParticleGeneratorInfo info) {});

		return particles;
	}
	/// Returns the total cross section as given by Pythia.
	/// If called before event generation has stopped, value is meaningless.
	double sigma() const {
		return pythia->info.sigmaGen();
	}
	/// Returns the total weight of all the generated events as given by Pythia.
	/// If called before event generation has stopped, value is meaningless.
	double total_weight() const {
		return pythia->info.weightSum();
	}
};

#endif // PARTICLE_GENERATOR_H
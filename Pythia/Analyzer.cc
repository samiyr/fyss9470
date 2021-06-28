#ifndef ANALYZER_H
#define ANALYZER_H

#include "Pythia8/Pythia.h"
#include "Histogram.cc"
#include "Helpers.cc"
#include "Constants.cc"
#include "ParticleGenerator.cc"
#include "AnalysisParameters.cc"

using namespace Pythia8;

class Analyzer {
public:
	AnalysisParameters parameters;

	std::vector<double> bins;

	ValueHistogram<double> histogram;

	double sigma_sps_1 = 0.0;
	double sigma_sps_2 = 0.0;
	double sigma_sps = 0.0;

	// Disgusting hack
	// 0.0 * 1.0 * 2.0 * 3.0 * 4.0 * 5.0 * 6.0
	//  |___SPS_1___|___SPS_2___|___SPS_T___|
	// ValueHistogram<double> sigmas = ValueHistogram<double>({RangedContainer<double>(0.0, 2.0, 0.0), RangedContainer<double>(2.0, 4.0, 0.0), RangedContainer<double>(4.0, 6.0, 0.0)});
	// ValueHistogram<double> sigma_sps_2 = ValueHistogram<double>({-1.0, 1.0});
	// ValueHistogram<double> sigma_sps_1 = ValueHistogram<double>({-1.0, 1.0});
	// ValueHistogram<double> sigma_sps = ValueHistogram<double>({-1.0, 1.0});

	Analyzer(AnalysisParameters params, std::vector<double> b) {
		parameters = params;
		bins = b;
		histogram = ValueHistogram<double>(bins);
	}

	template <typename F>
	void book(std::vector<ParticleContainer> *input, F lambda) {
		const auto N = input->size();

		for (std::vector<ParticleContainer>::size_type i = 0; i < N; i++) {
			const ParticleContainer particle1 = (*input)[i];
			const bool check11 = parameters.pT_small.in_range(particle1.pT) && parameters.y_small.in_range(particle1.y);
			const bool check12 = parameters.pT_large.in_range(particle1.pT) && parameters.y_large.in_range(particle1.y);

			// pT and/or y ranges are assumed to be non-overlapping, i.e. a particle passes at most one filter
			if (check11) {
				// #pragma omp atomic
				sigma_sps_1 += particle1.event_weight;
				// sigmas.fill(1.0, particle1.event_weight);
			} else if (check12) {
				// #pragma omp atomic
				sigma_sps_2 += particle1.event_weight;
				// sigmas.fill(3.0, particle1.event_weight);
			} else {
				continue;
			}

			const double phi1 = particle1.phi;
			for (std::vector<ParticleContainer>::size_type j = i + 1; j < N; j++) {
				const ParticleContainer particle2 = (*input)[j];
				const bool check21 = parameters.pT_small.in_range(particle2.pT) && parameters.y_small.in_range(particle2.y);
				const bool check22 = parameters.pT_large.in_range(particle2.pT) && parameters.y_large.in_range(particle2.y);
				if (!((check21 && !check11) || (check22 && !check12))) {
					continue;
				}
				const double phi2 = particle2.phi;

				const double delta_phi = abs(phi1 - phi2);
				const double value = min(delta_phi, 2 * M_PI - delta_phi);

				sigma_sps += particle1.event_weight;
				// cout << "filling pair: " << particle1.event_weight << "\n";
				// sigmas.fill(5.0, particle1.event_weight);

				lambda(value);
			}
		}
	}

	// void finalize(ParticleGenerator *particle_generator) {
	// 	// #pragma omp critical
	// 	// {
	// 		// const double factor = particle_generator->sigma() / particle_generator->total_weight();
	// 		// cout << "factor: " << factor << "\n";
	// 		// sigmas *= factor;
	// 		// sigma_sps_2 *= factor;
	// 	// }
		
	// }

	void book(std::vector<ParticleContainer> *input) {
		book(input, [this](double value) {
			histogram.fill(value);
		});
	}
};

#endif // CORRELATION_ANALYZER_H
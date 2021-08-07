#ifndef ANALYZER_H
#define ANALYZER_H

#include "Pythia8/Pythia.h"
#include "Histogram.cc"
#include "Helpers.cc"
#include "Constants.cc"
#include "ParticleGenerator.cc"
#include <cassert>

using namespace Pythia8;

/// Analyzer for azimuthal correlations.
class Analyzer {
public:
	/// Analyzer parameters.
	struct Parameters {
		/// Transverse momenta filter sets.
		OptionalRange<double> pT_small;
		OptionalRange<double> pT_large;
		/// Rapidity filter sets.
		OptionalRange<double> y_small;
		OptionalRange<double> y_large;
		/// Filename of the exported file. 
		/// If set to `nullopt`, no data is exported.
		std::optional<std::string> filename;
		/// Parameters for the DPS model. 
		double m;
		double sigma_eff;
		/// Constructs an empty parameter set with no cuts in transverse momenta or rapidity.
		Parameters() : 
		Parameters(OptionalRange<double>(), OptionalRange<double>(), OptionalRange<double>(), OptionalRange<double>()) {}
		/// Constructs a parameter set with the cuts in transverse momenta and rapidity given as OptionalRange.
		Parameters(
			OptionalRange<double> pT_1, 
			OptionalRange<double> pT_2, 
			OptionalRange<double> y_1, 
			OptionalRange<double> y_2, 
			std::optional<std::string> fn = std::nullopt,
			double _m = Defaults::m,
			double _sigma_eff = Defaults::sigma_eff)
			:pT_small(pT_1), 
			pT_large(pT_2), 
			y_small(y_1), 
			y_large(y_2), 
			filename(fn), 
			m(_m), 
			sigma_eff(_sigma_eff) {
		}
		/// Constructs a parameter set with the cuts in transverse momenta and rapidity given as individual endpoints.
		Parameters(
			std::optional<double> pT_1_l,
			std::optional<double> pT_1_u,
			std::optional<double> pT_2_l,
			std::optional<double> pT_2_u,
			std::optional<double> y_1_l,
			std::optional<double> y_1_u,
			std::optional<double> y_2_l,
			std::optional<double> y_2_u,
			std::optional<string> fn = std::nullopt,
			double _m = Defaults::m,
			double _sigma_eff = Defaults::sigma_eff)
			:pT_small(OptionalRange<double>(pT_1_l, pT_1_u)), 
			pT_large(OptionalRange<double>(pT_2_l, pT_2_u)),
			y_small(OptionalRange<double>(y_1_l, y_1_u)), 
			y_large(OptionalRange<double>(y_2_l, y_2_u)),
			filename(fn), 
			m(_m), 
			sigma_eff(_sigma_eff) {
		}
	};
	/// The analysis parameters.
	Analyzer::Parameters parameters;
	/// Histogram bins.
	std::vector<double> bins;
	/// Output histogram.
	ValueHistogram<double> histogram;
	/// Total weights of trigger and associated particles, as well as of the paired particles.
	double N_trigger = 0.0;
	double N_assoc = 0.0;
	double N_pair = 0.0;
	/// Constructs an analyzer with the given parameters and histogram bins.
	Analyzer(Analyzer::Parameters params, std::vector<double> b) 
	: parameters(params), bins(b), histogram(ValueHistogram<double>(bins)) {}
	/// Analyze the particles generated at an event.
	void book(std::vector<ParticleContainer> *input) {
		const auto N = input->size();

		for (std::vector<ParticleContainer>::size_type i = 0; i < N; i++) {
			const ParticleContainer particle1 = (*input)[i];
			const bool check11 = parameters.pT_small.in_range(particle1.pT) && parameters.y_small.in_range(particle1.y);
			const bool check12 = parameters.pT_large.in_range(particle1.pT) && parameters.y_large.in_range(particle1.y);

			// pT and/or y ranges are assumed to be non-overlapping, i.e. a particle passes at most one filter
			if (check11) {
				N_trigger += particle1.event_weight;
			} else if (check12) {
				N_assoc += particle1.event_weight;
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

				N_pair += particle1.event_weight;

				histogram.fill(value, particle1.event_weight);
			}
		}
	}

};

#endif // CORRELATION_ANALYZER_H
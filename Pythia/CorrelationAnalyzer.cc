#ifndef CORRELATION_ANALYZER_H
#define CORRELATION_ANALYZER_H

#include "Pythia8/Pythia.h"
#include "Histogram.cc"
#include "Helpers.cc"
#include "Constants.cc"

using namespace Pythia8;

struct CorrelationAnalyzerParameters {
	OptionalRange<double> pT_small;
	OptionalRange<double> pT_large;

	OptionalRange<double> y_small;
	OptionalRange<double> y_large;

	/// Filename of the exported file. 
	/// If set to `nullopt`, no data is exported.
	std::optional<std::string> filename;

	CorrelationAnalyzerParameters() {
		pT_small = OptionalRange<double>();
		pT_large = OptionalRange<double>();
		y_small = OptionalRange<double>();
		y_large = OptionalRange<double>();
		filename = std::nullopt;
	}

	CorrelationAnalyzerParameters(
		OptionalRange<double> pT_1, 
		OptionalRange<double> pT_2, 
		OptionalRange<double> y_1, 
		OptionalRange<double> y_2, 
		std::optional<std::string> fn = std::nullopt) {
		pT_small = pT_1;
		pT_large = pT_2;
		y_small = y_1;
		y_large = y_2;
		filename = fn;
	}

	CorrelationAnalyzerParameters(
		std::optional<double> pT_1_l,
		std::optional<double> pT_1_u,
		std::optional<double> pT_2_l,
		std::optional<double> pT_2_u,
		std::optional<double> y_1_l,
		std::optional<double> y_1_u,
		std::optional<double> y_2_l,
		std::optional<double> y_2_u,
		std::optional<string> fn = std::nullopt) {
		pT_small = OptionalRange<double>(pT_1_l, pT_1_u);
		pT_large = OptionalRange<double>(pT_2_l, pT_2_u);
		y_small = OptionalRange<double>(y_1_l, y_1_u);
		y_large = OptionalRange<double>(y_2_l, y_2_u);
		filename = fn;
	}
};

class CorrelationAnalyzer {
public:
	CorrelationAnalyzerParameters parameters;

	std::vector<double> bins;

	ValueHistogram<double> histogram;

	CorrelationAnalyzer(CorrelationAnalyzerParameters params, std::vector<double> b) {
		parameters = params;
		bins = b;
		histogram = ValueHistogram<double>(bins);
	}

	void book(std::vector<ParticleContainer> *input) {
		const auto N = input->size();

		for (std::vector<ParticleContainer>::size_type i = 0; i < N; i++) {
			const ParticleContainer particle1 = (*input)[i];
			const bool check11 = parameters.pT_small.in_range(particle1.pT) && parameters.y_small.in_range(particle1.y);
			const bool check12 = parameters.pT_large.in_range(particle1.pT) && parameters.y_large.in_range(particle1.y);
			if (!(check11 || check12)) {
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

				histogram.fill(value);
			}
		}
	}
};

#endif // CORRELATION_ANALYZER_H
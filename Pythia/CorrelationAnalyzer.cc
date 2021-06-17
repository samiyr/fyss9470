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
	std::vector<CorrelationAnalyzerParameters>::size_type run_index;
	std::vector<CorrelationAnalyzerParameters>::size_type run_count;

	std::vector<std::vector<ParticleContainer>> *particles;

	CorrelationAnalyzer(CorrelationAnalyzerParameters params, std::vector<std::vector<ParticleContainer>> *ps) {
		parameters = params;
		particles = ps;
	}

	ValueHistogram<double> analyze() const {
		ValueHistogram<double> hist(bins);

		const auto total_N = particles->size();
		int current_N = 0;
		const auto threshold = total_N / Defaults::status_threshold;

		for (auto &list : *particles) {
			const auto N = list.size();
			for (std::vector<ParticleContainer>::size_type i = 0; i < N; i++) {
				const Particle particle1 = list[i].particle;
				const bool check11 = parameters.pT_small.in_range(particle1.pT()) && parameters.y_small.in_range(particle1.y());
				const bool check12 = parameters.pT_large.in_range(particle1.pT()) && parameters.y_large.in_range(particle1.y());
				if (!(check11 || check12)) {
					continue;
				}
				const double phi1 = particle1.phi();
				for (std::vector<ParticleContainer>::size_type j = i + 1; j < N; j++) {
					const Particle particle2 = list[j].particle;
					const bool check21 = parameters.pT_small.in_range(particle2.pT()) && parameters.y_small.in_range(particle2.y());
					const bool check22 = parameters.pT_large.in_range(particle2.pT()) && parameters.y_large.in_range(particle2.y());
					if (!((check21 && !check11) || (check22 && !check12))) {
						continue;
					}
					const double phi2 = particle2.phi();

					const double delta_phi = abs(phi1 - phi2);
					const double value = min(delta_phi, 2 * M_PI - delta_phi);

					hist.fill(value);
				}
			}
			current_N++;

			if (current_N % threshold == 0) {
				cout << "analysis " << run_index + 1 << "/" << run_count << ": " << (double)current_N / total_N * Defaults::status_threshold << "%\n";
			}
		}

		return hist;		
	}
};

#endif // CORRELATION_ANALYZER_H
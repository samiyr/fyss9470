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
		/// Defaults to null. If a value is given, a nuclear histogram is generated by
		/// adding the effect of nuclear PDFs afterwards instead of during event generation.
		/// This assumes that nuclear PDFs weren't used during event generation.
		std::optional<Beam> beam;
		/// Working directory for the exported file.
		/// If set to `nullopt`, the working directory
		/// is taken from the Experiment class.
		std::optional<std::string> working_directory;
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
			double _sigma_eff = Defaults::sigma_eff,
			std::optional<Beam> _beam = std::nullopt,
			std::optional<std::string> _wd = std::nullopt)
			:pT_small(pT_1), 
			pT_large(pT_2), 
			y_small(y_1), 
			y_large(y_2), 
			filename(fn), 
			m(_m), 
			sigma_eff(_sigma_eff),
			beam(_beam),
			working_directory(_wd) {}
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
			double _sigma_eff = Defaults::sigma_eff,
			std::optional<Beam> _beam = std::nullopt,
			std::optional<std::string> _wd = std::nullopt)
			:pT_small(OptionalRange<double>(pT_1_l, pT_1_u)), 
			pT_large(OptionalRange<double>(pT_2_l, pT_2_u)),
			y_small(OptionalRange<double>(y_1_l, y_1_u)), 
			y_large(OptionalRange<double>(y_2_l, y_2_u)),
			filename(fn), 
			m(_m), 
			sigma_eff(_sigma_eff),
			beam(_beam),
			working_directory(_wd) {}
	};
	/// The analysis parameters.
	Analyzer::Parameters parameters;

	/// Histogram bins.
	std::vector<double> bins;

	/// Output histogram.
	ValueHistogram<double> histogram;

	/// Output histogram with post-hoc nuclear modifications.
	ValueHistogram<double> nuclear_histogram;

	/// Momentum fraction histograms.
	/// The "pre" histogram is filled once per event, regardless of cuts in pT or y.
	/// The "post" histogram is filled once per accepted pion pair.
	ValueHistogram<double> x1_pre_histogram;
	ValueHistogram<double> x2_pre_histogram;
	ValueHistogram<double> x1_post_histogram;
	ValueHistogram<double> x2_post_histogram;

	/// Total weights of trigger and associated particles, as well as of the paired particles.
	double N_trigger = 0.0;
	double N_assoc = 0.0;
	double N_pair = 0.0;

	/// Parton distribution function pointers.	
	PDFPtr protonPDF;
	std::optional<EPPS16 *> nuclearPDF;

	/// Constructs an analyzer with the given parameters and histogram bins.
	Analyzer(Analyzer::Parameters params, std::vector<double> b) 
	: parameters(params), 
	bins(b), 
	histogram(ValueHistogram<double>(bins)),
	nuclear_histogram(ValueHistogram<double>(bins)),
	x1_pre_histogram(ValueHistogram<double>(fixed_range(0.0, 1.0, 20))),
	x2_pre_histogram(ValueHistogram<double>(fixed_range(0.0, 1.0, 20))),
	x1_post_histogram(ValueHistogram<double>(fixed_range(0.0, 1.0, 20))),
	x2_post_histogram(ValueHistogram<double>(fixed_range(0.0, 1.0, 20))) {
		protonPDF = make_shared<LHAGrid1>(2212, "13", "/Users/samiyrjanheikki/pythia8306/share/Pythia8/pdfdata/");
		if (parameters.beam) {
			nuclearPDF = new EPPS16((*parameters.beam).nucleus.pdg_code(), 1, "/Users/samiyrjanheikki/pythia8306/share/Pythia8/pdfdata/", protonPDF); 
		}
	}

	/// Analyze the particles generated at an event.
	void book(std::vector<ParticleContainer> *input, Info *info) {
		const auto N = input->size();

   		const double x1 = info->x1();
   		const int id2 = info->id2pdf();
   		const double x2 = info->x2();
   		const double Q2 = info->Q2Fac();

   		double xf2 = 1.0;
   		double xfA = 1.0;

   		if (nuclearPDF) {
   			xf2 = protonPDF->xf(id2, x2, Q2);
   			xfA = nuclearPDF.value()->xf(id2, x2, Q2);
   		}

   		const double event_weight = info->weight();

   		x1_pre_histogram.fill(x1, event_weight);
   		x2_pre_histogram.fill(x2, event_weight);

		for (std::vector<ParticleContainer>::size_type i = 0; i < N; i++) {
			const ParticleContainer particle1 = (*input)[i];
			const bool check11 = parameters.pT_small.in_range(particle1.pT) && parameters.y_small.in_range(particle1.y);
			const bool check12 = parameters.pT_large.in_range(particle1.pT) && parameters.y_large.in_range(particle1.y);

			// pT and/or y ranges are assumed to be non-overlapping, i.e. a particle passes at most one filter
			if (check11) {
				N_trigger += event_weight;
			} else if (check12) {
				N_assoc += event_weight;
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

				N_pair += event_weight;

				histogram.fill(value, event_weight);
				if (nuclearPDF) {
					nuclear_histogram.fill(value, event_weight * xfA / xf2);
				}
				x1_post_histogram.fill(x1, event_weight);
   				x2_post_histogram.fill(x2, event_weight);
			}
		}
	}

};

#endif // CORRELATION_ANALYZER_H
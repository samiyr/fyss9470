#ifndef BEAM_H
#define BEAM_H

#include "Pythia8/Pythia.h"
#include <cassert>
#include "NumberListReader.cc"
#include "Constants.cc"

using namespace Pythia8;

/// A container of beam parameters.
struct Beam {
	/// Nuclear parton distribution function sets available.
	/// Usage of others than None require installation of the
	/// PDF grid files, see the section "Nuclear modifications of parton densities"
	/// in the Pythia online manual for "PDF Selection".
	enum class NuclearPDF : int {
		None = 0,
		EPS09LO = 1,
		EPS09NLO = 2,
		EPPS16NLO = 3
	};
	/// Representation of a nucleus.
	struct Nucleus {
		/// Atomic number Z, the number of protons.
		int atomic_number;
		/// Mass number A, the total number of protons and neutrons.
		int mass_number;
		///Total number of strange quarks.
		int strange_count;
		/// Isomer level.
		int isomer_level;
		/// Constructs a nucleus representation for the proton.
		Nucleus() : Nucleus(1, 1, 0, 0) {}
		/// Constructs a nucleus with the given parameters. B
		/// y default, the number of strange quarks and isomer levels are set to zero.
		Nucleus(int Z, int A, int ns = 0, int I = 0) : atomic_number(Z), mass_number(A), strange_count(ns), isomer_level(I) {}
		/// Returns the PDG code for the nucleus. See https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf.
		/// Al-27 = 1000130270
		/// Au-197 = 1000791970
		int pdg_code() const {
			int base = 1'000'000'000;
			int L = strange_count * 10'000'000;
			int Z = atomic_number * 10'000;
			int A = mass_number * 10;
			int I = isomer_level;
			return base + L + Z + A + I;
		}
	};
	/// The nucleus of the beam particle.
	Nucleus nucleus;
	/// The nuclear parton distribution functions used for the beam.
	NuclearPDF pdf;
	/// A flag determining whether to enable nuclear parton distribution functions.
	bool use_hard_npdf;
	/// Constructs the default beam, which is a proton with no nPDFs.
	Beam() : nucleus(Nucleus()), pdf(NuclearPDF::None), use_hard_npdf(false) {}
	/// Constructs a beam with the given nucleus, nPDF set. By default, nPDFs are disabled.
	Beam(Nucleus _nucleus, NuclearPDF _pdf = NuclearPDF::None, bool _use_hard_npdf = false) 
	: nucleus(_nucleus), pdf(_pdf), use_hard_npdf(_use_hard_npdf) {}
	/// Constructs a beam with the given atomic and mass numbers, and the nPDF set.
	/// By default, nPDFs are disabled.
	Beam(int Z, int A, NuclearPDF _pdf = NuclearPDF::None, bool _use_hard_npdf = false) 
	: nucleus(Nucleus(Z, A)), pdf(_pdf), use_hard_npdf(_use_hard_npdf) {}
	/// Apply the appropriate settings to the Pythia settings object.
	void apply_to(Settings &settings, string beam) {
		settings.flag("PDF:useHardNPDF" + beam, use_hard_npdf);
		settings.mode("PDF:nPDFSet" + beam, static_cast<int>(pdf));
		settings.mode("PDF:nPDFBeam" + beam, nucleus.pdg_code());
	}
	/// Define the textual format of the beam for outputting to streams.
	friend ostream& operator<<(ostream& os, Beam const & beam) {
		os << beam.nucleus.atomic_number << ", " << beam.nucleus.mass_number << ", " << beam.nucleus.pdg_code() << "; ";
		switch(beam.pdf) {
			case NuclearPDF::None:
				os << "None";
				break;
			case NuclearPDF::EPS09LO:
				os << "EPS09LO";
				break;
			case NuclearPDF::EPS09NLO:
				os << "EPS09NLO";
				break;
			case NuclearPDF::EPPS16NLO:
				os << "EPPS16NLO";
				break;
		}
		if (beam.use_hard_npdf) {
			os << " (nPDF on)";
		} else {
			os << " (nPDF off)";
		}
		return os;
    }

    std::optional<NumberListReader<int>> ncoll_list() {
    	if (nucleus.atomic_number == 13 && nucleus.mass_number == 27) {
    		return NumberListReader<int>(Constants::pAl_ncoll_list_file);
    	} else if (nucleus.atomic_number == 79 && nucleus.mass_number == 197) {
    		return NumberListReader<int>(Constants::pAu_ncoll_list_file);
    	} else {
    		return std::nullopt;
    	}
    }
};

/// Calculates the effective cross section for the given beam and with the given
/// effective cross section for p+p collisions. Currently, only nuclei with mass numbers
/// A = 27 (Al) and A = 197 (Au) are supported.

double calculate_sigma_eff(int B, double sigma_pp) {
	double geometric_integral;

	switch(B) {
		case 1:
			geometric_integral = 0.0;
			break;
		case 197:
			geometric_integral = 29.353;
			break;
		case 27:
			geometric_integral = 1.700;
			break;
		default:
			assert(false);
			break;
	}

	return (B * B * sigma_pp) / (B * B + (B - 1) * geometric_integral * sigma_pp);
}

double calculate_sigma_eff(Beam b, double sigma_pp) {
	return calculate_sigma_eff(b.nucleus.mass_number, sigma_pp);
}

#endif // BEAM_H
#ifndef BEAM_H
#define BEAM_H

#include "Pythia8/Pythia.h"
#include <cassert>

using namespace Pythia8;

struct Beam {
	enum class NuclearPDF : int {
		None = 0,
		EPS09LO = 1,
		EPS09NLO = 2,
		EPPS16NLO = 3
	};
	struct Nucleus {
		int atomic_number;
		int mass_number;
		int strange_count;
		int isomer_level;

		Nucleus() : Nucleus(1, 1, 0, 0) {}

		Nucleus(int Z, int A, int ns = 0, int I = 0) : atomic_number(Z), mass_number(A), strange_count(ns), isomer_level(I) {}

		int pdg_code() const {
			int base = 1'000'000'000;
			int L = strange_count * 10'000'000;
			int Z = atomic_number * 10'000;
			int A = mass_number * 10;
			int I = isomer_level;
			return base + L + Z + A + I;
		}
	};

	Nucleus nucleus;
	NuclearPDF pdf;
	bool use_hard_npdf;

	Beam() {
		nucleus = Nucleus();
		pdf = NuclearPDF::None;
		use_hard_npdf = false;
	}

	Beam(Nucleus _nucleus, NuclearPDF _pdf = NuclearPDF::None, bool _use_hard_npdf = false) : nucleus(_nucleus), pdf(_pdf), use_hard_npdf(_use_hard_npdf) {}

	Beam(int Z, int A, NuclearPDF _pdf = NuclearPDF::None, bool _use_hard_npdf = false) : nucleus(Nucleus(Z, A)), pdf(_pdf), use_hard_npdf(_use_hard_npdf) {}

	void apply_to(Settings &settings, string beam) {
		settings.flag("PDF:useHardNPDF" + beam, use_hard_npdf);
		settings.mode("PDF:nPDFSet" + beam, static_cast<int>(pdf));
		settings.mode("PDF:nPDFBeam" + beam, nucleus.pdg_code());
	}
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
};


double calculate_sigma_eff(Beam b, double sigma_pp) {
	double geometric_integral;
	const int B = b.nucleus.mass_number;

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

#endif // BEAM_H
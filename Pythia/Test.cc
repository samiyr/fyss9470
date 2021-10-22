#ifndef TEST_H
#define TEST_H

#include "Pythia8/Pythia.h"
#include "Histogram.cc"

using namespace Pythia8;

int main() {
	Pythia pythia;
	pythia.readString("Beams:idA = 2212");
	pythia.readString("Beams:idB = 2212");
	pythia.readString("Random:setSeed = on");
	pythia.readString("Beams:eCM = 200");
	pythia.readString("PhaseSpace:pTHatMin = 1.5");
	pythia.readString("Print:quiet = on");
	pythia.readString("Next:numberCount = 10000");
	pythia.readString("Random:seed = 1");
	pythia.readString("PartonLevel:MPI = off");
	pythia.readString("HardQCD:All = on");
	pythia.init();

	int N_trigger = 0;
	int N_assoc = 0;
	int N_pair = 0;
	double sigma_gen = 0;
	int events = 0;

	// histogram from 0 to pi with 20 bins
	ValueHistogram<double> histogram(fixed_range(0.0, M_PI, 20));

	for (int i = 0; i < 2'000'000; i++) {
		if (!pythia.next()) {
			continue;
		}

		std::vector<Particle> particles;
		Event event = pythia.event;
		Info info = pythia.info;

		// only the value at the end is needed, but this is simpler
		sigma_gen = info.sigmaGen();
		events = info.weightSum();

		// pick out pions
		for (int j = 0; j < event.size(); j++) {
			Particle particle = event[j];
			
			if (particle.id() == 111) {
				particles.push_back(particle);
			}
		}

		// pair the pions
		for (unsigned long i1 = 0; i1 < particles.size(); i1++) {
			Particle particle1 = particles[i1];
			bool check11 = (particle1.pT() >= 1.0 && particle1.pT() < 1.4) 
							&& (particle1.y() >= 2.6 && particle1.y() < 4.1);
			bool check12 = (particle1.pT() >= 1.4 && particle1.pT() < 2.0) 
							&& (particle1.y() >= 2.6 && particle1.y() < 4.1);

			if (check11) {
				// pion with smaller pT
				N_trigger += 1;
			} else if (check12) {
				// pion with larger pT
				N_assoc += 1;
			} else {
				continue;
			}

			double phi1 = particle1.phi();

			for (unsigned long i2 = i1 + 1; i2 < particles.size(); i2++) {
				Particle particle2 = particles[i2];
				bool check21 = (particle2.pT() >= 1.0 && particle2.pT() < 1.4) 
							&& (particle2.y() >= 2.6 && particle2.y() < 4.1);
				bool check22 = (particle2.pT() >= 1.4 && particle2.pT() < 2.0) 
							&& (particle2.y() >= 2.6 && particle2.y() < 4.1);

				if (!((check21 && !check11) || (check22 && !check12))) {
					continue;
				}

				double phi2 = particle2.phi();
				double delta_phi = abs(phi1 - phi2);
				double value = min(delta_phi, 2 * M_PI - delta_phi);
				N_pair += 1;
				histogram.fill(value);
			}
		}
	}

	cout << "N_trigger = " << N_trigger << "\n";
	cout << "N_assoc = " << N_assoc << "\n";
	cout << "N_pair = " << N_pair << "\n";
	cout << "N = " << events << "\n";
	cout << "sigma_gen = " << sigma_gen << "\n";

	cout << histogram << "\n";

	// delta_phi = pi / 20 = 0.157
	histogram *= 1.0 / (N_assoc * 0.157);
	histogram += (0.5 * sigma_gen * N_trigger) / (M_PI * 10 * events);

	cout << histogram << "\n";

	return 0;
}

#endif // TEST_H
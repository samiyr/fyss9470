#ifndef TEST2_H
#define TEST2_H

#include "Pythia8/Pythia.h"
#include "Histogram.cc"
#include "Beam.cc"
#include "NumberListReader.cc"

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
	pythia.readString("PartonLevel:MPI = on");
	pythia.readString("SoftQCD:nonDiffractive = on");

	pythia.init();

	int N_trigger = 0;
	int N_assoc = 0;
	int N_pair = 0;
	int events = 0;
	int nexts = 0;
	int failures = 0;
	int selected = 0;
	int rejected = 0;

	// histogram from 0 to pi with 20 bins
	auto points = fixed_range(0.0, M_PI, 20);
	ValueHistogram<double> histogrampp(points);
	ValueHistogram<double> histogramAl(points);
	ValueHistogram<double> histogramAu(points);

	ValueHistogram<double> x1histogram(fixed_range(0.0, 1.0, 20));
	ValueHistogram<double> x2histogram(fixed_range(0.0, 1.0, 20));

	NumberListReader<int> pAl_ncoll("./n_coll_pAl.txt");

	int retry_counter = 0;

	for (int i = 0; i < 10'000; i++) {
		events++;
		double x1_tot = 0.0;

		const int n_coll = 1; //pAl_ncoll(i);

		std::vector<Particle> particles;
		// collision loop
		for (int i_coll = 0; i_coll < n_coll; i_coll++) {
			nexts++;
			if (!pythia.next()) {
				failures++;
				continue;
			}

			Event event = pythia.event;
			Info info = pythia.info;

			// pick out pions
			for (int j = 0; j < event.size(); j++) {
				Particle particle = event[j];
				
				if (particle.id() == 111) {
					particles.push_back(particle);
				}
			}
			// calculate x1
			const int n_MPI = pythia.info.nMPI();
     		double x1_sum = 0.0;
     		
     		for (int i_MPI = 0; i_MPI < n_MPI; i_MPI++) {
       			double x1 = (*info.beamAPtr)[i_MPI].x();
       			x1_sum += x1;
     		}
     		x1_tot += x1_sum;
		}

		// check energy limit
		if (x1_tot > 1.0) {
			// cout << "x1 = " << x1_tot << " exceeded 1 in event " << i << " with ncoll = " << n_coll << "\n";
			retry_counter++;
			if (retry_counter > 0) {
				rejected++;
			} else {
				i--;
			}
			continue;
		}
		retry_counter = 0;
		selected++;
		// pair the pions after collision loop
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

				histogramAl.fill(value, 1);
			}
		}
	}

	cout << "N_trigger = " << N_trigger << "\n";
	cout << "N_assoc = " << N_assoc << "\n";
	cout << "N_pair = " << N_pair << "\n";
	cout << "N = " << events << "\n";
	cout << "nexts = " << nexts << "\n";
	cout << "generated events = " << nexts - failures << "\n";
	cout << "selected = " << selected << "\n";
	cout << "rejected = " << rejected << "\n";

	cout << histogramAl << "\n";

	histogramAl *= 1.0 / (N_assoc * 0.157);

	cout << histogramAl << "\n";

	return 0;
}

#endif // TEST2_H
#include "Pythia8/Pythia.h"
#include <string>
#include "PionGenerator.cc"
#include <boost/histogram.hpp>
#include <math.h>

using namespace Pythia8;
using namespace boost;

#define INCLUDE_DECAY 	true
#define Y_MIN	 		-0.35
#define Y_MAX 	 		0.35

void cross_section(double energy, int count) {
	cout << "Starting experiment with E = " << energy << ", N = " << count << "\n";
	cout << "Generating pions" << "\n";

	PionGenerator generator(energy, count);

	generator.include_decayed = INCLUDE_DECAY;
	generator.y_min = Y_MIN;
	generator.y_max = Y_MAX;

	generator.initialize();

	const std::vector<Particle> pions = generator.generate();
	cout << "Generated " << pions.size() << " pions" << "\n";

	const std::vector<double> pTs = find_pT(pions);
	const double min_pT = *std::min_element(pTs.begin(), pTs.end());
	const double max_pT = *std::max_element(pTs.begin(), pTs.end());
	const double mean_pT = mean(pTs);
	cout << "pT in range [" << min_pT << ", " << max_pT << "], mean = " << mean_pT << "\n";

	const double sigma = generator.sigma();
	cout << "sigma = " << sigma << "\n";

	auto axis = histogram::axis::variable({
		1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0
	}, "pT");
	auto hist = histogram::make_histogram(axis);
	hist.fill(pTs);

	cout << "Non-normalized pT histogram" << "\n";
	print_histogram(hist, axis);
	cout << "\n";

	std::vector<Bin<double>> normalized;
	for (histogram::axis::index_type i = 0; i < axis.size(); i++) {
		const int N = (int)hist.at(i);
		const double dy = generator.y_max - generator.y_min;
		const auto b = axis.bin(i);
		const double dpT = b.width();
		const double pT = b.center();
		const double dsigma = N * sigma / (2 * M_PI * count * dy * dpT * pT);
		Bin<double> bin(b.lower(), b.upper(), dsigma);
		normalized.push_back(bin);
	}

	cout << "Normalized pT histogram" << "\n";
	print_bins(normalized);
}

int main() {
	cross_section(200, 10000);
	
	return 0;
}


#include "Pythia8/Pythia.h"
#include <string>
#include "PionGenerator.cc"
#include "Histogram.cc"
#include "Helpers.cc"
#include <math.h>

using namespace Pythia8;

#define INCLUDE_DECAY 	true
#define Y_MIN	 		-0.35
#define Y_MAX 	 		0.35

std::vector<RangedContainer<double>> normalize(Histogram<double> hist, PionGenerator *generator, bool constant_pT) {
	const double sigma = generator->sigma();
	const int count = generator->event_count;
	std::vector<RangedContainer<double>> normalized;
	for (auto bin : hist.bins) {
		const int N = (int)bin.size();
		const double dy = generator->y_max - generator->y_min;
		const double dpT = bin.width();
		const double pT = bin.center();
		double dsigma;
		if (constant_pT) {
			dsigma = N * sigma / (2 * M_PI * count * dy * dpT * pT);
		} else {
			dsigma = std::accumulate(bin.contents.begin(), bin.contents.end(), 0.0, [dy, dpT, sigma, count](double sum, double pT) {
				double new_value = sigma / (2 * M_PI * count * dy * dpT * pT);
				return sum + new_value;
			});
		}
		RangedContainer<double> container(bin.lower, bin.upper, dsigma);
		normalized.push_back(container);
	}
	return normalized;
}

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

	Histogram<double> hist = Histogram<double>({
		1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0
	});

	hist.fill(pTs);

	cout << "Non-normalized pT histogram" << "\n";
	hist.print();
	cout << "\n";

	auto constant_normalized = normalize(hist, &generator, true);
	auto nonconstant_normalized = normalize(hist, &generator, false);

	cout << "Normalized constant pT histogram" << "\n";
	print_containers(constant_normalized);
	cout << "\n";
	cout << "Normalized nonconstant pT histogram" << "\n";
	print_containers(nonconstant_normalized);
	export_containers(constant_normalized, "constant_pT.csv");
	export_containers(nonconstant_normalized, "nonconstant_pT.csv");
}

int main() {
	cross_section(200, 50000);

	return 0;
}


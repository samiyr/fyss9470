#include "Pythia8/Pythia.h"
#include <string>
#include "ParticleGenerator.cc"
#include "Histogram.cc"
#include "Helpers.cc"
#include <math.h>

using namespace Pythia8;

#define INCLUDE_DECAY 	true
#define Y_MIN	 		-0.35
#define Y_MAX 	 		0.35

std::vector<RangedContainer<double>> normalize(Histogram<double> hist, int count, double sigma, bool constant_pT) {
	std::vector<RangedContainer<double>> normalized;
	for (auto bin : hist.bins) {
		const int N = (int)bin.size();
		const double dy = Y_MAX - Y_MIN;
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

std::vector<RangedContainer<double>> combine(std::vector<std::vector<RangedContainer<double>>> containers, std::vector<double> weights) {
	std::vector<RangedContainer<double>> result;
	const auto reference = containers.front();
	const auto N = reference.size();
	for (std::vector<RangedContainer<double>>::size_type i = 0; i < N; i++) {
		const double lower = reference[i].lower;
		const double upper = reference[i].upper;

		double value = 0;

		for (std::vector<RangedContainer<double>>::size_type j = 0; j < containers.size(); j++) {
			const double weight = weights[j];
			value += weight * containers[j][i].value;
		}

		RangedContainer<double> container(lower, upper, value);
		result.push_back(container);
	}
	return result;
}

void cross_section(double energy, int count, std::vector<double> bins, std::vector<double> pT_hat_bins) {
	std::vector<std::vector<RangedContainer<double>>> containers;
	std::vector<double> weights;

	for (std::vector<double>::size_type i = 0; i < pT_hat_bins.size() - 1; i++) {
		const double pT_hat_min = pT_hat_bins[i];
		const double pT_hat_max = pT_hat_bins[i + 1];

		ParticleGenerator generator(energy, count);

		generator.include_decayed = INCLUDE_DECAY;
		generator.y_min = Y_MIN;
		generator.y_max = Y_MAX;
		generator.pT_hat_min = pT_hat_min;
		generator.pT_hat_max = pT_hat_max;

		generator.initialize();

		const std::vector<Particle> pions = generator.generate();
		const double sigma = generator.sigma();

		std::vector<double> pTs = find_pTs(pions);

		Histogram<double> partial = Histogram<double>(bins);
		partial.fill(pTs);

		const auto partial_container = normalize(partial, count, sigma, false);
		containers.push_back(partial_container);

		weights.push_back(1);
	}

	const auto combined = combine(containers, weights);

	cout << "Normalized pT histogram" << "\n";
	print_containers(combined);
	cout << "\n";
	export_containers(combined, "pT_histogram.csv");
}

void azimuth_correlation(double energy, int count) {
	cout << "Starting experiment with E = " << energy << ", N = " << count << "\n";
	cout << "Generating pions" << "\n";

	ParticleGenerator generator(energy, count);

	generator.include_decayed = INCLUDE_DECAY;
	generator.y_min = Y_MIN;
	generator.y_max = Y_MAX;

	generator.initialize();

	const std::vector<Particle> pions = generator.generate();
	cout << "Generated " << pions.size() << " pions" << "\n";

	std::vector<double> deltas;
	for (std::vector<Particle>::size_type i = 0; i < pions.size(); i++) {
		const Particle p1 = pions[i];
		const double phi1 = p1.phi();
		for (std::vector<Particle>::size_type j = i + 1; j < pions.size(); j++) {
			const Particle p2 = pions[j];
			const double phi2 = p2.phi();
			const double delta_phi = abs(phi1 - phi2);
			deltas.push_back(min(delta_phi, 2 * M_PI - delta_phi));
		}
	}

	ofstream file;
	file.open("delta_phi.csv");
	file << std::setprecision(12);
	for (auto delta : deltas) {
		file << delta << "\n";
	}
	file.close();
}

int main() {
	const std::vector<double> bins = {
		1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0
	};
	const std::vector<double> pT_hat_bins = {
		2.0, 5.0, 10.0, 40.0, -1
	};
	cross_section(200, 10000, bins, pT_hat_bins);

	return 0;
}


#include "Pythia8/Pythia.h"
#include <string>
#include "ParticleGenerator.cc"
#include "Histogram.cc"
#include "Helpers.cc"
#include <math.h>
#include <execution>

using namespace Pythia8;

ValueHistogram<double> combine(std::vector<ValueHistogram<double>> containers, std::vector<double> weights) {
	auto reference = containers.front();
	const auto N = reference.size();
	ValueHistogram<double> result(N);
	for (std::vector<RangedContainer<double>>::size_type i = 0; i < N; i++) {
		double lower = reference[i].range.start;
		double upper = reference[i].range.end;

		double value = 0;

		for (std::vector<RangedContainer<double>>::size_type j = 0; j < containers.size(); j++) {
			const double weight = weights[j];
			value += weight * containers[j][i].value;
		}

		RangedContainer<double> container(lower, upper, value);
		result.containers.push_back(container);
	}
	return result;
}

void cross_section(double energy, int count, std::vector<double> bins, std::vector<double> pT_hat_bins, Range<double> y_range, bool include_decayed, bool use_biasing) {
	std::vector<ValueHistogram<double>> containers;
	std::vector<double> weights;

	for (std::vector<double>::size_type i = 0; i < pT_hat_bins.size() - 1; i++) {
		const double pT_hat_min = pT_hat_bins[i];
		const double pT_hat_max = pT_hat_bins[i + 1];

		ParticleGenerator generator(energy, count);

		generator.include_decayed = include_decayed;
		generator.y_range = y_range;
		generator.pT_hat_range = Range<double>(pT_hat_min, pT_hat_max);
		generator.use_biasing = use_biasing;

		generator.initialize();

		const std::vector<ParticleContainer> pions = generator.generate();
		const double sigma = generator.sigma();

		const std::vector<double> pTs = find_pTs(pions);
		const std::vector<double> event_weights = find_event_weights(pions);

		Histogram<double> partial = Histogram<double>(bins);
		partial.fill(pTs, event_weights);

		const auto partial_container = partial.normalize(generator.total_weight(), sigma, y_range, use_biasing);
		containers.push_back(partial_container);

		weights.push_back(1.0);
	}

	const auto combined = combine(containers, weights);

	cout << "Normalized pT histogram" << "\n";
	combined.print();
	cout << "\n";
	combined.export_histogram("pT_histogram.csv");
}

void azimuth_correlation(double energy, int count, Range<double> y_range, bool include_decayed) {
	cout << "Starting experiment with E = " << energy << ", N = " << count << "\n";
	cout << "Generating pions" << "\n";

	ParticleGenerator generator(energy, count);

	generator.include_decayed = include_decayed;
	generator.y_range = y_range;

	generator.initialize();

	const std::vector<ParticleContainer> pions = generator.generate();
	cout << "Generated " << pions.size() << " pions" << "\n";

	std::vector<double> deltas;
	for (std::vector<ParticleContainer>::size_type i = 0; i < pions.size(); i++) {
		const Particle p1 = pions[i].particle;
		const double phi1 = p1.phi();
		for (std::vector<Particle>::size_type j = i + 1; j < pions.size(); j++) {
			const Particle p2 = pions[j].particle;
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
	cross_section(
		200, 						// center-of-mass energy in GeV
		10000, 						// number of events
		bins, 						// pT histogram bins
		pT_hat_bins, 				// pT_hat bins
		Range<double>(-0.35, 0.35), // rapidity range
		true, 						// include decayed particles
		true						// use event weighting
	);

	return 0;
}


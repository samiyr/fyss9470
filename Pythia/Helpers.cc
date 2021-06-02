#include "Pythia8/Pythia.h"
#include <algorithm>
#include <string>
#include <utility>
#include <numeric>
#include <boost/histogram.hpp>
#include <iomanip>

using namespace Pythia8;
using namespace boost;

template <typename T>
bool contains(std::vector<T> vec, T element) {
	return std::find(vec.begin(), vec.end(), element) != vec.end();
}

std::vector<Particle> all_particles(std::vector<Event> events) {
	std::vector<Particle> particles;
	for (Event event : events) {
		const int particle_count = event.size();
		for (int i = 0; i < particle_count; ++i) {
			const Particle particle = event[i];
			particles.push_back(particle);
		}
	}
	return particles;
}

std::vector<double> find_azimuths(std::vector<Particle> particles) {
	std::vector<double> phis;
	for (Particle particle : particles) {
		phis.push_back(particle.phi());
	}
	return phis;
}
std::vector<double> find_pT(std::vector<Particle> particles) {
	std::vector<double> phis;
	for (Particle particle : particles) {
		phis.push_back(particle.pT());
	}
	return phis;
}
template <typename T>
T sum(std::vector<T> v) {
	T sum = std::accumulate(std::begin(v), std::end(v), 0.0);
	return sum;
}
template <typename T>
T mean(std::vector<T> v) {
	return sum(v) / v.size();
}

void print_event(Event event) {
	for (int i = 0; i < event.size(); ++i) {
		const Particle particle = event[i];
		cout << particle.name() << "\n";
	}
}

template <typename Histogram, typename Axis>
void print_histogram(Histogram hist, Axis axis) {
	for (histogram::axis::index_type i = 0; i < axis.size(); i++) {
		const int size = (int)hist.at(i);
		cout << "[" << axis.bin(i).lower() << ", " << axis.bin(i).upper() << "): " << size << "\n";
	}
}

template <typename T>
struct Bin {
	double start;
	double end;

	T value;

	Bin(double s, double e, T v) {
		start = s;
		end = e;
		value = v;
	}

	double width() {
		return end - start;
	}
};

template <typename T>
void print_with_precision(T value, int precision, bool newline = true) {
	cout << std::setprecision(precision) << value;
	if (newline) {
		cout << "\n";
	}
	cout << std::setprecision(3);
}

template <typename T>
void print_bins(std::vector<Bin<T>> bins) {
	for (auto bin : bins) {
		const T value = bin.value;
		cout << "[" << bin.start << ", " << bin.end << "): ";
		print_with_precision(value, 8);
	}
}
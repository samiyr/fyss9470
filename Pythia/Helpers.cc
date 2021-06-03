#ifndef HELPERS_H
#define HELPERS_H

#include "Pythia8/Pythia.h"
#include <algorithm>
#include <string>
#include <utility>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Pythia8;

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
std::vector<double> find_pTs(std::vector<Particle> particles) {
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


template <typename T>
void print_with_precision(T value, int precision, bool newline = true) {
	cout << std::setprecision(precision) << value;
	if (newline) {
		cout << "\n";
	}
	cout << std::setprecision(3);
}


template <typename T>
std::vector<T> range(T lower, T upper, T step = 1) {
	const size_t size = (upper - lower) / step;
	std::vector<T> v(size);

	T new_val = lower;
	std::generate(begin(v), end(v), [step, &new_val]() {
		T current = new_val;
		new_val += step;
		return current;
	});
	return v;
}

// https://stackoverflow.com/a/17299623
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v) {
        total_size += sub.size();
    }
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v) {
        result.insert(result.end(), sub.begin(), sub.end());
    }
    return result;
}

#endif // HELPERS_H
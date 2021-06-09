#ifndef HELPERS_H
#define HELPERS_H

#include "Pythia8/Pythia.h"
#include <numeric>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "ParticleContainer.cc"

using namespace Pythia8;

template <typename T>
bool contains(std::vector<T> vec, T element) {
	return std::find(vec.begin(), vec.end(), element) != vec.end();
}

template <typename T>
struct Range {
	T start;
	T end;
	Range() {
		start = T(0);
		end = T(0);
	}
	Range(T s, T e) {
		start = s;
		end = e;
	}
	T width() const {
		return end - start;
	}
	T center() const {
		return (start + end) / 2;
	}
	bool in_range(T v) const {
		return (v >= start && v < end);
	}
	std::string extent() const {
		return "[" + std::to_string(start) + ", " + std::to_string(end) + ")";
	}
};
template <typename T>
struct OptionalRange {
	std::optional<T> start;
	std::optional<T> end;
	OptionalRange() {
		start = std::nullopt;
		end = std::nullopt;
	}
	OptionalRange(std::optional<T> s, std::optional<T> e) {
		start = s;
		end = e;
	}
	std::optional<T> width() const {
		if (start && end) {
			return *end - *start;
		}
		return std::nullopt;
	}
	std::optional<T> center() const {
		if (start && end) {
			return (*start + *end) / 2;
		}
	}
	bool in_range(T v) const {
		if (start && v < *start) {
			return false;
		}
		if (end && v >= *end) {
			return false;
		}
		return true;
	}
	std::string extent() const {
		std::string lower;
		if (start) {
			lower = "[" + std::to_string(*start);
		} else {
			lower = "(-∞";
		}
		std::string middle = ", ";
		std::string upper;
		if (end) {
			upper = std::to_string(*end) + ")";
		} else {
			upper = "∞)";
		}
		return lower + middle + upper;
	}
};

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

std::vector<double> find_azimuths(std::vector<ParticleContainer> particles) {
	std::vector<double> phis;
	for (ParticleContainer particle : particles) {
		phis.push_back(particle.particle.phi());
	}
	return phis;
}
std::vector<double> find_pTs(std::vector<ParticleContainer> particles) {
	std::vector<double> phis;
	for (ParticleContainer particle : particles) {
		phis.push_back(particle.particle.pT());
	}
	return phis;
}
std::vector<double> find_pT_hats(std::vector<ParticleContainer> particles) {
	std::vector<double> phis;
	for (ParticleContainer particle : particles) {
		phis.push_back(particle.pT_hat);
	}
	return phis;
}
std::vector<double> find_event_weights(std::vector<ParticleContainer> particles) {
	std::vector<double> weights;
	for (ParticleContainer particle : particles) {
		weights.push_back(particle.event_weight);
	}
	return weights;
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
std::vector<T> create_range(T lower, T upper, T step = T(1)) {
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

template <typename T>
std::vector<T> fixed_range(T lower, T upper, int N) {
	const T step = (upper - lower) / N;
	std::vector<T> v(N + 1);

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

template <typename T>
std::vector<T> product(const std::vector<T> v1, const std::vector<T> v2) {
	std::vector<T> r;
	r.reserve(v1.size());
	std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(r), std::multiplies<T>());
	return r;
}

#endif // HELPERS_H
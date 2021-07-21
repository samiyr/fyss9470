#ifndef HELPERS_H
#define HELPERS_H

#include "Pythia8/Pythia.h"
#include <numeric>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "ParticleContainer.cc"
#include "Beam.cc"

using namespace Pythia8;

/**
 * 
 * Miscellaneous types
 * 
 */

/**
 * Types of normalization to apply in data analysis.
 */
enum class Normalization {
	/// No normalization
	None, 
	/// Normalize by integral
	Unity, 
	/// Normalize by event count
	Count,

	STARC
};
/**
 * Enum encapsulating the strategy with multiparton interactions.
 */
enum class MPIStrategy {
	/// No MPI
	Disabled, 
	/// Use Pythia's MPI model
	PythiaMPI, 
	/// Use an analytic DPS model
	DPS
};

enum class Process {
	HardQCD,
	SoftQCDNonDiffractive
};

template <typename T>
bool contains(std::vector<T> *vec, T element) {
	for (auto &v : *vec) {
		if (v == element) {
			return true;
		}
	}
	return false;
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
	bool complement;

	OptionalRange() {
		start = std::nullopt;
		end = std::nullopt;
		complement = false;
	}
	OptionalRange(std::optional<T> s, std::optional<T> e, bool c = false) {
		start = s;
		end = e;
		complement = c;
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
		const bool ans = !((start && v < *start) || (end && v >= *end));
		if (complement) {
			return !ans;
		} else {
			return ans;
		}
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
	static bool disjoint(OptionalRange<T> a, OptionalRange<T> b) {
		const auto start_A = a.start ? *a.start : -std::numeric_limits<T>::infinity();
		const auto end_A = a.end ? *a.end : std::numeric_limits<T>::infinity();
		const auto start_B = b.start ? *b.start : -std::numeric_limits<T>::infinity();
		const auto end_B = b.end ? *b.end : std::numeric_limits<T>::infinity();

		const bool comparison = end_A <= start_B || end_B <= start_A;
		const bool either_complemented = a.complement || b.complement;
		if (either_complemented) {
			return !comparison;
		} else {
			return comparison;
		}
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
		phis.push_back(particle.phi);
	}
	return phis;
}
std::vector<double> find_pTs(std::vector<ParticleContainer> particles) {
	std::vector<double> phis;
	for (ParticleContainer particle : particles) {
		phis.push_back(particle.pT);
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
template <typename T>
T variance(std::vector<T> v) {
	const T m = mean(v);
	T accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const T d) {
	    accum += (d - m) * (d - m);
	});
	return accum / (v.size() - 1);
}
template<typename T>
T standard_deviation(std::vector<T> v) {
	return sqrt(variance(v));
}
template <typename T>
T standard_error_of_mean(std::vector<T> v) {
	return sqrt(variance(v) / v.size());
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

std::string bool_to_string(bool in) {
	return in ? "true" : "false";
}


std::string to_string(Normalization norm) {
	switch(norm) {
		case Normalization::None:
			return "none";
			break;
		case Normalization::Unity:
			return "unity";
			break;
		case Normalization::Count:
			return "count";
			break;
		case Normalization::STARC:
			return "STARC";
			break;
	}
}
std::string to_string(Process in) {
	switch(in) {
		case Process::HardQCD:
			return "HardQCD";
			break;
		case Process::SoftQCDNonDiffractive:
			return "SoftQCDNonDiffractive";
			break;
	}
}
std::string to_string(MPIStrategy mpi) {
	switch(mpi) {
		case MPIStrategy::Disabled:
			return "disabled";
			break;
		case MPIStrategy::PythiaMPI:
			return "PythiaMPI";
			break;
		case MPIStrategy::DPS:
			return "DPS";
			break;
	}
}
#endif // HELPERS_H
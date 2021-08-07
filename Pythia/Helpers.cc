#ifndef HELPERS_H
#define HELPERS_H

#include "Pythia8/Pythia.h"
#include <numeric>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "ParticleContainer.cc"
#include "Beam.cc"
#include <optional>
#include <filesystem>

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
/**
 * Types of collision processes in Pythia.
 */
enum class Process {
	/// Enables HardQCD:all in Pythia
	HardQCD,
	/// Enables SoftQCD:nonDiffractive in Pythia
	SoftQCDNonDiffractive
};

/**
 *
 * List operations
 * 
 */

template <typename T>
/// Checks if a vector contains some element.
/// Stops immediately if an element is found.
bool contains(std::vector<T> *vec, T element) {
	for (auto &v : *vec) {
		if (v == element) {
			return true;
		}
	}
	return false;
}
/// Returns a vector of all particles in a given Pythia event.
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
/// Returns a vector of values corresponding to the azimuthal angles (phi) of the given particles.
std::vector<double> find_azimuths(std::vector<ParticleContainer> particles) {
	std::vector<double> phis;
	for (ParticleContainer particle : particles) {
		phis.push_back(particle.phi);
	}
	return phis;
}
/// Returns a vector of values corresponding to the transverse momenta (pT) of the given particles.
std::vector<double> find_pTs(std::vector<ParticleContainer> particles) {
	std::vector<double> pTs;
	for (ParticleContainer particle : particles) {
		pTs.push_back(particle.pT);
	}
	return pTs;
}
/// Returns a vector of values corresponding to the event weights associated with given particle containers.
std::vector<double> find_event_weights(std::vector<ParticleContainer> particles) {
	std::vector<double> weights;
	for (ParticleContainer particle : particles) {
		weights.push_back(particle.event_weight);
	}
	return weights;
}

template <typename T>
/// Calculates the sum of the values in the input vector.
T sum(std::vector<T> v) {
	T sum = std::accumulate(std::begin(v), std::end(v), T(0.0));
	return sum;
}
template <typename T>
/// Calculates the pointwise product of two vectors. Assumes that the vectors
/// have the same length. Violation of this assumption is not checked and has unknown consequences.
std::vector<T> product(const std::vector<T> v1, const std::vector<T> v2) {
	std::vector<T> r;
	r.reserve(v1.size());
	std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(r), std::multiplies<T>());
	return r;
}
template <typename T>
/// Calculates the arithmetic mean of the values in the input vector.
T mean(std::vector<T> v) {
	return sum(v) / v.size();
}
template <typename T>
/// Calculates the variance of the values in the input vector. This implementation
/// divides by N - 1.
T variance(std::vector<T> v) {
	const T m = mean(v);
	T accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const T d) {
	    accum += (d - m) * (d - m);
	});
	return accum / (v.size() - 1);
}
template<typename T>
/// Calculates the standard deviation of the values in the input vector.
T standard_deviation(std::vector<T> v) {
	return sqrt(variance(v));
}
template <typename T>
/// Calculates the standard error of the mean, given by standard_deviation / sqrt(N), 
/// of the values in the input vector.
T standard_error_of_mean(std::vector<T> v) {
	return sqrt(variance(v) / v.size());
}

template <typename T>
/// Flattens a 2D matrix to a 1D vector.
/// Taken from // https://stackoverflow.com/a/17299623
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

/**
 *
 * Range types
 * 
 */

template <typename T>
/**
 * Specifies a finite range [a, b).
 */
struct Range {
	/// Lower end of the range, a.
	T start;
	/// Upper end of the range, b.
	T end;
	/// Initializes an empty range, setting both endpoints to zero.
	Range() : start(T(0)), end(T(0)) { }
	/// Initializes a range with the given endpoints.
	/// Assumes that start <= end. This is not checked
	/// and a violation leads to unknown behaviour.
	Range(T s, T e) : start(s), end(e) { }
	/// Returns the width of the range.
	T width() const {
		return end - start;
	}
	/// Returns the center point of the range.
	T center() const {
		return (start + end) / 2;
	}
	/// Checks whether a given value is within the range.
	/// Returns true if value == start and false if value == end.
	bool in_range(T v) const {
		return (v >= start && v < end);
	}
	/// Returns a textual representation of the range.
	std::string extent() const {
		return "[" + std::to_string(start) + ", " + std::to_string(end) + ")";
	}
};
template <typename T>
/**
 * Specifies a possibly infinite range [a, b).
 */
struct OptionalRange {
	/// Lower end point of the range, a. 
	/// A nullopt value represents negative infinity,
	/// a range with no lower bound.
	std::optional<T> start;
	/// Upper end point of the range, b.
	/// A nullopt value represents positive infinity,
	/// a range with no upper bound.
	std::optional<T> end;
	/// Whether the range should be interpreted as
	/// the complement of [a, b), i.e. (-inf, a) union [b, inf).
	bool complement;
	/// Initializes the range as the real axis (-inf, inf), setting both endpoints to
	/// nullopt and complement to false.
	OptionalRange() : start(std::nullopt), end(std::nullopt), complement(false) {}
	/// Initializes the range with the given endpoints. Third argument, complement, defaults to false.
	OptionalRange(std::optional<T> s, std::optional<T> e, bool c = false) : start(s), end(e), complement(c) {}
	/// Returns the width of the range. Return value of nullopt indicates an infinite range.
	std::optional<T> width() const {
		if (start && end) {
			return *end - *start;
		}
		return std::nullopt;
	}
	/// Returns the center point of the range. Return value of nullopt indicates that no center point exists,
	/// i.e. the range is infinite.
	std::optional<T> center() const {
		if (start && end) {
			return (*start + *end) / 2;
		}
		return std::nullopt;
	}
	/// Checks whether a given value is within the range.
	/// Returns true if value == start and false if value == end.
	bool in_range(T v) const {
		const bool ans = !((start && v < *start) || (end && v >= *end));
		if (complement) {
			return !ans;
		} else {
			return ans;
		}
	}
	/// Returns a textual representation of the range.
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

template <typename T>
/// Creates a list of evenly spaced values from lower to upper in steps of step, 
/// which by default is 1.
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
/// Creates a list of N evenly spaced values from lower to upper.
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

/**
 *
 *
 * Printing
 * 
 * 
 */

/// Prints the Pythia event to stdout.
void print_event(Event event) {
	for (int i = 0; i < event.size(); ++i) {
		const Particle particle = event[i];
		cout << particle.name() << "\n";
	}
}

template <typename T>
/// Prints a value to stdout with the given precision and by default, appends a newline.
void print_with_precision(T value, int precision, bool newline = true) {
	cout << std::setprecision(precision) << value;
	if (newline) {
		cout << "\n";
	}
	cout << std::setprecision(3);
}

/**
 *
 *
 * Type conversions
 * 
 * 
 */

/// Returns the string 'true' or 'false' based on the input boolean.
std::string to_string(bool in) {
	return in ? "true" : "false";
}
/// Returns the name of the given normalization scheme.
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
	return "";
}
/// Returns the name of the given Pythia process type.
std::string to_string(Process in) {
	switch(in) {
		case Process::HardQCD:
			return "HardQCD";
			break;
		case Process::SoftQCDNonDiffractive:
			return "SoftQCDNonDiffractive";
			break;
	}
	return "";
}
/// Returns the name of the given MPI strategy.
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
	return "";
}
/// Given an optional string, returns the string if it exists and the empty string if not.
std::string convert_optional_string(std::optional<std::string> s) {
	return s ? *s : "";
}

/**
 *
 *
 * Filesystem
 * 
 * 
 */

/// Constructs a std::filesystem::path based on the parent directory, filename and extensions.
/// By default, creates the parent directory if it doesn't exists.
/// For example, the call construct_path("foo/bar", "file", "txt") gives a path representation of
/// 'foo/bar/file.txt'.
std::filesystem::path construct_path(
	std::optional<std::filesystem::path> parent_directory, 
	std::string filename, 
	std::string extension, 
	bool create_parent = true) {
	std::filesystem::path path;
	if (parent_directory) {
		path = *parent_directory / filename;

	} else {
		path = filename;
	}
	
	path.replace_extension(extension);

	if (create_parent) {
		const auto directory = path.parent_path();
		std::filesystem::create_directories(directory);
	}

	return path;
}

#endif // HELPERS_H
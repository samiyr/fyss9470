#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "Helpers.cc"
#include <iostream>

/**
 * Bin represents a single histogram bin. 
 * Unlike usual bins, Bin keeps track of the binned items. 
 * For a more traditional bin, see RangedContainer.
 * Stores content with generic type T.
 */
template <typename T>
struct Bin {
	/// Extent of the bin, [range.start, range.end).
	Range<double> range;

	/// Stored contents of the bin.
	std::vector<T> contents;
	/// Weights associated with the contents.
	std::vector<T> weights;
	/// Calculates the sum of the weights.
	T weight_sum() {
		return std::accumulate(weights.begin(), weights.end(), typename std::vector<T>::value_type(0.0));
	}
	/// Initializes a bin, with no content by default.
	Bin(double l, double u, std::vector<T> c = {}, std::vector<T> w = {}) : range(l, u) {
		contents = c;
		weights = w;
	}
	/// Returns the number of items stored.
	auto size() {
		return contents.size();
	}
	/// Provides access to the underlaying contents.
	T& operator[](int index) {
		return contents[index];
	}
};

/**
 * A traditional bin, with an extent and a value.
 */
template <typename T>
struct RangedContainer {
	/// Extent of the bin, [range.start, range.end).
	Range<double> range;
	/// Value associated with the bin.
	T value;
	/// Initializes a bin.
	RangedContainer(double l, double u, T v) : range(l, u) {
		value = v;
	}
};

/**
 * Represents a histogram of RangedContainer types. 
 * Bins only store a single value.
 */
template <typename T>
class ValueHistogram {
public:
	/// Bins stored in the histogram.
	std::vector<RangedContainer<T>> containers;
	/// Initializes a histogram, with no content by default.
	ValueHistogram(std::vector<RangedContainer<T>> c = {}) {
		containers = c;
	}
	/// Initializes a histogram and reserves capacity.
	ValueHistogram(int capacity) {
		containers.reserve(capacity);
	}

	ValueHistogram(std::vector<double> points) {
		std::vector<RangedContainer<T>> c;
		for (typename std::vector<T>::size_type i = 0; i < points.size() - 1; i++) {
			RangedContainer<T> container(points[i], points[i + 1], T(0));
			c.push_back(container);
		}
		containers = c;
	}
	/// Returns the number of bins.
	auto size() const {
		return containers.size();
	}
	/// Provides access to the underlaying bins.
	RangedContainer<T> &operator[](int index) {
		return containers[index];
	}
	/// Appends a value `v` to the appropriate bin. 
	/// If no bin is found, discards the value silently.
	void fill(double v) {
		for (auto &container : containers) {
			if (container.range.in_range(v)) {
				container.value += T(1);
				break;
			}
		}
	}
	/// Prints the histogram to stdout.
	void print() const {
		for (auto container : containers) {
			const T value = container.value;
			cout << container.range.extent() << ": ";
			print_with_precision(value, 8);
		}
	}

	void print_with_bars() const {
		T total = T(0);
		for (auto &container : containers) {
			total += container.value;
		}
		for (auto &container : containers) {
			const T value = container.value;
			const T percentage = ((double)value / (double)total) * 100;
			const int count = round(percentage);
			cout << container.range.extent() << ": ";
			for (int i = 0; i < count; i++) {
				cout << "#";
			}
			cout << "\n";
		}
	}
	/// Exports the histogram to a file `filename`.
	void export_histogram(std::string filename, int precision = 12) const {
		ofstream file;
		file.open(filename);
		file << std::setprecision(precision);
		for (auto container : containers) {
			const double center = container.range.center();
			const T value = container.value;
			file << center << "," << value << "\n";
		}
		file.close();
	}
};

/**
 * Represents a numbered histogram of Bin types.
 * Bins store the original contents.
 */
template <typename T>
class Histogram {
public:
	/// Bins stored in the histogram.
	std::vector<Bin<T>> bins;
	/// The points at which bins change, start and end points included.
	std::vector<double> axis;
	/// Returns the lower edge of the histogram.
	double lower_edge() const {
		return axis.front();
	}
	/// Returns the upper edge of the histogram.
	double upper_edge() const {
		return axis.back();
	}
	/// Initializes a histogram with `count` bins, extending from `lower` to `upper`.
	Histogram(int count, double lower, double upper) {
		const double step = (upper - lower) / count;
		const std::vector<double> points = create_range(lower, upper, step);
		Construct(points);
	}
	/// Initializes a histogram with explicit bin boundary points.
	Histogram(std::vector<double> points) {
		Construct(points);
	}
	/// Appends a value `v` with weight `w` to the appropriate bin. 
	/// If no bin is found, discards the value silently.
	void fill(T v, T w) {
		for (auto &bin : bins) {
			if (bin.range.in_range(v)) {
				bin.contents.push_back(v);
				bin.weights.push_back(w);
				break;
			}
		}
	}
	/// Appends a collection of values `v` with weights `w` to the appropriate bins.
	/// Calls `fill` for the individual values.
	void fill(std::vector<T> v, std::vector<T> w) {
		for (typename std::vector<T>::size_type i = 0; i < v.size(); i++) {
			fill(v[i], w[i]);
		}
	}
	/// Appends a collection of values `v` to the appropriate bins.
	/// Calls `fill` for the individual values.
	void fill(std::vector<T> v) {
		for (typename std::vector<T>::size_type i = 0; i < v.size(); i++) {
			fill(v[i], 1);
		}
	}
	/// Prints the contents of the bins to stdout.
	void print() const {
		std::for_each(std::begin(bins), std::end(bins), [](Bin<T>&bin) {
			std::cout << "[" << bin.lower << ", " << bin.upper << "): " << bin.size() << "\n";
		});
	}
	/// Normalizes the histogram to cross section and returns a `ValueHistogram`. 
	ValueHistogram<T> normalize(double total_weight, double sigma, OptionalRange<double> rapidity, bool use_weights) const {
		ValueHistogram<T> normalized;
		for (auto bin : bins) {
			const double y_min = rapidity.start.has_value() ? *rapidity.start : 0.0;
			const double y_max = rapidity.end.has_value() ? *rapidity.end : 2 * M_PI;
			const double dy = y_max - y_min;
			const double dpT = bin.range.width();
			// dsigma = 1 / 2π * N / N_ev * sigma / dy pT dpT
			double dsigma = 0.0;
			for (typename std::vector<T>::size_type i = 0; i < bin.contents.size(); i++) {
				const double pT = bin.contents[i];
				double weight_factor = 1.0;
				if (use_weights) {
					weight_factor = bin.weights[i];
				}
				dsigma += weight_factor / pT;
			}
			dsigma *= sigma / (2 * M_PI * total_weight * dy * dpT);
			RangedContainer<T> container(bin.range.start, bin.range.end, dsigma);
			normalized.containers.push_back(container);
		}
		return normalized;
	}

	ValueHistogram<T> export_to_values() const {
		std::vector<RangedContainer<T>> containers;
		ValueHistogram<T> hist;
		for (auto bin : bins) {
			const auto value = bin.size();
			RangedContainer<T> container(bin.range.start, bin.range.end, value);
			hist.containers.push_back(container);
		}
		return hist;
	}
private:
	void Construct(std::vector<double> points) {
		axis = points;
		std::vector<Bin<T>> b;
		for (std::vector<double>::size_type i = 0; i < points.size() - 1; i++) {
			Bin<T> bin(points[i], points[i + 1]);
			b.push_back(bin);
		}
		bins = b;
	}
};

template <typename T>
ValueHistogram<T> combine(std::vector<ValueHistogram<T>> containers) {
	auto reference = containers.front();
	const auto N = reference.size();
	ValueHistogram<T> result(N);
	for (typename std::vector<RangedContainer<T>>::size_type i = 0; i < N; i++) {
		const double lower = reference[i].range.start;
		const double upper = reference[i].range.end;

		T value = 0;

		for (typename std::vector<RangedContainer<T>>::size_type j = 0; j < containers.size(); j++) {
			value += containers[j][i].value;
		}

		const RangedContainer<T> container(lower, upper, value);
		result.containers.push_back(container);
	}
	return result;
}


#endif // HISTOGRAM_H
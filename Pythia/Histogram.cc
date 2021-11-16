#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "Helpers.cc"
#include "Around.cc"
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
	Bin(double l, double u, std::vector<T> c = {}, std::vector<T> w = {}) : range(l, u), contents(c), weights(w) {}
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
	Around<T> value;
	/// The number of hits in the bin. Starts from zero and is incremented by one every time the value changes.
	unsigned long long hit_count;
	/// Initializes a bin with the given lower and upper bounds and the value. By default, hit count is set to zero
	RangedContainer(double l, double u, T v, unsigned long long h = 0) : range(l, u), value(v), hit_count(h) {}
	/// Initializes a bin with the given lower and upper bounds and the value with uncertainty.
	/// By default, hit count is set to zero
	RangedContainer(double l, double u, Around<T> v, unsigned long long h = 0) : range(l, u), value(v), hit_count(h) {}
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
	ValueHistogram(std::vector<RangedContainer<T>> c = {}) : containers(c) {}
	/// Initializes a histogram and reserves capacity.
	ValueHistogram(int capacity) {
		containers.reserve(capacity);
	}
	/// Initializes a histogram with bin edges given by the input vector.
	/// For example, input of {0.0, 1.0, 2.0} corresponds to a histogram with
	/// two bins, the first [0.0, 1.0) and the second [1.0, 2.0).
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
	/// Returns a list of the widths of each bin.
	auto widths() const {
		std::vector<double> result;
		for (auto &container : containers) {
			result.push_back(container.range.width());
		}
		return result;
	}
	/// Appends a value `v` to the appropriate bin. 
	/// If no bin is found, discards the value silently.
	void fill(double v, T w = T(1)) {
		for (auto &container : containers) {
			if (container.range.in_range(v)) {
				container.value += w;
				container.hit_count += 1;
				break;
			}
		}
	}
	/// Prints the histogram to stdout.
	void print() const {
		for (auto container : containers) {
			const Around<T> value = container.value;
			cout << container.range.extent() << ": ";
			print_with_precision(value, 8);
		}
	}
	/// Defines a formatting for the histogram when outputting to streams.
	friend ostream& operator<<(ostream& os, ValueHistogram<T> const & hist) {
		os << "\n";
		T total = T(0);
		for (auto &container : hist.containers) {
			total += container.value.value;
		}
		for (auto &container : hist.containers) {
			const Around<T> value = container.value;
			const T percentage = ((double)value.value / (double)total) * 100;
			const int count = round(percentage);
			os << container.range.extent() << ": " << "(" << value << ") [" << container.hit_count << "]\t";
			for (int i = 0; i < count; i++) {
				os << "#";
			}
			os << "\n";
		}
		return os << "\n";
    }
    /// Prints the histogram to stdout.
	void print_with_bars() const {
		cout << *this;
	}
	/// Exports the histogram to a file filename. By default, precision is set to 12.
	void export_histogram(std::string filename, int precision = 12) const {
		ofstream file;
		file.open(filename);
		file << std::setprecision(precision);
		for (auto container : containers) {
			const double center = container.range.center();
			const Around<T> value = container.value;
			const T error = value.error ? *value.error : T(0);
			const auto hit_count = container.hit_count;
			file << center << "," << value.value << "," << error << "," << hit_count << "\n";
		}
		file.close();
	}
	/// Returns the total value of the histogram across all bins.
	double total() const {
		return std::accumulate(containers.begin(), containers.end(), 0.0, [](double partial, RangedContainer<T> current) {
			return partial + current.value.value;
		});
	}

	template<typename F>
	/// Applies a transformation to the histogram. The transformation, given by a lambda function,
	/// is called with each container and must return a new value based on the container.
	static ValueHistogram<T> transform(ValueHistogram<T> input, F lambda) {
		ValueHistogram<T> result;
		for (auto container : input.containers) {
			const auto new_value = lambda(container);
			RangedContainer<T> new_container(container.range.start, container.range.end, new_value, container.hit_count);
			result.containers.push_back(new_container);
		}
		return result;
	}
	/// Returns a histogram normalized to unity, i.e. the integral of the normalized histogram equals 1.
	ValueHistogram<double> normalize_to_unity() const {
		const double sum_total = total();
		const auto normalized = ValueHistogram<double>::transform(*this, [sum_total](RangedContainer<double> container) {
			const double factor = container.range.width() * sum_total;
			const Around<double> new_value = factor == 0 ? 0 : container.value / factor;
			return new_value;
		});
		return normalized;
	}
	/// Returns a histogram normalized with a procedure outlined in the STAR presentation. 
	/// Shouldn't be used.
	ValueHistogram<double> normalize_to_star_C(double N_assoc) const {
		ValueHistogram<double> normalized;
		for (auto container : containers) {
			const double delta_phi = container.range.width();
			const double factor = N_assoc * delta_phi;
			const Around<double> new_value = factor == 0 ? 0 : container.value / factor;
			RangedContainer<double> new_container(container.range.start, container.range.end, new_value, container.hit_count);
			normalized.containers.push_back(new_container);
		}
		return normalized;
	}
	ValueHistogram<double> divide_by_bin_width() const {
		const auto result = ValueHistogram<double>::transform(*this, [](RangedContainer<double> container) {
			const double factor = container.range.width();
			const Around<double> new_value = factor == 0 ? 0 : container.value / factor;
			return new_value;
		});
		return result;
	}
	template <typename V>
	/// Returns a histogram with the contents divided by a constant.
	ValueHistogram<double> normalize_by(V constant) const {
		return *this / constant;
	}
	/// Combines a list of histograms into one histogram by taking the sum across each histogram for all bins.
	static ValueHistogram<T> combine(std::vector<ValueHistogram<T>> _containers, bool calculate_error = false) {
		auto reference = _containers.front();
		const auto N = reference.size();
		ValueHistogram<T> result(N);
		for (typename std::vector<RangedContainer<T>>::size_type i = 0; i < N; i++) {
			const double lower = reference[i].range.start;
			const double upper = reference[i].range.end;

			std::vector<T> values;
			unsigned long long total_hits = 0;

			for (typename std::vector<RangedContainer<T>>::size_type j = 0; j < _containers.size(); j++) {
				const auto c = _containers[j][i];
				const auto v = c.value;
				values.push_back(v.value);
				total_hits += c.hit_count;
			}
			const auto value = sum(values);
			const auto error = standard_error_of_mean(values);

			const Around<T> c_value = calculate_error ? Around<T>(value, error) : Around<T>(value);

			const RangedContainer<T> container(lower, upper, c_value, total_hits);
			result.containers.push_back(container);
		}
		return result;
	}
	/// Calculates the estimated statistical error and erases all other errors. The error is given by
	/// 1 / sqrt(N), where N is the hit count of a bin. This corresponds with a relative error.
	static ValueHistogram<double> calculate_statistical_error(ValueHistogram<double> hist) {
		const auto transformed = ValueHistogram<double>::transform(hist, [](RangedContainer<double> container) {
			const auto value = container.value.value;
			const auto hit_count = container.hit_count;
			const auto error = 1.0 / sqrt((double)hit_count);
			return Around<double>(value, error);
		});
		return transformed;
	}
	/**
	 * Arithmetic operations for histograms.
	 */
	ValueHistogram<T>& operator+=(ValueHistogram<T> rhs) {
		*this = *this + rhs;
		return *this;
	}
	ValueHistogram<T> operator+(ValueHistogram<T>& rhs) {
		ValueHistogram<T> lhs = *this;
		return combine({lhs, rhs});
	}
	ValueHistogram<double>& operator+=(double constant) {
		for (typename std::vector<RangedContainer<double>>::size_type i = 0; i < size(); i++) {
			const auto container = containers[i];
			const Around<double> new_value = container.value + constant;
			RangedContainer<double> new_container(container.range.start, container.range.end, new_value, container.hit_count);
			containers[i] = new_container;
		}
		return *this;
	}
	ValueHistogram<double>& operator+=(Around<double> constant) {
		for (typename std::vector<RangedContainer<double>>::size_type i = 0; i < size(); i++) {
			const auto container = containers[i];
			const Around<double> new_value = container.value + constant;
			RangedContainer<double> new_container(container.range.start, container.range.end, new_value, container.hit_count);
			containers[i] = new_container;
		}
		return *this;
	}
	ValueHistogram<double> operator+(double constant) const {
		ValueHistogram<double> h = *this;
		h *= constant;
		return h;
	}
	ValueHistogram<double>& operator*=(double constant) {
		for (typename std::vector<RangedContainer<double>>::size_type i = 0; i < size(); i++) {
			const auto container = containers[i];
			const Around<double> new_value = container.value * constant;
			RangedContainer<double> new_container(container.range.start, container.range.end, new_value, container.hit_count);
			containers[i] = new_container;
		}
		return *this;	
	}
	ValueHistogram<double>& operator*=(Around<double> constant) {
		for (typename std::vector<RangedContainer<double>>::size_type i = 0; i < size(); i++) {
			const auto container = containers[i];
			const Around<double> new_value = container.value * constant;
			RangedContainer<double> new_container(container.range.start, container.range.end, new_value, container.hit_count);
			containers[i] = new_container;
		}
		return *this;
	}
	ValueHistogram<double> operator*(double constant) const {
		ValueHistogram<double> h = *this;
		h *= constant;
		return h;
	}

	ValueHistogram<double>& operator/=(double constant) {
		*this *= (1.0 / constant);
		return *this;
	}
	ValueHistogram<double> operator/(double constant) const {
		ValueHistogram<double> h = *this;
		h /= constant;
		return h;
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
			const auto size = bin.contents.size();
			for (typename std::vector<T>::size_type i = 0; i < size; i++) {
				const double pT = bin.contents[i];
				double weight_factor = 1.0;
				if (use_weights) {
					weight_factor = bin.weights[i];
				}
				dsigma += weight_factor / pT;
			}
			dsigma *= sigma / (2 * M_PI * total_weight * dy * dpT);
			RangedContainer<T> container(bin.range.start, bin.range.end, dsigma, size);
			normalized.containers.push_back(container);
		}
		return normalized;
	}
	/// Converts the histogram to a ValueHistogram by taking the number of elements as the size of the bin.
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
	/// Histogram constructor with bin edges given as a list.
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



#endif // HISTOGRAM_H
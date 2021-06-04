#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <algorithm>
#include <execution>
#include "Helpers.cc"
#include <iostream>

template <typename T>
struct Bin {
	Range<double> range;

	std::vector<T> contents;
	std::vector<T> weights;
	T weight_sum() {
		return std::accumulate(weights.begin(), weights.end(), 0.0);
	}

	Bin(double l, double u, std::vector<T> c = {}, std::vector<T> w = {}) : range(l, u) {
		contents = c;
		weights = w;
	}

	auto size() {
		return contents.size();
	}
	T& operator[](int index) {
		return contents[index];
	}
};

template <typename T>
struct RangedContainer {
	Range<double> range;

	T value;

	RangedContainer(double l, double u, T v) : range(l, u) {
		value = v;
	}
};

template <typename T>
class ValueHistogram {
public:
	std::vector<RangedContainer<T>> containers;

	ValueHistogram(std::vector<RangedContainer<T>> c = {}) {
		containers = c;
	}
	ValueHistogram(int capacity) {
		containers.reserve(capacity);
	}

	typename std::vector<RangedContainer<T>>::size_type size() const {
		return containers.size();
	}
	RangedContainer<T> &operator[](int index) {
		return containers[index];
	}
};

template <typename T>
class Histogram {
public:
	std::vector<Bin<T>> bins;
	std::vector<double> axis;

	double lower_edge() {
		return axis.front();
	}
	double upper_edge() {
		return axis.back();
	}

	Histogram(int count, double lower, double upper) {
		const double step = (upper - lower) / count;
		const std::vector<double> points = create_range(lower, upper, step);
		Construct(points);
	}
	Histogram(std::vector<double> points) {
		Construct(points);
	}

	void fill(T v, T w) {
		for (auto &bin : bins) {
			if (bin.range.in_range(v)) {
				bin.contents.push_back(v);
				bin.weights.push_back(w);
				break;
			}
		}
	}
	void fill(std::vector<T> v, std::vector<T> w) {
		for (typename std::vector<T>::size_type i = 0; i < v.size(); i++) {
			fill(v[i], w[i]);
		}
	}

	void print() {
		std::for_each(std::begin(bins), std::end(bins), [](Bin<T>&bin) {
			std::cout << "[" << bin.lower << ", " << bin.upper << "): " << bin.size() << "\n";
		});
	}

	ValueHistogram<T> normalize(int count, double sigma, Range<double> rapidity, bool use_weights = true) {
	std::vector<RangedContainer<T>> normalized;
	for (auto bin : bins) {
		const double dy = rapidity.width();
		const double dpT = bin.range.width();
		double dsigma = 0.0;
		double w_sum = 1.0;
		if (use_weights) {
			w_sum = bin.weight_sum();
		}
		for (typename std::vector<T>::size_type i = 0; i < bin.contents.size(); i++) {
			const double pT = bin.contents[i];
			double weight_factor = 1.0;
			if (use_weights) {
				const double w = bin.weights[i];
				weight_factor = w / w_sum;
			}
			dsigma += weight_factor * sigma / (2 * M_PI * count * dy * dpT * pT);
		}
		RangedContainer<T> container(bin.range.start, bin.range.end, dsigma);
		normalized.push_back(container);
	}
	return normalized;
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
void print_containers(ValueHistogram<T> hist) {
	for (auto container : hist.containers) {
		const T value = container.value;
		cout << container.range.extent() << ": ";
		print_with_precision(value, 8);
	}
}
template <typename T>
void export_containers(ValueHistogram<T> hist, std::string filename, int precision = 12) {
	ofstream file;
	file.open(filename);
	file << std::setprecision(precision);
	for (auto container : hist.containers) {
		const T center = container.range.center();
		const T value = container.value;
		file << center << "," << value << "\n";
	}
	file.close();
}

#endif // HISTOGRAM_H
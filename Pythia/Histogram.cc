#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <algorithm>
#include <execution>
#include "Helpers.cc"
#include <iostream>

template <typename T>
struct Bin {
	double lower;
	double upper;

	std::vector<T> contents;

	Bin(double l, double u, std::vector<T> c = {}) {
		lower = l;
		upper = u;
		contents = c;
	}

	auto size() {
		return contents.size();
	}
	double width() {
		return upper - lower;
	}
	double center() {
		return (lower + upper) / 2;
	}
	T& operator[](int index) {
		return contents[index];
	}
	bool in_range(T v) {
		return (v >= lower && v < upper);
	}
	std::string extent() {
		return "[" + std::to_string(lower) + ", " + std::to_string(upper) + ")";
	}
};

template <typename T>
struct RangedContainer {
	double lower;
	double upper;

	T value;

	RangedContainer(double l, double u, T v) {
		lower = l;
		upper = u;
		value = v;
	}

	double width() {
		return upper - lower;
	}
	double center() {
		return (lower + upper) / 2;
	}
	bool in_range(T v) {
		return (v >= lower && v < upper);
	}
	std::string extent() {
		return "[" + std::to_string(lower) + ", " + std::to_string(upper) + ")";
	}
};

template <typename T>
class Histogram {
public:
	std::vector<Bin<T>> bins;
	std::vector<double> axis;

	double lower_edge() {
		return axis.front().lower;
	}
	double upper_edge() {
		return axis.back().upper;
	}

	Histogram(int count, double lower, double upper) {
		const double step = (upper - lower) / count;
		const std::vector<double> points = range(lower, upper, step);
		Construct(points);
	}
	Histogram(std::vector<double> points) {
		Construct(points);
	}

	void fill(T v) {
		for (auto &bin : bins) {
			if (bin.in_range(v)) {
				bin.contents.push_back(v);
				break;
			}
		}
	}
	void fill(std::vector<T> v) {
		std::for_each(std::begin(v), std::end(v), [this](T&value) {
			fill(value);
		});
	}

	void print() {
		std::for_each(std::begin(bins), std::end(bins), [](Bin<T>&bin) {
			std::cout << "[" << bin.lower << ", " << bin.upper << "): " << bin.size() << "\n";
		});
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
void print_containers(std::vector<RangedContainer<T>> containers) {
	for (auto container : containers) {
		const T value = container.value;
		cout << container.extent() << ": ";
		print_with_precision(value, 8);
	}
}
template <typename T>
void export_containers(std::vector<RangedContainer<T>> containers, std::string filename, int precision = 12) {
	ofstream file;
	file.open(filename);
	file << std::setprecision(precision);
	for (auto container : containers) {
		const T center = container.center();
		const T value = container.value;
		file << center << "," << value << "\n";
	}
	file.close();
}

#endif // HISTOGRAM_H
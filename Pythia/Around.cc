#ifndef AROUND_H
#define AROUND_H

#include "Helpers.cc"
#include <cassert>
#include <optional>

template <typename T>
/// A wrapper for a value with an associated error. When used with arithmetic operations,
/// errors are combined using error propagation.
struct Around {
	/// The actual value.
	T value;
	/// Uncertainty associated with the value. A value of nullopt represents no uncertainty.
	std::optional<T> error;
	/// Constructs an Around object with the given value. Second optional argument specifies
	/// the associated uncertainty.
	Around(T v, std::optional<T> e = std::nullopt) : value(v), error(e) {}
	/// Constructs an Around object based on a list of values. The list must be non-empty.
	/// Violation of this assumption causes an assertion. The value is taken to be the mean
	/// of the values in the list and the uncertainty is given by the standard error the mean.
	Around(std::vector<T> v) {
		assert(v.size() != 0);
		value = mean(v);
		error = standard_error_of_mean(v);
	}
	/// Specifies the format for outputting to streams. For objects with uncertainty,
	/// the formatting is 'value ± uncertainty', while for objects with no (nullopt) uncertainty,
	/// the formatting 'value' is used.
	friend ostream& operator<<(ostream& os, Around<T> const & around) {
		if (around.error) {
			return os << around.value << " ± " << *around.error;
		} else {
			return os << around.value;
		}
    }
    /// Overloads the += operator, for expressions of Around += constant.
    Around<T>& operator+=(T c) {
		value += c;
		return *this;
	}
};

/**
 * Rules of arithmetic based on error propagation formulas.
 */

template <typename T>
const Around<T> operator+(const Around<T>& lhs, const T& rhs) {
	return Around<T>(lhs.value + rhs, *lhs.error);
}
template <typename T>
const Around<T> operator+(const T& lhs, const Around<T>& rhs) {
	return rhs + lhs;
}
template <typename T>
const Around<T> operator+(const Around<T>& lhs, const Around<T>& rhs) {
	const T x = lhs.value;
	const T dx = lhs.error ? *lhs.error : T(0);
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x + y, sqrt(dx * dx + dy * dy));
}

template <typename T>
const Around<T> operator-(const Around<T>& lhs, const Around<T>& rhs) {
	const T x = lhs.value;
	const T dx = lhs.error ? *lhs.error : T(0);
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x - y, sqrt(dx * dx + dy * dy));
}

template <typename T>
const Around<T> operator*(const Around<T>& lhs, const T& rhs) {
	if (lhs.error) {
		return Around<T>(lhs.value * rhs, *lhs.error * rhs);
	} else {
		return Around<T>(lhs.value * rhs);
	}
}
template <typename T>
const Around<T> operator*(const T& lhs, const Around<T>& rhs) {
	return rhs * lhs;
}
template <typename T>
const Around<T> operator*(const Around<T>& lhs, const Around<T>& rhs) {
	const T x = lhs.value;
	const T dx = lhs.error ? *lhs.error : T(0);
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x * y, sqrt(dx * dx * y * y + dy * dy * x * x));
}

template <typename T>
const Around<T> operator/(const Around<T>& lhs, const T& rhs) {
	if (lhs.error) {
		return Around<T>(lhs.value / rhs, *lhs.error / rhs);
	} else {
		return Around<T>(lhs.value / rhs);
	}
}
template <typename T>
const Around<T> operator/(const T& lhs, const Around<T>& rhs) {
	const T x = lhs;
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x / y, x * dy / (y * y));
}
template <typename T>
const Around<T> operator/(const Around<T>& lhs, const Around<T>& rhs) {
	const T x = lhs.value;
	const T dx = lhs.error ? *lhs.error : T(0);
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x / y, sqrt(dx * dx / (y * y) + x * x * dy * dy / (y * y * y * y)));
}

#endif // AROUND_H
#ifndef AROUND_H
#define AROUND_H

#include "Helpers.cc"

template <typename T>
struct Around {
	T value;
	std::optional<T> error;

	Around(T _value, std::optional<T> _error = std::nullopt) {
		value = _value;
		error = _error;
	}
	Around(std::vector<T> v) {
		assert(v.size() != 0);
		value = mean(v);
		error = standard_error_of_mean(v);
	}

	friend ostream& operator<<(ostream& os, Around<T> const & around) {
		if (around.error) {
			return os << around.value << " ± " << *around.error;
		} else {
			return os << around.value;
		}
    }
    Around<T>& operator+=(T c) {
		value += c;
		return *this;
	}
};

template <typename T>
Around<T> operator+(Around<T> lhs, const T& rhs) {
	return Around<T>(lhs.value + rhs, *lhs.error);
}
template <typename T>
Around<T> operator+(T lhs, const Around<T>& rhs) {
	return rhs + lhs;
}
template <typename T>
Around<T> operator+(const Around<T>& lhs, Around<T> rhs) {
	const T x = lhs.value;
	const T dx = lhs.error ? *lhs.error : T(0);
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x + y, sqrt(dx * dx + dy * dy));
}

template <typename T>
Around<T> operator*(Around<T> lhs, const T& rhs) {
	if (lhs.error) {
		return Around<T>(lhs.value * rhs, *lhs.error * rhs);
	} else {
		return Around<T>(lhs.value * rhs);
	}
}
template <typename T>
Around<T> operator*(T lhs, const Around<T>& rhs) {
	return rhs * lhs;
}
template <typename T>
Around<T> operator*(const Around<T>& lhs, Around<T> rhs) {
	const T x = lhs.value;
	const T dx = lhs.error ? *lhs.error : T(0);
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x * y, sqrt(dx * dx * y * y + dy * dy * x * x));
}
template <typename T>
Around<T> operator/(Around<T> lhs, const T& rhs) {
	if (lhs.error) {
		return Around<T>(lhs.value / rhs, *lhs.error / rhs);
	} else {
		return Around<T>(lhs.value / rhs);
	}
}
template <typename T>
Around<T> operator/(Around<T> lhs, const Around<T>& rhs) {
	const T x = lhs.value;
	const T dx = lhs.error ? *lhs.error : T(0);
	const T y = rhs.value;
	const T dy = rhs.error ? *rhs.error : T(0);
	return Around<T>(x / y, sqrt(dx * dx / (y * y) + x * x * dy * dy / (y * y * y * y)));
}

#endif // AROUND_H
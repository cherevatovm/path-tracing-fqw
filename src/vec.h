#ifndef VEC_H
#define VEC_H
#include <math.h>
#include "utility.h"

constexpr float eps = std::numeric_limits<float>::epsilon();

struct Vec {
	double x, y, z;

	Vec(double x_ = 0, double y_ = 0, double z_ = 0): x(x_), y(y_), z(z_) {}
	void print() const { printf("%f, %f, %f\n", x, y, z); }
	Vec operator+(const Vec& v) const { return Vec(x + v.x, y + v.y, z + v.z); }
	Vec& operator+=(const Vec& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	Vec& operator-=(const Vec& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	Vec operator-(const Vec& v) const { return Vec(x - v.x, y - v.y, z - v.z); }
	Vec operator-() const { return Vec(-x, -y, -z); }
	Vec operator*(double scalar) const { return Vec(x * scalar, y * scalar, z * scalar); }
	Vec& operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}
	Vec operator%(const Vec& v) const {
		return Vec(
			diff_of_products(y, v.z, z, v.y),
			diff_of_products(z, v.x, x, v.z),
			diff_of_products(x, v.y, y, v.x));
	}

	double operator[](int i) const { return (i == 0) ? x : (i == 1) ? y : z; }
	double dot_prod(const Vec& v) const { return x * v.x + y * v.y + z * v.z; }
	Vec mult(const Vec& v) const { return Vec(x * v.x, y * v.y, z * v.z); }
	Vec permute(const Vec& new_pos) const { return Vec((*this)[(int)new_pos[0]], (*this)[(int)new_pos[1]], (*this)[(int)new_pos[2]]); }
	double length_sqrd() const { return x * x + y * y + z * z; }
	double length() const { return sqrt(length_sqrd()); }

	Vec& norm() {
		double coef = 1 / length();
		x *= coef;
		y *= coef;
		z *= coef;
		return *this;
	}

	int min_comp_ind() {
		double min_val = std::min(x, std::min(y, z));
		if (std::abs(x - min_val) < eps)
			return 0;
		else if (std::abs(y - min_val) < eps)
			return 1;
		else
			return 2;
	}

	int max_comp_ind() {
		double max_val = std::max(x, std::max(y, z));
		if (std::abs(x - max_val) < eps)
			return 0;
		else if (std::abs(y - max_val) < eps)
			return 1;
		else
			return 2;
	}

	inline static Vec min(const Vec& v1, const Vec& v2) {
		return Vec(std::min(v1.x, v2.x),
			std::min(v1.y, v2.y),
			std::min(v1.z, v2.z));
	}

	inline static Vec max(const Vec& v1, const Vec& v2) {
		return Vec(std::max(v1.x, v2.x),
			std::max(v1.y, v2.y),
			std::max(v1.z, v2.z));
	}

	inline static Vec abs(const Vec& v) { return Vec(std::abs(v.x), std::abs(v.y), std::abs(v.z)); }

	void rotate_y(double alpha, const Vec& center) {
		double xc = x - center.x;
		double zc = z - center.z;
		double nx = xc * cos(alpha) - zc * sin(alpha);
		x = nx + center.x;
		double nz = xc * sin(alpha) + zc * cos(alpha);
		z = nz + center.z;
	}
};

void create_orthonorm_sys(const Vec& v1, Vec& v2, Vec& v3) {
	// Projection to y = 0 plane and normalized orthogonal vector construction
	if (std::abs(v1.x) > std::abs(v1.y)) {
		double inv_len = 1.0 / sqrt(v1.x * v1.x + v1.z * v1.z);
		v2 = Vec(-v1.z * inv_len, 0.0, v1.x * inv_len);
	}
	// Projection to x = 0 plane and normalized orthogonal vector construction
	else {
		double inv_len = 1.0 / sqrt(v1.y * v1.y + v1.z * v1.z);
		v2 = Vec(0.0, v1.z * inv_len, -v1.y * inv_len);
	}
	v3 = v1 % v2;
}

Vec sample_hemisphere(double u1, double u2) {
	const double r = sqrt(1.0 - u1 * u1);
	const double phi = 2 * M_PI * u2;
	return Vec(cos(phi) * r, sin(phi) * r, u1);
}

inline double cos_theta(const Vec& w) { return w.z; }
inline double abs_cos_theta(const Vec& w) { return std::abs(w.z); }
inline double cos_2_theta(const Vec& w) { return w.z * w.z; }
inline double sin_2_theta(const Vec& w) { return std::max(1 - cos_2_theta(w), 0.0); }
inline double sin_theta(const Vec& w) { return std::sqrt(sin_2_theta(w)); }
inline double tan_theta(const Vec& w) { return sin_theta(w) / cos_theta(w); }
inline double tan_2_theta(const Vec& w) { return sin_2_theta(w) / cos_2_theta(w); }

#endif

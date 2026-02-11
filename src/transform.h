#ifndef TRANSFORM_H
#define TRANSFORM_H
#include "square_matrix.h"
#include "vec.h"
#include "bounds.h"
#include <math.h>

class Transform {
private:
	SquareMatrix<4> matr, inv_matr;

public:
	Transform() = default;

	Transform(const SquareMatrix<4>& m, const SquareMatrix<4>& inv_m)
		: matr(m), inv_matr(inv_m) {
	}

	bool operator==(const Transform& t) const { return t.matr == matr; }

	Vec apply_transform(const Vec& v, bool is_point = true) const {
		double x_new = matr.contents[0][0] * v.x + matr.contents[0][1] * v.y + matr.contents[0][2] * v.z + matr.contents[0][3];
		double y_new = matr.contents[1][0] * v.x + matr.contents[1][1] * v.y + matr.contents[1][2] * v.z + matr.contents[1][3];
		double z_new = matr.contents[2][0] * v.x + matr.contents[2][1] * v.y + matr.contents[2][2] * v.z + matr.contents[2][3];
		double w_new = matr.contents[3][0] * v.x + matr.contents[3][1] * v.y + matr.contents[3][2] * v.z + matr.contents[3][3];
		Vec res(x_new, y_new, z_new);
		if (abs(w_new - 1) > eps && abs(w_new) > eps && is_point)
			res *= (1 / w_new);
		return res;
	}

	Vec transform_normal(const Vec& n) const {
		return Vec(inv_matr.contents[0][0] * n.x + inv_matr.contents[1][0] * n.y + inv_matr.contents[2][0] * n.z,
			inv_matr.contents[0][1] * n.x + inv_matr.contents[1][1] * n.y + inv_matr.contents[2][1] * n.z,
			inv_matr.contents[0][2] * n.x + inv_matr.contents[1][2] * n.y + inv_matr.contents[2][2] * n.z);
	}

	Bounds3 transform_bounds(const Bounds3& aabb) const {
		Bounds3 aabb_t;
		for (int i = 0; i < 8; ++i)
			aabb_t = Bounds3::find_union(aabb_t, apply_transform(aabb.get_corner(i)));
		return aabb_t;
	}

	const SquareMatrix<4>& get_matr() const { return matr; }
	const SquareMatrix<4>& get_inverse_matr() const { return inv_matr; }

	static Transform translate(const Vec& delta) {
		SquareMatrix<4> m = { { 1, 0, 0, delta.x },
							 { 0, 1, 0, delta.y },
							 { 0, 0, 1, delta.z },
							 { 0, 0, 0, 1 } };
		SquareMatrix<4> m_inv = { { 1, 0, 0, -delta.x },
								  { 0, 1, 0, -delta.y },
								  { 0, 0, 1, -delta.z },
								  { 0, 0, 0, 1 } };
		return Transform(m, m_inv);
	}

	static Transform scale(double x, double y, double z) {
		SquareMatrix<4> m = { { x, 0, 0, 0 },
							  { 0, y, 0, 0 },
							  { 0, 0, z, 0 },
							  { 0, 0, 0, 1 } };
		SquareMatrix<4> m_inv = { { 1 / x, 0, 0, 0 },
								  { 0, 1 / y, 0, 0 },
								  { 0, 0, 1 / z, 0 },
								  { 0, 0, 0, 1 } };
		return Transform(m, m_inv);
	}

	static Transform rotate_around_axis(double theta, int rotate_index) {
		double theta_rad = theta * (M_PI / 180.0);
		double sin_theta = std::sin(theta_rad);
		if (abs(sin_theta) < eps) sin_theta = 0;
		double cos_theta = std::cos(theta_rad);
		if (abs(cos_theta) < eps) cos_theta = 0;
		switch (rotate_index) {
		case 0:
		{
			SquareMatrix<4> m = { { 1, 0, 0, 0 },
								  { 0, cos_theta, -sin_theta, 0 },
								  { 0, sin_theta, cos_theta, 0 },
								  { 0, 0, 0, 1 } };
			return Transform(m, m.transpose());
		}
		case 1:
		{
			SquareMatrix<4> m = { { cos_theta, 0, sin_theta, 0 },
								  { 0, 1, 0, 0 },
								  { -sin_theta, 0, cos_theta, 0 },
								  { 0, 0, 0, 1 } };
			return Transform(m, m.transpose());
		}
		default:
		{
			SquareMatrix<4> m = { { cos_theta, -sin_theta, 0, 0 },
								  { sin_theta, cos_theta, 0, 0 },
								  { 0, 0, 1, 0 },
								  { 0, 0, 0, 1 } };
			return Transform(m, m.transpose());
		}
		}
	}

	static Transform rotate_around_point(const Vec& point, double theta, int rotate_index) {
		Transform translate_to = translate(-point);
		Transform rot_axis = rotate_around_axis(theta, rotate_index);
		SquareMatrix<4> rot_by_center = translate_to.get_inverse_matr() * rot_axis.get_matr() *
			translate_to.get_matr();
		SquareMatrix<4> rot_by_center_inv = translate_to.get_inverse_matr() * rot_axis.get_inverse_matr() *
			translate_to.get_matr();
		return Transform(rot_by_center, rot_by_center_inv);
	}

	static Transform scale_around_point(const Vec& point, double x, double y, double z) {
		Transform translate_to = translate(-point);
		Transform sc = scale(x, y, z);
		SquareMatrix<4> scale_by_center = translate_to.get_inverse_matr() * sc.get_matr() * translate_to.get_matr();
		SquareMatrix<4> scale_by_center_inv = translate_to.get_inverse_matr() * sc.get_inverse_matr() *
			translate_to.get_matr();
		return Transform(scale_by_center, scale_by_center_inv);
	}
};

#endif
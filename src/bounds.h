#ifndef AABB_H
#define AABB_H
#include <vector>
#include <limits>
#include "vec.h"

class Bounds3 {
private:
	Vec p_min, p_max;

public:
	Bounds3() {
		constexpr double min = std::numeric_limits<double>::lowest();
		constexpr double max = std::numeric_limits<double>::max();
		p_min = Vec(max, max, max);
		p_max = Vec(min, min, min);
	}

	Bounds3(const Vec& _min, const Vec& _max): p_min(_min), p_max(_max) {}

	bool is_empty() const {
		return p_min.x >= p_max.x || p_min.y >= p_max.y || p_min.z >= p_max.z;
	}
	
	bool is_degen() const {
		return p_min.x > p_max.x || p_min.y > p_max.y || p_min.z > p_max.z;
	}

	Vec operator[](int i) const { return (i == 0) ? p_min : p_max; }
	Vec& operator[](int i) { return (i == 0) ? p_min : p_max; }

	Vec centroid() const { return (p_min + p_max) * 0.5; }

	Vec get_corner(int i) const {
		return Vec((*this)[(i & 1)].x,
			(*this)[(i & 2) ? 1 : 0].y,
			(*this)[(i & 4) ? 1 : 0].z);
	}

	int largest_dim() const {
		Vec diag = p_max - p_min;
		if (diag.x > diag.y && diag.x > diag.z)
			return 0;
		else if (diag.y > diag.z)
			return 1;
		return 2;
	}

	bool intersect(const Ray& r, const Vec& reverse_dir, double t_max = LDBL_MAX, 
		double* hit_t0 = nullptr, double* hit_t1 = nullptr) const {
		double t0 = 0, t1 = t_max;

		for (int i = 0; i < 3; ++i) {
			double t_near = (p_min[i] - r.orig[i]) * reverse_dir[i];
			double t_far = (p_max[i] - r.orig[i]) * reverse_dir[i];

			if (t_near > t_far)
				std::swap(t_near, t_far);

			t0 = t_near > t0 ? t_near : t0;
			t1 = t_far < t1 ? t_far : t1;
			if (t0 > t1)
				return false;
		}

		if (hit_t0)
			*hit_t0 = t0;
		if (hit_t1)
			*hit_t1 = t1;

		return true;
	}

	static inline Bounds3 find_union(const Bounds3& aabb, const Vec& p) {
		return Bounds3(Vec::min(aabb.p_min, p),
			Vec::max(aabb.p_max, p));
	}

	static inline Bounds3 find_union(const Bounds3& aabb1, const Bounds3& aabb2) {
		return Bounds3(Vec::min(aabb1.p_min, aabb2.p_min), 
			Vec::max(aabb1.p_max, aabb2.p_max));
	}

	static inline Bounds3 find_box_intersect(const Bounds3& aabb1, const Bounds3& aabb2) {
		return Bounds3(Vec::max(aabb1.p_min, aabb2.p_min),
			Vec::min(aabb1.p_max, aabb2.p_max));
	}
};

#endif
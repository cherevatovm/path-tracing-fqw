#ifndef INTERSECTION_H
#define INTERSECTION_H
#include "shape.h"

class Shape;

class Intersection {
public:
	double t;
	const Shape* object = nullptr;
	Vec hit_point, baryc_coords;

	Intersection() = default;
	Intersection(double t_, const Shape* object_ = nullptr, const Vec& baryc_coords_ = Vec()) :
		t(t_), object(object_), baryc_coords(baryc_coords_) {
	}

	void calc_inters_point(const Ray& r) { hit_point = r.orig + r.dir * t; }
	bool operator<(const Intersection& inters) const { return t < inters.t; }
};

#endif
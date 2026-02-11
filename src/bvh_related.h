#ifndef BVH_H
#define BVH_H
#include "shape.h"
#include "bounds.h"

class BoundingInfo {
public:
	Shape* obj;	
	Bounds3 aabb;
	Vec centroid;

	BoundingInfo(Shape* obj_) :
		obj(obj_), aabb(obj->get_bounds()), centroid(aabb.centroid()) {};

	BoundingInfo(Shape* obj_, const Bounds3& aabb_) : 
		obj(obj_), aabb(aabb_), centroid(aabb.centroid()) {};
};

class BVH_Node {
public:
	Bounds3 aabb;	
	uint8_t obj_count, split_axis;

	union {
		uint32_t first_obj_ind;	// For leaf nodes: index of the first obj in a leaf
		uint32_t second_child_ind;	// For interior nodes
	};
};

class CentroidComparator {
private:
	size_t dim;

public:
	CentroidComparator(size_t dim_) : dim(dim_) {};

	bool operator()(const BoundingInfo& first, const BoundingInfo& second) const {
		return first.centroid[dim] < second.centroid[dim];
	}
};

#endif
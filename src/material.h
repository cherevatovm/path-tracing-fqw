#ifndef MATERIAL_H
#define MATERIAL_H
#include "vec.h"

enum Refl_type { DIFF, SPEC, REFR, ROUGH, ROUGH_DIEL };

struct Material {
	Vec color;
	Refl_type brdf;
	double refr_ind;
	double roughness;
	double metal_prop;

	Material(const Vec& color_ = Vec(0.7, 0.7, 0.7), Refl_type brdf_ = DIFF,
		double refr_ind_ = 1 + 1e-4, double roughness_ = 1e-4, double metal_prop_ = 1e-4) :
		brdf(brdf_), color(color_), 
		refr_ind(refr_ind_), roughness(roughness_), metal_prop(metal_prop_) {}
};

Material diffuse = {
	Vec(0.7, 0.7, 0.7),
	DIFF
};
Material mirror = {
	Vec(0.85, 0.85, 0.85),
	SPEC
};
Material glass = {
	Vec(0.85, 0.85, 0.85),
	REFR,
	1.5
};
Material frosted_glass = {
	Vec(1, 1, 1),
	ROUGH_DIEL,
	1.5,
	0.5
};
Material ruby = {
	Vec(0.87, 0, 0.13),
	ROUGH_DIEL,
	1.77
};
Material gold = {
	Vec(1.0, 0.86, 0.57),
	ROUGH,
	0.2773,
	0.1,
	0.5
};
Material iron = {
	Vec(0.77, 0.78, 0.78),
	ROUGH,
	2.9304,
	0.3,
	0.5
};

std::vector<Material> materials = { diffuse, mirror, glass, frosted_glass, ruby, gold, iron };

#endif

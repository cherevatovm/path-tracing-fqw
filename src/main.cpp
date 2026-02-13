#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define _USE_MATH_DEFINES

#include <GL/glew.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <locale.h>
#include <cstdint>
#include <ctime>
#include <chrono>
#include <random>
#include <vector>

#include "vec.h"
#include "ray.h"
#include "shape.h"
#include "material.h"
#include "microfacets.h"
#include "bvh_related.h"
#include "my_random.h"
#include "utility.h"

#include "ImFileDialog.h"
#include "stb_image_write.h"

using namespace std::chrono;

static void glfw_error_callback(int error, const char* description) {
	fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

int width = 1024, height = 768, samples;
Vec* image;
std::string ui_text = "";
std::vector<int> model_choices = { 0, 0 };
std::vector<int> material_choices = { 0, 0 };
std::vector<double> roughness_choices = { 0.0, 0.0 };

class Scene {
private:
	std::vector<BoundingInfo> all_bounding_data;
	std::vector<BVH_Node> bv_hierarchy;
	std::vector<Shape*> bounded_objects;
	std::vector<Shape*> unbounded_objects;
	std::vector<Sphere*> light_sources;
	std::vector<Mesh*> meshes;

	Intersection intersect_scene(const Ray& r) {
		Intersection closest_inters(LDBL_MAX);
		for (Shape* obj : unbounded_objects) {
			Intersection inters = obj->intersect(r);
			if (inters.object && inters < closest_inters)
				closest_inters = inters;
		}

		if (!bv_hierarchy.empty()) {
			uint32_t stack[64];
			uint32_t cur_stack_size = 0;
			uint32_t node_num = 0;
			Vec reverse_dir(1 / r.dir.x, 1 / r.dir.y, 1 / r.dir.z);
			while (true) {			
				BVH_Node& node = bv_hierarchy[node_num];
				if (node.aabb.intersect(r, reverse_dir, closest_inters.t)) {
					if (node.obj_count > 0) {	
						// Leaf
						for (int i = 0; i < node.obj_count; ++i) {
							Intersection inters = bounded_objects[node.first_obj_ind + i]->intersect(r);
							if (inters.object && inters < closest_inters)
								closest_inters = inters;
						}
						if (cur_stack_size == 0)
							break;
						node_num = stack[--cur_stack_size];
					}
					else {	
						// Interior node
						if (r.dir[node.split_axis] < 0) {
							stack[cur_stack_size++] = node_num + 1;
							node_num = node.second_child_ind;
						}
						else {
							stack[cur_stack_size++] = node.second_child_ind;
							node_num = node_num + 1;
						}
					}
				}
				else {
					if (cur_stack_size == 0)
						break;
					node_num = stack[--cur_stack_size];
				}
			}
		}

		return closest_inters;
	}

	void split_bounds(const std::vector<BoundingInfo>::iterator begin, 
		const std::vector<BoundingInfo>::iterator end, uint8_t max_leaf_obj_num) {
		bv_hierarchy.push_back(BVH_Node());
		auto cur_node = bv_hierarchy.end();
		--cur_node;

		uint32_t obj_cnt = end - begin;
		if (obj_cnt <= max_leaf_obj_num) {
			// Make cur_node a leaf
			cur_node->obj_count = obj_cnt;
			for (auto i = begin; i != end; ++i)
				cur_node->aabb = Bounds3::find_union(cur_node->aabb, i->aabb);
			cur_node->first_obj_ind = begin - all_bounding_data.begin();
		}
		else {
			// Make cur_node an interior node
			cur_node->obj_count = 0;
			Bounds3 centroid_aabb;
			for (auto i = begin; i != end; ++i)
				centroid_aabb = Bounds3::find_union(centroid_aabb, i->centroid);

			cur_node->split_axis = centroid_aabb.largest_dim();
			auto median = begin + (end - begin) / 2;
			std::nth_element(begin, median, end, CentroidComparator(cur_node->split_axis));

			// Recursively call for left and right child
			size_t first_child_ind = bv_hierarchy.size();
			split_bounds(begin, median, max_leaf_obj_num);
			cur_node->second_child_ind = bv_hierarchy.size();
			split_bounds(median, end, max_leaf_obj_num);
			
			cur_node->aabb = Bounds3::find_union(bv_hierarchy[first_child_ind].aabb, bv_hierarchy[cur_node->second_child_ind].aabb);
		}
	}

	Vec path_tracing(const Ray& r, int depth, unsigned short* Xi, int E = 1) {
		#ifdef _DEBUG
			const int max_depth = 80;
		#else
			const int max_depth = 1000;
		#endif
		if (depth > max_depth) return Vec();
		Intersection inters = intersect_scene(r);
		if (!inters.object) return Vec(); // Return black if there is no intersection 

		inters.calc_inters_point(r);
		Vec n = inters.object->get_normal(inters);
		Vec nl = n.dot_prod(r.dir) < 0 ? n : -n;
		Vec color = inters.object->material.color;
		double rr_prob = std::max(color.x, std::max(color.y, color.z));

		// Russian Roulette for path termination
		if (++depth > 5 || !rr_prob) {
			if (erand48(Xi) < rr_prob)
				color = color * (1 / rr_prob);
			else
				return inters.object->emis * E;
		}

		// Diffuse BRDF
		if (inters.object->material.brdf == DIFF) {
			double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi);
			double r2_sqrt = sqrt(r2);

			Vec w = nl, u, v;
			create_orthonorm_sys(w, u, v);
			Vec d = (u * cos(r1) * r2_sqrt + v * sin(r1) * r2_sqrt + w * sqrt(1 - r2)).norm();

			// For hemisphere sampling
			/*
			Vec samp_dir_ = sample_hemisphere(erand48(Xi), erand48(Xi));
			Vec d;
			d.x = Vec(u.x, v.x, w.x).dot_prod(samp_dir_);
			d.y = Vec(u.y, v.y, w.y).dot_prod(samp_dir_);
			d.z = Vec(u.z, v.z, w.z).dot_prod(samp_dir_);
			d.norm();
			*/

			Vec e;
			for (int i = 0; i < light_sources.size(); ++i) {
				// Create orthonormal coord system and sample direction by solid angle
				Vec sw = (light_sources[i]->center - inters.hit_point).norm(), su, sv;
				create_orthonorm_sys(sw, su, sv);

				double cos_a_max = sqrt(1 - light_sources[i]->radius * light_sources[i]->radius /
					(inters.hit_point - light_sources[i]->center).dot_prod(inters.hit_point - light_sources[i]->center));
				double eps1 = erand48(Xi), eps2 = erand48(Xi);
				double cos_a = 1 - eps1 + eps1 * cos_a_max;
				double sin_a = sqrt(1 - cos_a * cos_a);
				double phi = 2 * M_PI * eps2;

				// For hemisphere sampling
				/*
				Vec rand = sample_hemisphere(erand48(Xi), erand48(Xi));
				Vec samp_dir;
				samp_dir.x = Vec(su.x, sv.x, sw.x).dot_prod(rand);
				samp_dir.y = Vec(su.y, sv.y, sw.y).dot_prod(rand);
				samp_dir.z = Vec(su.z, sv.z, sw.z).dot_prod(rand);
				samp_dir.norm();
				*/

				Vec samp_dir = (su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a).norm();

				// shoot shadow rays
				if (intersect_scene(Ray(inters.hit_point, samp_dir)).object == light_sources[i]) {
					double omega = 2 * M_PI * (1 - cos_a_max);
					e = e + color.mult(light_sources[i]->emis * samp_dir.dot_prod(nl) * omega) * (1 / M_PI); // 1/PI for BRDF
				}
			}

			return inters.object->emis * E + e + color.mult(path_tracing(Ray(inters.hit_point, d), depth, Xi, 0));
		}
		// Specular BRDF
		else if (inters.object->material.brdf == SPEC)
			return inters.object->emis + color.mult(path_tracing(Ray(inters.hit_point, r.dir - n * n.dot_prod(r.dir) * 2), depth, Xi)); // Angle of incidence == angle of reflection
		else if (inters.object->material.brdf == ROUGH) {
			double roughness = inters.object->material.roughness;
			double refr_ind = inters.object->material.refr_ind;

			Vec wo = -r.dir;
			Vec wm = sample_wm(Xi, nl, roughness);
			Vec wi = wm * wm.dot_prod(wo) * 2 - wo;

			if (nl.dot_prod(wi) <= 0) return inters.object->emis;

			Vec diffuse = color * (1 - inters.object->material.metal_prop);
			double Fr = fresnel_schlick(std::max(wm.dot_prod(wo), 0.0),
				((refr_ind - 1) * (refr_ind - 1)) / ((refr_ind + 1) * (refr_ind + 1)));
			Vec specular = path_tracing(Ray(inters.hit_point, wi), depth, Xi) *
				brdf_div_by_pdf(nl, wo, wi, wm, roughness, Fr);

			return inters.object->emis + diffuse.mult(specular);
		}
		else if (inters.object->material.brdf == ROUGH_DIEL) {
			double roughness = inters.object->material.roughness;
			double refr_ind1 = 1, refr_ind2 = inters.object->material.refr_ind;
			double F0 = ((refr_ind2 - refr_ind1) * (refr_ind2 - refr_ind1)) /
				((refr_ind2 + refr_ind1) * (refr_ind2 + refr_ind1));

			Vec wo = -r.dir;
			Vec wm = sample_wm(Xi, nl, roughness);
			Vec wi_refl = wm * wm.dot_prod(wo) * 2 - wo;
			
			bool into = wm.dot_prod(n) > 0;
			double refr_ratio = into ? refr_ind1 / refr_ind2 : refr_ind2 / refr_ind1;
			double cos_incid_angle = r.dir.dot_prod(wm);
			double cos2t = 1 - refr_ratio * refr_ratio * (1 - cos_incid_angle * cos_incid_angle);
			double Fr = fresnel_schlick(std::max(wo.dot_prod(wm), 0.0), F0);

			if (cos2t < 0)
				return inters.object->emis + color.mult(path_tracing(Ray(inters.hit_point, wi_refl), depth, Xi));
			
			Vec wi_refr = (r.dir * refr_ratio - wm * ((into ? 1 : -1) *
				(cos_incid_angle * refr_ratio + sqrt(cos2t)))).norm();

			double Tr = 1 - Fr;
			double prob = 0.25 + 0.5 * Fr;
			double refl_prob = Fr / prob, trans_prob = Tr / (1 - prob);

			if (depth > 2) {
				if (erand48(Xi) < prob) {
					if (nl.dot_prod(wi_refl) <= 0) return inters.object->emis;
					return inters.object->emis + color.mult(path_tracing(Ray(inters.hit_point, wi_refl), depth, Xi)) * refl_prob *
						brdf_div_by_pdf(nl, wo, wi_refl, wm, roughness, Fr);
				}
				else {
					return inters.object->emis + color.mult(path_tracing(Ray(inters.hit_point, wi_refr), depth, Xi)) * trans_prob *
						btdf_div_by_pdf(wo, wi_refr, wm, roughness, Tr, refr_ratio);
				}
			}
			else {
				return inters.object->emis + color.mult(path_tracing(Ray(inters.hit_point, wi_refl), depth, Xi) * brdf_div_by_pdf(nl, wo, wi_refl, wm, roughness, Fr) +
					path_tracing(Ray(inters.hit_point, wi_refr), depth, Xi) * btdf_div_by_pdf(wo, wi_refr, wm, roughness, Tr, refr_ratio));
			}
		}
		
		// Refractive BRDF
		Ray refl_ray(inters.hit_point, r.dir - n * n.dot_prod(r.dir) * 2);
		bool into = n.dot_prod(nl) > 0; // Check if a ray goes in or out 
		double refr_ind1 = 1, refr_ind2 = inters.object->material.refr_ind;
		double refr_ratio = into ? refr_ind1 / refr_ind2 : refr_ind2 / refr_ind1;
		double cos_incid_angle = r.dir.dot_prod(nl);
		double cos2t = 1 - refr_ratio * refr_ratio * (1 - cos_incid_angle * cos_incid_angle);

		if (cos2t < 0)
			return inters.object->emis + color.mult(path_tracing(refl_ray, depth, Xi)); // Total internal reflection 

		Vec t_dir = (r.dir * refr_ratio - n * ((into ? 1 : -1) * (cos_incid_angle * refr_ratio + sqrt(cos2t)))).norm();
		double a = refr_ind2 - refr_ind1;
		double b = refr_ind2 + refr_ind1;

		double F0 = a * a / (b * b);
		double Fr = fresnel_schlick((into ? -cos_incid_angle : t_dir.dot_prod(n)), F0);
		double Tr = 1 - Fr;

		double prob = 0.25 + 0.5 * Fr;
		double refl_prob = Fr / prob, trans_prob = Tr / (1 - prob);

		Vec result;
		if (depth > 2) {
			// Russian roulette
			if (erand48(Xi) < prob)
				result = path_tracing(refl_ray, depth, Xi) * refl_prob;
			else
				result = path_tracing(Ray(inters.hit_point, t_dir), depth, Xi) * trans_prob;
		}
		else
			result = path_tracing(refl_ray, depth, Xi) * Fr + path_tracing(Ray(inters.hit_point, t_dir), depth, Xi) * Tr;

		return inters.object->emis + color.mult(result);
	}

public:
	Scene() {
		light_sources = {
			new Sphere(5.5, Vec(50, 81.6 - 15.5, 90), Material(Vec(), DIFF), Vec(1, 1, 1) * 40),
			//new Sphere(2.75, Vec(70, 81.6 - 30.5, 81.6), Vec(0, 1, 1) * 20,  Vec(), DIFF),
			//new Sphere(1.375, Vec(20, 81.6 - 20, 81.6), Vec(1, 1, 0) * 40,  Vec(), DIFF)
		};
		
		// radius, position, material
		unbounded_objects = {
			new Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Material(Vec(0.75, 0.25, 0.25), DIFF)), // Left wall
			new Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Material(Vec(0.25, 0.25, 0.75), DIFF)), // Right wall
			new Sphere(1e5, Vec(50, 40.8, 1e5), Material(Vec(0.25, 0.75, 0.25), DIFF)), // Back wall
			new Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Material(Vec(0.75, 0.75, 0.75), DIFF)), // Front wall
			new Sphere(1e5, Vec(50, 1e5, 81.6), Material(Vec(0.75, 0.75, 0.75), DIFF)), // Floor
			new Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Material(Vec(0.75, 0.75, 0.75), DIFF)), // Ceiling
		};
	}

	~Scene() {
		for (Mesh* m : meshes)
			delete m;
		for (Shape* o : unbounded_objects)
			delete o;
		for (Shape* o : bounded_objects)
			delete o;
	}

	void add_object(Shape* obj) {
		Bounds3 aabb = obj->get_bounds();
		if (aabb.is_empty())
			unbounded_objects.push_back(obj);
		else
			all_bounding_data.push_back(BoundingInfo(obj, aabb));
	}

	void fill() {
		for (int i = 0; i < material_choices.size(); ++i) {
			Mesh* m = new Mesh(materials[material_choices[i]]);
			m->material.roughness = roughness_choices[i];
			meshes.push_back(m);

			int res = -1;
			switch (model_choices[i]) {
			case -1:
				res = 0;
				break;
			case 0: {
				res = 0;
				Sphere* sph = new Sphere(16.5, Vec(33, 16.5, 65), materials[material_choices[i]]);
				sph->material.roughness = roughness_choices[i];
				if (i != 0)
					sph->center = Vec(73, 16.5, 95);
				add_object(sph);
				break;
			}
			case 1:
				res = m->load_model("3d_models/cube.obj");
				if (res == 0) {
					m->scale(15, 15, 15);
					if (i == 0) {
						m->rotate(30, 1);
						m->translate(Vec(30, 15, 80));
					}
					else {
						m->rotate(-30, 1);
						m->translate(Vec(75, 15, 90));
					}
				}
				break;
			case 2:
				res = m->load_model("3d_models/pinetree.obj");
				if (res == 0) {
					m->scale(0.2, 0.2, 0.2, false);
					if (i == 0)
						m->translate(Vec(100, 20, 80));
					else {
						m->rotate(30, 1);
						m->translate(Vec(142, 20, 95));
					}
				}
				break;
			case 3:
				res = m->load_model("3d_models/dog.obj");
				if (res == 0) {
					if (i == 0) {
						m->rotate(20, 1);
						m->translate(Vec(30, 0, 60));
					}
					else {
						m->rotate(-20, 1);
						m->translate(Vec(70, 0, 70));
					}
				}
				break;
			case 4:
				res = m->load_model("3d_models/skull.obj");
				if (res == 0) {
					m->scale(5, 5, 5);
					if (i == 0) {
						m->rotate(30, 1);
						m->translate(Vec(30, 25, 80));
					}
					else {
						m->rotate(-15, 1);
						m->translate(Vec(65, 25, 100));
					}
				}
				break;
			case 5:
				res = m->load_model("3d_models/heart.obj");
				if (res == 0) {
					m->scale(2, 2, 2);
					m->rotate(-90, 0);
					if (i == 0) {
						m->rotate(15, 1);
						m->translate(Vec(35, 25, 60));
					}
					else {
						m->rotate(-25, 1);
						m->translate(Vec(72, 23, 85));
					}
				}
				break;
			case 6:
				res = m->load_model("3d_models/gear.obj");
				if (res == 0) {
					m->scale(4, 4, 4);
					if (i == 0) {
						m->rotate(30, 1);
						m->translate(Vec(32, 18, 75));
					}
					else {
						m->rotate(-35, 1);
						m->translate(Vec(72, 18, 100));
					}
				}
				break;
			case 7:
				res = m->load_model("3d_models/GLman02.obj");
				if (res == 0) {
					m->scale(0.4, 0.4, 0.4);
					if (i == 0) {
						m->rotate(30, 1);
						m->translate(Vec(28, -62, 85));
					}
					else {
						m->rotate(-45, 1);
						m->translate(Vec(75, -62, 95));
					}
				}
				break;
			case 8:
				res = m->load_model("3d_models/weird_sphere.obj");
				if (res == 0) {				
					m->scale(25, 25, 25);
					if (i == 0) {
						m->rotate(-45, 1);
						m->translate(Vec(32, 30, 80));
					}
					else {
						m->rotate(45, 1);
						m->translate(Vec(73, 30, 105));
					}
				}
				break;
			}
			if (res != 0) {
				Sphere* sph = new Sphere(16.5, Vec(33, 16.5, 65), materials[material_choices[i]]);
				sph->material.roughness = roughness_choices[i];
				if (i != 0)
					sph->center = Vec(73, 16.5, 95);
				add_object(sph);
			}

			for (auto* f : m->faces)
				add_object(f);
		}
		
		for (Sphere* l_s : light_sources)
			add_object(l_s);
	}

	void rebuild_bvh(uint8_t max_leaf_obj_num = 1) {
		if (max_leaf_obj_num <= 0)
			max_leaf_obj_num = 1;

		// Depth-first representation of bvh stored in bv_hierarchy vec
		bv_hierarchy.clear();
		bv_hierarchy.reserve(all_bounding_data.size() * 4);
		split_bounds(all_bounding_data.begin(), all_bounding_data.end(), max_leaf_obj_num);

		bounded_objects.clear();
		bounded_objects.reserve(all_bounding_data.size());
		for (BoundingInfo& b_inf : all_bounding_data)
			bounded_objects.push_back(b_inf.obj);
	}

	void render_with_supersampling() {
		Vec result;
		Ray camera(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());
		Vec camera_x = Vec(width * 0.5135 / height),
			camera_y = (camera_x % camera.dir).norm() * 0.5135;
		const int bar_width = 50; // Progress bar width
		//const double sigma = 1.5; // Standard deviation
#pragma omp parallel for schedule(dynamic, 1) private(result)       // Use OpenMP 
		for (int y = 0; y < height; y++) {                       // Go through image rows 
#pragma omp critical
			{
				float percent = 100.0f * y / (height - 1);
				int filled = static_cast<int>(percent * bar_width / 100.0f);
				filled = std::max(0, std::min(bar_width, filled));
				std::string bar(filled, '=');
				bar.append(bar_width - filled, ' ');
				fprintf(stderr, "\rRendering (%d samples per pixel): [%s] %5.2f%%", samples * 4, bar.c_str(), percent);
				fflush(stderr);
			}
			for (unsigned short x = 0, Xi[3] = { 0, 0, y * y * y }; x < width; x++) { // Go through image cols 
				// 2x2 subpixel rows 
				for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++)
					// 2x2 subpixel cols 
					for (int sx = 0; sx < 2; sx++, result = Vec()) {
						for (int s = 0; s < samples; s++) {
							// For Gaussian filter
							/*
							double r1 = erand48(Xi), r2 = erand48(Xi);
							double z1, z2;
							box_muller(r1, r2, z1, z2);
							double offset_x = z1 * sigma, offset_y = z2 * sigma;
							offset_x = std::clamp(offset_x, -1.0, 1.0);
							offset_y = std::clamp(offset_y, -1.0, 1.0);
							*/

							double r1 = 2 * erand48(Xi), offset_x = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
							double r2 = 2 * erand48(Xi), offset_y = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
							
							Vec d = camera_x * (((sx + 0.5 + offset_x) / 2 + x) / width - 0.5) +
								camera_y * (((sy + 0.5 + offset_y) / 2 + y) / height - 0.5) + camera.dir;
							result = result + path_tracing(Ray(camera.orig + d * 140, d.norm()), 0, Xi) * (1.0 / samples);
						}
						image[i] = image[i] + Vec(clamp01(result.x), clamp01(result.y), clamp01(result.z)) * 0.25;
					}
			}
		}

	}

	void render_wout_supersampling(int samples) {
		Vec result;
		Ray camera(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());
		Vec camera_x = Vec(width * 0.5135 / height),
			camera_y = (camera_x % camera.dir).norm() * 0.5135;
		const int bar_width = 50; // Progress bar width
		const double sigma = 1.5; // Standard deviation
//#pragma omp parallel for schedule(dynamic, 1) private(result)       // Use OpenMP 
		for (int y = 0; y < height; y++) {                       // Go through image rows 
//#pragma omp critical
			{
				float percent = 100.0f * y / (height - 1);
				int filled = static_cast<int>(percent * bar_width / 100.0f);
				filled = std::max(0, std::min(bar_width, filled));
				std::string bar(filled, '=');
				bar.append(bar_width - filled, ' ');
				fprintf(stderr, "\rRendering (%d samples per pixel): [%s] %5.2f%%", samples, bar.c_str(), percent);
				fflush(stderr);
			}
			for (unsigned short x = 0, Xi[3] = { 0, 0, y * y * y }; x < width; x++) { // Go through image cols 
				int i = (height - y - 1) * width + x; // Index in the image array
				result = Vec(); // Reset result for the current pixel

				for (int s = 0; s < samples; s++) {
					double r1 = erand48(Xi), r2 = erand48(Xi);
					double z1, z2;
					box_muller(r1, r2, z1, z2);
					double offset_x = z1 * sigma, offset_y = z2 * sigma;
					offset_x = std::clamp(offset_x, -1.0, 1.0);
					offset_y = std::clamp(offset_y, -1.0, 1.0);

					// For tent filter
					//double r1 = 2 * erand48(Xi), offset_x = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
					//double r2 = 2 * erand48(Xi), offset_y = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

					Vec d = camera_x * ((x + offset_x + 0.5) / width - 0.5) +
						camera_y * ((y + offset_y + 0.5) / height - 0.5) + camera.dir;
					result = result + path_tracing(Ray(camera.orig + d * 140, d.norm()), 0, Xi) * (1.0 / samples);
				}
				image[i] = image[i] + Vec(clamp01(result.x), clamp01(result.y), clamp01(result.z));
			}
		}
	}
};

GLuint create_texture() {
	auto* pixels = new unsigned char[width * height * 3];
	for (int i = 0; i < width * height; ++i) {
		pixels[i * 3] = to_int(image[i].x);
		pixels[i * 3 + 1] = to_int(image[i].y);
		pixels[i * 3 + 2] = to_int(image[i].z);
	}

	GLuint texture_id;
	glGenTextures(1, &texture_id);
	glBindTexture(GL_TEXTURE_2D, texture_id);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	delete[] pixels;

	return texture_id;
}

void user_interaction() {
	std::cout << "\n          +-----------------+" << std::endl;
	std::cout << "         /                 /|" << std::endl;
	std::cout << "        /                 / |" << std::endl;
	std::cout << "       +-----------------+  |" << std::endl;
	std::cout << "       |                 |  |" << std::endl;
	std::cout << "       |                 |  |" << std::endl;
	std::cout << "       |                 |  |" << std::endl;
	std::cout << "       |    ?       ?    |  |" << std::endl;
	std::cout << "       |  Model1  Model2 | /" << std::endl;
	std::cout << "       |                 |/" << std::endl;
	std::cout << "       +-----------------+\n" << std::endl;

	for (int i = 0; i < model_choices.size(); ++i) {
		std::cout << "Choose model #" << i + 1 << ":\n1 - Sphere\n2 - Cube\n3 - Pinetree\n4 - Dog\n5 - Skull\n6 - Heart\n7 - Gear\n8 - Person statue\n9 - Techno-sphere" << std::endl;
		std::cin >> model_choices[i];
		if (model_choices[i] >= 0 && model_choices[i] <= 9)
			--model_choices[i];
		else {
			std::cout << "Entered incorrect model ID, a sphere will be used instead." << std::endl;
			model_choices[i] = 0;
		}

		std::cout << "\nChoose a material for the model #" << i + 1 << ":\n1 - Diffuse\n2 - Mirror\n3 - Glass\n4 - Frosted glass\n5 - Ruby\n6 - Gold\n7 - Iron" << std::endl;
		std::cin >> material_choices[i];
		if (material_choices[i] >= 1 && material_choices[i] <= 7)
			--material_choices[i];
		else {
			std::cout << "Entered incorrect choice, the model will be diffuse.\n" << std::endl;
			material_choices[i] = 0;
		}
		std::cout << std::endl;

		if (material_choices[i] >= 3 && material_choices[i] <= 6) {
			std::cout << "Enter roughness of the material (belongs to [0..1] range): ";
			std::cin >> roughness_choices[i];
			roughness_choices[i] = std::clamp(roughness_choices[i], 1e-4, 1.0);
			std::cout << std::endl;
		}
	}
}

int main() {
	setlocale(LC_ALL, "RUSSIAN");
	while (true) {
		std::cout << "Enter image resolution (in the format 'width height', both values must be >= 100):\n";
		std::cin >> width;
		std::cin >> height;
		if (width >= 100 && height >= 100)
			break;
		else
			std::cout << "Entered incorrect resolution value.\n" << std::endl;
	}
	while (true) {
		std::cout << "Enter number of samples per pixel (must be >= 4): ";
		std::cin >> samples;
		if (samples >= 4) {
			samples /= 4;
			break;
		}
		else
			std::cout << "Entered incorrect samples per pixel value.\n" << std::endl;
	}
	user_interaction();

	image = new Vec[width * height];
	Scene sc;
	sc.fill();
	sc.rebuild_bvh(2);
	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	sc.render_with_supersampling();
	//sc.render_wout_supersampling(samples * 4);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	printf("\nTime elapsed: %d ms;\n", std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit())
		return -1;
	int win_width = width < 1024 ? 1024 : width + 15;
	int win_height = height < 768 ? 768 : height + 84;
	GLFWwindow* window = glfwCreateWindow(win_width, win_height, "Path Tracing", NULL, NULL);
	if (window == NULL)
		return -1;
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		printf("Failed to initialize GLEW\n");
		return 0;
	}

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 130");

	//ImFileDialog requires you to set the CreateTexture and DeleteTexture
	
	ifd::FileDialog::Instance().CreateTexture = [](uint8_t* data, int w, int h, char fmt) -> void* {
		GLuint tex;

		glGenTextures(1, &tex);
		glBindTexture(GL_TEXTURE_2D, tex);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, (fmt == 0) ? GL_BGRA : GL_RGBA, GL_UNSIGNED_BYTE, data);
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);

		return reinterpret_cast<void*>(static_cast<uintptr_t>(tex));
	};
	ifd::FileDialog::Instance().DeleteTexture = [](void* tex) {
		GLuint texID = static_cast<GLuint>(reinterpret_cast<uintptr_t>(tex));
		glDeleteTextures(1, &texID);
	};

	GLuint tex_id = create_texture();
	// Save the result to png image
	auto* res_image = new unsigned char[width * height * 3];
	for (int i = 0; i < width * height; i++) {
		res_image[i * 3] = to_int(image[i].x);
		res_image[i * 3 + 1] = to_int(image[i].y);
		res_image[i * 3 + 2] = to_int(image[i].z);
	}

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::SetNextWindowPos({ 0, 0 });
		ImGui::SetNextWindowSize({ io.DisplaySize.x, io.DisplaySize.y });
		ImGui::Begin("Rendered image", (bool*)0, ImGuiWindowFlags_NoResize |
			ImGuiWindowFlags_NoCollapse);
		ImGui::SetCursorPos(ImVec2(8, 27));
		if (ImGui::Button("Save image", ImVec2(150, 40))) {
			ifd::FileDialog::Instance().Save("ImageSaveDialog", "Save an image", "*.png;*.jpg {.png,.jpg}");
		}
		ImGui::SetCursorPos(ImVec2(165, 53));
		ImGui::Text("%s", ui_text.c_str());
		ImGui::SetCursorPos(ImVec2(8, 75));
		ImGui::Image((void*)(intptr_t)tex_id, ImVec2(width, height));
		ImGui::End();

		// File dialog
		if (ifd::FileDialog::Instance().IsDone("ImageSaveDialog")) {
			ui_text = "";
			if (ifd::FileDialog::Instance().HasResult()) {
				std::string res = ifd::FileDialog::Instance().GetResult().string();
				std::string ext = res.substr(res.length() - 3);
				ui_text = "Image has been saved.";
				if (ext == "png")
					stbi_write_png(res.c_str(), width, height, 3, res_image, width * 3);
				else if (ext == "jpg")
					stbi_write_jpg(res.c_str(), width, height, 3, res_image, 100);
				else
					ui_text = "";
			}
			ifd::FileDialog::Instance().Close();
		}

		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		glfwSwapBuffers(window);
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	glfwDestroyWindow(window);
	glfwTerminate();
	std::remove("imgui.ini");

	return 0;
}
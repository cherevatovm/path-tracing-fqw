#ifndef SHAPE_H
#define SHAPE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include "vec.h"
#include "ray.h"
#include "material.h"
#include "bounds.h"
#include "transform.h"
#include "utility.h"
#include "intersection.h"

using uint = unsigned int;

class Shape {
public:
	Material material;
	Vec emis;

	Shape() = default;

	Shape(Material material_, Vec emis_ = Vec()) :
		material(material_), emis(emis_) {};

	virtual Intersection intersect(const Ray& r) const = 0;
	
	virtual Bounds3 get_bounds() const = 0;

	virtual Vec get_normal(const Intersection& inters) const = 0;
};

class Sphere : public Shape {
public:
	Vec center;
	double radius;

	Sphere() = default;

	Sphere(double radius_, Vec center_, 
		Material material_, Vec emis_ = Vec()) : 
		Shape(material_, emis_),
		radius(radius_), center(center_) {};

	Sphere& operator=(const Sphere& sph) {
		center = sph.center;
		radius = sph.radius;
		material = sph.material;
		emis = sph.emis;
		return *this;
	}

	// Returns distance or 0 if there is no hit 
	Intersection intersect(const Ray& r) const override {
		Vec op = r.orig - center; // Solve t^2 * d.d + 2 * t * (o - p).d + (o - p).(o - p) - R^2 = 0
		double b = 2 * (op).dot_prod(r.dir);
		double c = op.dot_prod(op) - radius * radius;
		double discr = b * b - 4 * c;

		if (discr < 0) return Intersection();
		else discr = sqrt(discr);

		double t0 = -b + discr;
		double t1 = -b - discr;

		return (t1 > eps) ? Intersection(0.5 * t1, this) :
			((t0 > eps) ? Intersection(0.5 * t0, this) : Intersection());
	}

	Vec get_normal(const Intersection& inters) const override {
		return (inters.hit_point - center).norm();
	}

	Bounds3 get_bounds() const override {
		return Bounds3(Vec(center.x - radius, center.y - radius, center.z - radius), 
			Vec(center.x + radius, center.y + radius, center.z + radius));
	}
};

class Mesh {
private:
	class MeshTriangle : public Shape {
	public:
		const Mesh* mesh;
		uint vert0_ind, vert1_ind, vert2_ind;
		uint vert0_norm_ind, vert1_norm_ind, vert2_norm_ind;
		Vec surface_normal;
		bool is_flat;

		MeshTriangle() = default;

		MeshTriangle(const Mesh* mesh_, Material material_) : 
			Shape(material_, Vec()), mesh(mesh_) {}

		MeshTriangle(const Mesh* mesh_,
			uint v1_ind_, uint v2_ind_, uint v3_ind_,
			Material material_, Vec emis_ = Vec(),
			uint v1_n_ind_ = UINT32_MAX, uint v2_n_ind_ = UINT32_MAX, uint v3_n_ind_ = UINT32_MAX) :
			Shape(material_, emis_), mesh(mesh_),
			vert0_ind(v1_ind_), vert1_ind(v2_ind_), vert2_ind(v3_ind_),
			vert0_norm_ind(v1_n_ind_), vert1_norm_ind(v2_n_ind_), vert2_norm_ind(v3_n_ind_) 
		{
			surface_normal = ((mesh->vertices[vert1_ind] - mesh->vertices[vert0_ind]) %
				(mesh->vertices[vert2_ind] - mesh->vertices[vert0_ind])).norm();
		}

		MeshTriangle& operator=(const MeshTriangle& t) {
			mesh = t.mesh;
			vert0_ind = t.vert0_ind;
			vert1_ind = t.vert1_ind;
			vert2_ind = t.vert2_ind;
			vert0_norm_ind = t.vert0_norm_ind;
			vert1_norm_ind = t.vert1_norm_ind;
			vert2_norm_ind = t.vert2_norm_ind;
			material = t.material;
			emis = t.emis;
			return *this;
		}

		void calc_surface_normal() {
			surface_normal = ((mesh->vertices[vert1_ind] - mesh->vertices[vert0_ind]) %
				(mesh->vertices[vert2_ind] - mesh->vertices[vert0_ind])).norm();
		}

		void set_vert_ind(uint i, uint new_val) { 
			i == 0 ? vert0_ind = new_val : (i == 1 ? vert1_ind = new_val : vert2_ind = new_val);
		}
		
		void set_norm_ind(uint i, uint new_val) { 
			i == 0 ? vert0_norm_ind = new_val : (i == 1 ? vert1_norm_ind = new_val : vert2_norm_ind = new_val); 
		}

		void check_if_flat() {
			Vec normal0 = mesh->vert_normals[vert0_norm_ind];
			normal0.norm();
			Vec normal1 = mesh->vert_normals[vert1_norm_ind];
			normal1.norm();
			Vec normal2 = mesh->vert_normals[vert2_norm_ind];
			normal2.norm();
			is_flat = (normal0 - normal1).length() < eps && (normal0 - normal2).length() < eps;
		}

		Intersection intersect(const Ray& r) const override {
			Vec e0 = mesh->vertices[vert1_ind] - mesh->vertices[vert0_ind];
			Vec e1 = mesh->vertices[vert2_ind] - mesh->vertices[vert0_ind];
			Vec ray_cr_e1 = r.dir % e1;
			double det = e0.dot_prod(ray_cr_e1);

			if (abs(det) < eps)
				return {};

			double inv_det = 1.0 / det;
			Vec s = r.orig - mesh->vertices[vert0_ind];
			double v = inv_det * s.dot_prod(ray_cr_e1);

			if (v < -eps || v > 1 + eps)
				return {};

			Vec s_cross_e1 = s % e0;
			double w = inv_det * r.dir.dot_prod(s_cross_e1);

			if (w < -eps || v + w > 1 + eps)
				return {};

			double t = inv_det * e1.dot_prod(s_cross_e1);
			if (t > eps)
				return Intersection(t, this, Vec(1 - v - w, v, w));
			
			return {};
		}

		Vec get_normal(const Intersection& inters) const override {
			if (vert0_norm_ind == UINT32_MAX || is_flat)
				return surface_normal;
			else {
				Vec hp_normal = mesh->vert_normals[vert0_norm_ind] * inters.baryc_coords.x +
					mesh->vert_normals[vert1_norm_ind] * inters.baryc_coords.y +
					mesh->vert_normals[vert2_norm_ind] * inters.baryc_coords.z;
				if (hp_normal.length() > eps)
					hp_normal.norm();
				else
					hp_normal = surface_normal;
				return hp_normal;
			}
		}

		Bounds3 get_bounds() const override {
			const Vec& v0 = mesh->vertices[vert0_ind];
			const Vec& v1 = mesh->vertices[vert1_ind];
			const Vec& v2 = mesh->vertices[vert2_ind];

			Vec p_min = Vec::min(v0, Vec::min(v1, v2));
			Vec p_max = Vec::max(v0, Vec::max(v1, v2));

			p_min -= Vec(eps, eps, eps);
			p_max += Vec(eps, eps, eps);

			return Bounds3(p_min, p_max);
		}

		Vec compute_barycentric_coords(const Vec& hit_point) const {
			Vec v0 = mesh->vertices[vert1_ind] - mesh->vertices[vert0_ind];
			Vec v1 = mesh->vertices[vert2_ind] - mesh->vertices[vert0_ind];
			Vec v2 = hit_point - mesh->vertices[vert0_ind];

			double d00 = v0.dot_prod(v0);
			double d01 = v0.dot_prod(v1);
			double d11 = v1.dot_prod(v1);
			double d20 = v2.dot_prod(v0);
			double d21 = v2.dot_prod(v1);

			double denom = d00 * d11 - d01 * d01;
			if (abs(denom) < eps)
				return Vec(-1, -1, -1);
			double v = (d11 * d20 - d01 * d21) / denom;
			double w = (d00 * d21 - d01 * d20) / denom;
			double u = 1.0 - v - w;
			Vec res1(u, v, w);

			return Vec(u, v, w);
		}
	};

	Vec centroid() const {
		double x_sum = 0, y_sum = 0, z_sum = 0;
		for (const auto& vertex : vertices) {
			x_sum += vertex.x;
			y_sum += vertex.y;
			z_sum += vertex.z;
		}
		size_t count = vertices.size();
		return { x_sum / count, y_sum / count, z_sum / count };
	}

	void affine_transformation(const Transform& tr) {
		for (Vec& vert : vertices)
			vert = tr.apply_transform(vert);
		for (Vec& norm : vert_normals)
			norm = tr.transform_normal(norm);
		for (MeshTriangle* face : faces)
			face->calc_surface_normal();
	}

	void calc_vertex_normals() {
		vert_normals.clear();
		for (uint i = 0; i < vertices.size(); ++i) {
			std::list<Vec> triangle_normals;
			for (MeshTriangle* face : faces) {
				if (face->vert0_ind == i || face->vert1_ind == i || face->vert2_ind == i)
					triangle_normals.push_back(face->surface_normal);
			}

			Vec vert_norm;
			for (const Vec& normal : triangle_normals)
				vert_norm += normal;
			vert_norm *= (1.0 / triangle_normals.size());
			vert_normals.push_back(vert_norm);
		}
		for (MeshTriangle* face : faces) {
			face->vert0_norm_ind = face->vert0_ind;
			face->vert1_norm_ind = face->vert1_ind;
			face->vert2_norm_ind = face->vert2_ind;
		}
	}

	void process_face(MeshTriangle* face, const std::vector<std::string>& face_data) {
		uint vertex_ind, norm_ind, tex_coord_ind;
		for (int i = 0; i < 3; ++i) {
			std::istringstream iss(face_data[i]);

			iss >> vertex_ind;
			face->set_vert_ind(i, --vertex_ind);
			char ch1 = iss.peek();
			if (ch1 == '/') {
				iss.ignore();
				ch1 = iss.peek();
				if (ch1 == '/') {
					iss.ignore();
					iss >> norm_ind;
					face->set_norm_ind(i, --norm_ind);
				}
				else if (isdigit(ch1)) {
					iss >> tex_coord_ind;
					//res_vert.tex_coords = vert_tex_coords[--tex_coord_ind];
					ch1 = iss.peek();
					if (ch1 == '/') {
						iss.ignore();
						iss >> norm_ind;
						face->set_norm_ind(i, --norm_ind);
					}
				}
			}
		}
		if (vertex_ind > 1000000 || norm_ind > 1000000)
			std:: cout << "pizdec" << std::endl;
		face->calc_surface_normal();
	}

	void clear() {
		for (int i = 0; i < faces.size(); ++i)
			delete faces[i];
		faces.clear();
		vertices.clear();
		vert_normals.clear();
	}

public:
	std::vector<MeshTriangle*> faces;
	std::vector<Vec> vertices;
	std::vector<Vec> vert_normals;
	Material material;

	Mesh(Material material_) : material(material_) {};

	Mesh(const std::vector<Vec>& vertices_, const std::vector<Vec>& vert_normals_,
		const std::vector<uint>& indices_, Material material_) :
		vertices(vertices_), vert_normals(vert_normals_), material(material_) {
		bool has_normals = vert_normals_.size() > 0;
		uint inc = has_normals ? 6 : 3;
		for (int i = 0; i < indices_.size(); i += inc) {
			if (has_normals) {
				faces.push_back(new MeshTriangle(this, indices_[i], indices_[i + 2], indices_[i + 4],
					material_, Vec(), indices_[i + 1], indices_[i + 3], indices_[i + 5]));
			}
			else
				faces.push_back(new MeshTriangle(this, indices_[i], indices_[i + 1], indices_[i + 2], material_, Vec()));
		}
		if (!has_normals)
			calc_vertex_normals();
	};

	void rotate(double theta, int rotate_index, bool by_center = true) { 
		if (by_center)
			affine_transformation(Transform::rotate_around_point(centroid(), theta, rotate_index)); 
		else
			affine_transformation(Transform::rotate_around_axis(theta, rotate_index));
	}
	
	void scale(double x, double y, double z, bool by_center = true) {
		if (by_center)
			affine_transformation(Transform::scale_around_point(centroid(), x, y, z));
		else
			affine_transformation(Transform::scale(x, y, z));
	}

	void translate(const Vec& delta) {
		affine_transformation(Transform::translate(delta));
	}

	int load_model(const std::string& file_name) {
		std::ifstream file(file_name);
		if (!file.is_open()) {
			std::cerr << "Failed to open file: " << file_name << std::endl;
			return -1;
		}
		clear();
		std::string line;
		while (std::getline(file, line)) {
			std::istringstream iss(line);
			std::string type;
			iss >> type;
			if (type.empty() || type[0] == '#') {
				continue;
			}
			if (type == "v") {
				float x, y, z;
				iss >> x >> y >> z;
				vertices.emplace_back(x, y, z);
			}
			else if (type == "vn") {
				float x, y, z;
				iss >> x >> y >> z;
				vert_normals.emplace_back(x, y, z);
			}
			else if (type == "f") {
				auto spl = split(line);
				for (int i = 3; i < spl.size(); ++i) {
					MeshTriangle* cur_triangle = new MeshTriangle(this, material);
					faces.push_back(cur_triangle);
					process_face(cur_triangle, std::vector<std::string>{ spl[1], spl[i - 1], spl[i] });
				}
			}
		}
		file.close();
		if (vert_normals.size() > 0)
			for (MeshTriangle* face : faces) face->check_if_flat();
		return 0;
	}

	static Mesh* create_hexahedron(Vec center_, double radius_, 
		double alpha_, Material material_) {
		Vec p1{ center_.x + radius_, center_.y + radius_, center_.z + radius_ };
		Vec p7{ center_.x - radius_, center_.y - radius_, center_.z - radius_ };
		Vec p2{ p1.x, p1.y, p7.z };
		Vec p3{ p7.x, p1.y, p7.z };
		Vec p4{ p7.x, p1.y, p1.z };
		Vec p5{ p1.x, p7.y, p1.z };
		Vec p6{ p1.x, p7.y, p7.z };
		Vec p8{ p7.x, p7.y, p1.z };

		p1.rotate_y(alpha_, center_);
		p2.rotate_y(alpha_, center_);
		p3.rotate_y(alpha_, center_);
		p4.rotate_y(alpha_, center_);
		p5.rotate_y(alpha_, center_);
		p6.rotate_y(alpha_, center_);
		p7.rotate_y(alpha_, center_);
		p8.rotate_y(alpha_, center_);

		return new Mesh(std::vector<Vec>{ p1, p2, p3, p4, p5, p6, p7, p8 },
			std::vector<Vec>(),
			std::vector<uint>{ 0, 1, 2, 2, 3, 0,
			7, 6, 5, 5, 4, 7, 1, 5, 6, 6, 2,
			1, 0, 4, 5, 5, 1, 0, 7, 6, 2,
			2, 3, 7, 0, 3, 7, 7, 4, 0 }, material_);
	}
};

#endif
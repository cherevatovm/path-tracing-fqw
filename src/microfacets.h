#ifndef MICROFACETS_H
#define MICROFACETS_H
#define _USE_MATH_DEFINES
#include "vec.h"

Vec sample_wm(unsigned short* Xi, const Vec& n, double roughness) {
    double r1 = 2 * M_PI * erand48(Xi);
    double r2 = erand48(Xi);
    double alpha = roughness * roughness;
    double theta = atan(alpha * sqrt(r2) / sqrt(1 - r2));
    double phi = r1;

    Vec local_dir = Vec(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    );

    // Transform to world coords
    Vec w = n, u, v;
    create_orthonorm_sys(w, u, v);
    return (u * local_dir.x + v * local_dir.y + w * local_dir.z).norm();
}

inline double fresnel_schlick(double cos_theta, double F0) { return F0 + (1 - F0) * std::pow(1 - cos_theta, 5); }

inline double lambda(const Vec& w, double roughness) {
    double tan2theta = tan_2_theta(w);
    if (std::isinf(tan2theta)) return 0;
    return (std::sqrt(1 + roughness * roughness * tan2theta) - 1) / 2;
}

inline double G1(const Vec& w, double roughness) { return 1 / (1 + lambda(w, roughness)); }
inline double G(const Vec& wo, const Vec& wi, double roughness) { return 1 / (1 + lambda(wo, roughness) + lambda(wi, roughness)); }

inline double D(const Vec& wm, double roughness) {
    double tan2theta = tan_2_theta(wm);
    if (std::isinf(tan2theta)) return 0;
    double cos4theta = cos_2_theta(wm) * cos_2_theta(wm);
    double alpha2 = roughness * roughness;
    double e = alpha2 + tan2theta;
    return alpha2 / (M_PI * cos4theta * e * e);
}

inline double D(const Vec& w, const Vec& wm, double roughness) { 
    return G1(w, roughness) / abs_cos_theta(w) * D(wm, roughness) * std::abs(w.dot_prod(wm)); }

double beckmann_distribution(const Vec& n, const Vec& h, double roughness) {
    double n_dot_h = std::max(n.dot_prod(h), 0.0);
    double alpha = roughness * roughness;
    double alpha2 = alpha * alpha;
    double n_dot_h_sqr = n_dot_h * n_dot_h;
    double tan2 = (1 - n_dot_h_sqr) / n_dot_h_sqr;
    return std::exp(-tan2 / alpha2) / (M_PI * alpha2 * n_dot_h_sqr * n_dot_h_sqr);
}

double ggx_distribution(const Vec& n, const Vec& h, double roughness) {
    double alpha = roughness * roughness;
    double alpha2 = alpha * alpha;
    double n_dot_h = std::max(n.dot_prod(h), 0.0);
    double denom = n_dot_h * n_dot_h * (alpha2 - 1) + 1;
    return alpha2 / (M_PI * denom * denom);
}

double geom_cook_torrance(const Vec& n, const Vec& v, const Vec& l, const Vec& h, double roughness) {
    double n_dot_v = std::max(n.dot_prod(v), 0.0);
    double n_dot_l = std::max(n.dot_prod(l), 0.0);
    double n_dot_h = std::max(n.dot_prod(h), 0.0);
    double v_dot_h = std::max(v.dot_prod(h), 0.0);
    double G1 = (2 * n_dot_h * n_dot_v) / v_dot_h;
    double G2 = (2 * n_dot_h * n_dot_l) / v_dot_h;
    return std::min(1.0, std::min(G1, G2));
}

inline double geom_schlick_ggx(double n_dot_v, double roughness) {
    double k = (roughness + 1) * (roughness + 1) / 8;
    return n_dot_v / (n_dot_v * (1 - k) + k + 1e-6); // +1e-6 to prevent divison by 0
}

double geom_smith(const Vec& n, const Vec& v, const Vec& l, double roughness) {
    double n_dot_v = std::max(n.dot_prod(v), 0.0);
    double n_dot_l = std::max(n.dot_prod(l), 0.0);
    double ggx1 = geom_schlick_ggx(n_dot_v, roughness);
    double ggx2 = geom_schlick_ggx(n_dot_l, roughness);
    return ggx1 * ggx2;
}

double brdf_div_by_pdf(const Vec& n, const Vec& wo, const Vec& wi,
    const Vec& wm, double roughness, double Fr_) {
    double D_ = D(wm, roughness);
    double G_ = geom_smith(n, wo, wi, roughness);

    double brdf = (D_ * G_ * Fr_) / (4.0 * abs(n.dot_prod(wo)) * abs(n.dot_prod(wi)));
    brdf *= std::abs(n.dot_prod(wi));
    double pdf = (D_ * std::max(n.dot_prod(wm), 0.0)) / (4 * std::max(wm.dot_prod(wo), 1e-6));

    return brdf / pdf;
}

double btdf_div_by_pdf(const Vec& wo, const Vec& wi,
    const Vec& wm, double roughness, double Tr_, double refr_ratio) {
    double D_ = D(wm, roughness);
    double G_ = G(wo, wi, roughness);

    double denom = std::pow(wi.dot_prod(wm) + wo.dot_prod(wm) / refr_ratio, 2);
    double dwm_dwi = std::abs(wi.dot_prod(wm)) / denom;
    double pdf = D(wo, wm, roughness) * dwm_dwi;
    double btdf = (D_ * G_ * Tr_) * std::abs(wi.dot_prod(wm) * wo.dot_prod(wm) / denom);

    return btdf / pdf;
}

#endif

#pragma once

#include "torrey.h"
#include "vector.h"
#include "scene.h"

struct BRDFEvalRecord {
    Vector3 value;
    Real pdf;
};

inline BRDFEvalRecord eval_brdf(const Scene &scene,
                                const Intersection &isect,
                                const Vector3 &n,
                                const Vector3 &wi,
                                const Vector3 &wo,
                                bool debug = false) {
    const Material &mat = scene.materials[isect.material_id];

    if (auto *diffuse = std::get_if<Diffuse>(&mat)) {
        Vector3 value = eval(diffuse->reflectance, isect.uv) *
                        (max(dot(wo, n), Real(0)) / c_PI);
        Real pdf = (max(dot(wo, n), Real(0)) / c_PI);
        return {value, pdf};
    } else if (auto *mirror = std::get_if<Mirror>(&mat)) {
        // pure specular reflection is handled somewhere else
        return {Vector3{0, 0, 0}, 0};
    } else if (auto *plastic = std::get_if<Plastic>(&mat)) {
        Real F0s = square((plastic->ior - 1) / (plastic->ior + 1));
        Vector3 F0{F0s, F0s, F0s};
        Vector3 F = schlick_fresnel(F0, dot(n, wi));
        // pure specular reflection is handled somewhere else
        Vector3 value = 
            (1 - F) * eval(plastic->reflectance, isect.uv) *
            (max(dot(wo, n), Real(0)) / c_PI);
        Real pdf = (1 - F[0]) * (max(dot(wo, n), Real(0)) / c_PI);
        return {value, pdf};
    } else if (auto *phong = std::get_if<Phong>(&mat)) {
        Vector3 r = -wi + 2 * dot(wi, n) * n;
        Real r_dot_wo = dot(r, wo);
        Real n_dot_wo = dot(n, wo);
        if (r_dot_wo > 0 && n_dot_wo > 0) {
            Real phong_response = 
                ((phong->exponent + 1) / (2 * c_PI)) *
                pow(r_dot_wo, phong->exponent);
            Vector3 value = 
                eval(phong->reflectance, isect.uv) * 
                phong_response;
            Real pdf = ((phong->exponent + 1) / (2 * c_PI)) *
                pow(r_dot_wo, phong->exponent);
            return {value, pdf};
        } else {
            // actual PDF may not be zero but it doesn't matter
            return {Vector3{0, 0, 0}, 0};
        }
    } else if (auto *blinn = std::get_if<BlinnPhong>(&mat)) {
        Vector3 h = normalize(wi + wo);
        Real n_dot_wo = dot(n, wo);
        Real h_dot_wo = dot(h, wo);
        Real n_dot_h = dot(n, h);
        if (n_dot_wo > 0 && h_dot_wo > 0 && n_dot_h > 0) {
            Real e = blinn->exponent;
            Real normalization = (e + 2) / (4 * c_PI * (2 - pow(2, -e/2)));
            Real blinnphong_response = 
                normalization * pow(n_dot_h, e);
            Vector3 F0 = eval(blinn->reflectance, isect.uv);
            Vector3 F = schlick_fresnel(F0, dot(h, wo));
            Vector3 value = F * blinnphong_response;
            Real pdf = (e + 1) * pow(n_dot_h, e) / 
                       (2 * c_PI * 4 * h_dot_wo);
            return {value, pdf};
        } else {
            // actual PDF may not be zero but it doesn't matter
            return {Vector3{0, 0, 0}, 0};
        }
    } else if (auto *blinn_mic = std::get_if<BlinnPhongMicrofacet>(&mat)) {
        Vector3 h = normalize(wi + wo);
        Real n_dot_wo = dot(n, wo);
        if (n_dot_wo > 0) {
            Real e = blinn_mic->exponent;
            Vector3 F0 = eval(blinn_mic->reflectance, isect.uv);
            Vector3 F = schlick_fresnel(F0, dot(h, wo));
            Real D = ((e + 2) / (2 * c_PI)) * pow(dot(h, n), e);
            Real G = smithG1(wi, h, n, e) *
                     smithG1(wo, h, n, e);
            Vector3 value = F * D * G / (4 * dot(n, wi));
            Real pdf = (e + 1) * pow(dot(n, h), e) / 
                       (2 * c_PI * 4 * dot(wo, h));
            return {value, pdf};
        } else {
            // actual PDF may not be zero but it doesn't matter
            return {Vector3{0, 0, 0}, 0};
        }
    }
    assert(false);
    return {Vector3{0, 0, 0}, 0};
}

struct SamplingRecord {
    Vector3 wo;
    bool is_pure_specular = false;
    // Undefined if is_pure_specular = false
    Vector3 weight = Vector3{0, 0, 0}; 
};

inline SamplingRecord sample_brdf(
        const Scene &scene,
        const Intersection &isect,
        const Vector3 &n,
        const Vector3 &wi,
        pcg32_state &rng) {
    const Material &mat = scene.materials[isect.material_id];
    if (auto *diffuse = std::get_if<Diffuse>(&mat)) {
        Frame frame(n);
        Vector2 u{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Vector3 wo = to_world(frame, sample_cos_hemisphere(u));
        return {wo};
    } else if (auto *mirror = std::get_if<Mirror>(&mat)) {
        Vector3 wo = -wi + 2 * dot(wi, n) * n;
        Vector3 F0 = eval(mirror->reflectance, isect.uv);
        Vector3 F = schlick_fresnel(F0, dot(n, wo));
        return {wo, true /* is_pure_specular */, F};
    } else if (auto *plastic = std::get_if<Plastic>(&mat)) {
        Real F0s = square((plastic->ior - 1) / (plastic->ior + 1));
        Vector3 F0{F0s, F0s, F0s};
        Vector3 F = schlick_fresnel(F0, dot(n, wi));
        if (next_pcg32_real<Real>(rng) <= F[0]) {
            // Sample mirror reflection with probability F
            Vector3 wo = -wi + 2 * dot(wi, n) * n;
            // Weight is 1 since F cancels out
            return {wo, true, Vector3{1, 1, 1}};
        } else {
            // Sample diffuse with probability (1-F)
            Frame frame(n);
            Vector2 u{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Vector3 wo = to_world(frame, sample_cos_hemisphere(u));
            return {wo};
        }
    } else if (auto *phong = std::get_if<Phong>(&mat)) {
        Vector2 u{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Vector3 wo_local = sample_cos_n_hemisphere(u, phong->exponent);
        Frame frame(-wi + 2 * dot(wi, n) * n);
        Vector3 wo_world = to_world(frame, wo_local);
        return {wo_world};
    } else if (auto *blinn = std::get_if<BlinnPhong>(&mat)) {
        Vector2 u{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Vector3 h_local = sample_cos_n_hemisphere(u, blinn->exponent);
        Frame frame(n);
        Vector3 h_world = to_world(frame, h_local);
        Vector3 wo = -wi + 2 * dot(wi, h_world) * h_world;
        return {wo};
    } else if (auto *blinn_mic = std::get_if<BlinnPhongMicrofacet>(&mat)) {
        Vector2 u{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Vector3 h_local = sample_cos_n_hemisphere(u, blinn_mic->exponent);
        Frame frame(n);
        Vector3 h_world = to_world(frame, h_local);
        Vector3 wo = -wi + 2 * dot(wi, h_world) * h_world;
        return {wo};
    }
    assert(false);
    return {};
}

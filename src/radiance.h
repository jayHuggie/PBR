#pragma once

#include "torrey.h"
#include "vector.h"
#include "scene.h"
#include "sample.h"

inline Vector3 radiance_hw4 (const Scene &scene, ray input_ray, pcg32_state &rng, int depth, bool debug = false) {
    if (depth <= 0) {
        return Vector3{0, 0, 0};
    }
    if (auto hit_isect = intersect(scene, input_ray)) {
        Vector3 p = hit_isect->position;
        Vector3 n = hit_isect->shading_normal;
        Vector3 L = Vector3{0, 0, 0};
        if (hit_isect->area_light_id != -1) {
            const DiffuseAreaLight &light =
                std::get<DiffuseAreaLight>(scene.lights[hit_isect->area_light_id]);
            if (dot(-input_ray.dir, n) > 0) {
                L += light.radiance;
            }
        }
        if (dot(-input_ray.dir, n) < 0) {
            n = -n;
        }

        SamplingRecord sampling = sample_brdf(
            scene, *hit_isect, n, -input_ray.dir, rng);
        ray refl_ray{p, sampling.wo, Real(1e-4), infinity<Real>()};
        if (sampling.is_pure_specular) {
            if (max(sampling.weight) > 0) {
                L += sampling.weight * 
                    radiance_hw4(scene, refl_ray, rng, depth - 1);
            }
        } else {
            BRDFEvalRecord brdf = eval_brdf(
                scene, *hit_isect, n, -input_ray.dir, sampling.wo, debug);
            if (max(brdf.value) > 0 && brdf.pdf > 0) {
                L += (brdf.value / brdf.pdf) *
                    radiance_hw4(scene, refl_ray, rng, depth - 1);
            }
        }
        return L;
    }
    return scene.background_color;
}


inline Vector3 radiance_mis(const Scene &scene, ray input_ray, pcg32_state &rng, int depth) {
    if (depth <= 0) {
        return Vector3{0, 0, 0};
    }
    if (auto hit_isect = intersect(scene, input_ray)) {
        Vector3 p = hit_isect->position;
        Vector3 n = hit_isect->shading_normal;
        Vector3 L = Vector3{0, 0, 0};
        if (hit_isect->area_light_id != -1) {
            const DiffuseAreaLight &light =
                std::get<DiffuseAreaLight>(scene.lights[hit_isect->area_light_id]);
            if (dot(-input_ray.dir, n) > 0) {
                L += light.radiance;
            }
        }
        if (dot(-input_ray.dir, n) < 0) {
            n = -n;
        }

        SamplingRecord sampling;
        if (scene.lights.size() > 0 && next_pcg32_real<Real>(rng) < Real(0.5)) {
            // Sample light sources
            // Uniformly pick a light
            int light_id = std::clamp(
                int(next_pcg32_real<Real>(rng) * scene.lights.size()),
                0, (int)scene.lights.size() - 1);
            const Light &light = scene.lights[light_id];
            if (auto *area_light = std::get_if<DiffuseAreaLight>(&light)) {
                const Shape &shape = scene.shapes[area_light->shape_id];
                Vector2 u{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                sampling.wo = sample_dir_to_shape(shape, p, u);
            }
        } else {
            // Sample BRDF
            sampling = sample_brdf(scene, *hit_isect, n, -input_ray.dir, rng);
        }

        ray refl_ray{p, sampling.wo, Real(1e-4), infinity<Real>()};
        if (sampling.is_pure_specular) {
            if (scene.lights.size() > 0) {
                // sample_brdf_PDF
                sampling.weight *= Real(2);
            }
            L += sampling.weight * 
                radiance_mis(scene, refl_ray, rng, depth - 1);
        } else {
            Real pdf = 0;
            BRDFEvalRecord brdf = eval_brdf(
                scene, *hit_isect, n, -input_ray.dir, sampling.wo);
            if (max(brdf.value) > 0) {
                if (scene.lights.size() > 0) {
                    pdf += Real(0.5) * brdf.pdf;
                } else {
                    pdf += brdf.pdf;
                }
                for (const Light &light : scene.lights) {
                    if (auto *area_light = std::get_if<DiffuseAreaLight>(&light)) {
                        const Shape &shape = scene.shapes[area_light->shape_id];
                        pdf += Real(0.5) * pdf_sample_dir_to_shape(
                            shape, p, sampling.wo) / scene.lights.size();
                    }
                }
                if (max(brdf.value) > 0 && pdf > 0) {
                    L += (brdf.value / pdf) *
                        radiance_mis(scene, refl_ray, rng, depth - 1);
                }
            }
        }
        return L;
    }
    return scene.background_color;
}


// HW 3
inline Vector3 radiance_new(const Scene &scene, ray input_ray, pcg32_state &rng) {
    
    // Emission
    // add emission of the intersected shape
    if (auto hit_isect = intersect(scene, input_ray)) {
        Vector3 p = hit_isect->position;
        Vector3 n = hit_isect->shading_normal;
        Vector3 L = Vector3{0, 0, 0};
        if (hit_isect->area_light_id != -1) {
            const DiffuseAreaLight &light =
                std::get<DiffuseAreaLight>(scene.lights[hit_isect->area_light_id]);
            if (dot(-input_ray.dir, n) > 0) {
                L += light.radiance;
            }
        }
        if (dot(-input_ray.dir, n) < 0) {
            n = -n;
        }
        const Material &mat = scene.materials[hit_isect->material_id];

        // Direct lighting
        // Currently we only evaluate the diffuse part.
        Vector3 refl = Vector3{0, 0, 0};
        if (auto *diffuse = std::get_if<Diffuse>(&mat)) {
            refl = eval(diffuse->reflectance, hit_isect->uv);
        } else if (auto *plastic = std::get_if<Plastic>(&mat)) {
            Real F0s = square((plastic->ior - 1) / (plastic->ior + 1));
            Vector3 F0{F0s, F0s, F0s};
            Vector3 F = schlick_fresnel(F0, dot(n, -input_ray.dir));
            refl = (1 - F) * eval(plastic->reflectance, hit_isect->uv);
        } else {
            refl = Vector3{0, 0, 0};
        }

        if (max(refl) > Real(0)) {
            // Loop over the lights
            for (const Light &light : scene.lights) {
                if (auto *point_light = std::get_if<PointLight>(&light)) {
                    Vector3 l = point_light-> position - p;
                    ray shadow_ray {p, normalize(l), Real(1e-4), (1 - Real(1e-4)) * length(l)};
                    if (!occluded(scene, shadow_ray)) {
                        L += (max(dot(n, normalize(l)), Real(0)) / c_PI) *
                             (point_light->intensity / length_squared(l)) *
                             refl;
                    }
                } else if (auto *area_light = std::get_if<DiffuseAreaLight>(&light)) {
                    const Shape &shape = scene.shapes[area_light->shape_id];
                    Vector2 u{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                    PointAndNormal pn = sample_on_shape(shape, u);
                    Vector3 l = pn.position - p;
                    Real cos_l = -dot(pn.normal, normalize(l));
                    ray shadow_ray{p, normalize(l), Real(1e-4), (1 - Real(1e-4)) * length(l)};
                    if (cos_l > 0 && !occluded(scene, shadow_ray)) {
                        Real pdf = pdf_sample_on_shape(shape, pn.position);
                        L += (max(dot(n, normalize(l)), Real(0)) / c_PI) *
                             (area_light->radiance / length_squared(l)) *
                             cos_l * refl / pdf;
                    }
                }
            }
        }

        // Scattering
        if (auto *mirror = std::get_if<Mirror>(&mat)) {
            ray refl_ray{p, input_ray.dir - 2 * dot(input_ray.dir, n) * n, Real(1e-4), infinity<Real>()};
            Vector3 F0 = eval(mirror->reflectance, hit_isect->uv);
            Vector3 F = schlick_fresnel(F0, dot(n, refl_ray.dir));
            L += F * radiance_new(scene, refl_ray, rng);
        } else if (auto *plastic = std::get_if<Plastic>(&mat)) {
            ray refl_ray{p, input_ray.dir - 2 * dot(input_ray.dir, n) * n, Real(1e-4), infinity<Real>()};
            Real F0s = square((plastic-> ior - 1) / (plastic-> ior + 1));
            Vector3 F0{F0s, F0s, F0s};
            Vector3 F = schlick_fresnel(F0, dot(n, -input_ray.dir));
            L += F * radiance_new(scene, refl_ray, rng);
        }

        return L;
    }
    return scene.background_color;
}


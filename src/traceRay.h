#pragma once
#include "vector.h"
#include "ray.h"
#include "scene.h"

double hit_multiple_spheres(const Vector3& center, double radius, const ray& r) {

    // ray's orgin - sphere's orgin (move the sphere back to the origin)
    Vector3 oc =  r.origin() - center;

    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius * radius;
    auto discriminant = b * b - 4.0 * a * c;

    double t0 = (-b + sqrt(discriminant) ) / (2.0 * a);
    double closestT = (-b - sqrt(discriminant) ) / (2.0 * a);
     
    if (discriminant < 0.0)
        return -1.0;

   return closestT;
}

Vector3 trace_ray(const ray& r, const Scene& scene) {
    const Sphere* closestSphere = nullptr;
    double hitDistance = std::numeric_limits<Real>::max();
    double EPSILON = 0.0001;

    // Iterate over all spheres in the scene
    for (const auto& shape : scene.shapes) {
        if (std::holds_alternative<Sphere>(shape)) {
            const auto& sphere = std::get<Sphere>(shape);
            double t = hit_multiple_spheres(sphere.position, sphere.radius, r);
            if (t > 0 && t < hitDistance) {
                hitDistance = t;
                closestSphere = &sphere;
            }
        }
    }

    if (closestSphere != nullptr) {
        Vector3 shading_point = r.origin() + hitDistance * r.direction();
        Vector3 normal = normalize(shading_point - closestSphere->position);
        const auto& material = scene.materials[closestSphere->material_id];

        Vector3 hit_color(0.0, 0.0, 0.0);

        if (std::holds_alternative<Mirror>(material)) {
            const auto& mirror = std::get<Mirror>(material);

            // Calculate reflection direction and create a reflected ray
            Vector3 reflection_dir = r.direction() - 2.0 * dot(r.direction(), normal) * normal;
            ray reflected_ray(shading_point + EPSILON * reflection_dir, reflection_dir);

            // Call trace_ray() recursively
            if (std::holds_alternative<Vector3>(mirror.reflectance)) {
                return std::get<Vector3>(mirror.reflectance) * trace_ray(reflected_ray, scene);
            }
        } else if (std::holds_alternative<Diffuse>(material)) {
            const auto& diffuse = std::get<Diffuse>(material);

            for (const auto& light : scene.lights) {
                if (std::holds_alternative<PointLight>(light)) {
                    const auto& point_light = std::get<PointLight>(light);
                    Vector3 light_dir = normalize(point_light.position - shading_point);
                    double visibility = 1.0;
                    ray shadow_ray(shading_point + EPSILON * light_dir, light_dir);
                    // Calculate the distance between shading point and light source
                    double d = distance(shading_point, point_light.position);

                    for (const auto& s : scene.shapes) {
                        if (std::holds_alternative<Sphere>(s)) {
                            const auto& sphere = std::get<Sphere>(s);
                            if (&sphere == closestSphere) continue; // Skip the sphere we hit
                            double t = hit_multiple_spheres(sphere.position, sphere.radius, shadow_ray);
                            if (t > 0 && t < distance(shadow_ray.origin(), point_light.position)) {
                                // There is an intersection between the shadow ray and other spheres,
                                // hence the point is in shadow
                                visibility = 0.0;
                                break;
                            }
                        }
                    }

                    // Calculate diffuse color
                    if (std::holds_alternative<Vector3>(diffuse.reflectance)) {
                        Vector3 diffuse_color = std::get<Vector3>(diffuse.reflectance) * point_light.intensity *
                                                std::max(0.0, dot(normal, light_dir));
                        diffuse_color /= (c_PI * d * d); // Divide by pi * d^2

                        hit_color += visibility * diffuse_color;
                    }
                }
            }
            return hit_color;
        }
    } else {
        return Vector3(0.5, 0.5, 0.5); // Default background color
    }
}

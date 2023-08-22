#include "hw1.h"
#include "hw1_scenes.h"
#include "ray.h"
#include "vector.h"
#include "pcg.h"
#include "parallel.h"
#include <random>

#include <iostream>

using namespace hw1;

Image3 hw_1_1(const std::vector<std::string> &/*params*/) {
    // Homework 1.1: generate camera rays and output the ray directions
    // The camera is positioned at (0, 0, 0), facing towards (0, 0, -1),
    // with an up vector (0, 1, 0) and a vertical field of view of 90 degree.

    Image3 img(640 /* width */, 480 /* height */);
    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            
            // Map pixel coordinates to image plane coordinates
            Real pixel_x = (x + Real(0.5)) / img.width;
            Real pixel_y = (y + Real(0.5)) / img.height;
            
            // Calculate camera ray direction
            Real view_x = 2.0 * pixel_x - 1.0;
            Real view_y = 2.0 * pixel_y - 1.0;
            view_x *= aspect_ratio;

            // Calculate camera ray direction
            Vector3 origin = Vector3 (0.0,0.0,0.0);
            Vector3 ray_dir = 
                Vector3 (view_x, // Map to [-aspect_ratio, aspect_ratio] range
                -view_y, // Map to [-1, 1] range and flip vertically
                -1.0);               // Fixed Z direction for camera facing (0, 0, -1)
            

            Vector3 normal = normalize(ray_dir); // Normalize the ray direction to unit length

            img(x, y) = normal;
        }
    }
    return img;
}


double hit_sphere(const Vector3& center, double radius, const ray& r) {
    // ray's orgin - sphere's orgin (move the sphere back to the origin)
    Vector3 oc =  r.origin() - center;

    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b * b - 4.0 * a * c;

    double t0 = (-b + sqrt(discriminant) ) / (2.0 * a);
    double closestT = (-b - sqrt(discriminant) ) / (2.0 * a);
     
    if (discriminant < 0.0)
        return -1.0;

    return closestT;
}



Vector3 ray_color(const ray& r) {
    Vector3 sphereCenter = Vector3(0.0,0.0,-2.0);
    
    auto t = hit_sphere(sphereCenter, 1.0, r);

    if (t > 0.0) {
        Vector3 intersection_point = r.at(t);
        
        // normalize (hit point - center of sphere) to get normal
        Vector3 N = normalize(intersection_point - sphereCenter);
        return 0.5* Vector3(N.x + 1, N.y + 1, N.z + 1);
    }

    return Vector3(0.5, 0.5, 0.5);
}

Image3 hw_1_2(const std::vector<std::string> &/*params*/) {
    // Homework 1.2: intersect the rays generated from hw_1_1
    // with a unit sphere located at (0, 0, -2)

    Image3 img(640 /* width */, 480 /* height */);
    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            
            // Map pixel coordinates to image plane coordinates
            Real pixel_x = (x + Real(0.5)) / img.width;
            Real pixel_y = (y + Real(0.5)) / img.height;
            
            // Calculate camera ray direction
            Real view_x = 2.0 * pixel_x - 1.0;
            Real view_y = 2.0 * pixel_y - 1.0;
            view_x *= aspect_ratio;

            // Calculate camera ray direction
            Vector3 origin = Vector3 (0.0,0.0,0.0);
            Vector3 ray_dir = 
                Vector3 (view_x, // Map to [-aspect_ratio, aspect_ratio] range
                -view_y, // Map to [-1, 1] range and flip vertically
                -1.0);               // Fixed Z direction for camera facing (0, 0, -1)
            

            ray r (origin, ray_dir);
            Vector3 col = ray_color(r);
            img(x,y) = col;
            
        }
    }



    return img;
}

Image3 hw_1_3(const std::vector<std::string> &params) {
    // Homework 1.3: add camera control to hw_1_2. 
    // We will use a look at transform:
    // The inputs are "lookfrom" (camera position),
    //                "lookat" (target),
    //                and the up vector
    // and the vertical field of view (in degrees).
    // If the user did not specify, fall back to the default
    // values below.
    // If you use the default values, it should render
    // the same image as hw_1_2.

    Vector3 lookfrom = Vector3{0, 0,  0};
    Vector3 lookat   = Vector3{0, 0, -2};
    Vector3 up       = Vector3{0, 1,  0};
    Real    vfov     = 90;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-lookfrom") {
            Real x = std::stof(params[++i]);
            Real y = std::stof(params[++i]);
            Real z = std::stof(params[++i]);
            lookfrom = Vector3{x, y, z};
        } else if (params[i] == "-lookat") {
            Real x = std::stof(params[++i]);
            Real y = std::stof(params[++i]);
            Real z = std::stof(params[++i]);
            lookat = Vector3{x, y, z};
        } else if (params[i] == "-up") {
            Real x = std::stof(params[++i]);
            Real y = std::stof(params[++i]);
            Real z = std::stof(params[++i]);
            up = Vector3{x, y, z};
        } else if (params[i] == "-vfov") {
            vfov = std::stof(params[++i]);
        }
    }

    // avoid unused warnings
    UNUSED(lookfrom);
    UNUSED(lookat);
    UNUSED(up);
    UNUSED(vfov);


    Image3 img(640 /* width */, 480 /* height */);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;


    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            
            // Map pixel coordinates to image plane coordinates
            Real pixel_x = (x + Real(0.5)) / img.width;
            Real pixel_y = (y + Real(0.5)) / img.height;

            // Calculate camera ray direction
            Vector3 direction = (2.0 * pixel_x - 1.0) * half_width * u +
                                (1.0 - 2.0 * pixel_y) * half_height * v -
                                w;

            ray r(lookfrom, normalize(direction));
            Vector3 col = ray_color(r);
            img(x, y) = col;
        }
    }



    return img;
}


double hit_multiple_spheres1(const Vector3& center, double radius, const ray& r) {

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


Image3 hw_1_4(const std::vector<std::string> &params) {
    // Homework 1.4: render the scenes defined in hw1_scenes.h
    // output their diffuse color directly.
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = std::stoi(params[0]);
    UNUSED(scene_id); // avoid unused warning
    // Your scene is hw1_scenes[scene_id]

    Scene scene = hw1_scenes[scene_id];

    // Extract camera information from the scene
    Camera camera = scene.camera;
    Vector3 lookfrom = camera.lookfrom;
    Vector3 lookat = camera.lookat;
    Vector3 up = camera.up;
    Real vfov = camera.vfov;

    std::vector<Sphere> spheres = scene.shapes;
    std::vector<Material> materials = scene.materials;
    
    int num_spheres = spheres.size();


    Image3 img(640 /* width */, 480 /* height */);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;


    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            
            // Map pixel coordinates to image plane coordinates
            Real pixel_x = (x + Real(0.5)) / img.width;
            Real pixel_y = (y + Real(0.5)) / img.height;

            // Calculate camera ray direction
            Vector3 direction = (2.0 * pixel_x - 1.0) * half_width * u +
                                (1.0 - 2.0 * pixel_y) * half_height * v -
                                w;

            ray r(lookfrom, normalize(direction));

            Vector3 hit_color(0.5, 0.5, 0.5); // Default color if no intersection

            const Sphere* closestSphere = nullptr;
            double hitDistance = std::numeric_limits<Real>::max();

            // Iterate over all spheres in the scene
            for (int i = 0; i < num_spheres; i++) {
                const auto& sphere = spheres[i];
                double t = hit_multiple_spheres1(sphere.center, sphere.radius, r);
                if (t > 0 && t < hitDistance) {
                    hitDistance = t;
                    closestSphere = &sphere;
                }
            }

            if (closestSphere != nullptr) {
                hit_color = scene.materials[closestSphere->material_id].color;
            }

            img(x, y) = hit_color;
        }
    }
    return img;
}

Image3 hw_1_5(const std::vector<std::string> &params) {
    // Homework 1.5: render the scenes defined in hw1_scenes.h,
    // light them using the point lights in the scene.
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = std::stoi(params[0]);
    UNUSED(scene_id); // avoid unused warning
    // Your scene is hw1_scenes[scene_id]

    Scene scene = hw1_scenes[scene_id];
    double EPSILON = 0.0001;

    // Extract camera information from the scene
    Camera camera = scene.camera;
    Vector3 lookfrom = camera.lookfrom;
    Vector3 lookat = camera.lookat;
    Vector3 up = camera.up;
    Real vfov = camera.vfov;

    std::vector<Sphere> spheres = scene.shapes;
    std::vector<Material> materials = scene.materials;
    
    int num_spheres = spheres.size();


    Image3 img(640 /* width */, 480 /* height */);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;


    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            
            // Map pixel coordinates to image plane coordinates
            Real pixel_x = (x + Real(0.5)) / img.width;
            Real pixel_y = (y + Real(0.5)) / img.height;

            // Calculate camera ray direction
            Vector3 direction = (2.0 * pixel_x - 1.0) * half_width * u +
                                (1.0 - 2.0 * pixel_y) * half_height * v -
                                w;

            ray r(lookfrom, normalize(direction));

            Vector3 hit_color(0.5, 0.5, 0.5); // Default color if no intersection

            const Sphere* closestSphere = nullptr;
            double hitDistance = std::numeric_limits<Real>::max();

            // Iterate over all spheres in the scene
            for (int i = 0; i < num_spheres; i++) {
                const auto& sphere = spheres[i];
                double t = hit_multiple_spheres1(sphere.center, sphere.radius, r);
                if (t > 0 && t < hitDistance) {
                    hitDistance = t;
                    closestSphere = &sphere;
                }
            }

            if (closestSphere != nullptr) {
                Vector3 shading_point = r.origin() + hitDistance * r.direction();
                Vector3 normal = normalize(shading_point - closestSphere->center);
                Material material = materials[closestSphere->material_id];

                Vector3 hit_color(0.0, 0.0, 0.0);

                for (const auto& light : scene.lights) {
                    Vector3 light_dir = normalize(light.position - shading_point);
                    double visibility = 1.0;
                    ray shadow_ray(shading_point + EPSILON * light_dir, light_dir);
                    // Calculate the distance between shading point and light source
                    double d = distance(shading_point, light.position);

                    for (int i = 0; i < num_spheres; i++) {
                        if (&spheres[i] == closestSphere) continue; // Skip the sphere we hit
                        double t = hit_sphere(spheres[i].center, spheres[i].radius, shadow_ray);
                        if (t > 0 && t < distance(shadow_ray.origin(), light.position)) {
                            // There is an intersection between the shadow ray and other spheres,
                            // hence the point is in shadow
                            visibility = 0.0;
                            break;
                        }
                    }

                    // Calculate diffuse color
                    Vector3 diffuse_color = material.color * light.intensity *
                                            std::max(0.0, dot(normal, light_dir));
                    diffuse_color /= (c_PI * d * d); // Divide by pi * d^2

                    hit_color += visibility * diffuse_color;
                }

                img(x, y) = hit_color;

            } else {
                img(x, y) = Vector3(0.5, 0.5, 0.5);
            }
        }     
    }

    return img;
}

Image3 hw_1_6(const std::vector<std::string> &params) {
    // Homework 1.6: add antialiasing to homework 1.5
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = 0;
    int spp = 64;

    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            scene_id = std::stoi(params[i]);
        }
    }

    Scene scene = hw1_scenes[scene_id];
    double EPSILON = 0.0001;

    // Extract camera information from the scene
    Camera camera = scene.camera;
    Vector3 lookfrom = camera.lookfrom;
    Vector3 lookat = camera.lookat;
    Vector3 up = camera.up;
    Real vfov = camera.vfov;

    std::vector<Sphere> spheres = scene.shapes;
    std::vector<Material> materials = scene.materials;
    
    int num_spheres = spheres.size();


    Image3 img(160 /* width */, 120 /* height */);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;

    // Create PCG random number generator with a fixed seed
    pcg32_state rng;
    rng.state = 42u;
    rng.inc = 54u;

    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            Vector3 hit_color(0.0, 0.0, 0.0); // Accumulated color for antialiasing

            for (int s = 0; s < spp; s++) {
                // Map pixel coordinates to image plane coordinates with jittering
                Real jitter_x = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
                Real jitter_y = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
            
            // Map pixel coordinates to image plane coordinates
            Real pixel_x = (x + Real(0.5) + jitter_x) / img.width;
            Real pixel_y = (y + Real(0.5) + jitter_y) / img.height;

            // Calculate camera ray direction
            Vector3 direction = (2.0 * pixel_x - 1.0) * half_width * u +
                                (1.0 - 2.0 * pixel_y) * half_height * v -
                                w;

            ray r(lookfrom, normalize(direction));

            Vector3 hit_color(0.5, 0.5, 0.5); // Default color if no intersection

            const Sphere* closestSphere = nullptr;
            double hitDistance = std::numeric_limits<Real>::max();

            // Iterate over all spheres in the scene
            for (int i = 0; i < num_spheres; i++) {
                const auto& sphere = spheres[i];
                double t = hit_multiple_spheres1(sphere.center, sphere.radius, r);
                if (t > 0 && t < hitDistance) {
                    hitDistance = t;
                    closestSphere = &sphere;
                }
            }

            if (closestSphere != nullptr) {
                Vector3 shading_point = r.origin() + hitDistance * r.direction();
                Vector3 normal = normalize(shading_point - closestSphere->center);
                Material material = materials[closestSphere->material_id];

                Vector3 hit_color(0.0, 0.0, 0.0);

                for (const auto& light : scene.lights) {
                    Vector3 light_dir = normalize(light.position - shading_point);
                    double visibility = 1.0;
                    ray shadow_ray(shading_point + EPSILON * light_dir, light_dir);
                    // Calculate the distance between shading point and light source
                    double d = distance(shading_point, light.position);

                    for (int i = 0; i < num_spheres; i++) {
                        if (&spheres[i] == closestSphere) continue; // Skip the sphere we hit
                        double t = hit_sphere(spheres[i].center, spheres[i].radius, shadow_ray);
                        if (t > 0 && t < distance(shadow_ray.origin(), light.position)) {
                            // There is an intersection between the shadow ray and other spheres,
                            // hence the point is in shadow
                            visibility = 0.0;
                            break;
                        }
                    }

                    // Calculate diffuse color
                    Vector3 diffuse_color = material.color * light.intensity *
                                            std::max(0.0, dot(normal, light_dir));
                    diffuse_color /= (c_PI * d * d); // Divide by pi * d^2

                    hit_color += visibility * diffuse_color;
                }

                img(x, y) += hit_color / static_cast<double>(spp); // Add to accumulated color

            } else {
                img(x, y) += Vector3(0.5, 0.5, 0.5) / static_cast<double>(spp); // Add to accumulated color
            }
            //hit_color = hit_color / static_cast<double>(spp);
            }
        }     
    }

    return img;
}

Vector3 trace_ray(const ray& r, const Scene& scene, const std::vector<Material>& materials, const std::vector<Sphere>& spheres) {
    const Sphere* closestSphere = nullptr;
    double hitDistance = std::numeric_limits<Real>::max();
    double EPSILON = 0.0001;

    // Iterate over all spheres in the scene
    for (int i = 0; i < spheres.size(); i++) {
        const auto& sphere = spheres[i];
        double t = hit_multiple_spheres1(sphere.center, sphere.radius, r);
        if (t > 0 && t < hitDistance) {
            hitDistance = t;
            closestSphere = &sphere;
        }
    }

    if (closestSphere != nullptr) {
        Vector3 shading_point = r.origin() + hitDistance * r.direction();
        Vector3 normal = normalize(shading_point - closestSphere->center);
        Material material = materials[closestSphere->material_id];

        Vector3 hit_color(0.0, 0.0, 0.0);

        if (material.type == MaterialType::Mirror) {
            // Calculate reflection direction and create a reflected ray
            Vector3 reflection_dir = r.direction() - 2.0 * dot(r.direction(), normal) * normal;
            ray reflected_ray(shading_point + EPSILON * reflection_dir, reflection_dir);

            // Call trace_ray() recursively
            return material.color * trace_ray(reflected_ray, scene, materials, spheres);
        } else {
            for (const auto& light : scene.lights) {
                Vector3 light_dir = normalize(light.position - shading_point);
                double visibility = 1.0;
                ray shadow_ray(shading_point + EPSILON * light_dir, light_dir);
                // Calculate the distance between shading point and light source
                double d = distance(shading_point, light.position);

                for (int i = 0; i < spheres.size(); i++) {
                    if (&spheres[i] == closestSphere) continue; // Skip the sphere we hit
                    double t = hit_sphere(spheres[i].center, spheres[i].radius, shadow_ray);
                    if (t > 0 && t < distance(shadow_ray.origin(), light.position)) {
                        // There is an intersection between the shadow ray and other spheres,
                        // hence the point is in shadow
                        visibility = 0.0;
                        break;
                    }
                }

                // Calculate diffuse color
                Vector3 diffuse_color = material.color * light.intensity *
                                        std::max(0.0, dot(normal, light_dir));
                diffuse_color /= (c_PI * d * d); // Divide by pi * d^2

                hit_color += visibility * diffuse_color;
            }

            return hit_color;
        }
    } else {
        return Vector3(0.5, 0.5, 0.5); // Default background color
    }
}

Image3 hw_1_7(const std::vector<std::string> &params) {
    // Homework 1.7: add mirror materials to homework 1.6
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = 0;
    int spp = 64;

    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            scene_id = std::stoi(params[i]);
        }
    }

   Scene scene = hw1_scenes[scene_id];
    double EPSILON = 0.0001;

    // Extract camera information from the scene
    Camera camera = scene.camera;
    Vector3 lookfrom = camera.lookfrom;
    Vector3 lookat = camera.lookat;
    Vector3 up = camera.up;
    Real vfov = camera.vfov;

    std::vector<Sphere> spheres = scene.shapes;
    std::vector<Material> materials = scene.materials;
    
    int num_spheres = spheres.size();


    Image3 img(640 /* width */, 480 /* height */);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;

    // Create PCG random number generator with a fixed seed
    pcg32_state rng;
    rng.state = 42u;
    rng.inc = 54u;

    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            Vector3 hit_color(0.0, 0.0, 0.0); // Accumulated color for antialiasing

            for (int s = 0; s < spp; s++) {
                // Map pixel coordinates to image plane coordinates with jittering
                Real jitter_x = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
                Real jitter_y = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
            
            // Map pixel coordinates to image plane coordinates
            Real pixel_x = (x + Real(0.5) + jitter_x) / img.width;
            Real pixel_y = (y + Real(0.5) + jitter_y) / img.height;

            // Calculate camera ray direction
            Vector3 direction = (2.0 * pixel_x - 1.0) * half_width * u +
                                (1.0 - 2.0 * pixel_y) * half_height * v -
                                w;

            ray r(lookfrom, normalize(direction));


            Vector3 sample_color = trace_ray(r, scene, materials, spheres);
            img(x, y) += sample_color / static_cast<double>(spp);

            }
        }     
    }

    return img;
}

Image3 hw_1_8(const std::vector<std::string> &params) {
    // Homework 1.8: parallelize HW 1.7
     if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = 0;
    int spp = 64;

    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            scene_id = std::stoi(params[i]);
        }
    }

   Scene scene = hw1_scenes[scene_id];
    double EPSILON = 0.0001;

    // Extract camera information from the scene
    Camera camera = scene.camera;
    Vector3 lookfrom = camera.lookfrom;
    Vector3 lookat = camera.lookat;
    Vector3 up = camera.up;
    Real vfov = camera.vfov;

    std::vector<Sphere> spheres = scene.shapes;
    std::vector<Material> materials = scene.materials;
    
    int num_spheres = spheres.size();

    Image3 img(1280 /* width */, 960 /* height */);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;

    constexpr int tile_size = 16;
    int num_tiles_x = (img.width + tile_size - 1) / tile_size;
    int num_tiles_y = (img.height + tile_size - 1) / tile_size;

    parallel_for([&](const Vector2i &tile) {
        pcg32_state rng;
        rng.state = 42u + tile[1] * num_tiles_x + tile[0];
        rng.inc = 54u;

        int x0 = tile[0] * tile_size;
        int x1 = std::min(x0 + tile_size, img.width);
        int y0 = tile[1] * tile_size;
        int y1 = std::min(y0 + tile_size, img.height);

        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Vector3 hit_color(0.0, 0.0, 0.0); // Accumulated color for antialiasing

                for (int s = 0; s < spp; s++) {
                    // Map pixel coordinates to image plane coordinates with jittering
                    Real jitter_x = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
                    Real jitter_y = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
                
                    // Map pixel coordinates to image plane coordinates
                    Real pixel_x = (x + Real(0.5) + jitter_x) / img.width;
                    Real pixel_y = (y + Real(0.5) + jitter_y) / img.height;

                    // Calculate camera ray direction
                    Vector3 direction = (2.0 * pixel_x - 1.0) * half_width * u +
                                        (1.0 - 2.0 * pixel_y) * half_height * v -
                                        w;

                    ray r(lookfrom, normalize(direction));


                    Vector3 sample_color = trace_ray(r, scene, materials, spheres);
                    img(x, y) += sample_color / static_cast<double>(spp);

                }
            }     
        }
    }, Vector2i(num_tiles_x, num_tiles_y));
    return img;
}

// Helper function to check if a sphere overlaps with any existing sphere
bool is_overlapping(const Vector3& center, Real radius, const std::vector<Sphere>& spheres) {
    for (const auto& sphere : spheres) {
        Real min_distance = sphere.radius + radius;
        if (length(sphere.center - center) < min_distance) {
            return true;
        }
    }
    return false;
}

// Same as trace_ray, but returns a different background color if not hit (for hw1_9)
Vector3 trace_ray_2(const ray& r, const Scene& scene, const std::vector<Material>& materials, const std::vector<Sphere>& spheres) {
    const Sphere* closestSphere = nullptr;
    double hitDistance = std::numeric_limits<Real>::max();
    double EPSILON = 0.0001;

    // Iterate over all spheres in the scene
    for (int i = 0; i < spheres.size(); i++) {
        const auto& sphere = spheres[i];
        double t = hit_multiple_spheres1(sphere.center, sphere.radius, r);
        if (t > 0 && t < hitDistance) {
            hitDistance = t;
            closestSphere = &sphere;
        }
    }

    if (closestSphere != nullptr) {
        Vector3 shading_point = r.origin() + hitDistance * r.direction();
        Vector3 normal = normalize(shading_point - closestSphere->center);
        Material material = materials[closestSphere->material_id];

        Vector3 hit_color(0.0, 0.0, 0.0);

        if (material.type == MaterialType::Mirror) {
            // Calculate reflection direction and create a reflected ray
            Vector3 reflection_dir = r.direction() - 2.0 * dot(r.direction(), normal) * normal;
            ray reflected_ray(shading_point + EPSILON * reflection_dir, reflection_dir);

            // Call trace_ray() recursively
            return material.color * trace_ray(reflected_ray, scene, materials, spheres);
        } else {
            for (const auto& light : scene.lights) {
                Vector3 light_dir = normalize(light.position - shading_point);
                double visibility = 1.0;
                ray shadow_ray(shading_point + EPSILON * light_dir, light_dir);
                // Calculate the distance between shading point and light source
                double d = distance(shading_point, light.position);

                for (int i = 0; i < spheres.size(); i++) {
                    if (&spheres[i] == closestSphere) continue; // Skip the sphere we hit
                    double t = hit_sphere(spheres[i].center, spheres[i].radius, shadow_ray);
                    if (t > 0 && t < distance(shadow_ray.origin(), light.position)) {
                        // There is an intersection between the shadow ray and other spheres,
                        // hence the point is in shadow
                        visibility = 0.0;
                        break;
                    }
                }

                // Calculate diffuse color
                Vector3 diffuse_color = material.color * light.intensity *
                                        std::max(0.0, dot(normal, light_dir));
                diffuse_color /= (c_PI * d * d); // Divide by pi * d^2

                hit_color += visibility * diffuse_color;
            }

            return hit_color;
        }
    } else {
        return Vector3(0.78, 0.83, 0.98); // Default background color
    }
}

Vector3 random_in_unit_disk(pcg32_state &rng) {
    Vector3 p;
    do {
        p = 2.0 * Vector3(static_cast<Real>(next_pcg32(rng)) / float(UINT32_MAX), static_cast<Real>(next_pcg32(rng) / float(UINT32_MAX)), 0.0) - Vector3(1.0, 1.0, 0.0);
    } while (dot(p, p) >= 1.0);
    return p;
}


// Design my own scene
// type "./torrey -hw 1_9 5" to run

// Big two Spheres in the middle with mirror material
// Randomly generated spheres with 30+ different materials
// The center for the random spheres is randomly generated,
// but the radius is fixed to 0.5, and I made sure that the spheres
// are not overlapping with each other

// I also implemented Defocus Blur.
// I haven't yet found a perfect value for the scene
// because I was out of time. However, I believe that the implementation
// is mostly successful.

Image3 hw_1_9(const std::vector<std::string> &params) {
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = 5;
    int spp = 64;

    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            scene_id = std::stoi(params[i]);
        }
    }

    Scene scene = hw1_scene_5;
    std::vector<Sphere>& spheres = scene.shapes;
    std::vector<Material> materials = scene.materials;

    // Set up the random number generator
    std::mt19937 rng(42); // Seed with a fixed value for reproducibility
    std::uniform_real_distribution<Real> choose_mat(0.0, 1.0);
    std::uniform_real_distribution<Real> center_x_dist(-10.0, 1.0);
    std::uniform_real_distribution<Real> center_z_dist(-10.0, 1.0);
    std::uniform_int_distribution<int> choose_color(3, 32);

    std::uniform_real_distribution<Real> rgb_dist(0.0, 1.0); // Range of RGB values

    constexpr int num_spheres = 20;
    constexpr Real sphere_radius = 0.2;

    spheres.reserve(spheres.size() + num_spheres*num_spheres);



    for (int x = 0; x < num_spheres; x++) {
        for (int z = 0; z < num_spheres; z++) {
            auto choose_mat_val = choose_mat(rng);
            Vector3 center;
            do {
                center = Vector3(
                    x + 1.1 * center_x_dist(rng),
                    0.2,
                    z + 1.1 * center_z_dist(rng)
                );
            } while (is_overlapping(center, sphere_radius, spheres));

                spheres.push_back({center, sphere_radius, choose_color(rng)});

        }
        
    }
    double EPSILON = 0.0001;

    // Extract camera information from the scene
    Camera camera = scene.camera;
    Vector3 lookfrom = camera.lookfrom;
    Vector3 lookat = camera.lookat;
    Vector3 up = camera.up;
    Real vfov = camera.vfov;
    
    Real aperture = 0.045;
    Real focus_dist = 6.0;
    //Real lens_radius = 0.1; // Controls the amount of defocus blur
    Real lens_radius = aperture / 2.0;

    Image3 img(1080 /* width */, 1350 /* height */);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;

    constexpr int tile_size = 16;
    int num_tiles_x = (img.width + tile_size - 1) / tile_size;
    int num_tiles_y = (img.height + tile_size - 1) / tile_size;

    parallel_for([&](const Vector2i &tile) {
        pcg32_state rng;
        rng.state = 42u + tile[1] * num_tiles_x + tile[0];
        rng.inc = 54u;

        int x0 = tile[0] * tile_size;
        int x1 = std::min(x0 + tile_size, img.width);
        int y0 = tile[1] * tile_size;
        int y1 = std::min(y0 + tile_size, img.height);

        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Vector3 hit_color(0.0, 0.0, 0.0); // Accumulated color for antialiasing

                for (int s = 0; s < spp; s++) {
                    // Map pixel coordinates to image plane coordinates with jittering
                    Real jitter_x = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
                    Real jitter_y = (next_pcg32(rng) / float(UINT32_MAX)) - 0.5;
                
                    // Map pixel coordinates to image plane coordinates
                    Real pixel_x = (x + Real(0.5) + jitter_x) / img.width;
                    Real pixel_y = (y + Real(0.5) + jitter_y) / img.height;

                    Vector3 rd = lens_radius * random_in_unit_disk(rng);
                    Vector3 offset = u * rd.x + v * rd.y;
                    Vector3 new_lookfrom = lookfrom + offset;
                    Vector3 direction = (2.0 * pixel_x - 1.0) * half_width * u * focus_dist +
                                        (1.0 - 2.0 * pixel_y) * half_height * v * focus_dist -
                                        w * focus_dist - offset;
                    ray r(new_lookfrom, normalize(direction));


                    Vector3 sample_color = trace_ray_2(r, scene, materials, spheres);
                    img(x, y) += sample_color / static_cast<double>(spp);

                }
            }     
        }
    }, Vector2i(num_tiles_x, num_tiles_y));
    return img;
}
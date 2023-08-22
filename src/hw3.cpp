#include "hw3.h"
#include "parse_scene.h"
#include "timer.h"
#include "ray.h"
#include "vector.h"
#include "pcg.h"
#include <random>
#include "scene.h"
#include "parse_obj.h"
#include "progressreporter.h"
#include "aabb.h"
#include "bvh.h"
#include "parallel.h"
#include "camera.h"
#include "radiance.h"

static Vector3 radiance_hw3(const Scene &scene, ray input_ray, bool use_shading_normal) {
    Vector3 color = scene.background_color;
    if (auto hit_isect = intersect(scene, input_ray)) {
        Vector3 p = hit_isect->position;
        Vector3 n = use_shading_normal ?
            hit_isect->shading_normal : hit_isect->geometric_normal;
        if (dot(-input_ray.dir, n) < 0) {
            n = -n;
        }
        const Material &mat = scene.materials[hit_isect->material_id];
        if (auto *diffuse = std::get_if<Diffuse>(&mat)) {
            Vector3 L = Vector3{0, 0, 0};
            // Loop over the lights
            for (const Light &light : scene.lights) {
                // Assume point lights for now.
                const PointLight &point_light = std::get<PointLight>(light);
                Vector3 l = point_light.position - p;
                ray shadow_ray{p, normalize(l), Real(1e-4), (1 - Real(1e-4)) * length(l)};
                if (!occluded(scene, shadow_ray)) {
                    L += (max(dot(n, normalize(l)), Real(0)) / c_PI) *
                         (point_light.intensity / length_squared(l)) *
                         eval(diffuse->reflectance, hit_isect->uv);
                }
            }
            color = L;
        } else if (auto *mirror = std::get_if<Mirror>(&mat)) {
            ray refl_ray{p, input_ray.dir - 2 * dot(input_ray.dir, n) * n, Real(1e-4), infinity<Real>()};
            return eval(mirror->reflectance, hit_isect->uv) *
                radiance_hw3(scene, refl_ray, use_shading_normal);
        }
    }
    return color;
}

Image3 hw_3_1(const std::vector<std::string> &params) {
    // Homework 3.1: image textures
    if (params.size() < 1) {
        return Image3(0, 0);
    }

    Timer timer;
    tick(timer);
    ParsedScene parsed_scene = parse_scene(params[0]);
    std::cout << "Scene parsing done. Took " << tick(timer) << " seconds." << std::endl;

    tick(timer);
    Scene scene(parsed_scene);
    std::cout << "Scene construction done. Took " << tick(timer) << " seconds." << std::endl;


    int spp = 16;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        }
    }

    // Create PCG random number generator with a fixed seed
    pcg32_state rng = init_pcg32();

    Vector3 lookfrom = scene.camera.lookfrom;
    Vector3 lookat   = scene.camera.lookat;
    Vector3 up       = scene.camera.up;
    Real    vfov     = scene.camera.vfov;

    // Create an output image with the correct dimensions
    Image3 img(scene.width, scene.height);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;

    ProgressReporter reporter(img.width * img.height);

    tick(timer);
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
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

                ray r(lookfrom, normalize(direction), Real(0), infinity<Real>());
                Vector3 col;
                
                col = radiance_hw3(scene, r, false  /*use_shading_normal*/);
            
                img(x, y) += col / Real(spp); // Add to accumulated color
            }
                    reporter.update(1);
        }

    }


    reporter.done();
    std::cout << "Rendering done. Took " << tick(timer) << " seconds." << std::endl;
    return img;
}

Image3 hw_3_2(const std::vector<std::string> &params) {
    // Homework 3.2: shading normals
    if (params.size() < 1) {
        return Image3(0, 0);
    }
    Timer timer;
    tick(timer);
    ParsedScene parsed_scene = parse_scene(params[0]);
    std::cout << "Scene parsing done. Took " << tick(timer) << " seconds." << std::endl;

    tick(timer);
    Scene scene(parsed_scene);
    std::cout << "Scene construction done. Took " << tick(timer) << " seconds." << std::endl;


    int spp = 16;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        }
    }

    // Create PCG random number generator with a fixed seed
    pcg32_state rng = init_pcg32();

    Vector3 lookfrom = scene.camera.lookfrom;
    Vector3 lookat   = scene.camera.lookat;
    Vector3 up       = scene.camera.up;
    Real    vfov     = scene.camera.vfov;

    // Create an output image with the correct dimensions
    Image3 img(scene.width, scene.height);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;

    ProgressReporter reporter(img.width * img.height);

    tick(timer);
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
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

                ray r(lookfrom, normalize(direction), Real(0), infinity<Real>());
                Vector3 col;
                
                col = radiance_hw3(scene, r,
                        true /*use_shading_normal*/);
            
                img(x, y) += col / Real(spp); // Add to accumulated color
            }
                    reporter.update(1);
        }

    }


    reporter.done();
    std::cout << "Rendering done. Took " << tick(timer) << " seconds." << std::endl;
    return img;
}

Image3 hw_3_3(const std::vector<std::string> &params) {
if (params.size() < 1) {
        return Image3(0, 0);
    }

    Timer timer;
    tick(timer);
    ParsedScene parsed_scene = parse_scene(params[0]);
    std::cout << "Scene parsing done. Took " << tick(timer) << " seconds." << std::endl;

    tick(timer);
    Scene scene(parsed_scene);
    std::cout << "Scene construction done. Took " << tick(timer) << " seconds." << std::endl;
    int spp = scene.samples_per_pixel;

    Image3 img(scene.width, scene.height);

    int w = img.width;
    int h = img.height;

    CameraRayData cam_ray_data = compute_camera_ray_data(
        scene.camera, w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    tick(timer);
    parallel_for([&](const Vector2i &tile) {
        // Use a different rng stream for each thread.
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                for (int s = 0; s < spp; s++) {
                    Real u, v;
                    u = (x + next_pcg32_real<Real>(rng)) / w;
                    v = (y + next_pcg32_real<Real>(rng)) / h;

                    ray ray = generate_primary_ray(cam_ray_data, u, v);
                    img(x, y) += radiance_new(scene, ray, rng) / Real(spp);
                }
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    std::cout << "Rendering done. Took " << tick(timer) << " seconds." << std::endl;

    return img;
}

///// since I failed to complete 3_4, I tried to make the custom scene first
Image3 hw_3_4(const std::vector<std::string> &params) {
    // Homework 3.4: area lights
    if (params.size() < 1) {
        return Image3(0, 0);
    }

    Timer timer;
    tick(timer);
    ParsedScene parsed_scene = parse_scene(params[0]);
    std::cout << "Scene parsing done. Took " << tick(timer) << " seconds." << std::endl;

    tick(timer);
    Scene scene(parsed_scene);
    std::cout << "Scene construction done. Took " << tick(timer) << " seconds." << std::endl;


    int spp = 16;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        }
    }

    // Create PCG random number generator with a fixed seed
    pcg32_state rng = init_pcg32();

    Vector3 lookfrom = scene.camera.lookfrom;
    Vector3 lookat   = scene.camera.lookat;
    Vector3 up       = scene.camera.up;
    Real    vfov     = scene.camera.vfov;

    // Create an output image with the correct dimensions
    Image3 img(scene.width, scene.height);

    float aspect_ratio = static_cast<float>(img.width) / static_cast<float>(img.height);

    Vector3 w = normalize(lookfrom - lookat);
    Vector3 u = normalize(cross(up, w));
    Vector3 v = cross(w, u);

    // Convert vfov to radians
    Real theta = radians(vfov);

    // Calculate image plane height and width
    Real half_height = tan(theta / 2);
    Real half_width = aspect_ratio * half_height;

    ProgressReporter reporter(img.width * img.height);

    tick(timer);
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
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

                ray r(lookfrom, normalize(direction), Real(0), infinity<Real>());
                Vector3 col;
                
                col = radiance_new(scene, r, rng);
            
                img(x, y) += col / Real(spp); // Add to accumulated color
            }
                    reporter.update(1);
        }

    }


    reporter.done();
    std::cout << "Rendering done. Took " << tick(timer) << " seconds." << std::endl;
    return img;
}

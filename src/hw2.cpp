#include "hw2.h"
#include "parse_scene.h"
#include "print_scene.h"
#include "timer.h"
#include "ray.h"
#include "vector.h"
#include "pcg.h"
#include <random>
#include "scene.h"
//#include "traceRay.h"
#include "parse_obj.h"
#include "progressreporter.h"
//#include "radiance.h"
#include "aabb.h"
#include "bvh.h"
#include "camera.h"



Vector3 ray_color_triangle(const ray& r, const Vector3& p0, const Vector3& p1, const Vector3& p2) {
    
    Real t, u, v;
    if (intersect_triangle(r, p0, p1, p2, t, u, v)) {
        // Define vertex colors
        Vector3 col_top = Vector3{1, 0, 0};    // Red
        Vector3 col_left = Vector3{0, 1, 0};   // Green
        Vector3 col_right = Vector3{0, 0, 1};  // Blue

        // Interpolate colors based on barycentric coordinates
        Vector3 col = col_top * (1 - u - v) + col_left * u + col_right * v;

        return col;
    }

    // Background color
    return Vector3{0.5, 0.5, 0.5};
}

static Vector3 radiance_old(const Scene &scene, ray input_ray) {
    Vector3 color = Vector3{0.5, 0.5, 0.5};
    if (auto hit_isect = intersect(scene, input_ray)) {
        Vector3 p = hit_isect->position;
        Vector3 n = hit_isect->geometric_normal;
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
                    ConstantTexture c = std::get<ConstantTexture>(diffuse->reflectance);
                    L += (max(dot(n, normalize(l)), Real(0)) / c_PI) *
                         (point_light.intensity / length_squared(l)) *
                         c.color;
                }
            }
            color = L;
        } else if (auto *mirror = std::get_if<Mirror>(&mat)) {
            ConstantTexture c = std::get<ConstantTexture>(mirror->reflectance);
            ray refl_ray{p, input_ray.dir - 2 * dot(input_ray.dir, n) * n, Real(1e-4), infinity<Real>()};
            return c.color * radiance_old(scene, refl_ray);
        }
    }
    return color;
}

Image3 hw_2_1(const std::vector<std::string> &params) {
    // Homework 2.1: render a single triangle and outputs
    // its barycentric coordinates.
    // We will use the following camera parameter
    // lookfrom = (0, 0,  0)
    // lookat   = (0, 0, -1)
    // up       = (0, 1,  0)
    // vfov     = 45
    // and we will parse the triangle vertices from params
    // The three vertices are stored in v0, v1, and v2 below.


    std::vector<float> tri_params;
    int spp = 16;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            tri_params.push_back(std::stof(params[i]));
        }
    }

    if (tri_params.size() < 9) {
        // Not enough parameters to parse the triangle vertices.
        return Image3(0, 0);
    }

    // Create PCG random number generator with a fixed seed
    pcg32_state rng = init_pcg32();

    Vector3 lookfrom = Vector3{0, 0,  0};
    Vector3 lookat   = Vector3{0, 0, -1};
    Vector3 up       = Vector3{0, 1,  0};
    Real    vfov     = 45;
 
    Image3 img(640 /* width */, 480 /* height */);

    Vector3 p0{tri_params[0], tri_params[1], tri_params[2]};
    Vector3 p1{tri_params[3], tri_params[4], tri_params[5]};
    Vector3 p2{tri_params[6], tri_params[7], tri_params[8]};

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
                Vector3 col = ray_color_triangle(r, p0, p1, p2);
                img(x, y) += col / Real(spp); // Add to accumulated color
            }
        }
    }



    return img;
}

Vector3 ray_color_triangle_mesh(const ray& r, const std::vector<Vector3>& positions, const std::vector<Vector3i>& indices) {
    Real closest_t = std::numeric_limits<Real>::max();
    Vector3 col;
    bool intersected = false;

    for (const auto& tri_idx : indices) {
        Vector3 p0 = positions[tri_idx[0]];
        Vector3 p1 = positions[tri_idx[1]];
        Vector3 p2 = positions[tri_idx[2]];

        Real t, u, v;
        if (intersect_triangle(r, p0, p1, p2, t, u, v)) {
            if (t < closest_t) {
                closest_t = t;
                col = Vector3{1 - u - v, u, v};
                intersected = true;
            }
        }
    }

    if (intersected) {
        return col;
    }

    // Background color
    return Vector3{0.5, 0.5, 0.5};
}

Image3 hw_2_2(const std::vector<std::string> &params) {
    // Homework 2.2: render a triangle mesh.
    // We will use the same camera parameter:
    // lookfrom = (0, 0,  0)
    // lookat   = (0, 0, -1)
    // up       = (0, 1,  0)
    // vfov     = 45
    // and we will use a fixed triangle mesh: a tetrahedron!

    int spp = 16;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        }
    }

    // Create PCG random number generator with a fixed seed
    pcg32_state rng = init_pcg32();

    Vector3 lookfrom = Vector3{0, 0,  0};
    Vector3 lookat   = Vector3{0, 0, -1};
    Vector3 up       = Vector3{0, 1,  0};
    Real    vfov     = 45;
 
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


    std::vector<Vector3> positions = {
        Vector3{ 0.0,  0.5, -2.0},
        Vector3{ 0.0, -0.3, -1.0},
        Vector3{ 1.0, -0.5, -3.0},
        Vector3{-1.0, -0.5, -3.0}
    };
    std::vector<Vector3i> indices = {
        Vector3i{0, 1, 2},
        Vector3i{0, 3, 1},
        Vector3i{0, 2, 3},
        Vector3i{1, 2, 3}
    };

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

                ray r(lookfrom, normalize(direction));
                Vector3 col = ray_color_triangle_mesh(r, positions, indices);
                img(x, y) += col / Real(spp); // Add to accumulated color
            }
        }
    }
    
    return img;
}


Image3 hw_2_3(const std::vector<std::string> &params) {
    // Homework 2.3: render a scene file provided by our parser.
    if (params.size() < 1) {
        return Image3(0, 0);
    }

    int spp = 16;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        }
    }

    Timer timer;
    tick(timer);
    ParsedScene parsed_scene = parse_scene(params[0]);
    std::cout << "Scene parsing done. Took " << tick(timer) << " seconds." << std::endl;

    tick(timer);
    Scene scene(parsed_scene);
    std::cout << "Scene construction done. Took " << tick(timer) << " seconds." << std::endl;


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

    // Check if there are triangle meshes in the scene
    bool has_triangle_meshes = false;
    for (const auto& shape : scene.shapes) {
        if (std::holds_alternative<Triangle>(shape)) {
            has_triangle_meshes = true;
            break;
        }
    }

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
                
                col = radiance_old(scene, r);
            
                img(x, y) += col / Real(spp); // Add to accumulated color
            }
            reporter.update(1);
        }
    }

    reporter.done();
    std::cout << "Rendering done. Took " << tick(timer) << " seconds." << std::endl;

    return img;
}

Image3 hw_2_4(const std::vector<std::string> &params) {
    // Homework 2.4: render the AABBs of the scene.
    if (params.size() < 1) {
        return Image3(0, 0);
    }

    int spp = 16;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        }
    }

    Timer timer;
    tick(timer);
    ParsedScene parsed_scene = parse_scene(params[0]);
    std::cout << "Scene parsing done. Took " << tick(timer) << " seconds." << std::endl;

    tick(timer);
    Scene scene(parsed_scene);
    std::cout << "Scene construction done. Took " << tick(timer) << " seconds." << std::endl;


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


    std::vector<AABB> bboxes(scene.shapes.size());
    for (int i = 0; i < (int)bboxes.size(); i++) {
        if (auto *sph = std::get_if<Sphere>(&scene.shapes[i])) {
            Vector3 p_min = sph->position - sph->radius;
            Vector3 p_max = sph->position + sph->radius;
            bboxes[i] = AABB{p_min, p_max};
        } else if (auto *tri = std::get_if<Triangle>(&scene.shapes[i])) {
            const TriangleMesh *mesh = tri->mesh;
            Vector3i index = mesh->indices[tri->face_index];
            Vector3 p0 = mesh->positions[index[0]];
            Vector3 p1 = mesh->positions[index[1]];
            Vector3 p2 = mesh->positions[index[2]];
            Vector3 p_min = min(min(p0, p1), p2);
            Vector3 p_max = max(max(p0, p1), p2);
            bboxes[i] = AABB{p_min, p_max};
        }
    }

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
                    bool check_hit = false;
                    for (const AABB &bbox : bboxes) {
                        if (hit(bbox, r)) {
                            check_hit = true;
                            break;
                        }
                    }
                    if (check_hit) {
                        img(x, y) += Vector3{1, 1, 1} / Real(spp);
                    } else {
                        img(x, y) += Vector3{0.5, 0.5, 0.5} / Real(spp);
                    }
            }

        }
    }

    return img;
}

Image3 hw_2_5(const std::vector<std::string> &params) {
    // Homework 2.5: rendering with BVHs
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
                
                col = radiance_old(scene, r);
            
                img(x, y) += col / Real(spp); // Add to accumulated color
            }
                    reporter.update(1);
        }
    }


    reporter.done();
    std::cout << "Rendering done. Took " << tick(timer) << " seconds." << std::endl;
    return img;
}


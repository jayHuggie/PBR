#pragma once
#include "vector.h"
#include "matrix.h"
#include "torrey.h"
#include <variant>
#include "aabb.h"
#include <optional>
#include "bvh.h"
#include "image.h"
#include <cmath>
#include "pcg.h"
#include <random>
#include "texture.h"
#include "camera.h"
#include "shape.h"
#include "material.h"
#include "light.h"

struct Scene {
    Scene(const ParsedScene &scene);
    Camera camera;
    int width, height;
    std::vector<Shape> shapes;
    std::vector<Material> materials;
    std::vector<Light> lights;

    Vector3 background_color;
    int samples_per_pixel;
    // For the Triangle in the shapes to reference to.
    std::vector<TriangleMesh> meshes;

    std::vector<BVHNode> bvh_nodes;
    int bvh_root_id;

    bool hit(const ray& r, Real t_min, Real t_max, Intersection& rec) const;
};

std::optional<Intersection> intersect(const Scene &scene, ray ray);
inline bool occluded(const Scene &scene, const ray &ray);


inline Scene::Scene(const ParsedScene &scene) :
        camera(from_parsed_camera(scene.camera)),
        width(scene.camera.width),
        height(scene.camera.height),
        background_color(scene.background_color),
        samples_per_pixel(scene.samples_per_pixel) {
    // Extract triangle meshes from the parsed scene.
    int tri_mesh_count = 0;
    for (const ParsedShape &parsed_shape : scene.shapes) {
        if (std::get_if<ParsedTriangleMesh>(&parsed_shape)) {
            tri_mesh_count++;
        }
    }
    meshes.resize(tri_mesh_count);

    // Extract the shapes
    tri_mesh_count = 0;
    for (int i = 0; i < (int)scene.shapes.size(); i++) {
        const ParsedShape &parsed_shape = scene.shapes[i];
        if (auto *sph = std::get_if<ParsedSphere>(&parsed_shape)) {
            if (sph->area_light_id >= 0) {
                const ParsedDiffuseAreaLight &light =
                    std::get<ParsedDiffuseAreaLight>(
                        scene.lights[sph->area_light_id]);
                lights.push_back(DiffuseAreaLight{
                    (int)shapes.size(), light.radiance});
            }
            shapes.push_back(
                Sphere{sph->position,
                    sph->radius,
                    sph->material_id,
                    sph->area_light_id});
        } else if (auto *parsed_mesh = std::get_if<ParsedTriangleMesh>(&parsed_shape)) {
            meshes[tri_mesh_count] = TriangleMesh{
                parsed_mesh->positions,
                parsed_mesh->indices,
                parsed_mesh->normals,
                parsed_mesh->uvs,
                parsed_mesh->material_id,
                parsed_mesh->area_light_id};
            // Extract all the individual triangles
            for (int face_index = 0; face_index < (int)parsed_mesh->indices.size(); face_index++) {
                if (parsed_mesh->area_light_id >= 0) {
                    const ParsedDiffuseAreaLight &light =
                        std::get<ParsedDiffuseAreaLight>(
                            scene.lights[parsed_mesh->area_light_id]);
                    lights.push_back(DiffuseAreaLight{
                        (int)shapes.size(), light.radiance});
                }
                shapes.push_back(Triangle{face_index, &meshes[tri_mesh_count]});
            }
            tri_mesh_count++;
        } else {
            // Unhandled case
            assert(false);
        }
    }

    // Copy the materials
    for (const ParsedMaterial &parsed_mat : scene.materials) {
        if (auto *diffuse = std::get_if<ParsedDiffuse>(&parsed_mat)) {
            materials.push_back(Diffuse{
                parsed_color_to_texture(diffuse->reflectance)});
        } else if (auto *mirror = std::get_if<ParsedMirror>(&parsed_mat)) {
            materials.push_back(Mirror{
                parsed_color_to_texture(mirror->reflectance)});
        } else if (auto *plastic = std::get_if<ParsedPlastic>(&parsed_mat)) {
            materials.push_back(Plastic{
                plastic->eta,
                parsed_color_to_texture(plastic->reflectance)});
        } else if (auto *phong = std::get_if<ParsedPhong>(&parsed_mat)) {
            materials.push_back(Phong{
                parsed_color_to_texture(phong->reflectance),
                phong->exponent});
        } else if (auto *blinnphong = std::get_if<ParsedBlinnPhong>(&parsed_mat)) {
            materials.push_back(BlinnPhong{
                parsed_color_to_texture(blinnphong->reflectance),
                blinnphong->exponent});
        } else if (auto *blinn_mic = std::get_if<ParsedBlinnPhongMicrofacet>(&parsed_mat)) {
            materials.push_back(BlinnPhongMicrofacet{
                parsed_color_to_texture(blinn_mic->reflectance),
                blinn_mic->exponent});
        } else {
            // Unhandled case
            assert(false);
        }
    }

    // Copy the lights
    for (const ParsedLight &parsed_light : scene.lights) {
        if (auto *point_light = std::get_if<ParsedPointLight>(&parsed_light)) {
            lights.push_back(PointLight{point_light->intensity, point_light->position});
        } 
        /* 
        else if (auto *area_light = std::get_if<ParsedDiffuseAreaLight>(&parsed_light)) {
            lights.push_back(DiffuseAreaLight{area_light->shape_id, area_light->radiance});
        }*/
    }


    // Build BVH
    std::vector<BBoxWithID> bboxes(shapes.size());
    for (int i = 0; i < (int)bboxes.size(); i++) {
        if (auto *sph = std::get_if<Sphere>(&shapes[i])) {
            Vector3 p_min = sph-> position - sph->radius;
            Vector3 p_max = sph->position + sph->radius;
            bboxes[i] = {AABB{p_min, p_max}, i};
        } else if (auto *tri = std::get_if<Triangle>(&shapes[i])) {
            const TriangleMesh *mesh = tri->mesh;
            Vector3i index = mesh->indices[tri->face_index];

            //printf("index: %d, %d, %d\n", index[0], index[1], index[2]);
            Vector3 p0 = mesh->positions[index[0]];
            //printf("p0: %f, %f, %f\n", p0[0], p0[1], p0[2]);
            Vector3 p1 = mesh->positions[index[1]];
            //printf("p1: %f, %f, %f\n", p1[0], p1[1], p1[2]);
            Vector3 p2 = mesh->positions[index[2]];
            //printf("p2: %f, %f, %f\n", p2[0], p2[21], p2[2]);
            Vector3 p_min = min(min(p0, p1), p2);
            //printf("p_min: %f, %f, %f\n", p_min[0], p_min[1], p_min[2]);
            Vector3 p_max = max(max(p0, p1), p2);
            //printf(" p_max: %f, %f, %f\n", p_max[0], p_max[1], p_max[2]);
            bboxes[i] = {AABB{p_min, p_max}, i};
            

        }
    }
    bvh_root_id = construct_bvh(bboxes, bvh_nodes);
}


inline std::optional<Intersection> intersect(const Scene &scene, const BVHNode &node, ray ray) {
    if (node.primitive_id != -1) {
        return intersect(scene.shapes[node.primitive_id], ray);
    }
    const BVHNode &left = scene.bvh_nodes[node.left_node_id];
    const BVHNode &right = scene.bvh_nodes[node.right_node_id];
    std::optional<Intersection> isect_left;
    if (hit(left.box, ray)) {
        //printf("Step 2!\n");
        isect_left = intersect(scene, left, ray);
        if (isect_left) {
            ray.tfar = isect_left->distance;
        }
    }
    if (hit(right.box, ray)) {
        // Since we've already set ray.tfar to the left node
        // if we still hit something on the right, it's closer
        // and we should return that.
        if (auto isect_right = intersect(scene, right, ray)) {
            return isect_right;
        }
    }
    return isect_left;
}

inline bool occluded(const Scene &scene, const BVHNode &node, ray ray) {
    if (node.primitive_id != -1) {
        return occluded(scene.shapes[node.primitive_id], ray);
    }
    const BVHNode &left = scene.bvh_nodes[node.left_node_id];
    const BVHNode &right = scene.bvh_nodes[node.right_node_id];
    if (hit(left.box, ray)) {
        if (occluded(scene, left, ray)) {
            return true;
        }
    }
    if (hit(right.box, ray)) {
        if (occluded(scene, right, ray)) {
            return true;
        }
    }
    return false;
}

inline std::optional<Intersection> intersect(const Scene &scene, ray ray) {
    return intersect(scene, scene.bvh_nodes[scene.bvh_root_id], ray);
}

inline bool occluded(const Scene &scene, const ray &ray) {

    return occluded(scene, scene.bvh_nodes[scene.bvh_root_id], ray);
}

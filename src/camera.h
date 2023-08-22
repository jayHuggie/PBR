#pragma once

#include "parse_scene.h"

struct ray;

struct Camera {
    Vector3 lookfrom;
    Vector3 lookat;
    Vector3 up;
    Real vfov;
};

struct CameraRayData {
    Vector3 origin;
    Vector3 top_left_corner;
    Vector3 horizontal;
    Vector3 vertical;
};

inline Camera from_parsed_camera(const ParsedCamera &cam) {
    return Camera{cam.lookfrom, cam.lookat, cam.up, cam.vfov};
}

inline CameraRayData compute_camera_ray_data(const Camera &cam, int width, int height) {
    Real aspect_ratio = Real(width) / Real(height);
    Real viewport_height = 2.0 * tan(radians(cam.vfov / 2));
    Real viewport_width = aspect_ratio * viewport_height;

    Vector3 cam_dir = normalize(cam.lookat - cam.lookfrom);
    Vector3 right = normalize(cross(cam_dir, cam.up));
    Vector3 new_up = cross(right, cam_dir);

    Vector3 origin = cam.lookfrom;
    Vector3 horizontal = viewport_width * right;
    Vector3 vertical = viewport_height * new_up;
    Vector3 top_left_corner =
        origin - horizontal / Real(2) + vertical / Real(2) + cam_dir;
    return CameraRayData{origin, top_left_corner, horizontal, vertical};
}

inline ray generate_primary_ray(const CameraRayData &cam_ray_data, Real u, Real v) {
    Vector3 dir = normalize(cam_ray_data.top_left_corner +
                            u * cam_ray_data.horizontal -
                            v * cam_ray_data.vertical - cam_ray_data.origin);
    return ray{cam_ray_data.origin, dir, Real(0), infinity<Real>()};
}
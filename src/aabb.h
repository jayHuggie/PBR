#pragma once
#include "torrey.h"
#include "parse_scene.h"
#include "vector.h"
#include "ray.h"

struct AABB {
    Vector3 p_min = Vector3{infinity<Real>(),
                            infinity<Real>(),
                            infinity<Real>()};
    Vector3 p_max = Vector3{-infinity<Real>(),
                            -infinity<Real>(),
                            -infinity<Real>()};
};

inline bool hit(const AABB &bbox, ray ray) {
	// https://raytracing.github.io/books/RayTracingTheNextWeek.html#boundingvolumehierarchies/anoptimizedaabbhitmethod
    for (int i = 0; i < 3; i++) {
        Real inv_dir = Real(1) / ray.dir[i];
        Real t0 = (bbox.p_min[i] - ray.orig[i]) * inv_dir;
        Real t1 = (bbox.p_max[i] - ray.orig[i]) * inv_dir;
        if (inv_dir < 0) {
            std::swap(t0, t1);
        }
        ray.tnear = t0 > ray.tnear ? t0 : ray.tnear;
        ray.tfar = t1 < ray.tfar ? t1 : ray.tfar;
        if (ray.tfar < ray.tnear) {
            return false;
        }
    }
    return true;
}
inline int largest_axis(const AABB &box) {
    Vector3 extent = box.p_max - box.p_min;
    if (extent.x > extent.y && extent.x > extent.z) {
        return 0;
    } else if (extent.y > extent.x && extent.y > extent.z) {
        return 1;
    } else { // z is the largest
        return 2;
    }
}

inline AABB surrounding_box(const AABB& box1, const AABB& box2) {
    Vector3 p_min = Vector3{
        std::min(box1.p_min.x, box2.p_min.x),
        std::min(box1.p_min.y, box2.p_min.y),
        std::min(box1.p_min.z, box2.p_min.z)};
    Vector3 p_max = Vector3{
        std::max(box1.p_max.x, box2.p_max.x),
        std::max(box1.p_max.y, box2.p_max.y),
        std::max(box1.p_max.z, box2.p_max.z)};
    return AABB{p_min, p_max};
}

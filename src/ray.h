#pragma once

#include "vector.h"

struct ray {
    public:
        ray() {}
        ray(const Vector3& origin, const Vector3& direction)
            : orig(origin), dir(direction)
        {}

        ray(const Vector3& origin, const Vector3& direction, Real tnear, Real tfar)
            : orig(origin), dir(direction), tnear(tnear), tfar(tfar)
        {}

        Vector3 origin() const  { return orig; }
        Vector3 direction() const { return dir; }

        Vector3 at(double t) const {
            return orig + t*dir;
        }

    public:
        Vector3 orig;
        Vector3 dir;
        Real tnear, tfar;
};
#pragma once

#include "torrey.h"
#include "texture.h"
#include "vector.h"
#include <variant>
#include "shape.h"
#include "frame.h"
#include "pcg.h"

struct Scene;

struct Diffuse {
    Texture reflectance;
};
struct Mirror {
    Texture reflectance;
};

struct Plastic {
    Real ior;
    Texture reflectance;
};

struct Phong {
    Texture reflectance;  // Ks
    Real exponent;  // alpha
};

struct BlinnPhong {
    Texture reflectance; // Ks
    Real exponent; // alpha
};

struct BlinnPhongMicrofacet {
    Texture reflectance; // Ks
    Real exponent; // alpha
};

using Material = std::variant<Diffuse, Mirror, Plastic, 
                              Phong, BlinnPhong, BlinnPhongMicrofacet>;

enum class MaterialType {
    Diffuse,
    Mirror,
    Plastic,
    Phong,
    BlinnPhong,
    BlinnPhongMicrofacet
};

inline Vector3 schlick_fresnel(const Vector3 &F0, Real cos_theta) {
    return F0 + (1 - F0) * pow(1 - cos_theta, 5);
}

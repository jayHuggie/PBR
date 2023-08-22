#pragma once

#include "torrey.h"
#include "vector.h"
#include <optional>
#include <variant>
#include <vector>
#include "ray.h"

struct Sphere {
    Vector3 position;
    Real radius;
    int material_id = -1;
    int area_light_id = -1;

    void get_sphere_uv(const Vector3& point, double& u, double& v) const {
        Vector3 p = (point - position) / radius;  // Transform to unit sphere at origin
        auto theta = acos(p.y);
        auto phi = atan2(-p.z, p.x) + c_PI;

        u = phi / (2.0 * c_PI);   
        v = theta / c_PI;   

    }
};

struct TriangleMesh {
    std::vector<Vector3> positions;
    std::vector<Vector3i> indices;
    std::vector<Vector3> normals;
    std::vector<Vector2> uvs;
    int material_id = -1;
    int area_light_id = -1;
};

struct Triangle {
    int face_index;
    const TriangleMesh *mesh;
};

using Shape = std::variant<Sphere, Triangle, TriangleMesh>;

struct Intersection {
    Vector3 position;
    Vector3 geometric_normal;
    Vector3 shading_normal;
    Real distance;
    Vector2 uv;
    int material_id;
    int area_light_id;
};

// Möller-Trumbore intersection algorithm
inline bool intersect_triangle(const ray& r, const Vector3& v0, const Vector3& v1, const Vector3& v2, Real& t, Real& u, Real& v) {
    Vector3 edge1 = v1 - v0;
    Vector3 edge2 = v2 - v0;

    Vector3 h = cross(r.direction(), edge2);
    Real a = dot(edge1, h);

    const Real EPSILON = 1e-5;
    if (std::abs(a) < EPSILON)
        return false;

    Real f = 1.0 / a;
    Vector3 s = r.origin() - v0;
    u = f * dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    Vector3 q = cross(s, edge1);
    v = f * dot(r.direction(), q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    t = f * dot(edge2, q);

    return t > EPSILON;
}

inline std::optional<Vector3>
    intersect_triangle(const ray &ray,
                       const Vector3 &p0,
                       const Vector3 &p1,
                       const Vector3 &p2) {
    Vector3 e1 = p1 - p0;
    Vector3 e2 = p2 - p0;
    Vector3 s1 = cross(ray.dir, e2);
    Real divisor = dot(s1, e1);
    if (divisor == 0) {
        return {};
    }
    Real inv_divisor = 1 / divisor;
    Vector3 s = ray.orig - p0;
    Real u = dot(s, s1) * inv_divisor;
    Vector3 s2 = cross(s, e1);
    Real v = dot(ray.dir, s2) * inv_divisor;
    Real t = dot(e2, s2) * inv_divisor;

    if (t > ray.tnear && t < ray.tfar && u >= 0 && v >= 0 && u + v <= 1) {
        return Vector3{u, v, t};
    }

    return {};
}


// For sphere
inline bool solve_quadratic(Real a, Real b, Real c, Real *t0, Real *t1) {
    // Degenerated case
    if (a == 0) {
        if (b == 0) {
            return false;
        }
        *t0 = *t1 = -c / b;
        return true;
    }

    Real discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return false;
    }
    Real root_discriminant = sqrt(discriminant);
    if (b >= 0) {
        *t0 = (- b - root_discriminant) / (2 * a);
        *t1 = 2 * c / (- b - root_discriminant);
    } else {
        *t0 = 2 * c / (- b + root_discriminant);
        *t1 = (- b + root_discriminant) / (2 * a);
    }
    return true;
}


inline std::optional<Intersection> intersect(const Sphere &sph, const ray &ray) {

    Vector3 oc = ray.orig - sph.position;
    Real A = dot(ray.dir, ray.dir);
    Real B = 2 * dot(ray.dir, oc);
    Real C = dot(oc, oc) - sph.radius * sph.radius;
    Real t0, t1;
    if (!solve_quadratic(A, B, C, &t0, &t1)) {
        // No intersection
        return {};
    }
    if (t0 > t1) {
        std::swap(t0, t1);
    }
    Real t = t0;
    if (t0 >= ray.tnear && t0 < ray.tfar) {
        t = t0;
    }
    if (t1 >= ray.tnear && t1 < ray.tfar && t < ray.tnear) {
        t = t1;
    }

    if (t >= ray.tnear && t < ray.tfar) {
        // Record the intersection
        Vector3 p = ray.orig + t * ray.dir;
        Vector3 n = normalize(p - sph.position);
        Real theta = acos(n.y);
        Real phi = atan2(-n.z, n.x) + c_PI;
        Vector2 uv{phi / (2 * c_PI), theta / c_PI};
        return Intersection{p, // position
                            n, // geometric normal
                            n, // shading normal
                            t, // distance
                            uv, // UV
                            sph.material_id,
                            sph.area_light_id};
    }
    return {};
}


inline std::optional<Intersection> intersect(const Triangle &tri, const ray &ray) {
	const TriangleMesh &mesh = *tri.mesh;
	Vector3i index = mesh.indices[tri.face_index];
	Vector3 p0 = mesh.positions[index[0]];
	Vector3 p1 = mesh.positions[index[1]];
	Vector3 p2 = mesh.positions[index[2]];
	if (auto uvt_ = intersect_triangle(ray, p0, p1, p2)) {
		Vector2 b = Vector2{uvt_->x, uvt_->y};
		Real t = uvt_->z;
		Vector3 p = (1 - b[0] - b[1]) * p0 + b[0] * p1 + b[1] * p2;
		Vector3 geometric_normal = normalize(cross(p1 - p0, p2 - p0));
        Vector2 uv = b;
        if (mesh.uvs.size() > 0) {
            Vector2 uv0 = mesh.uvs[index[0]];
            Vector2 uv1 = mesh.uvs[index[1]];
            Vector2 uv2 = mesh.uvs[index[2]];
            uv = (1 - b[0] - b[1]) * uv0 + b[0] * uv1 + b[1] * uv2;
        }

        Vector3 shading_normal = geometric_normal;
        if (mesh.normals.size() > 0) {
            Vector3 n0 = mesh.normals[index[0]];
            Vector3 n1 = mesh.normals[index[1]];
            Vector3 n2 = mesh.normals[index[2]];
            shading_normal = normalize((1 - b[0] - b[1]) * n0 + b[0] * n1 + b[1] * n2);
        }
		return Intersection{p, // position
                            geometric_normal,
                            shading_normal,
                            t, // distance
                            uv,
                            mesh.material_id,
                            mesh.area_light_id};
	}
	return {};
}

inline std::optional<Intersection> intersect(const Shape &shape, const ray &ray) {
	if (auto *sph = std::get_if<Sphere>(&shape)) {
		return intersect(*sph, ray);
	} else if(auto *tri = std::get_if<Triangle>(&shape)) {
		return intersect(*tri, ray);
	} else {
		assert(false);
		return {};
	}
}


inline bool occluded(const Shape &shape, const ray &ray) {
    return bool(intersect(shape, ray));
}

/////// HW 3-4 

struct PointAndNormal {
	Vector3 position;
	Vector3 normal;
};

inline PointAndNormal sample_on_shape(const Shape &shape,
                               const Vector2 &u) {
    if (auto *sph = std::get_if<Sphere>(&shape)) {
        Real theta = acos(std::clamp(1 - 2 * u[0], Real(-1), Real(1)));
        Real phi = 2 * c_PI * u[1];
        Vector3 n{
            cos(phi) * sin(theta),
            sin(phi) * sin(theta),
            cos(theta)
        };
        Vector3 p = sph->radius * n + sph-> position;
        return {p, n};
    } else if (auto *tri = std::get_if<Triangle>(&shape)) {
        Real b1 = 1 - sqrt(max(u[0], Real(0)));
        Real b2 = u[1] * sqrt(max(u[0], Real(0)));
        const TriangleMesh *mesh = tri->mesh;
        Vector3i index = mesh->indices[tri->face_index];
        Vector3 p0 = mesh->positions[index[0]];
        Vector3 p1 = mesh->positions[index[1]];
        Vector3 p2 = mesh->positions[index[2]];
        Vector3 p = (1 - b1 - b2) * p0 + b1 * p1 + b2 * p2;
        Vector3 geometric_normal = normalize(cross(p1 - p0, p2 - p0));
        Vector3 shading_normal = geometric_normal;
        if (mesh->normals.size() > 0) {
            Vector3 n0 = mesh->normals[index[0]];
            Vector3 n1 = mesh->normals[index[1]];
            Vector3 n2 = mesh->normals[index[2]];
            shading_normal = normalize((1 - b1 - b2) * n0 + b1 * n1 + b2 * n2);
        }
        if (dot(geometric_normal, shading_normal) < 0) {
            geometric_normal = -geometric_normal;
        }
        return {p, geometric_normal};
    } else {
        assert(false);
        return {};
    }
}

inline Real pdf_sample_on_shape(const Shape &shape,
                         const Vector3 &p) {
    if (auto *sph = std::get_if<Sphere>(&shape)) {
        return Real(1) / (4 * c_PI * sph->radius * sph->radius);
    } else if (auto *tri = std::get_if<Triangle>(&shape)) {
        const TriangleMesh *mesh = tri->mesh;
        Vector3i index = mesh->indices[tri->face_index];
        Vector3 p0 = mesh->positions[index[0]];
        Vector3 p1 = mesh->positions[index[1]];
        Vector3 p2 = mesh->positions[index[2]];
        return Real(1) / (Real(0.5) * length(cross(p1 - p0, p2 - p0)));
    } else {
        assert(false);
        return Real(0);
    }
}

inline Vector3 sample_dir_to_shape(const Shape &shape,
                            const Vector3 &p,
                            const Vector2 &u) {
    Vector3 lp = sample_on_shape(shape, u).position;
    return normalize(lp - p);
}

inline Real pdf_sample_dir_to_shape(const Shape &shape,
                             const Vector3 &p,
                             const Vector3 &dir) {
    Real pdf = 0;
    ray input_ray{p, dir, Real(1e-4), infinity<Real>()};
    while (auto isect = intersect(shape, input_ray)) {
        // Convert to solid angle measure
        pdf += pdf_sample_on_shape(shape, isect->position) *
               distance_squared(p, isect->position) /
               fabs(dot(dir, isect->geometric_normal));
        input_ray = ray{isect->position, dir, Real(1e-4), infinity<Real>()};
    }
    return pdf;
}

inline Vector3 sample_cos_hemisphere(const Vector2 &u) {
    Real phi = c_TWOPI * u[0];
    Real tmp = sqrt(std::clamp(1 - u[1], Real(0), Real(1)));
    return Vector3{
        cos(phi) * tmp, sin(phi) * tmp,
        sqrt(std::clamp(u[1], Real(0), Real(1)))
    };
}

// For Phong/Blinn
inline Vector3 sample_cos_n_hemisphere(const Vector2 &u, Real exponent) {
    Real phi = c_TWOPI * u[0];
    Real cos_theta = pow(u[1], Real(1)/(exponent + 1));
    Real sin_theta = sqrt(std::clamp(1 - cos_theta * cos_theta, Real(0), Real(1)));
    return Vector3{
        cos(phi) * sin_theta,
        sin(phi) * sin_theta,
        cos_theta
    };
}

inline Real smithG1(const Vector3 &v,
                    const Vector3 &m,
                    const Vector3 &n,
                    Real e) {
    if (dot(v, m) <= 0) {
        return 0;
    }

    // tan^2 + 1 = 1/cos^2
    Real cos_theta = dot(v, n);
    Real tan_theta = sqrt(1/square(cos_theta) - 1);
    Real a = sqrt(Real(0.5) * e + 1) / tan_theta;
    if (a < Real(1.6)) {
        return (Real(3.535) * a + Real(2.181) * square(a))
             / (Real(1.0) + Real(2.276) * a + Real(2.577) * square(a));
    } else {
        return 1;
    }
}

#ifndef VECTOR_H
#define VECTOR_H

// Raytracer headers.
#include "utils.h"

// Standard headers.
#include <math.h>
#include <iostream>


class vec3
{
  public:

    // Constructors.
    vec3() : v{0.0, 0.0, 0.0} {}
    vec3(double x) : v{x, x, x} {} 
    vec3(double x, double y, double z) : v{x, y, z} {}

    // Getters.
    double x() const { return v[0]; }
    double y() const { return v[1]; }
    double z() const { return v[2]; }

    // Functions.

    vec3 operator-() const { return vec3(-v[0], -v[1], -v[2]); }

    vec3& operator+=(const vec3 &c)
    {
        v[0] += c.x();
        v[1] += c.y();
        v[2] += c.z();
        return *this;
    }

    vec3& operator*=(const double c)
    {
        v[0] *= c;
        v[1] *= c;
        v[2] *= c;
        return *this;
    }

    void print_vec()
    {
        std::cerr<<"x "<<v[0]<<" "<<"y "<<v[1]<<" "<<"z "<<v[2]<<" "<<std::endl;
    }

    void print_vec() const
    {
        std::cerr<<"x "<<v[0]<<" "<<"y "<<v[1]<<" "<<"z "<<v[2]<<" "<<std::endl;
    }

    inline double length() const 
    {
        return std::sqrt(length_squared());
    }

    inline double length_squared() const
    {
        return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    }

    inline static vec3 random_vector()
    {
        return vec3(
            random_double(), 
            random_double(), 
            random_double());
    }

    inline static vec3 random_vector(double min, double max)
    {
        return vec3(
            random_double(min, max), 
            random_double(min, max), 
            random_double(min, max));
    }

    // Attributes.
    double v[3];
    
};

// Functions.

inline vec3 operator*(double t, const vec3 &v)
{
    return vec3(t * v.x(), t * v.y(), t * v.z());
}

inline vec3 operator*(const vec3 &v, double t)
{
    return vec3(t * v.x(), t * v.y(), t * v.z());
}

inline vec3 operator*(const vec3 &u, const vec3 &v)
{
    return vec3(u.x() * v.x(), u.y() * v.y(), u.z() * v.z());
}

inline vec3 operator+(const vec3 &u, const vec3  &v)
{
    return vec3(v.x() + u.x(), v.y() + u.y(), v.z() + u.z());
}

inline vec3 operator-(const vec3 &u, const vec3  &v)
{
    return vec3(u.x() - v.x(), u.y() - v.y(), u.z() - v.z());
}

inline vec3 operator/(vec3 v, double t)
{
    return (1/t) * v;
}

inline vec3 normalize(vec3 v)
{
    return v / v.length();
}

inline double dot(const vec3 &u, const vec3 &v)
{
    return u.x() * v.x() + u.y() * v.y() + u.z() * v.z(); 
}

inline vec3 cross(const vec3 &u, const vec3& v)
{
    return vec3(u.y() * v.z() - u.z() * v.y(),
                u.z() * v.x() - u.x() * v.z(),
                u.x() * v.y() - u.y() * v.x());
}

inline vec3 floor(vec3 p)
{
    return vec3(
        floor(p.x()),
        floor(p.y()),
        floor(p.z()));
}

// Create random point in unit sphere
// using rejection sampling.
// Poor sampling for indirect diffuse lambertian.
vec3 random_point_in_unit_sphere()
{
    while(true)
    {
        auto p = vec3::random_vector(-1, 1);
        if(p.length_squared() >= 1)
            continue; // reject
        return p;
    }
}

// Better (more uniform) sampling for indirect diffuse lambertian.
vec3 random_unit_vector_on_sphere()
{
    return normalize(random_point_in_unit_sphere());
}

// Uniform hemisphere sampling for indirect diffuse.
vec3 uniform_hemisphere_sampler(const vec3& normal)
{
    vec3 vector_in_unit_sphere = random_point_in_unit_sphere();
    if(dot(vector_in_unit_sphere, normal) > 0.0)
        return vector_in_unit_sphere; // in the same hemisphere as normal
    else
        return -vector_in_unit_sphere; // flip in normal hemisphere
    
}

// Mirror reflection.
vec3 reflect(const vec3& v, const vec3& n)
{
    return v - 2 * dot(v, n) * n;
}

// Perfect refraction (Look in book for derivation).
vec3 refract(const vec3& v, const vec3& n, double etai_over_etat)
{
    auto cos_theta = fmin(dot(-v, n), 1.0);
    vec3 refracted_perpendicular = etai_over_etat * (v + cos_theta * n);
    vec3 refracted_parallel = -sqrt(1 - fabs(refracted_perpendicular.length_squared())) * n;
    return refracted_perpendicular + refracted_parallel;
}


// Aliases.
using point3 = vec3; // 3D point
using color = vec3; // RGB color

#endif

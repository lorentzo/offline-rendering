
#ifndef RAY_H
#define RAY_H

// Raytracer include header files.
#include "vector.h"

// STD headers
#include<iostream>

class ray
{
  public:
    
    // Constructors.
    ray() {}
    ray(const point3& origin, const vec3& direction)
        : orig(origin), dir(normalize(direction))
    {}

    // Getters.
    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }

    // Functions.
    point3 at(double t) const
    {
        return orig + t * dir;
    }

    void print_ray() const
    {
        std::cerr<<"o.x "<<orig.x()<<" o.y "<<orig.y()<<" o.z "<<orig.z()<<std::endl;
        std::cerr<<"d.x "<<dir.x()<<" d.y "<<dir.y()<<" d.z "<<dir.z()<<std::endl;
    }

    // Attributes.
    point3 orig;
    vec3 dir;
};

#endif

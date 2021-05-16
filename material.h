
#ifndef MATERIAL_H
#define MATERIAL_H

#include "ray.h"
#include "vector.h"

// Materials for objects in the scene.
enum MaterialType {diffuse, specular, glossy, volume, nonReflective_diffuse_emission};

/*
class Material
{
  public:
    virtual bool scatter(
        const ray& outgoing,
        const IntersectionContext& intersection_context,
        color& attenuation,
        ray& incoming) const = 0;
};
*/


#endif

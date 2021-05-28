
#ifndef MATERIAL_H
#define MATERIAL_H

#include "ray.h"
#include "vector.h"
#include "intersection.h"

// Materials for objects in the scene.
enum MaterialType 
{
  diffuse, 
  specular, 
  glossy, 
  volume, 
  nonReflective_diffuse_emission,
  constant_medium};

/*
// Abstraction for all materials.
class Material
{
  public:
    virtual bool scatter(
        const ray& outgoing,
        const IntersectionContext& intersection_context,
        color& attenuation,
        ray& incoming) const = 0;
};

class ConstantMediumMaterial : Material
{
  public:
    ConstantMediumMaterial(color _medium_color) : medium_color(_medium_color) {}

    virtual bool scatter(
        const ray& outgoing, 
        const IntersectionContext& intersection_context,
        color& attenuation,
        ray& scattered) const override 
        {
            scattered = ray(intersection_context.p, random_unit_vector_on_sphere());
            attenuation = medium_color;
            return true;
        }

  private:
    color medium_color; // attenuation
};
*/


#endif

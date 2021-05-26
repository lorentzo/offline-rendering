
#ifndef INTERSECTION_H
#define INTERSECTION_H

// Raytracer header files.
#include "vector.h"
#include "ray.h"
#include "material.h"

// Standard header files.
#include <vector>
#include <memory> // shared_ptr
#include <math.h>

struct IntersectionContext
{
    point3 p;
    vec3 normal;
    double t;
    double t2;
    ray r;
    bool front_face;

    MaterialType material_type;

    // Design decision: we want normal always to be pointing
    // against the ray.
    // Outward normal is one pointing out of object's
    // face e.g. sphere: center - intersection_point.
    inline void set_face_normal(const ray& r, const vec3& outward_normal)
    {
        // Angle between outward normal and ray determines
        // from which side the ray is.
        front_face = dot(r.direction(), outward_normal) < 0.0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

// Abstraction for all objects.
class Intersectable
{
  public:
    virtual bool intersect(
        const ray& r, 
        double t_min, 
        double t_max,
        IntersectionContext& intersection_context) const = 0;
};

// List of intersectable objects.
class IntersectableList : public Intersectable
{
  public:
    // Constructors.
    IntersectableList() {}
    IntersectableList(std::shared_ptr<Intersectable> object) { add(object); }

    // Functions.
    void clear() { objects.clear(); }
    void add(std::shared_ptr<Intersectable> object) { objects.push_back(object); }
    virtual bool intersect(
        const ray& r, 
        double t_min, 
        double t_max,
        IntersectionContext& intersection_context) const override;

    // Attributes.
    std::vector<std::shared_ptr<Intersectable>> objects;
};

bool IntersectableList::intersect(
    const ray& r, 
    double t_min, 
    double t_max,
    IntersectionContext& intersection_context) const
{
    bool hit_anything = false;
    auto closest_hit = t_max;
    IntersectionContext closest_intersection_context;

    for(const auto& object : objects)
    {
        if(object->intersect(r, t_min, t_max, closest_intersection_context)){
            hit_anything = true;
            closest_hit = closest_intersection_context.t;
            intersection_context = closest_intersection_context;
            break;
        }
    }
    return hit_anything;
}

// Sphere object.
class Sphere : public Intersectable
{
  public:
    // Constructors.
    Sphere() {}
    Sphere(
        point3 _center, 
        double _radius,
        MaterialType _material_type) 
      : center(_center)
      , radius(_radius) 
      , material_type(_material_type) {};

    virtual bool intersect(
        const ray& r, 
        double t_min, 
        double t_max,
        IntersectionContext& intersection_context) const override;

    void print_info()
    {
        std::cerr<<"c.x "<<center.x()<<" c.y "<<center.y()<<" c.z "<<center.z()<<std::endl;
        std::cerr<<"r "<<radius<<std::endl;
    }

    // Attributes.
    point3 center;
    double radius;
    MaterialType material_type;
};

// Sphere intersection.
bool Sphere::intersect(
    const ray& r, 
    double t_min, 
    double t_max,
    IntersectionContext& intersection_context) const
{
    // Quadratic equation of t describing ray-sphere intersection.
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction()); // or r.direction().length_squared()
    auto b = 2.0 * dot(r.direction(), oc);
    auto c = dot(oc, oc) - radius * radius;
    auto discriminant = b * b - 4 * a * c; // tells if there are real solutions
    if (discriminant < 0)
        return false;

    // Find nearest solution (t, nearest intersection) in range [t_min, t_max].
    // TODO other t can be used for volumentrics.
    
    auto sqrt_discriminant = sqrt(discriminant);
    /*
    auto root = (-b - sqrt_discriminant) / (2 * a); // set to root1
    if(root < t_min || root > t_max)
    {
        root = (-b + sqrt_discriminant) / (2 * a); // set to root2
        if(root < t_min || root > t_max)
            return false;
    }
    */
    auto t1 = (-b - sqrt_discriminant) / (2 * a);
    auto t2 = (-b + sqrt_discriminant) / (2 * a);
    if(t2 < t_min || t1 > t_max) // not visible
        return false;

    intersection_context.t = t1;
    intersection_context.t2 = t2;
    intersection_context.p = r.at(intersection_context.t);
    vec3 outward_normal = normalize(intersection_context.p - center);
    intersection_context.set_face_normal(r, outward_normal);
    intersection_context.r = r; // can be set independent of intersection solution.
    intersection_context.material_type = material_type;
    return true;
}

// AABB.
class AABB 
{
  public:
    AABB() {}
    AABB(const point3& min, const point3& max) {minimum = min; maximum = max;}

    point3 min() const { return minimum; }
    point3 max() const { return maximum; }

    bool hit(const ray& r, double t_min, double t_max) const {
        for (int a = 0; a < 3; a++) {
            auto t0 = fmin((minimum[a] - r.origin()[a]) / r.direction()[a],
                            (maximum[a] - r.origin()[a]) / r.direction()[a]);
            auto t1 = fmax((minimum[a] - r.origin()[a]) / r.direction()[a],
                            (maximum[a] - r.origin()[a]) / r.direction()[a]);
            t_min = fmax(t0, t_min);
            t_max = fmin(t1, t_max);
            if (t_max <= t_min)
                return false;
        }
        return true;
    }

  private:
    point3 minimum;
    point3 maximum;
}

#endif

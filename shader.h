
#ifndef SHADER_H
#define SHADER_H

// Raytracer include.
#include "vector.h"
#include "ray.h"
#include "environment.h"
#include "intersection.h"
#include "utils.h"
#include "noise.h"

// Standard headers
#include <iostream>

//
// NPR: normals.
//

color shade_normal(IntersectionContext intersection_context)
{
    vec3 normal = intersection_context.normal;
    color normal_rgb = 0.5 * color(normal.x() + 1, normal.y() + 1, normal.z() + 1);
    return normal_rgb;
}

//
// Environment shading.
//
color environment_shading_engine(ray r, EnvironmentType environment_type, NishitaSky nishita_sky)
{
    if (environment_type == EnvironmentType::black)
        return black_env();
    if (environment_type == EnvironmentType::white)
        return white_env();
    if (environment_type == EnvironmentType::gradient)
        return gradient_env(r);
    if (environment_type == EnvironmentType::sky)
        return sky_env(r);
    if (environment_type == EnvironmentType::nishita)
        return nishita_sky.compute_color(r, 0, infinity);
    else
        return color(1.0, 0.0, 1.0); // color code for non-existing environment
}

// 
// Recursive generic shading for objects (surfaces and volumes).
// Calculates scattered ray and attenuation for intersected object.
//
color recursive_shading_engine(
    const ray& r, 
    const IntersectableList& scene,
    IntersectionContext& intersection_context,
    int depth,
    EnvironmentType environment_type,
    NishitaSky nishita_sky)
{
    // Stop recursive calls after certain depth.
    if(depth <= 0)
        return color(0.0, 0.0, 0.0);

    if(scene.intersect(r, 0.001, infinity, intersection_context)) //[0.001, inf] due to shadow acne problem
    {
        if(intersection_context.material_type == MaterialType::diffuse)
        {            
            // Find random point in unit sphere on the side of surface as the intersection (next ray origin).
            point3 scattering_target = 
                intersection_context.p + intersection_context.normal + random_unit_vector_on_sphere(); // or random_unit_vector_on_sphere()

            // Calculate scattering direction (new ray) based on that point.
            ray scattering_ray = ray(intersection_context.p, scattering_target - intersection_context.p);

            // Attenuate (absorb) when hitting and continue tracing.
            return 0.5 * recursive_shading_engine(scattering_ray, scene, intersection_context, depth-1, environment_type, nishita_sky);
        }
        else if(intersection_context.material_type == MaterialType::specular)
        {
            // Mirror reflection of ray on intersection point.
            point3 scattering_direction = 
                reflect(normalize(r.direction()), intersection_context.normal);

            // Calculate scattering direction from intersection point to reflected direction.
            ray scattering_ray = ray(intersection_context.p, scattering_direction);

            // Attenuate (absorb) when hitting and continue tracing.
            return 0.5 * recursive_shading_engine(scattering_ray, scene, intersection_context, depth-1, environment_type, nishita_sky);
        }
        else if(intersection_context.material_type == MaterialType::nonReflective_diffuse_emission)
        {
            // No reflection (recursive call).
            // Diffusive: same for all points on surface and incoming/viewing directions.
            // White color.
            float intensity = 5.0;
            return intensity * color(1.0, 1.0, 1.0);
        }
        else if(intersection_context.material_type == MaterialType::volume)
        {
            // Noise.
            //std::cerr<<intersection_context.t2<<" "<<intersection_context.t<<std::endl;
            float dt = 0.2;
            float n_segments = (intersection_context.t2 - intersection_context.t) / dt;
            color c = color(0.0);
            for(float t = intersection_context.t; t < intersection_context.t2; t += dt)
            {
                point3 p = intersection_context.p + intersection_context.r.direction() * t;
                c += worley3(p, 1.0);
            }
            //c.print_vec();
            return c / n_segments;
        }
        else
        {
            return color(1.0, 0.0, 1.0); // color code for non-existing material
        }
        
    }
    else
    {
        return environment_shading_engine(r, environment_type, nishita_sky);
    }
}

/*
//
// Diffuse 1.
//

color shade_indirect_diffuse1(
    const ray& r, 
    const IntersectableList& scene,
    IntersectionContext& intersection_context,
    int depth,
    EnvironmentType environment_type)
{
    // Stop recursive calls after certain depth.
    if(depth <= 0)
        return color(0.0, 0.0, 0.0);

    if(scene.intersect(r, 0.001, infinity, intersection_context)) //[0.001, inf] due to shadow acne problem
    {
        // Find random point in unit sphere on the side of surface as the intersection (next ray origin).
        point3 scattering_target = 
            intersection_context.p + intersection_context.normal + random_unit_vector_on_sphere(); // or random_unit_vector_on_sphere()

        // Calculate scattering direction (new ray) based on that point.
        ray scattering_ray = ray(intersection_context.p, scattering_target - intersection_context.p);

        // Attenuate (absorb) when hitting and continue tracing.
        return 0.5 * shade_indirect_diffuse1(scattering_ray, scene, intersection_context, depth-1, environment_type);
    }
    else
    {
        return environment_shading_engine(r, environment_type);
    }
}

//
// Diffuse2.
//

color shade_indirect_diffuse2(
    const ray& r, 
    const IntersectableList& scene,
    IntersectionContext& intersection_context,
    int depth,
    EnvironmentType environment_type)
{
    // Stop recursive calls after certain depth.
    if(depth <= 0)
        return color(0.0, 0.0, 0.0);

    if(scene.intersect(r, 0.001, infinity, intersection_context)) //[0.001, inf] due to shadow acne problem
    {
        // Find random point in unit sphere on the side of surface as the intersection (next ray origin).
        point3 scattering_target = 
            intersection_context.p + uniform_hemisphere_sampler(intersection_context.normal); // or random_unit_vector_on_sphere()

        // Calculate scattering direction (new ray) based on that point.
        ray scattering_ray = ray(intersection_context.p, scattering_target - intersection_context.p);

        // Attenuate (absorb) when hitting and continue tracing.
        return 0.5 * shade_indirect_diffuse2(scattering_ray, scene, intersection_context, depth-1, environment_type);
    }
    else
    {
        return environment_shading_engine(r, environment_type);
    }
}

//
// Specular.
//

color shade_indirect_specular(
    const ray& r, 
    const IntersectableList& scene,
    IntersectionContext& intersection_context,
    int depth,
    EnvironmentType environment_type)
{
    // Stop recursive calls after certain depth.
    if(depth <= 0)
        return color(0.0, 0.0, 0.0);

    if(scene.intersect(r, 0.001, infinity, intersection_context)) //[0.001, inf] due to shadow acne problem
    {
        // Mirror reflection of ray on intersection point.
        point3 scattering_direction = 
            reflect(normalize(r.direction()), intersection_context.normal);

        // Calculate scattering direction from intersection point to reflected direction.
        ray scattering_ray = ray(intersection_context.p, scattering_direction);

        // Attenuate (absorb) when hitting and continue tracing.
        return 0.5 * shade_indirect_specular(scattering_ray, scene, intersection_context, depth-1, environment_type);
    }
    else
    {
        return environment_shading_engine(r, environment_type);
    }
}

//
// Specular glossy.
//

color shade_indirect_glossy(
    const ray& r, 
    const IntersectableList& scene,
    IntersectionContext& intersection_context,
    int depth,
    double fuzz,
    EnvironmentType environment_type)
{
    // Stop recursive calls after certain depth.
    if(depth <= 0)
        return color(0.0, 0.0, 0.0);

    if(scene.intersect(r, 0.001, infinity, intersection_context)) //[0.001, inf] due to shadow acne problem
    {
        // Mirror reflection of ray on intersection point.
        point3 scattering_direction = 
            reflect(normalize(r.direction()), intersection_context.normal);

        // Calculate scattering direction from intersection point to reflected direction.
        // Reflected direction is bit fuzzed to appear glossy.
        ray scattering_ray = ray(intersection_context.p, scattering_direction + fuzz * random_point_in_unit_sphere());

        // Attenuate (absorb) when hitting and continue tracing.
        return 0.5 * shade_indirect_glossy(scattering_ray, scene, intersection_context, depth-1, fuzz, environment_type);
    }
    else
    {
        return environment_shading_engine(r, environment_type);
    }
}

//
// Refraction.
//

// Helper function for refraction: 
// Schlick approx: glass reflectivity changes with angle.
double reflectance(double cosine, double refraction_ior_ratio)
{
    auto r0 = (1- refraction_ior_ratio) / (1 + refraction_ior_ratio);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}

color shade_indirect_refraction(
    const ray& r, 
    const IntersectableList& scene,
    IntersectionContext& intersection_context,
    int depth,
    double IOR,
    EnvironmentType environment_type)
{
    // Stop recursive calls after certain depth.
    if(depth <= 0)
        return color(0.0, 0.0, 0.0);

    if(scene.intersect(r, 0.001, infinity, intersection_context)) //[0.001, inf] due to shadow acne problem
    {
        // Set refraction IOR ratio.
        double refraction_ratio;
        if(intersection_context.front_face)
            refraction_ratio = (1.0 / IOR);
        else
            refraction_ratio = IOR;

        // Calculate scattering direction. Check for total internal reflection.
        vec3 incoming_dir = normalize(r.direction());
        vec3 scattering_direction;
        double cos_theta = fmin(1.0, dot(-incoming_dir, intersection_context.normal));
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        bool cannot_refract = refraction_ratio * sin_theta > 1.0; // no solution to snell
        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            scattering_direction = reflect(incoming_dir, intersection_context.normal);
        else
            scattering_direction = refract(incoming_dir, intersection_context.normal, refraction_ratio);

        // Calculate scattering ray from intersection point to scattering direction.
        ray scattering_ray = ray(intersection_context.p, scattering_direction);

        // (No) attenuation (absorb) when hitting and continue tracing.
        return 1.0 * shade_indirect_refraction(scattering_ray, scene, intersection_context, depth-1, IOR, environment_type);
    }
    else
    {
        return environment_shading_engine(r, environment_type);
    }
}
*/


#endif

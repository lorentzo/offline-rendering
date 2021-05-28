
#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

// Raytracer headers
#include "intersection.h"
#include "material.h"
#include "vector.h"
#include "ray.h"
#include "utils.h"

// STD headers
#include <iostream>

enum EnvironmentType 
{
    black, 
    white, 
    sky, 
    gradient, 
    nishita
};

// Non emissive. Constant black.
color black_env()
{
    return color(0.0, 0.0, 0.0);
}

// Max emissive. Constant white.
color white_env()
{
    return color(0.3, 0.3, 0.3);
}

// Linearly blend white and blue depending on film height y after scaling [-1, 1].
// t = 0: blue, t = 0: white.
// lerp: blended_val = (1 - t) * white + t * blue
color sky_env(const ray& r)
{
    vec3 unit_direction = normalize(r.direction()); // get film point height [-1, 1]
    auto t = 0.5 * (unit_direction.y() + 1.0); // scale to [0, 1]
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0); // use for blending
}

// Idea:
// vec3( double(i) / (image_width-1), double(j) / (image_height-1), c).
color gradient_env(const ray& r)
{
    vec3 unit_direction = normalize(r.direction());
    auto x = 0.5 * (unit_direction.x() + 1.0);
    auto y = 0.5 * (unit_direction.y() + 1.0);
    return color(x, y, 0.55);
}

// Based on: https://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds/simulating-sky
class NishitaSky
{
  public:
    NishitaSky(
        float _earth_radius,
        float _atmosphere_radius,
        float _view_samples,
        float _light_samples,
        vec3 _sun_direction,
        color _sun_intensity,
        float _Hr,
        float _Hm)
    {
        earth_radius = _earth_radius;
        atmosphere_radius = _atmosphere_radius;
        atmosphere = Sphere(point3(0.0, 0.0, 0.0), atmosphere_radius, MaterialType::volume);
        view_samples = _view_samples;
        light_samples = _light_samples;
        sun_direction = _sun_direction;
        sun_intensity = _sun_intensity;
        Hr = _Hr;
        Hm = _Hm;
    }

    color compute_color( // currently only for Earth-Sky perspective
        const ray& view_ray,
        float tmin,
        float tmax);

    float earth_radius;
    float atmosphere_radius;
    Sphere atmosphere;
    float view_samples;
    float light_samples;
    vec3 sun_direction;
    color sun_intensity;
    float Hr;
    float Hm;
    const vec3 beta_rayleigh = vec3(3.8e-6f, 13.5e-6f, 33.1e-6f);
    const vec3 beta_mie = vec3(21e-6f, 21e-6f, 21e-6f);
    const float anisotropic_mie = 0.76;
};

color NishitaSky::compute_color(
        const ray& view_ray,
        float tmin,
        float tmax) 
{
    // Compute view ray and atmosphere intersection.
    float t0view, t1view, tview_curr;
    IntersectionContext view_intersection_context;
    if(!atmosphere.intersect(view_ray, tmin, tmax, view_intersection_context))
        return color(0.0, 0.0, 0.0);
    //view_ray.print_ray();
    //atmosphere.print_info();
    
    t0view = tmin;
    t1view = view_intersection_context.t;
    tview_curr = t0view;
    //std::cerr<<"t0 "<<t0view<<" t1 "<<t1view<<std::endl;

    // Calculate viewing integration segment length.
    float view_segment = (t1view - t0view) / view_samples;
    //std::cerr<<"view_segment "<<view_segment<<std::endl;

    // Phase functions.
    float mu = dot(view_ray.direction(), sun_direction);
    float phaseR = (3.0 / (16.0 * pi)) * (1 + mu * mu);
    float phaseM = 
        (3.0 / (8.0 * pi)) * 
        (((1.0 - anisotropic_mie * anisotropic_mie) * (1.0 + mu * mu)) /
        ((2.0 + anisotropic_mie * anisotropic_mie) * pow((1.0 + anisotropic_mie * anisotropic_mie - 2 * anisotropic_mie * mu), 1.5)));
    //std::cerr << phaseR << " " << phaseM << std::endl;
    // Final color.
    color final_color_rayleigh;
    color final_color_mie;

    // Prepare transmittance, phase function and scattering coeff variables.
    vec3 transmittance_view_rayleigh;
    vec3 transmittance_view_mie;
    vec3 scattering_view_rayleigh;
    vec3 scattering_view_mie;

    // For each viewing sample.
    for (int i = 0; i < view_samples; ++i)
    {
        // Calculate viewing sample position and height.
        // NOTE: Prenosim origin (0.0, 0.0, 0.0), ali promatrac je Er + 2!
        vec3 view_sample_position = 
            (view_ray.origin()) + (tview_curr + view_segment / 2.0) * view_ray.direction();
        float view_sample_height = view_sample_position.length() - earth_radius;
        if (view_sample_height < 0.0)
            std::cerr<<"view sample height negative!" << std::endl;
        //view_ray.dir.print_vec();
        //view_sample_position.print_vec();
        //std::cerr<<view_sample_height<<std::endl;

        // Calculate transmittance, phase function and scattering coeff contrib for view.
        transmittance_view_rayleigh += vec3(
            exp(-beta_rayleigh.x() * exp(-view_sample_height / Hr) * view_segment),
            exp(-beta_rayleigh.y() * exp(-view_sample_height / Hr) * view_segment),
            exp(-beta_rayleigh.z() * exp(-view_sample_height / Hr) * view_segment));
        //transmittance_view_rayleigh.print_vec();
        transmittance_view_mie += vec3(
            exp(-beta_mie.x() * exp(-view_sample_height / Hm) * view_segment),
            exp(-beta_mie.y() * exp(-view_sample_height / Hm) * view_segment),
            exp(-beta_mie.z() * exp(-view_sample_height / Hm) * view_segment));
        //transmittance_view_mie.print_vec();
        scattering_view_rayleigh += vec3(
            beta_rayleigh.x() * exp(-view_sample_height / Hr),
            beta_rayleigh.y() * exp(-view_sample_height / Hr),
            beta_rayleigh.z() * exp(-view_sample_height / Hr));
        //std::cerr<<"br "<<beta_rayleigh.x()<<" vsh "<<view_sample_height<<" hr "<<Hr<<std::endl;
        //scattering_view_rayleigh.print_vec();

        scattering_view_mie += vec3(
            beta_mie.x() * exp(-view_sample_height / Hm),
            beta_mie.y() * exp(-view_sample_height / Hm),
            beta_mie.z() * exp(-view_sample_height / Hm));
        //if(isnan(beta_mie.x() * exp(-view_sample_height / Hm)))
        //    std::cerr<<"nan sca"<<std::endl;
        //if(isnan(scattering_view_mie.x()))
        //    std::cerr<<"nan vec"<<std::endl;
        //std::cerr<<"bm "<<beta_mie.x()<<" vsh "<<view_sample_height<<" hm "<<Hm<<" b*exp "<<beta_mie.x()*exp(-view_sample_height / Hm)<<std::endl;
        //scattering_view_mie.print_vec();
        
        // Find intersection with atmosphere in sun direction.
        float t0light, t1light, tlight_curr;
        ray view_sample_sun_ray = ray(view_sample_position, sun_direction);
        IntersectionContext light_intersection_context;
        atmosphere.intersect(view_sample_sun_ray, 0.0, infinity, light_intersection_context); // is should always intersect
        t0light = 0.0;
        t1light = light_intersection_context.t;
        tlight_curr = t1light;

        // Calculate light integration segment length.
        float light_segment = (t1light - t0light) / light_samples;
        //std::cerr<<t1light<<std::endl;

        // Prepare transmittance, phase function, scattering coeff and sun in current view variables.
        vec3 transmittance_light_rayleigh;
        vec3 transmittance_light_mie;
        vec3 scattering_light_rayleigh;
        vec3 scattering_light_mie;
        color sun_color_rayleigh; // sun color for current view sample
        color sun_color_mie; // sun color for current view sample

        // For each light sample.
        for(int j = 0; j < light_samples; ++j)
        {
            // Calculate light sample position and height.
            vec3 light_sample_position =
                view_sample_position + (tlight_curr + light_segment / 2.0) * sun_direction;
            float light_sample_height = light_sample_position.length() - earth_radius;
            if (light_sample_height < 0.0)
                std::cerr<<"light sample height negative!" << std::endl;

            // Calculate transmittance, phase function and scattering coeff.
            transmittance_light_rayleigh += vec3(
                exp(-beta_rayleigh.x() * exp(-light_sample_height / Hr) * light_segment),
                exp(-beta_rayleigh.y() * exp(-light_sample_height / Hr) * light_segment),
                exp(-beta_rayleigh.z() * exp(-light_sample_height / Hr) * light_segment));
            //transmittance_light_rayleigh.print_vec();
            transmittance_light_mie += vec3(
                exp(-beta_mie.x() * exp(-light_sample_height / Hm) * light_segment),
                exp(-beta_mie.y() * exp(-light_sample_height / Hm) * light_segment),
                exp(-beta_mie.z() * exp(-light_sample_height / Hm) * light_segment));
            //transmittance_light_mie.print_vec();
            scattering_light_rayleigh += vec3(
                beta_rayleigh.x() * exp(-light_sample_height / Hr),
                beta_rayleigh.y() * exp(-light_sample_height / Hr),
                beta_rayleigh.z() * exp(-light_sample_height / Hr));
            //scattering_light_rayleigh.print_vec();
            scattering_view_mie += vec3(
                beta_mie.x() * exp(-light_sample_height / Hm),
                beta_mie.y() * exp(-light_sample_height / Hm),
                beta_mie.z() * exp(-light_sample_height / Hm));

            // Sum light sample contribution.
            sun_color_rayleigh += 
                transmittance_light_rayleigh *
                scattering_light_rayleigh *
                phaseR *
                sun_intensity *
                light_segment;

            sun_color_mie += 
                transmittance_light_mie *
                scattering_light_mie *
                phaseM *
                sun_intensity *
                light_segment;

            // Increase current light.
            tlight_curr += light_segment;
        }

        //sun_color_rayleigh.print_vec();

        // Sum view sample contribution for Rayleigh and Mie scattering.
        final_color_rayleigh += 
            transmittance_view_rayleigh * 
            scattering_view_rayleigh * 
            phaseR *
            sun_color_rayleigh *
            view_segment;

        final_color_mie += 
            transmittance_view_mie * 
            scattering_view_mie * 
            phaseM *
            sun_color_mie *
            view_segment;

        // Increase current view.
        tview_curr += view_segment;
    }

    final_color_rayleigh.print_vec();
    
    // Return final color with sun intensity.
    return final_color_rayleigh + final_color_mie;
}

#endif

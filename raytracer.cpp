
// Raytracer include files.
#include "vector.h"
#include "color.h"
#include "ray.h"
#include "environment.h"
#include "shader.h"
#include "intersection.h"
#include "material.h"
#include "camera.h"
#include "utils.h"


// Standard headers.
#include <iostream>
#include <memory> // shared_ptr


int main()
{
    // Image properties.
    // Using squared pixels.
    // Aspect ratio of viewport/film is equal to image aspect ratio.
    const auto aspect_ratio = 16.0 / 9.0; // w / h
    const int image_width = 800;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 1;
    const int max_depth = 2;

    // Camera properties.
    point3 look_from = point3(0.0, 0.0, 0.0); // TODO
    point3 look_at = point3(0.0, 1.0, 0.0);
    vec3 vector_up = vec3(1, 0, 0); // up != at!
    Camera camera(look_from, look_at, vector_up, 30.0, aspect_ratio);

    // Scene definition.
    IntersectableList scene;
    //auto sphere_earth = std::make_shared<Sphere>(point3(0.0, 0.0, 0.0), 6360e3, MaterialType::diffuse);
    //auto sphere2_atmosphere = std::make_shared<Sphere>(point3(0.0, 0.0, 0.0), 6420e3, MaterialType::diffuse);
    auto sphere = std::make_shared<Sphere>(point3(0.0, 10.0, 0.0), 2.0, MaterialType::volume);
    //auto sphere_light = std::make_shared<Sphere>(point3(0.3, 0.3, -1), 0.3, MaterialType::nonReflective_diffuse_emission);
    //scene.add(sphere_earth);
    //scene.add(sphere2_atmosphere);
    scene.add(sphere);
    //scene.add(sphere_light);
    NishitaSky nishita_sky(6360e3, 6420e3, 4, 2, vec3(1.0, 0.0, 0.0), vec3(20), 7994, 1200);
    

    // Rendering.
    // Use PPM, P3 image frame buffer. ASCII file.

    // PPM file header: 
    // Type:       P3\n
    // Resolution: image_width' 'image_height\n
    // Max color:  255\n
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    // RGB triplets: 
    // R G B ... R G B 
    // ...
    // R G B ... R G B 
    for (int j = image_height-1; j >= 0; --j)
    {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color(0.0, 0.0, 0.0);
            for(int s = 0; s < samples_per_pixel; ++s)
            {
                // Current point on film calculated from current pixel sample.
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);

                // Construct primary ray using eye-film point
                ray r = camera.construct_camera_ray(u, v);
                
                // Intersect with primary (camera) ray and shade.
                double t_min = 0.0;
                double t_max = 10000.0;
                IntersectionContext intersection_context;
                if(scene.intersect(r, 0.001, infinity, intersection_context))
                {   
                    // Intersects an object that can have different shading functions.
                    // Note that if recursive shading is used then environment can also be
                    // intersected and calculated here.
                    pixel_color += 
                        recursive_shading_engine(
                            r, 
                            scene, 
                            intersection_context,
                            max_depth, 
                            EnvironmentType::sky,
                            nishita_sky);
                    //std::cerr<<"Nema objekata, dakle ovo se ne bi trebalo ispisivati"<<std::endl;
                }
                else
                {
                    // Intersects environment.
                    // NPR Idea: environment can be different then environment used in object shading.
                    pixel_color += environment_shading_engine(r, EnvironmentType::sky, nishita_sky);
                }
            }                   
            // Write pixel value to image frame buffer.
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }
    std::cerr << "\nDone.\n";
}


#ifndef COLOR_H
#define COLOR_H

// Raytracer include headers.
#include "vector.h"
#include "utils.h"

// standard include headers.
#include <iostream>

void write_color(std::ostream &out, color pixel_color, int samples_per_pixel)
{
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Multiple samples per pixel are used and 
    // contributed in full color - perform scaling by 
    // dividing by number of samples.
    //
    // Perform gamma correction. Required of image viewers.
    // gamma 2 = pow(color, 1/gamma) = sqrt(color)
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Map [0,1] range to [0, 255] for PPM.
    // ART IDEA: write only one chanell (r, r, r)
    //std::cerr << g << " " << clamp(g, 0.0, 0.999) << std::endl;
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';

}

#endif

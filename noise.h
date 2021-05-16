
#ifndef NOISE_H
#define NOISE_H

// Renderer headers.
#include "vector.h"
#include "utils.h"

// STD hearder.
#include "math.h"

// Noise. https://gist.github.com/patriciogonzalezvivo/670c22f3966e662d2f83
// R(1) -> R(1)
float noise1(float n)
{
    return fract(std::sin(n) * 43758.5453123);
}

// R(3) -> R(3)
point3 noise3(point3 p)
{
    return point3(
        noise1(p.x()),
        noise1(p.y()),
        noise1(p.z()));
}

// Worley (Voronoi) noise.
// R(3) -> R(1)
float worley3(point3 sample, float density)
{
    point3 cell_center = floor(sample * density) + vec3(0.5, 0.5, 0.5);
    float F1 = float_infinity; // utils

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                point3 curr_cell_center = cell_center + vec3(i, j, k);
                point3 curr_feature = curr_cell_center + noise3(curr_cell_center);
                vec3 delta = curr_feature - sample;
                float dist = delta.length();
                if(dist < F1) F1 = dist;
            }
        }
    }
    return F1;
}

#endif

#ifndef UTILS_H
#define UTILS_H

#include <limits>
#include <cstdlib>
#include <math.h>
#include <iostream>

// Constants.
const double infinity = std::numeric_limits<double>::infinity();
const float float_infinity = std::numeric_limits<float>::infinity();
const double pi = 3.1415926535897932385;

// Functions.
inline double degrees_to_radians(double degrees)
{
    return degrees * pi / 180.0;
}

// Return random [0, 1>.
inline double random_double()
{
    return rand() / (RAND_MAX + 1.0);
}

// Return random [min, max].
inline double random_double(double min, double max)
{
    return (max - min) * random_double() + min;
}

inline double clamp(double x, double min, double max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

// Fraction part.
float fract(float a)
{
    return a - std::floor(a);
}

/*
inline double fmin(double a, double b)
{
    if(a < b) return a;
    else return b;
}

inline double fabs(double a)
{
    if (a < 0.0) return -a;
    else return a;
}
*/

#endif

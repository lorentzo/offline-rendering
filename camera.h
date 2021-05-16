
#ifndef CAMERA_H
#define CAMERA_H

// Raytracer include files.
#include "vector.h"
#include "ray.h"

class Camera
{
  public:
    // Constructor.
    Camera(
      point3 look_from,
      point3 look_at,
      vec3 vector_up,
      double _vertical_fov_deg,
      double _aspect_ratio)
    {
        // Viewport/film dimensions.
        auto theta = degrees_to_radians(_vertical_fov_deg);
        auto h = tan(theta / 2.0);
        auto viewport_height = 2.0 * h;
        auto viewport_width = _aspect_ratio * viewport_height;

        // Focal length.
        //auto focal_length = 1.0;

        // Position.
        auto w = normalize(look_from - look_at);
        auto u = normalize(cross(vector_up, w));
        auto v = cross(w, u);

        // Eye to film.
        eye = look_from; // right hand coordinate system
        horizontal = viewport_width * u; //  vec3(viewport_width, 0.0, 0.0)
        vertical = viewport_height * v; // vec3(0.0, viewport_height, 0.0)
        lower_left_corner = eye - horizontal / 2 - vertical / 2 - w;// - vec3(0.0, 0.0, focal_length);
    }

    // Functions.
    ray construct_camera_ray(double u, double v) const
    {
        return ray(eye, lower_left_corner + u * horizontal + v * vertical - eye);
    }

  private:
    point3 eye;
    point3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
};

#endif

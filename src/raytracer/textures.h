#include "raytracer.h"

struct uv
{
    float u;
    float v;
};

uv get_uv(const Hit &hit);
float get_scale_factor(const RayTraceScene &scene, const Hit &hit, const glm::vec3 &norm, const glm::vec3 &d,
                       const uv &uv_map, const int samples);
glm::vec3 get_texture(const Hit &hit, const FilterType filter_type, uv uv_map, const float scale_factor);
glm::vec3 get_bump_normal(const Hit &hit, const FilterType filter_type, uv uv_map, float bump_scale);

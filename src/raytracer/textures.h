#include "raytracer.h"

struct uv
{
    float u;
    float v;
};

uv get_uv(const Hit &hit);
float get_scale_factor(const RayTraceScene &scene, const Hit &hit, const glm::vec3 &norm, const glm::vec3 &d,
                       const uv &uv_map, const int samples);
glm::vec3 get_texture(const RenderShapeData *shape, const TextureFilterType filter_type, const float u, const float v,
                      const float scale_factor);

#include "raytracer.h"
#include "raytracescene.h"

glm::vec3 shadePixel(const RayTraceScene &scene, const Hit &hit, const glm::vec3 &world_camera,
                     unsigned int &last_ray_id, std::vector<unsigned int> &last_visited);

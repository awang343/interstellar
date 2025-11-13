#include "raytracer.h"
#include "raytracescene.h"

std::optional<Hit> checkIntersection(const RayTraceScene &scene, const glm::vec3 &dir);

std::optional<Hit> traverseKDTree(const RayTraceScene &scene, const glm::vec4 &start,
                                  const glm::vec3 &dir, bool cam_ray, unsigned int &ray_id,
                                  std::vector<unsigned int> &last_visited);

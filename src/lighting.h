#pragma once
#include "src/renderer.h"
#include "utils/scenefileparser.h"
#include <glm/glm.hpp>
// Both camera_direction and sphere_point should be in object space
glm::vec3 shadePixel(const Hit &hit, const ImageData &texture, const BumpMap &bump_map,
const std::vector<SceneLightData> &lights);

glm::vec3 shadePixel(const Hit &hit, const ImageData &texture, const BumpMap &bump_map,
const std::vector<SceneLightData> &lights, PrimitiveType objectType, const glm::vec3 &surfaceNormal);
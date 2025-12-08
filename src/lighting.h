#pragma once
#include "src/renderer.h"
#include <glm/glm.hpp>

// Enum of the types of virtual lights that might be in the scene
enum class LightType
{
    LIGHT_POINT,
    LIGHT_DIRECTIONAL,
    LIGHT_SPOT,
};

struct SceneLightData
{
    LightType type;

    glm::vec4 color;
    glm::vec3 function; // Attenuation function

    glm::vec4 pos;
    glm::vec3 dir;

    float penumbra;
    float angle;
};

struct Hit
{
    glm::vec4 point;
    glm::vec3 to_camera;
    SphereData &sphere;
};

struct BumpMap
{
    glm::vec2 *gradients;
    int width;
    int height;
};

// Both camera_direction and sphere_point should be in object space
glm::vec3 shadePixel(const Hit &hit, const ImageData &texture, const BumpMap &bump_map,
                     const std::vector<SceneLightData> &lights);

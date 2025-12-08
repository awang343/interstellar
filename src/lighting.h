#pragma once
#include <glm/glm.hpp>
#include "src/renderer.h"

// Enum of the types of virtual lights that might be in the scene
enum class LightType {
    LIGHT_POINT,
    LIGHT_DIRECTIONAL,
    LIGHT_SPOT,
};

struct SceneLightData {
    int id;
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


// Both camera_direction and sphere_point should be in object space
glm::vec3 shadePixel(const glm::vec3 &sphere_point, const glm::vec3 &dir_to_camera);

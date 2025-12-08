#include <glm/glm.hpp>

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

    glm::vec4 pos; // Position with CTM applied (Not applicable to directional lights)
    glm::vec4 dir; // Direction with CTM applied (Not applicable to point lights)

    float penumbra; // Only applicable to spot lights, in RADIANS
    float angle;    // Only applicable to spot lights, in RADIANS
};

struct Hit
{
    glm::vec3 point;
    glm::vec3 normal;
};


// Both camera_direction and sphere_point should be in object space
glm::vec3 shadePixel(const glm::vec3 &sphere_point, const glm::vec3 &dir_to_camera);

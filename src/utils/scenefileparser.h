#pragma once

#include "utils/rgba.h"
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <glm/glm.hpp>

// a struct containing info of an image
struct ImageData {
    int width  = 0;
    int height = 0;
    std::vector<RGBA> pixels;  // size = width * height
};

struct BumpMap
{
    glm::vec2 *gradients;
    int width;
    int height;
};

struct SceneInfo {
    QString upperTexturePath;
    QString lowerTexturePath;
    QString primitiveTexturePath;
    QString outputPath;

    float rho;
    float a;
    float M;

    int outWidth;
    int outHeight;

    float viewPlaneWidthAngle;  // in degrees

    float dt;
    float cameraDistance;
    float cameraTheta = M_PI/2.0f;
    float cameraPhi = 0.0f;
};

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


bool loadSceneInfoFromJson(const QString &configPath, SceneInfo &scene);

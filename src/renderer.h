#pragma once
#include <glm/glm.hpp>
#include "raytrace.h"

#include "raytrace.h"
#include "utils/rgba.h"
#include "utils/scenefileparser.h"

struct HitObjectData {
    glm::vec3 center;
    float size;
    float l;
};

struct Hit
{
    glm::vec4 point;
    glm::vec3 to_camera;
    HitObjectData &obj_data;
};

// this function samples the celestial sphere texture gased on given theta and phi
RGBA sampleCelestial(const ImageData &img, float theta, float phi);

// this function traces the rays using the equatorial symmetry method
// see https://www.youtube.com/watch?v=PVO8nvb1o2w for the explanation of this method
void render(RGBA *framebuffer,
            int outWidth,
            int outHeight,
            const ImageData &sphereUpper,
            const ImageData &sphereLower,
            const ImageData &primitiveTexture,
            float fovW,
            WormholeParams wp,
            float dt,
            float cameraDistance,
            HitObjectData sphereData,
            std::vector<SceneLightData> lights);

void render(RGBA *framebuffer,
            int outWidth,
            int outHeight,
            const ImageData &sphereUpper,
            const ImageData &sphereLower,
            const ImageData &primitiveTexture,
            float fovW,
            WormholeParams wp,
            float dt,
            float cameraDistance,
            float cameraTheta,
            float cameraPhi,
            HitObjectData sphereData,
            std::vector<SceneLightData> lights);


bool renderSingleImage(
    QImage outputImage, QString outputPath,
    RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
    const ImageData &sphereLower, const ImageData &primitiveTexture, float fovW,
    WormholeParams wp, float dt, float cameraDistance, HitObjectData sphereData, std::vector<SceneLightData> lights);


bool renderPath(
    QImage outputImage, QString outputPath,
    RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
    const ImageData &sphereLower, const ImageData &primitiveTexture, float fovW,
    WormholeParams wp, float dt, float cameraDistance, HitObjectData sphereData, std::vector<SceneLightData> lights,
    std::vector<glm::vec4> &keyFrames, int numPhotos);

bool renderFrames(
    QImage &outputImage,
    FrameData &frameData,
    RGBA *framebuffer,
    int outWidth,
    int outHeight,
    const ImageData &sphereUpper,
    const ImageData &sphereLower,
    float fovW,
    float dt,
    std::vector<SceneLightData> &lights);

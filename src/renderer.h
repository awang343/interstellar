#pragma once
#include <glm/glm.hpp>
#include "raytrace.h"

#include "raytrace.h"
#include "utils/rgba.h"
#include "utils/scenefileparser.h"

struct SphereData {
    glm::vec3 center;
    float radius;
    float l;
};

struct CubeData {
    glm::vec3 center;
    glm::vec3 halfExtents;
    float l;  
    float radius;  
};

enum class PrimitiveType {
    SPHERE,
    CUBE
};


struct Hit
{
    glm::vec4 point;
    glm::vec3 to_camera;
    PrimitiveType type;  // Which primitive was hit
    
    // Use a union since only one will be valid at a time
    union {
        SphereData sphere;
        CubeData cube;
    };
    
    // Constructor for sphere hits
    Hit(glm::vec4 p, glm::vec3 cam, const SphereData& s) 
        : point(p), to_camera(cam), type(PrimitiveType::SPHERE), sphere(s) {}
    
    // Constructor for cube hits
    Hit(glm::vec4 p, glm::vec3 cam, const CubeData& c) 
        : point(p), to_camera(cam), type(PrimitiveType::CUBE), cube(c) {}
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
            SphereData sphereData,
            CubeData cubeData,
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
            SphereData sphereData,
            CubeData cubeData,
            std::vector<SceneLightData> lights);


bool renderSingleImage(
    QImage outputImage, QString outputPath,
    RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
    const ImageData &sphereLower, const ImageData &primitiveTexture, float fovW,
    WormholeParams wp, float dt, float cameraDistance, SphereData sphereData, CubeData cubeData, std::vector<SceneLightData> lights);


bool renderPath(
    QImage outputImage, QString outputPath,
    RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
    const ImageData &sphereLower, const ImageData &primitiveTexture, float fovW,
    WormholeParams wp, float dt, float cameraDistance, SphereData sphereData, CubeData cubeData,std::vector<SceneLightData> lights,
    std::vector<glm::vec4> &keyFrames, int numPhotos);

#pragma once
#include <glm/glm.hpp>
#include "raytrace.h"

#include "raytrace.h"
#include "utils/rgba.h"
#include "utils/scenefileparser.h"

struct SphereData {
    glm::vec3 center;
    double radius;
    double l;
};

// this function samples the celestial sphere texture gased on given theta and phi
RGBA sampleCelestial(const ImageData &img, double theta, double phi);

// this function traces the rays using the equatorial symmetry method
// see https://www.youtube.com/watch?v=PVO8nvb1o2w for the explanation of this method
void renderEquatorial(RGBA *framebuffer,
                      int outWidth,
                      int outHeight,
                      const ImageData &sphereUpper,
                      const ImageData &sphereLower,
                      double fovW,
                      WormholeParams wp,
                      double dt,
                      double cameraDistance,
                      SphereData sphereData);


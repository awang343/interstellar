#pragma once

#include <QtCore>
#include <math.h>
#include <vector>
#include <glm/glm.hpp>
#include "utils/structs.h"

class FrameData
{
public:
    FrameData();

    glm::vec3 cameraAtT(float t);
    std::vector<Sphere> objectsAtT(float t);

    void loadDefaultObject();
    void loadDefaultCamera(float cameraDistance);

    void setTimeRange();
    float tMin;
    float tMax;

    int numPhotos;
    QString outputPath;
    std::vector<glm::vec4> cameraPoints;
    std::vector<Sphere> objects;
};



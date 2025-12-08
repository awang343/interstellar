#pragma once

#include "utils/rgba.h"
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>

// a struct containing info of an image
struct ImageData {
    int width  = 0;
    int height = 0;
    std::vector<RGBA> pixels;  // size = width * height
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


bool loadSceneInfoFromJson(const QString &configPath, SceneInfo &scene);

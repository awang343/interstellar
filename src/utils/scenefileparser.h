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

struct sceneInfo {
    QString upperTexturePath;
    QString lowerTexturePath;
    QString outputPath;

    double rho;
    double a;
    double M;

    int outWidth;
    int outHeight;

    double viewPlaneWidthAngle;  // in degrees

    double dt;
    double cameraDistance;
    double cameraTheta = M_PI/2.0f;
    double cameraPhi = 0.0f;
};


bool loadSceneInfoFromJson(const QString &configPath, sceneInfo &scene);

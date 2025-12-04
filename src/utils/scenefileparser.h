#pragma once

#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>

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
};


bool loadSceneInfoFromJson(const QString &configPath, sceneInfo &scene);

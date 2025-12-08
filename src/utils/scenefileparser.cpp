#include "scenefileparser.h"
#include <iostream>

bool loadSceneInfoFromJson(const QString &path, sceneInfo &outScene)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        std::cerr << "Could not open JSON config: "
                  << path.toStdString() << "\n";
        return false;
    }

    QByteArray raw = file.readAll();
    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(raw, &err);

    if (err.error != QJsonParseError::NoError) {
        std::cerr << "JSON parse error: "
                  << err.errorString().toStdString() << "\n";
        return false;
    }

    if (!doc.isObject()) {
        std::cerr << "Config file must contain a JSON object.\n";
        return false;
    }

    QJsonObject obj = doc.object();

    auto getString = [&](const QString &key, QString &dst) {
        if (!obj.contains(key) || !obj[key].isString()) {
            std::cerr << "Missing or invalid string field: "
                      << key.toStdString() << "\n";
            return false;
        }
        dst = obj[key].toString();
        return true;
    };

    auto getFloat = [&](const QString &key, float &dst) {
        if (!obj.contains(key) || !obj[key].isDouble()) {
            std::cerr << "Missing or invalid float field: "
                      << key.toStdString() << "\n";
            return false;
        }
        dst = static_cast<float>(obj[key].toDouble());
        return true;
    };

    auto getInt = [&](const QString &key, int &dst) {
        if (!obj.contains(key) || !obj[key].isDouble()) {
            std::cerr << "Missing or invalid int field: "
                      << key.toStdString() << "\n";
            return false;
        }
        dst = obj[key].toInt();
        return true;
    };

    // Read paths
    if (!getString("upperTexture",  outScene.upperTexturePath)) return false;
    if (!getString("lowerTexture",  outScene.lowerTexturePath)) return false;
    if (!getString("outputImage",   outScene.outputPath))       return false;

    // Wormhole parameters
    if (!getFloat("rho", outScene.rho)) return false;
    if (!getFloat("a",   outScene.a))   return false;
    if (!getFloat("M",   outScene.M))   return false;

    // Resolution
    if (!getInt("outWidth",  outScene.outWidth))  return false;
    if (!getInt("outHeight", outScene.outHeight)) return false;

    // FOV (in degrees), convert to radians
    float fovDeg;
    if (!getFloat("viewPlaneWidthAngle", fovDeg)) return false;
    outScene.viewPlaneWidthAngle = fovDeg * M_PI / 180.0;

    // NEW fields
    if (!getFloat("dt", outScene.dt)) return false;
    if (!getFloat("cameraDistance", outScene.cameraDistance)) return false;

    return true;
}



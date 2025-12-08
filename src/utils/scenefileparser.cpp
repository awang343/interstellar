#include "scenefileparser.h"
#include <iostream>
#include <QFileInfo>
#include <QDir>
#include <qjsonarray.h>

bool loadSceneInfoFromJson(const QString &path, SceneInfo &outScene)
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

    // ----------------------
    //   Helper lambdas
    // ----------------------

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

    auto getStringOptional = [&](const QString &key, QString &dst) {
        if (!obj.contains(key) || !obj[key].isString()) {
            obj[key] = "";
        }
        dst = obj[key].toString();
        return true;
    };

    auto getIntOptional = [&](const QString &key, int &dst) {
        if (!obj.contains(key) || !obj[key].isDouble()) {
            dst = 1;
            return false;
        }
        dst = obj[key].toInt();
        return true;
    };

    auto getBoolOptional = [&](const QString &key, bool &dst) {
        if (!obj.contains(key) || !obj[key].isBool()) {
            dst = false;
            return false;
        }
        dst = obj[key].toBool();
        return true;
    };

    auto getPathsOptional = [&](const QString &key, std::vector<glm::vec4> &dst) {
        dst.clear();

        if (!obj.contains(key)) {
            std::cerr << "Optional key '" << key.toStdString()
            << "' not found. Using empty path list.\n";
            return false;
        }

        if (!obj[key].isArray()) {
            std::cerr << "Field '" << key.toStdString()
            << "' must be an array.\n";
            return false;
        }

        QJsonArray outer = obj[key].toArray();
        dst.reserve(outer.size());

        for (int i = 0; i < outer.size(); i++) {
            if (!outer[i].isArray()) {
                std::cerr << "pathPoints[" << i << "] is not an array.\n";
                return false;
            }

            QJsonArray inner = outer[i].toArray();

            if (inner.size() != 4) {
                std::cerr << "pathPoints[" << i
                          << "] must contain exactly 4 numbers.\n";
                return false;
            }

            glm::vec4 v;
            for (int j = 0; j < 4; j++) {
                if (!inner[j].isDouble()) {
                    std::cerr << "pathPoints[" << i << "][" << j
                              << "] must be a number.\n";
                    return false;
                }
            }

            v.x = float(inner[0].toDouble());
            v.y = float(inner[1].toDouble());
            v.z = float(inner[2].toDouble());
            v.w = float(inner[3].toDouble());

            dst.push_back(v);
        }

        return true;
    };


    // ----------------------
    //   Texture paths
    // ----------------------

    if (!getString("upperTexture",     outScene.upperTexturePath))     return false;
    if (!getString("lowerTexture",     outScene.lowerTexturePath))     return false;
    if (!getString("primitiveTexture", outScene.primitiveTexturePath)) return false;  // NEW
    if (!getString("outputImage",      outScene.outputPath))           return false;

    // ----------------------
    //   Wormhole parameters
    // ----------------------

    if (!getFloat("rho", outScene.rho)) return false;
    if (!getFloat("a",   outScene.a))   return false;
    if (!getFloat("M",   outScene.M))   return false;

    // ----------------------
    //   Output resolution
    // ----------------------

    if (!getInt("outWidth",  outScene.outWidth))  return false;
    if (!getInt("outHeight", outScene.outHeight)) return false;

    // ----------------------
    //   FOV (degrees → radians)
    // ----------------------

    float fovDeg;
    if (!getFloat("viewPlaneWidthAngle", fovDeg)) return false;
    outScene.viewPlaneWidthAngle = fovDeg * M_PI / 180.0f;

    // ----------------------
    //   New fields
    // ----------------------

    if (!getFloat("dt", outScene.dt)) return false;
    if (!getFloat("cameraDistance", outScene.cameraDistance)) return false;

    // ----------------------
    //   Camera paths (optional)
    // ----------------------

    if (getBoolOptional("useCameraPath", outScene.usePaths) && outScene.usePaths){
        std::cout << "Using camera paths" << std::endl;
        getPathsOptional("pathPoints", outScene.paths);

        QFileInfo fi(outScene.outputPath);
        QDir outDir = fi.dir();
        if (!outDir.exists()) {
            outDir.mkpath(".");
        }

        getIntOptional("numPhotos", outScene.numPhotos);
        getStringOptional("outputFolder", outScene.outputPath);
    } else {
        std::cout << "Not using camera paths" << std::endl;
    }


    return true;
}



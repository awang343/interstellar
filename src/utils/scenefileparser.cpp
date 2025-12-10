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

    auto getObjectsOptional = [&](const QString &key, std::vector<Object> &dst) {
        dst.clear();

        if (!obj.contains(key)) {
            std::cerr << "Optional key '" << key.toStdString()
            << "' not found. Using empty objects list.\n";
            return false;
        }

        if (!obj[key].isArray()) {
            std::cerr << "Field '" << key.toStdString()
            << "' must be an array.\n";
            return false;
        }

        QJsonArray objectsArray = obj[key].toArray();
        dst.reserve(objectsArray.size());

        for (int i = 0; i < objectsArray.size(); ++i) {
            if (!objectsArray[i].isObject()) {
                std::cerr << "objects[" << i << "] is not an object.\n";
                return false;
            }

            QJsonObject objEntry = objectsArray[i].toObject();
            Object o;

            // --- type (optional, defaults to Sphere) ---
            if (objEntry.contains("type") && objEntry["type"].isString()) {
                QString typeStr = objEntry["type"].toString().toLower();
                if (typeStr == "cube") {
                    o.type = PrimitiveType::Cube;
                } else {
                    o.type = PrimitiveType::Sphere;
                }
            } else {
                o.type = PrimitiveType::Sphere;
            }

            // --- radius (for spheres) ---
            if (o.type == PrimitiveType::Sphere) {
                if (objEntry.contains("radius") && objEntry["radius"].isDouble()) {
                    o.radius = static_cast<float>(objEntry["radius"].toDouble());
                } else {
                    o.radius = 0.5f;
                }
            }

            // --- side (for cubes) ---
            if (o.type == PrimitiveType::Cube) {
                if (objEntry.contains("side") && objEntry["side"].isDouble()) {
                    o.side = static_cast<float>(objEntry["side"].toDouble());
                } else {
                    o.side = 1.0f;
                }
            }

            // --- textureFile ---
            if (!objEntry.contains("textureFile") || !objEntry["textureFile"].isString()) {
                std::cerr << "objects[" << i << "].textureFile is missing or not a string.\n";
                return false;
            }

            QString texPath = objEntry["textureFile"].toString();

            if (!loadImageToStruct(texPath, o.textureFile)) {
                std::cerr << "issue loading " << texPath.toStdString() << "\n";
                return false;
            }

            // --- bumpMapFile ---
            if (!objEntry.contains("bumpMapFile") || !objEntry["bumpMapFile"].isString()) {
                std::cerr << "objects[" << i << "].bumpMapFile is missing or not a string.\n";
                return false;
            }

            QString bumpPath = objEntry["bumpMapFile"].toString();

            if (!loadBumpMapToStruct(bumpPath, o.bumpMapFile)) {
                std::cerr << "issue loading " << bumpPath.toStdString() << "\n";
                return false;
            }

            // --- objectPoints ---
            if (!objEntry.contains("objectPoints") || !objEntry["objectPoints"].isArray()) {
                std::cerr << "objects[" << i << "].objectPoints is missing or not an array.\n";
                return false;
            }

            QJsonArray pointsArray = objEntry["objectPoints"].toArray();
            std::vector<std::vector<float>> points;
            points.reserve(pointsArray.size());

            for (int p = 0; p < pointsArray.size(); ++p) {
                if (!pointsArray[p].isArray()) {
                    std::cerr << "objects[" << i << "].objectPoints[" << p
                              << "] is not an array.\n";
                    return false;
                }

                QJsonArray inner = pointsArray[p].toArray();
                std::vector<float> row;
                row.reserve(inner.size());

                for (int j = 0; j < inner.size(); ++j) {
                    if (!inner[j].isDouble()) {
                        std::cerr << "objects[" << i << "].objectPoints[" << p
                                  << "][" << j << "] must be a number.\n";
                        return false;
                    }
                    row.push_back(static_cast<float>(inner[j].toDouble()));
                }

                points.push_back(std::move(row));
            }

            o.points = std::move(points);

            dst.push_back(std::move(o));
        }

        return true;
    };

    // ----------------------
    //   Texture paths
    // ----------------------

    if (!getString("upperTexture",     outScene.upperTexturePath))     return false;
    if (!getString("lowerTexture",     outScene.lowerTexturePath))     return false;
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
    //  FrameData (optional but always created)
    // ----------------------

    getIntOptional("numPhotos", outScene.frameData.numPhotos);

    if (outScene.frameData.numPhotos > 1) {
        std::cout << "numPhotos > 1" << std::endl;
        // output will have multiple frames, ensure there exists an output folder
        if (!getString("outputFolder", outScene.frameData.outputPath)) {
            std::cerr << "Multiple frames to be rendered but no outputFolder specified\n";
        }

        // ensures there exists an output folder at specified
        QDir outDir(outScene.frameData.outputPath);
        if (!outDir.exists()) {
            outDir.mkpath(".");
        }



    } else {
        // if numPhotos not specified, set to 1
        std::cout << "numPhotos = 1" << std::endl;
        outScene.frameData.numPhotos = 1;

        // set output to outputImage
        outScene.frameData.outputPath = outScene.outputPath;
    }

    // ----------------------
    //  Objects (optional, default to green sphere at (0,2,0) with r=0.5)
    // ----------------------

    if (!getObjectsOptional("objects", outScene.frameData.objects)) {
        std::cout << "Either no objects specified or issue with parsing. Using default." << std::endl;

        // if no objects in use, load basic object
        outScene.frameData.loadDefaultObject();
    } else {
        std::cout << "Objects loaded" << std::endl;
    }

    // ----------------------
    //   Camera paths (optional)
    // ----------------------

    if (!getPathsOptional("cameraPoints", outScene.frameData.cameraPoints)) {
        std::cout << "Either no camera points specified or issue with parsing. Using cameraDistance." << std::endl;
         outScene.frameData.loadDefaultCamera(outScene.cameraDistance);
    } else {
        std::cout << "Camera points loaded" << std::endl;
    }

    // Set tmin, tmax for FrameData
    outScene.frameData.setTimeRange();

    return true;
}



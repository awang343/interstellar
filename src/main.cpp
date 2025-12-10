#include <QCoreApplication>
#include <QCommandLineParser>
#include <QImage>
#include <QtCore>

#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "src/utils/rgba.h"
#include "src/raytrace.h"
#include "utils/scenefileparser.h"
#include "renderer.h"

// the main function for the program
int main(int argc, char *argv[]) {
    QCoreApplication app(argc, argv);

    QCommandLineParser parser;
    parser.setApplicationDescription("Relativistic wormhole raytracer (Dneg metric)");
    parser.addHelpOption();

    // Now we take a single positional argument: the JSON config file
    parser.addPositionalArgument("config",
                                 "Path to the JSON configuration file.");

    parser.process(app);

    const QStringList args = parser.positionalArguments();
    if (args.size() != 1) {
        std::cerr << "Usage: raytracer <config.json>\n";
        return 1;
    }

    const QString configPath = args[0];

    // 1. Load scene configuration from JSON
    SceneInfo scene;
    if (!loadSceneInfoFromJson(configPath, scene)) {
        return 1; // errors already printed by loader
    }

    // 2. Load input textures
    ImageData sphereUpper;
    ImageData sphereLower;
    ImageData primitiveTexture;

    if (!loadImageToStruct(scene.upperTexturePath, sphereUpper)) return 1;
    if (!loadImageToStruct(scene.lowerTexturePath, sphereLower)) return 1;

    // 3. Allocate output QImage using configured resolution
    QImage outputImage(scene.outWidth, scene.outHeight, QImage::Format_RGBA8888);
    if (outputImage.isNull()) {
        std::cerr << "Failed to allocate output QImage\n";
        return 1;
    }

    RGBA *framebuffer = reinterpret_cast<RGBA*>(outputImage.bits());

    WormholeParams wp{scene.rho, scene.a, scene.M};

    // make point light
    SceneLightData light{LightType::LIGHT_POINT, glm::vec4(1.0), glm::vec3(0.2, 0.05, 0.0), glm::vec4(3.0, 2.0, 3.0, -1.0), glm::vec3(-1.0, 0.0, -1.0), 1.0, 1.0};

    bool ok = renderFrames(
        outputImage,
        scene.frameData,
        framebuffer,
        scene.outWidth,
        scene.outHeight,
        sphereUpper,
        sphereLower,
        scene.viewPlaneWidthAngle,
        wp,
        scene.dt,
        std::vector<SceneLightData>{light});

    if (!ok) return 1;

    return 0;
}


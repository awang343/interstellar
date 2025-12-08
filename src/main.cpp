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
#include "lighting.h"


// this function loads the two celestial sphere textures into image data
bool loadImageToStruct(const QString &path, ImageData &out)
{
    QImage img(path);
    if (img.isNull()) {
        std::cerr << "Failed to load image: " << path.toStdString() << "\n";
        return false;
    }

    QImage converted = img.convertToFormat(QImage::Format_RGBA8888);

    out.width  = converted.width();
    out.height = converted.height();
    out.pixels.resize(out.width * out.height);

    for (int y = 0; y < out.height; ++y) {
        const uchar *line = converted.constScanLine(y);
        for (int x = 0; x < out.width; ++x) {
            const uchar *src = line + 4 * x;
            RGBA &dst = out.pixels[y * out.width + x];
            dst.r = src[0];
            dst.g = src[1];
            dst.b = src[2];
            dst.a = src[3];
        }
    }
    return true;
}

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
    if (!loadImageToStruct(scene.primitiveTexturePath, primitiveTexture)) return 1;

    // 3. Allocate output QImage using configured resolution
    QImage outputImage(scene.outWidth, scene.outHeight, QImage::Format_RGBA8888);
    if (outputImage.isNull()) {
        std::cerr << "Failed to allocate output QImage\n";
        return 1;
    }

    RGBA *framebuffer = reinterpret_cast<RGBA*>(outputImage.bits());

    WormholeParams wp{scene.rho, scene.a, scene.M};

    // make a sphere
    glm::vec3 spherePos(3.0, 0.0, 0.0);
    SphereData sphereData{spherePos, 0.5, length(spherePos)};

    // 4. Render using FOV from config (scene.viewPlaneWidthAngle is in radians)
    render(framebuffer,
           scene.outWidth,
           scene.outHeight,
           sphereUpper,
           sphereLower,
           primitiveTexture,
           scene.viewPlaneWidthAngle, // in radians
           wp,
           scene.dt,
           scene.cameraDistance,
           sphereData);

    // 5. Save output
    bool ok = outputImage.save(scene.outputPath);
    if (!ok) ok = outputImage.save(scene.outputPath, "PNG");

    if (!ok) {
        std::cerr << "Failed to save output image: "
                  << scene.outputPath.toStdString() << "\n";
        return 1;
    }

    std::cout << "Saved rendered image to "
              << scene.outputPath.toStdString() << "\n";
    return 0;
}


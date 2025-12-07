#include <QCoreApplication>
#include <QCommandLineParser>
#include <QImage>
#include <QtCore>

#include <vector>
#include <iostream>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "src/utils/rgba.h"
#include "src/raytrace.h"
#include "utils/scenefileparser.h"

// a struct containing info of an image
struct ImageData {
    int width  = 0;
    int height = 0;
    std::vector<RGBA> pixels;  // size = width * height
};


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


// this function samples the celestial sphere texture gased on given theta and phi
RGBA sampleCelestial(const ImageData &img, double theta, double phi)
{
    if (img.width == 0 || img.height == 0) {
        return {0, 0, 0, 255};
    }

    // Normalize φ to [0, 2π)
    double twoPi = 2.0 * M_PI;
    double phiNorm = std::fmod(phi, twoPi);
    if (phiNorm < 0.0) phiNorm += twoPi;

    // θ in [0, π] ideally; clamp
    double thetaClamped = std::max(0.0, std::min(M_PI, theta));

    // longitude-latitude map: u = φ / (2π), v = θ / π
    double u = phiNorm / twoPi;
    double v = thetaClamped / M_PI;

    int x = static_cast<int>(u * (img.width  - 1));
    int y = static_cast<int>(v * (img.height - 1));

    x = std::max(0, std::min(img.width  - 1, x));
    y = std::max(0, std::min(img.height - 1, y));

    return img.pixels[y * img.width + x];
}

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
    double cameraDistance) {

    // Camera setup
    const double l_c     = cameraDistance;
    const double theta_c = M_PI / 2.0; // equatorial
    const double phi_c   = 0; // will eventually make changes to this

    // Aspect ratio and vertical FOV
    double aspect   = static_cast<double>(outHeight) / static_cast<double>(outWidth);

    // Loop over pixels
    double planeWidth = tan(fovW * 0.5);
    double planeHeight = planeWidth * aspect;
    double lengthDiagonal = sqrt(planeWidth*planeWidth + planeHeight*planeHeight);

    // Integration constants (-0.01 seems good enough, but very slow...)
    const double tMin = -200.0;
    const double lMax = 20.0;

    // compute the number of rays required to trace and store along the equator
    const int precMultiplier = 2; // multiplier of the precision metric
    int numRays = (outWidth + outHeight) * precMultiplier;
    std::vector<RayState> rays(numRays+1);
    std::vector<glm::vec4> rayPositions((numRays+1) * tMin);

    // trace through the angles of these rays on the equatorial plane
    for (int i = 0; i < numRays+1; i++) {
        // Pinhole camera direction
        double px = 1.0;
        double py = lengthDiagonal * static_cast<double>(i) / static_cast<double>(numRays+1);
        double pz = 0.0;

        double normP = sqrt(px*px + py*py + pz*pz);
        px /= normP;
        py /= normP;
        pz /= normP;

        // Direction of incoming ray is -p
        double n_l     = -px;
        double n_phi   = -py;
        double n_theta = pz;

        // Initial state
        RayState ray;
        ray.l     = l_c;
        ray.theta = theta_c;
        ray.phi   = phi_c;

        double r_c = r_of_l(ray.l, wp);

        ray.p_l     = n_l;
        ray.p_theta = r_c * n_theta;
        ray.b       = r_c * sin(ray.theta) * n_phi;

        // numerical integration
        double t = 0.0;
        while (t > tMin) {
            rk4Step(ray, wp, dt);
            if (abs(ray.l) > lMax) {
                break;
            }
            // compute and store the position coordinates (x,y,z,l) in world space


            t += dt;
        }

        // add the ray to the ray state
        rays[i] = ray;
    }

    // sample the texture color based on the precomputed ray directions
    for (int j = 0; j < outHeight; ++j) {
        for (int i = 0; i < outWidth; ++i) {
            int idx = j * outWidth + i;

            // Normalized device coordinates in [-1,1]
            double ndcW = 2.0 * ((i + 0.5) / static_cast<double>(outWidth))  - 1.0;
            double ndcH = 2.0 * ((j + 0.5) / static_cast<double>(outHeight)) - 1.0;

            // ray direction
            glm::vec3 rayDir(1.0, planeWidth * (ndcW), planeHeight * (-ndcH));
            glm::vec2 normalizedYZ = normalize(rayDir.yz());

            // length at the yz-direction
            double lengthYZ = length(rayDir.yz());
            // get the precomputed ray
            int rayInd = static_cast<int>(ceil(static_cast<double>(numRays) * (lengthYZ / lengthDiagonal)));
            RayState equatorialRay = rays[rayInd];

            // normalize the angles

            // first, convert the final ray direction into a unit vector
            glm::vec4 finalEuclidean(cos(equatorialRay.phi), sin(equatorialRay.phi), 0.0, 0.0);
            // next, apply the transformation with respect to the tilting of the camera
            glm::vec3 cameraDirection(cos(phi_c), sin(phi_c), 0.0);
            // the angle to rotate about the camera direction
            float angle = atan2(normalizedYZ.y, normalizedYZ.x);
            glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), angle, glm::normalize(cameraDirection));
            glm::vec4 rotatedEuclidean = rotationMatrix * finalEuclidean;

            float thetaFinal = acos(rotatedEuclidean[2] / length(rotatedEuclidean));
            float phiFinal = atan2(rotatedEuclidean[1], rotatedEuclidean[0]);

            // break into cases where l<0 and l>0
            RGBA color = {0, 0, 0, 255};

            if (abs(equatorialRay.l) > lMax) {
                if (equatorialRay.l > 0.0) {
                    color = sampleCelestial(sphereUpper, thetaFinal, phiFinal);
                } else {
                    color = sampleCelestial(sphereLower, thetaFinal, phiFinal);
                }
            } else {
                color = {0, 0, 0, 255};
            }

            framebuffer[idx] = color;
        }
    }
}

// this function traces the rays
void render(
    RGBA *framebuffer,
    int outWidth,
    int outHeight,
    const ImageData &sphereUpper,
    const ImageData &sphereLower,
    double fovW,
    WormholeParams wp,
    double dt,
    double cameraDistance)
{

    // Camera setup
    const double l_c     = cameraDistance;
    const double theta_c = M_PI / 2.0; // equatorial
    const double phi_c   = 0;

    // Aspect ratio and vertical FOV
    double aspect   = static_cast<double>(outHeight) / static_cast<double>(outWidth);

    // Loop over pixels
    double planeWidth = tan(fovW * 0.5);
    double planeHeight = planeWidth * aspect;

    // Integration constants (-0.01 seems good enough, but very slow...)
    const double tMin = -200.0;
    const double lMax = 20.0;

    for (int j = 0; j < outHeight; ++j) {
        for (int i = 0; i < outWidth; ++i) {
            int idx = j * outWidth + i;

            if (idx % 10000 == 0) {
                std::cout << "ind: " << idx << std::endl;
            }

            // Normalized device coordinates in [-1,1]
            double ndcW = 2.0 * ((i + 0.5) / static_cast<double>(outWidth))  - 1.0;
            double ndcH = 2.0 * ((j + 0.5) / static_cast<double>(outHeight)) - 1.0;

            // Pinhole camera direction
            double px = 1.0;
            double py = planeWidth * (ndcW);
            double pz = planeHeight * (-ndcH);

            double normP = sqrt(px*px + py*py + pz*pz);
            px /= normP;
            py /= normP;
            pz /= normP;

            // Direction of incoming ray is -p
            double n_l     = -px;
            double n_phi   = -py;
            double n_theta = pz;

            // Initial state
            RayState ray;
            ray.l     = l_c;
            ray.theta = theta_c;
            ray.phi   = phi_c;

            double r_c = r_of_l(ray.l, wp);

            ray.p_l     = n_l;
            ray.p_theta = r_c * n_theta;
            ray.b       = r_c * sin(ray.theta) * n_phi;

            // numerical integration
            double t = 0.0;
            while (t > tMin) {
                rk4Step(ray, wp, dt);
                if (abs(ray.l) > lMax) {
                    break;
                }
                t += dt;
            }

            RGBA color = {0, 0, 0, 255};

            if (abs(ray.l) > lMax) {
                if (ray.l > 0.0) {
                    color = sampleCelestial(sphereUpper, ray.theta, ray.phi);
                } else {
                    color = sampleCelestial(sphereLower, ray.theta, ray.phi);
                }
            } else {
                color = {0, 0, 0, 255};
            }

            framebuffer[idx] = color;
        }
    }
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
    sceneInfo scene;
    if (!loadSceneInfoFromJson(configPath, scene)) {
        return 1; // errors already printed by loader
    }

    // 2. Load input textures
    ImageData sphereUpper;
    ImageData sphereLower;

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

    // 4. Render using FOV from config (scene.viewPlaneWidthAngle is in radians)
    renderEquatorial(framebuffer,
           scene.outWidth,
           scene.outHeight,
           sphereUpper,
           sphereLower,
           scene.viewPlaneWidthAngle, // in radians
           wp,
           scene.dt,
           scene.cameraDistance);

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


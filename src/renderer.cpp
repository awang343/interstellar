#include "renderer.h"
#include "camerapath.h"
#include "lighting.h"

#include <algorithm>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <qimage.h>
#include <vector>

// this function samples the celestial sphere texture gased on given theta and phi
RGBA sampleCelestial(const ImageData &img, float theta, float phi)
{
    if (img.width == 0 || img.height == 0)
    {
        return {0, 0, 0, 255};
    }

    // Normalize φ to [0, 2π)
    float twoPi = 2.0 * M_PI;
    float phiNorm = std::fmod(phi, twoPi);
    if (phiNorm < 0.0)
        phiNorm += twoPi;

    // θ in [0, π] ideally; clamp
    float thetaClamped = std::max(0.0, std::fmin(M_PI, theta));

    // longitude-latitude map: u = φ / (2π), v = θ / π
    float u = phiNorm / twoPi;
    float v = thetaClamped / M_PI;

    float x_float = u * (img.width - 1);
    float y_float = v * (img.height - 1);

    int x0 = static_cast<int>(std::floor(x_float));
    int x1 = x0 + 1;
    int y0 = static_cast<int>(std::floor(y_float));
    int y1 = y0 + 1;

    x0 = std::max(0, std::min(img.width - 1, x0));
    x1 = std::max(0, std::min(img.width - 1, x1));
    y0 = std::max(0, std::min(img.height - 1, y0));
    y1 = std::max(0, std::min(img.height - 1, y1));

    float alpha_x = x_float - std::floor(x_float);
    float alpha_y = y_float - std::floor(y_float);

    RGBA c00 = img.pixels[y0 * img.width + x0]; // top-left
    RGBA c10 = img.pixels[y0 * img.width + x1]; // top-right
    RGBA c01 = img.pixels[y1 * img.width + x0]; // bottom-left
    RGBA c11 = img.pixels[y1 * img.width + x1]; // bottom-right

    // Bilinear interpolation
    // horizontal interpolation
    float r_top = c00.r * (1.0f - alpha_x) + c10.r * alpha_x;
    float g_top = c00.g * (1.0f - alpha_x) + c10.g * alpha_x;
    float b_top = c00.b * (1.0f - alpha_x) + c10.b * alpha_x;
    float a_top = c00.a * (1.0f - alpha_x) + c10.a * alpha_x;

    float r_bottom = c01.r * (1.0f - alpha_x) + c11.r * alpha_x;
    float g_bottom = c01.g * (1.0f - alpha_x) + c11.g * alpha_x;
    float b_bottom = c01.b * (1.0f - alpha_x) + c11.b * alpha_x;
    float a_bottom = c01.a * (1.0f - alpha_x) + c11.a * alpha_x;

    // vertical interpolation
    uint8_t r_final = static_cast<uint8_t>(r_top * (1.0f - alpha_y) + r_bottom * alpha_y);
    uint8_t g_final = static_cast<uint8_t>(g_top * (1.0f - alpha_y) + g_bottom * alpha_y);
    uint8_t b_final = static_cast<uint8_t>(b_top * (1.0f - alpha_y) + b_bottom * alpha_y);
    uint8_t a_final = static_cast<uint8_t>(a_top * (1.0f - alpha_y) + a_bottom * alpha_y);

    return {r_final, g_final, b_final, a_final};
}

bool intersectCube(const glm::vec3 &rayOrigin, const glm::vec3 &rayDir,
                   const glm::vec3 &cubeCenter, float halfSide,
                   float &tNear, float &tFar)
{
    glm::vec3 tMin = (cubeCenter - glm::vec3(halfSide) - rayOrigin) / rayDir;
    glm::vec3 tMax = (cubeCenter + glm::vec3(halfSide) - rayOrigin) / rayDir;

    glm::vec3 t1 = glm::min(tMin, tMax);
    glm::vec3 t2 = glm::max(tMin, tMax);

    tNear = glm::max(glm::max(t1.x, t1.y), t1.z);
    tFar = glm::min(glm::min(t2.x, t2.y), t2.z);

    return tFar >= tNear && tFar > 0;
}

glm::vec3 getCubeNormal(const glm::vec3 &point, const glm::vec3 &cubeCenter, float halfSide)
{
    glm::vec3 local = point - cubeCenter;
    glm::vec3 absLocal = glm::abs(local);

    if (absLocal.x > absLocal.y && absLocal.x > absLocal.z)
    {
        return glm::vec3(glm::sign(local.x), 0, 0);
    }
    else if (absLocal.y > absLocal.z)
    {
        return glm::vec3(0, glm::sign(local.y), 0);
    }
    else
    {
        return glm::vec3(0, 0, glm::sign(local.z));
    }
}

// this function traces the rays using the equatorial symmetry method
// see https://www.youtube.com/watch?v=PVO8nvb1o2w for the explanation of this method
void render(RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
            const ImageData &sphereLower, const std::vector<Object> &objects, float fovW,
            WormholeParams wp, float dt, float cameraDistance, float cameraTheta, float cameraPhi,
            std::vector<SceneLightData> lights)
{

    // Camera setup
    const float l_c = cameraDistance;
    const float theta_c = cameraTheta; // equatorial
    const float phi_c = cameraPhi;     // will eventually make changes to this

    // Aspect ratio and vertical FOV
    float aspect = static_cast<float>(outHeight) / static_cast<float>(outWidth);

    // Loop over pixels
    float planeWidth = tan(fovW * 0.5);
    float planeHeight = planeWidth * aspect;
    float lengthDiagonal = sqrt(planeWidth * planeWidth + planeHeight * planeHeight);

    // Integration constants (-0.01 seems good enough, but very slow...)
    const float tMin = -200.0;
    const float lMax = 20.0;

    // compute the number of rays required to trace and store along the equator
    const int precMultiplier = 2; // multiplier of the precision metric
    int numRays = (outWidth + outHeight) * precMultiplier;
    std::vector<RayState> rays(numRays + 1);

    int numRayPositions = tMin * static_cast<int>(1 / dt);
    std::vector<glm::vec4> rayPositions((numRays + 1) * numRayPositions);
    std::vector<int> rayPositionCounts(numRays + 1); // the total number of positions stored

    // trace through the angles of these rays on the equatorial plane
    for (int i = 0; i < numRays + 1; i++)
    {
        // Pinhole camera direction
        float px = 1.0;
        float py = lengthDiagonal * static_cast<float>(i) / static_cast<float>(numRays + 1);
        float pz = 0.0;

        float normP = sqrt(px * px + py * py + pz * pz);
        px /= normP;
        py /= normP;
        pz /= normP;

        // Direction of incoming ray is -p
        float n_l = -px;
        float n_phi = -py;
        float n_theta = pz;

        // Initial state
        RayState ray;
        ray.l = l_c;
        ray.theta = theta_c;
        ray.phi = phi_c;

        float r_c = r_of_l(ray.l, wp);
        ray.p_l = n_l;
        ray.p_theta = r_c * n_theta;
        ray.b = r_c * sin(ray.theta) * n_phi;

        // numerical integration
        float t = 0.0;
        int counter = 0;
        while (t > tMin)
        {
            rk4Step(ray, wp, dt);
            if (abs(ray.l) > lMax)
            {
                break;
            }

            // compute and store the position coordinates (x,y,z,l) in world space
            if (ray.l > 0)
            {
                glm::vec4 pos(ray.l * sin(ray.theta) * cos(ray.phi),
                              ray.l * sin(ray.theta) * sin(ray.phi), ray.l * cos(ray.theta), ray.l);
                rayPositions[i * numRayPositions + counter] = pos;
            }
            else
            {
                glm::vec4 pos(-ray.l * sin(ray.theta) * cos(ray.phi),
                              -ray.l * sin(ray.theta) * sin(ray.phi), -ray.l * cos(ray.theta),
                              ray.l);
                rayPositions[i * numRayPositions + counter] = pos;
            }

            t += dt;
            counter += 1;
        }

        // add the ray to the ray state
        rays[i] = ray;
        // add the total count
        rayPositionCounts[i] = counter;
    }

    const int AA_SAMPLES = 9;

    // sample the texture color based on the precomputed ray directions
    for (int j = 0; j < outHeight; ++j)
    {
        for (int i = 0; i < outWidth; ++i)
        {
            int idx = j * outWidth + i;

            float r_accum = 0.0f;
            float g_accum = 0.0f;
            float b_accum = 0.0f;

            for (int sampleIdx = 0; sampleIdx < AA_SAMPLES; ++sampleIdx)
            {

                float offsetX = 0.0f;
                float offsetY = 0.0f;

                int gridX = sampleIdx % 3;
                int gridY = sampleIdx / 3;

                offsetX = (gridX - 1) * 0.33f;
                offsetY = (gridY - 1) * 0.33f;

                if (sampleIdx == 1)
                {
                    offsetX = 0.25f;
                }
                else if (sampleIdx == 2)
                {
                    offsetY = 0.25f;
                }
                else if (sampleIdx == 3)
                {
                    offsetX = 0.25f;
                    offsetY = 0.25f;
                }

                // Normalized device coordinates in [-1,1]
                float ndcW = 2.0 * ((i + 0.5 + offsetX) / static_cast<float>(outWidth)) - 1.0;
                float ndcH = 2.0 * ((j + 0.5 + offsetY) / static_cast<float>(outHeight)) - 1.0;

                // ray direction
                glm::vec3 rayDir(1.0, planeWidth * (ndcW), planeHeight * (-ndcH));
                glm::vec2 normalizedYZ = normalize(rayDir.yz());

                // length at the yz-direction
                float lengthYZ = length(rayDir.yz());
                // get the precomputed ray
                int rayInd =
                    static_cast<int>(ceil(static_cast<float>(numRays) * (lengthYZ / lengthDiagonal)));
                rayInd = std::min(rayInd, numRays);
                RayState equatorialRay = rays[rayInd];

                // normalize the angles
                // first, convert the final ray direction into a unit vector
                glm::vec4 finalEuclidean(cos(equatorialRay.phi), sin(equatorialRay.phi), 0.0, 0.0);
                // next, apply the transformation with respect to the tilting of the camera
                glm::vec3 cameraDirection(cos(phi_c), sin(phi_c), 0.0);
                // the angle to rotate about the camera direction
                float angle = atan2(normalizedYZ.y, normalizedYZ.x);
                glm::mat4 rotationMatrix =
                    glm::rotate(glm::mat4(1.0f), angle, glm::normalize(cameraDirection));
                glm::vec4 rotatedEuclidean = rotationMatrix * finalEuclidean;

                float thetaFinal = acos(rotatedEuclidean[2] / length(rotatedEuclidean));
                float phiFinal = atan2(rotatedEuclidean[1], rotatedEuclidean[0]);

                // sample the color
                RGBA color = {0, 0, 0, 255};

                // if it intersects the sphere, display the color on the sphere
                bool intersectsPrimitive = false;
                float closestDistance = std::numeric_limits<float>::max();
                RGBA finalColor = {0, 0, 0, 255};

                for (const Object &currentObject : objects)
                {
                    glm::vec3 objectPos(currentObject.points[0][0], currentObject.points[0][1], currentObject.points[0][2]);
                    SphereData objectData;
                    float cubeHalfSide = 0.0f;
                    float lMinObject, lMaxObject;

                    if (currentObject.type == PrimitiveType::Sphere)
                    {
                        objectData = SphereData(objectPos, currentObject.points[0][3], -length(objectPos));
                        lMinObject = objectData.l - objectData.radius;
                        lMaxObject = objectData.l + objectData.radius;
                    }
                    else
                    {
                        cubeHalfSide = currentObject.side * 0.5f;
                        objectData = SphereData(objectPos, cubeHalfSide * 1.732f, -length(objectPos));
                        lMinObject = objectData.l - cubeHalfSide * 1.732f;
                        lMaxObject = objectData.l + cubeHalfSide * 1.732f;
                    }

                    for (int k = 1; k < rayPositionCounts[rayInd]; k++)
                    {
                        glm::vec4 pos = rayPositions[k + rayInd * numRayPositions];
                        if (pos[3] >= lMinObject && pos[3] <= lMaxObject)
                        {
                            glm::vec3 posEuclidean(pos);
                            posEuclidean = glm::vec3(rotationMatrix * glm::vec4(posEuclidean, 1.0));

                            bool hitObject = false;
                            glm::vec3 normal;
                            glm::vec3 hitPoint;
                            float hitDistance;

                            if (currentObject.type == PrimitiveType::Sphere)
                            {
                                float dist = length(posEuclidean - objectData.center);
                                if (dist - objectData.radius < 0.0)
                                {
                                    hitObject = true;
                                    hitPoint = normalize(posEuclidean - objectData.center) * objectData.radius + objectData.center;
                                    normal = normalize(hitPoint - objectData.center);
                                    hitDistance = length(hitPoint - glm::vec3(rayPositions[k - 1 + rayInd * numRayPositions]));
                                }
                            }
                            else
                            {
                                glm::vec3 prevPos = glm::vec3(rayPositions[k - 1 + rayInd * numRayPositions]);
                                prevPos = glm::vec3(rotationMatrix * glm::vec4(prevPos, 1.0));
                                glm::vec3 rayDirection = normalize(posEuclidean - prevPos);

                                float tNear, tFar;
                                if (intersectCube(prevPos, rayDirection, objectData.center, cubeHalfSide, tNear, tFar))
                                {
                                    hitObject = true;
                                    hitPoint = prevPos + rayDirection * tNear;
                                    normal = getCubeNormal(hitPoint, objectData.center, cubeHalfSide);
                                    hitDistance = tNear;
                                }
                            }

                            if (hitObject && hitDistance < closestDistance)
                            {
                                closestDistance = hitDistance;
                                intersectsPrimitive = true;

                                // sample the texture color
                                glm::vec4 posIntersection;
                                if (pos[3] > 0)
                                {
                                    posIntersection = glm::vec4(hitPoint, 1.0);
                                }
                                else
                                {
                                    posIntersection = glm::vec4(hitPoint, -1.0);
                                }

                                glm::vec3 dirToCamera = normalize(
                                    glm::vec3(rayPositions[k - 1 + rayInd * numRayPositions] - pos));

                                Hit hit = {posIntersection, dirToCamera, objectData};
                                glm::vec3 color_vec =
                                    shadePixel(hit, currentObject.textureFile, BumpMap{nullptr, 0, 0},
                                               lights, currentObject.type, normal) *
                                    255.f;
                                finalColor = RGBA{static_cast<std::uint8_t>(std::min(255.f, color_vec.x)),
                                                  static_cast<std::uint8_t>(std::min(255.f, color_vec.y)),
                                                  static_cast<std::uint8_t>(std::min(255.f, color_vec.z)), 255};

                                break;
                            }
                        }
                    }
                }

                if (intersectsPrimitive)
                {
                    color = finalColor;
                }

                // if not, display the skybox texture
                if (!intersectsPrimitive && abs(equatorialRay.l) > lMax)
                {
                    if (equatorialRay.l > 0.0)
                    {
                        color = sampleCelestial(sphereUpper, thetaFinal, phiFinal);
                    }
                    else
                    {
                        color = sampleCelestial(sphereLower, thetaFinal, phiFinal);
                    }
                }

                r_accum += color.r;
                g_accum += color.g;
                b_accum += color.b;
            }

            RGBA finalPixelColor;
            finalPixelColor.r = static_cast<uint8_t>(r_accum / AA_SAMPLES);
            finalPixelColor.g = static_cast<uint8_t>(g_accum / AA_SAMPLES);
            finalPixelColor.b = static_cast<uint8_t>(b_accum / AA_SAMPLES);
            finalPixelColor.a = 255;

            framebuffer[idx] = finalPixelColor;
        }
    }
}

void render(RGBA *framebuffer,
            int outWidth,
            int outHeight,
            const ImageData &sphereUpper,
            const ImageData &sphereLower,
            const ImageData &primitiveTexture,
            float fovW,
            WormholeParams wp,
            float dt,
            float cameraDistance,
            SphereData sphereData,
            std::vector<SceneLightData> lights)
{
    // defaults if not specified
    float cameraTheta = M_PI / 2.0f;
    float cameraPhi = 0.f;
    Object obj;
    obj.type = PrimitiveType::Sphere;
    obj.radius = sphereData.radius;
    obj.textureFile = primitiveTexture;
    obj.points.push_back({sphereData.center.x, sphereData.center.y, sphereData.center.z, sphereData.radius, 0.0f});
    std::vector<Object> objects = {obj};
    render(framebuffer, outWidth, outHeight, sphereUpper, sphereLower, objects, fovW, wp, dt, cameraDistance, cameraTheta, cameraPhi, lights);
}

// // this function traces the rays
// void render(
//     RGBA *framebuffer,
//     int outWidth,
//     int outHeight,
//     const ImageData &sphereUpper,
//     const ImageData &sphereLower,
//     float fovW,
//     WormholeParams wp,
//     float dt,
//     float cameraDistance)
// {

//     // Camera setup
//     const float l_c     = cameraDistance;
//     const float theta_c = M_PI / 2.0; // equatorial
//     const float phi_c   = 0;

//     // Aspect ratio and vertical FOV
//     float aspect   = static_cast<float>(outHeight) / static_cast<float>(outWidth);

//     // Loop over pixels
//     float planeWidth = tan(fovW * 0.5);
//     float planeHeight = planeWidth * aspect;

//     // Integration constants (-0.01 seems good enough, but very slow...)
//     const float tMin = -200.0;
//     const float lMax = 20.0;

//     for (int j = 0; j < outHeight; ++j) {
//         for (int i = 0; i < outWidth; ++i) {
//             int idx = j * outWidth + i;

//             if (idx % 10000 == 0) {
//                 std::cout << "ind: " << idx << std::endl;
//             }

//             // Normalized device coordinates in [-1,1]
//             float ndcW = 2.0 * ((i + 0.5) / static_cast<float>(outWidth))  - 1.0;
//             float ndcH = 2.0 * ((j + 0.5) / static_cast<float>(outHeight)) - 1.0;

//             // Pinhole camera direction
//             float px = 1.0;
//             float py = planeWidth * (ndcW);
//             float pz = planeHeight * (-ndcH);

//             float normP = sqrt(px*px + py*py + pz*pz);
//             px /= normP;
//             py /= normP;
//             pz /= normP;

//             // Direction of incoming ray is -p
//             float n_l     = -px;
//             float n_phi   = -py;
//             float n_theta = pz;

//             // Initial state
//             RayState ray;
//             ray.l     = l_c;
//             ray.theta = theta_c;
//             ray.phi   = phi_c;

//             float r_c = r_of_l(ray.l, wp);

//             ray.p_l     = n_l;
//             ray.p_theta = r_c * n_theta;
//             ray.b       = r_c * sin(ray.theta) * n_phi;

//             // numerical integration
//             float t = 0.0;
//             while (t > tMin) {
//                 rk4Step(ray, wp, dt);
//                 if (abs(ray.l) > lMax) {
//                     break;
//                 }
//                 t += dt;
//             }

//             RGBA color = {0, 0, 0, 255};

//             if (abs(ray.l) > lMax) {
//                 if (ray.l > 0.0) {
//                     color = sampleCelestial(sphereUpper, ray.theta, ray.phi);
//                 } else {
//                     color = sampleCelestial(sphereLower, ray.theta, ray.phi);
//                 }
//             } else {
//                 color = {0, 0, 0, 255};
//             }

//             framebuffer[idx] = color;
//         }
//     }
// }

bool renderSingleImage(
    QImage outputImage, QString outputPath,
    RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
    const ImageData &sphereLower, const ImageData &primitiveTexture, float fovW,
    WormholeParams wp, float dt, float cameraDistance, SphereData sphereData, std::vector<SceneLightData> lights)
{

    render(framebuffer,
           outWidth,
           outHeight,
           sphereUpper,
           sphereLower,
           primitiveTexture,
           fovW, // in radians
           wp,
           dt,
           cameraDistance,
           sphereData,
           lights);

    // 5. Save output
    bool ok = outputImage.save(outputPath);
    if (!ok)
        ok = outputImage.save(outputPath, "PNG");

    if (!ok)
    {
        std::cerr << "Failed to save output image: "
                  << outputPath.toStdString() << "\n";
        return 1;
    }

    std::cout << "Saved rendered image to "
              << outputPath.toStdString() << "\n";
    return 0;
}

bool renderPath(
    QImage outputImage, QString outputPath,
    RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
    const ImageData &sphereLower, const ImageData &primitiveTexture, float fovW,
    WormholeParams wp, float dt, float cameraDistance, SphereData sphereData, std::vector<SceneLightData> lights,
    std::vector<glm::vec4> &keyframes, int numPhotos)
{

    // create paths
    cameraPath paths;
    paths.buildFromKeyframes(keyframes, numPhotos, pathtype::PATHTYPE_LINEAR);

    Object obj;
    obj.type = PrimitiveType::Sphere;
    obj.radius = sphereData.radius;
    obj.textureFile = primitiveTexture;
    obj.points.push_back({sphereData.center.x, sphereData.center.y, sphereData.center.z, sphereData.radius, 0.0f});
    std::vector<Object> objects = {obj};

    // generate each photo
    int photoNum = 0;
    for (auto cameraPoint : paths.getPoints())
    {
        QString output = outputPath + "/" + QString::number(photoNum) + ".png";
        photoNum += 1;
        float cameraTheta = cameraPoint.theta;
        float cameraDistance = cameraPoint.r;
        float cameraPhi = cameraPoint.phi;
        render(framebuffer, outWidth, outHeight, sphereUpper, sphereLower, objects, fovW, wp, dt, cameraDistance, cameraTheta, cameraPhi, lights);

        // load photo
        bool ok = outputImage.save(output);
        if (!ok)
            ok = outputImage.save(output, "PNG");

        if (!ok)
        {
            std::cerr << "Failed to save output image: "
                      << output.toStdString() << "\n";
            return 1;
        }
    }
    return 0;
}

bool renderFrames(QImage outputImage,
                  FrameData frameData,
                  RGBA *framebuffer,
                  int outWidth,
                  int outHeight,
                  const ImageData &sphereUpper,
                  const ImageData &sphereLower,
                  float fovW,
                  WormholeParams wp,
                  float dt,
                  std::vector<SceneLightData> lights)
{

    float tMin = frameData.tMin;
    float tMax = frameData.tMax;
    float tEnd = frameData.numPhotos;

    for (int i = 0; i < frameData.numPhotos; ++i)
    {
        float u = float(i) / tEnd;
        float t = tMin + u * (tMax - tMin);

        glm::vec3 camCoords = frameData.cameraAtT(t);
        auto objs = frameData.objectsAtT(t);

        render(framebuffer,
               outWidth,
               outHeight,
               sphereUpper,
               sphereLower,
               objs,
               fovW,
               wp,
               dt,
               camCoords.x,
               camCoords.y,
               camCoords.z,
               lights);

        // 3. Decide output path
        QString outFilePath;
        if (frameData.numPhotos == 1)
        {
            // Single frame: frameData.outputPath is a file path, e.g. "output/saturn.png"
            outFilePath = frameData.outputPath;
        }
        else
        {
            // Multiple frames: frameData.outputPath is a folder, e.g. "output/saturn_wide_path"
            QDir outDir(frameData.outputPath);

            // Use zero-padded filenames like 000.png, 001.png, ...
            QString fileName = QString("%1.png").arg(i, 3, 10, QChar('0'));
            outFilePath = outDir.filePath(fileName);
        }

        // 4. Save
        if (!outputImage.save(outFilePath))
        {
            qWarning() << "Failed to save frame to" << outFilePath;
            return false;
        }
    }
    return true;
}
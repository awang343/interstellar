#include "renderer.h"
#include "lighting.h"

#include <algorithm>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>
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

    int x = static_cast<int>(u * (img.width - 1));
    int y = static_cast<int>(v * (img.height - 1));

    x = std::max(0, std::min(img.width - 1, x));
    y = std::max(0, std::min(img.height - 1, y));

    return img.pixels[y * img.width + x];
}

// this function traces the rays using the equatorial symmetry method
// see https://www.youtube.com/watch?v=PVO8nvb1o2w for the explanation of this method
void render(RGBA *framebuffer, int outWidth, int outHeight, const ImageData &sphereUpper,
            const ImageData &sphereLower, const ImageData &primitiveTexture, float fovW,
            WormholeParams wp, float dt, float cameraDistance, SphereData sphereData)
{

    // Camera setup
    const float l_c = cameraDistance;
    const float theta_c = M_PI / 2.0; // equatorial
    const float phi_c = 0;            // will eventually make changes to this

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

    // the l interval that contains the sphere
    float lMinSphere = sphereData.l - sphereData.radius;
    float lMaxSphere = sphereData.l + sphereData.radius;

    // sample the texture color based on the precomputed ray directions
    for (int j = 0; j < outHeight; ++j)
    {
        for (int i = 0; i < outWidth; ++i)
        {
            int idx = j * outWidth + i;

            // Normalized device coordinates in [-1,1]
            float ndcW = 2.0 * ((i + 0.5) / static_cast<float>(outWidth)) - 1.0;
            float ndcH = 2.0 * ((j + 0.5) / static_cast<float>(outHeight)) - 1.0;

            // ray direction
            glm::vec3 rayDir(1.0, planeWidth * (ndcW), planeHeight * (-ndcH));
            glm::vec2 normalizedYZ = normalize(rayDir.yz());

            // length at the yz-direction
            float lengthYZ = length(rayDir.yz());
            // get the precomputed ray
            int rayInd =
                static_cast<int>(ceil(static_cast<float>(numRays) * (lengthYZ / lengthDiagonal)));
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
            bool intersectsSphere = false;

            for (int k = 1; k < rayPositionCounts[rayInd]; k++)
            {
                if (intersectsSphere)
                {
                    break;
                }
                glm::vec4 pos = rayPositions[k + rayInd * numRayPositions];
                if (pos[3] >= lMinSphere && pos[3] <= lMaxSphere)
                {
                    glm::vec3 posEuclidean(pos);
                    posEuclidean = glm::vec3(rotationMatrix * glm::vec4(posEuclidean, 1.0));
                    if (length(posEuclidean - sphereData.center) - sphereData.radius < 0.0)
                    {
                        intersectsSphere = true;

                        // sample the texture color
                        glm::vec4 posIntersection;
                        if (pos[3] > 0)
                        {
                            posIntersection = glm::vec4(
                                normalize(posEuclidean - sphereData.center) * sphereData.radius +
                                    sphereData.center,
                                1.0);
                        }
                        else
                        {
                            posIntersection = glm::vec4(
                                normalize(posEuclidean - sphereData.center) * sphereData.radius +
                                    sphereData.center,
                                -1.0);
                        }

                        glm::vec3 dirToCamera = normalize(
                            glm::vec3(rayPositions[k - 1 + rayInd * numRayPositions] - pos));

                        Hit hit = {posIntersection, dirToCamera, sphereData};
                        glm::vec3 color_vec =
                            shadePixel(hit, primitiveTexture, BumpMap{nullptr, 0, 0},
                                       std::vector<SceneLightData>{}) *
                            255.f;
                        RGBA finalColor = {static_cast<std::uint8_t>(std::min(255.f, color_vec.x)),
                                           static_cast<std::uint8_t>(std::min(255.f, color_vec.y)),
                                           static_cast<std::uint8_t>(std::min(255.f, color_vec.z)),
                                           255};
                    }
                }
            }

            // if not, display the skybox texture
            if (!intersectsSphere && abs(equatorialRay.l) > lMax)
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

            framebuffer[idx] = color;
        }
    }
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

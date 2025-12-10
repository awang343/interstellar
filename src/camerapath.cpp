#include "camerapath.h"

cameraPath::cameraPath() {}

float cameraPath::lerpAngle(float a0, float a1, float t) {
    // bring difference into [-pi, pi]
    float d = std::fmod(a1 - a0 + M_PI, 2.0f * M_PI);
    if (d < 0.0f)
        d += 2.0f * M_PI;
    d -= M_PI;

    return a0 + t * d;
}

float cameraPath::lerp(float p0, float p1, float t) {
    return p0 + t*(p1-p0);
}

void cameraPath::addPath(float start_r, float start_theta, float start_phi,
                         float end_r, float end_theta, float end_phi,
                         pathtype path, int numPoints)
{
    if (numPoints < 2) numPoints = 2;

    switch (path) {
    case PATHTYPE_LINEAR:
    default:
        for (int i = 0; i < numPoints; ++i) {
            float t = static_cast<float>(i) / (numPoints - 1);

            cameradata c;
            c.r     = lerp(start_r, end_r, t);
            c.theta = lerpAngle(start_theta, end_theta, t);
            c.phi   = lerpAngle(start_phi,   end_phi,   t);

            points.push_back(c);
        }
        break;

    case PATHTYPE_QUADRATIC:
        // placeholder: simple ease-in-out quadratic; adjust as needed
        for (int i = 0; i < numPoints; ++i) {
            float t_lin = static_cast<float>(i) / (numPoints - 1);
            float t = t_lin * t_lin; // quadratic easing (not symmetric)

            cameradata c;
            c.r     = start_r   + t * (end_r   - start_r);
            c.theta = lerpAngle(start_theta, end_theta, t);
            c.phi   = lerpAngle(start_phi,   end_phi,   t);

            points.push_back(c);
        }
        break;
    }
}


void cameraPath::buildFromKeyframes(const std::vector<glm::vec4> &keyframes,
                                    int numPhotos,
                                    pathtype path)
{
    points.clear();

    if (keyframes.size() < 2 || numPhotos <= 0) {
        return;
    }

    // Assume keyframes are sorted by t (the .w component)
    float tStart = keyframes.front().w;
    float tEnd   = keyframes.back().w;

    // Special case: only one photo → just take first keyframe
    if (numPhotos == 1) {
        cameradata c;
        c.r     = keyframes.front().x;
        c.theta = keyframes.front().y;
        c.phi   = keyframes.front().z;
        points.push_back(c);
        return;
    }

    // Global time spacing between photos
    float totalT = tEnd - tStart;
    float dt = totalT / float(numPhotos - 1);

    int seg = 0; // index of current segment [seg, seg+1]

    for (int i = 0; i < numPhotos; ++i) {
        float t = tStart + dt * float(i);

        // Clamp last sample exactly to tEnd
        if (i == numPhotos - 1) {
            t = tEnd;
        }

        // Move seg until t is within [keyframes[seg].w, keyframes[seg+1].w]
        while (seg + 1 < int(keyframes.size()) - 1 &&
               t > keyframes[seg + 1].w) {
            ++seg;
        }

        const glm::vec4 &p0 = keyframes[seg];
        const glm::vec4 &p1 = keyframes[seg + 1];

        float t0 = p0.w;
        float t1 = p1.w;
        float segLen = t1 - t0;

        float u = 0.0f;
        if (segLen > 0.0f) {
            u = (t - t0) / segLen;
            u = glm::clamp(u, 0.0f, 1.0f);
        }

        // Apply easing if quadratic
        switch (path) {
        case PATHTYPE_QUADRATIC:
            u = u * u; // simple ease-in example
            break;
        case PATHTYPE_LINEAR:
        default:
            break;
        }

        // Interpolate spherical parameters
        float r     = lerp(p0.x, p1.x, u);
        float theta = lerpAngle(p0.y, p1.y, u);
        float phi   = lerpAngle(p0.z, p1.z, u);

        cameradata c;
        c.r     = r;
        c.theta = theta;
        c.phi   = phi;
        points.push_back(c);
    }
}

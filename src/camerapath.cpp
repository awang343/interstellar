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

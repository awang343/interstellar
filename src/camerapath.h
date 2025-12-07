#ifndef CAMERAPATH_H
#define CAMERAPATH_H

#include <math.h>
#include <vector>

struct cameradata {
    float r     = -5.0f;        // camera distance
    float theta = M_PI / 2.0f;  // equatorial angle
    float phi   = 0.0f;         // other angle
};

enum pathtype {
    PATHTYPE_LINEAR,
    PATHTYPE_QUADRATIC
};

class cameraPath
{
public:
    cameraPath();
    std::vector<cameradata> getPoints() const {return points;}

    void addPath(float start_r, float start_theta, float start_phi,
                 float end_r, float end_theta, float end_phi,
                 pathtype path, int numPoints);


private:
    std::vector<cameradata> points;
    float lerpAngle(float a0, float a1, float t);
    float lerp(float p0, float p1, float t);
};

#endif // CAMERAPATH_H

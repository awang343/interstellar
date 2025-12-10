#include "framedata.h"
#include <qimage.h>
#include <algorithm>
#include <cmath>

FrameData::FrameData() {}

namespace
{

    // Angle interpolation with wrap-around at 2π
    inline float lerpAngle(float a0, float a1, float t)
    {
        float d = std::fmod(a1 - a0 + float(M_PI), 2.0f * float(M_PI));
        if (d < 0.0f)
            d += 2.0f * float(M_PI);
        d -= float(M_PI);
        return a0 + t * d;
    }

    // Find segment index such that t is in [times[i], times[i+1]]
    template <typename GetTime>
    int findSegment(float t, int count, GetTime getTime)
    {
        int seg = 0;
        while (seg + 1 < count - 1 && t > getTime(seg + 1))
        {
            ++seg;
        }
        return seg;
    }

} // namespace

void FrameData::loadDefaultObject()
{
    // remove all objects
    objects.clear();

    // create object
    Object o;
    o.type = PrimitiveType::Sphere;
    // set object's only position to be at (0,2,0)
    o.points.push_back({0.0, 2.0, 0.0, 0.5, 0});
    o.radius = 0.5;
    // set texture file, default to bricks
    loadImageToStruct("texture/bricks.jpg", o.textureFile);

    // set object's list to be only this object
    objects.push_back(o);
}

void FrameData::loadDefaultCamera(float cameraDistance)
{
    cameraPoints.clear();
    cameraPoints.push_back(glm::vec4(cameraDistance, M_PI / 2.0f, 0.f, 0.f));
}

glm::vec3 FrameData::cameraAtT(float t)
{
    if (cameraPoints.empty())
    {
        // Nothing specified: fall back to something sensible
        // (you can change this to your default camera)
        return glm::vec3(0.f, 0.f, 0.f);
    }

    // Assume cameraPoints sorted by .w (time)
    const glm::vec4 &first = cameraPoints.front();
    const glm::vec4 &last = cameraPoints.back();

    if (t <= first.w)
    {
        return glm::vec3(first.x, first.y, first.z);
    }
    if (t >= last.w)
    {
        return glm::vec3(last.x, last.y, last.z);
    }

    // Find segment containing t
    int n = static_cast<int>(cameraPoints.size());
    int seg = findSegment(t, n, [&](int i)
                          { return cameraPoints[i].w; });

    const glm::vec4 &p0 = cameraPoints[seg];
    const glm::vec4 &p1 = cameraPoints[seg + 1];

    float t0 = p0.w;
    float t1 = p1.w;
    float u = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0f;
    u = std::clamp(u, 0.0f, 1.0f);

    float r = std::lerp(p0.x, p1.x, u);

    float theta = lerpAngle(p0.y, p1.y, u);
    float phi = lerpAngle(p0.z, p1.z, u);

    // You’re using spherical coords elsewhere, so return (r, theta, phi)
    return glm::vec3(r, theta, phi);
}

std::vector<Object> FrameData::objectsAtT(float t)
{
    std::vector<Object> result;
    result.reserve(objects.size());

    for (const Object &o : objects)
    {
        if (o.points.empty())
        {
            continue;
        }

        const auto &pts = o.points;
        int m = static_cast<int>(pts.size());

        float tMin = (pts[0].size() >= 5) ? pts[0][4] : 0.f;
        float tMax = (pts[m - 1].size() >= 5) ? pts[m - 1][4] : tMin;

        std::vector<float> interp;

        if (t <= tMin)
        {
            interp = pts.front();
        }
        else if (t >= tMax)
        {
            interp = pts.back();
        }
        else
        {

            // Find segment [k, k+1]
            int seg = 0;
            while (seg + 1 < m - 1 && t > pts[seg + 1][4])
            {
                ++seg;
            }

            const auto &p0 = pts[seg];
            const auto &p1 = pts[seg + 1];

            if (p0.size() != p1.size() || p0.size() < 5)
            {
                continue;
            }

            float t0 = p0[4];
            float t1 = p1[4];
            float u = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0f;
            if (u < 0.f)
                u = 0.f;
            if (u > 1.f)
                u = 1.f;

            interp.resize(p0.size());
            for (size_t j = 0; j < p0.size(); ++j)
            {
                if (j == 4)
                {
                    // time slot
                    interp[j] = t;
                }
                else
                {
                    interp[j] = std::lerp(p0[j], p1[j], u);
                }
            }
        }

        Object out = o;
        out.points.clear();
        out.points.push_back(std::move(interp));

        result.push_back(std::move(out));
    }

    return result;
}

void FrameData::setTimeRange()
{
    float currentMin = cameraPoints.front().w;
    float currentMax = cameraPoints.back().w;

    for (auto const &obj : objects)
    {
        currentMin = fmin(currentMin, obj.points.front()[4]);
        currentMax = fmax(currentMax, obj.points.back()[4]);
    }

    tMin = currentMin;
    tMax = currentMax;
}

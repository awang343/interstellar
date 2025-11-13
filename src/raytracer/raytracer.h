#pragma once

#include "utils/ini_utils.h"
#include "utils/rgba.h"
#include "utils/sceneparser.h"

#include <glm/glm.hpp>

#define RAY_TRACE_MAX_DEPTH 4
#define RAY_TRACE_DEFAULT_SPP 64

// A forward declaration for the RaytraceScene class

class RayTraceScene;

// A class representing a ray-tracer

struct Hit
{
    float t;
    glm::vec3 point;
    glm::vec3 normal;
    const RenderShapeData *shape;
};

class RayTracer
{
  public:
    struct Config
    {
        bool enableShadow = false;
        bool enableReflection = false;
        bool enableRefraction = false;
        bool enableTextureMap = false;
        TextureFilterType textureFilterType = TextureFilterType::Nearest;
        bool enableParallelism = false;
        bool enableSuperSample = false;
        bool enableAcceleration = false;
        bool enableDepthOfField = false;
        int maxRecursiveDepth = RAY_TRACE_MAX_DEPTH;
        int samplesPerPixel = RAY_TRACE_DEFAULT_SPP;
        SuperSamplerPattern superSamplerPattern = SuperSamplerPattern::Grid;
        bool onlyRenderNormals = false;
        bool enableMipMapping = false;
    };

  public:
    RayTracer(Config config);

    // Renders the scene synchronously.
    // The ray-tracer will render the scene and fill imageData in-place.
    // @param imageData The pointer to the imageData to be filled.
    // @param scene The scene to be rendered.
    void sampler(std::vector<std::pair<float, float>> &samples);

    void render(RGBA *imageData, RayTraceScene &scene);
    void render_pixel(RGBA *imageData, RayTraceScene &scene, const int &i, const int &j, unsigned int &last_ray_id,
                      std::vector<unsigned int> &last_visited, std::vector<std::pair<float, float>> &samples);
    glm::vec3 shadePixel(const RayTraceScene &scene, const Hit &hit, const glm::vec3 &world_camera,
                         const glm::vec3 &world_dir, unsigned int &last_ray_id, std::vector<unsigned int> &last_visited,
                         int depth = 0, float passthrough_scale_factor = -1.f);

  private:
    const Config m_config;

    int scene_width;
    int scene_height;
};

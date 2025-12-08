#include "lighting.h"
#include "src/utils/rgba.h"
#include <glm/glm.hpp>
#include <vector>

struct BumpMap
{
    glm::vec2 *gradients;
    int width;
    int height;
};

struct ImageData
{
    int width = 0;
    int height = 0;
    std::vector<RGBA> pixels; // size = width * height
};

struct uv
{
    float u;
    float v;
};

enum class FilterType
{
    Nearest = 0,
    Bilinear = 1,
    Trilinear = 2,
};

uv get_uv(const Hit &hit);
glm::vec3 get_texture(const ImageData &texture, const Hit &hit, const FilterType filter_type, const uv &uv_map);
glm::vec3 get_bump_normal(const BumpMap &bump_map, const Hit &hit, const FilterType filter_type, const uv &uv_map,
                          float bump_scale);

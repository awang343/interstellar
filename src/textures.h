#include "src/utils/scenefileparser.h"
#include "lighting.h"
#include <glm/glm.hpp>

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
glm::vec3 get_texture(const ImageData &texture, const FilterType filter_type, const uv &uv_map);
glm::vec3 get_bump_normal(const BumpMap &bump_map, const FilterType filter_type, const uv &uv_map,
                          float bump_scale, glm::vec3 normal);

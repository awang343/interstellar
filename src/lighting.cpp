#include "lighting.h"
#include "textures.h"
#include <iostream>

static glm::vec3 calc_light_vec(const glm::vec3 &point, const SceneLightData &light)
{
    // Point is the intersection point
    switch (light.type)
    {
    case LightType::LIGHT_DIRECTIONAL:
        return -light.dir;
        break;
    case LightType::LIGHT_POINT:
        return glm::vec3(light.pos) - point;
        break;
    case LightType::LIGHT_SPOT:
        return glm::vec3(light.pos) - point;
        break;
    default:
        exit(1);
        break;
    }
}

static inline float calc_falloff(const float &inner, const float &outer, const float &angle)
{
    const float interm = (angle - inner) / (outer - inner);
    return -2 * interm * interm * interm + 3 * interm * interm;
}

glm::vec3 shadePixel(const Hit &hit, const ImageData &texture, const BumpMap &bump_map,
                     const std::vector<SceneLightData> &lights)
{
    static const float ka = 0.f; // HARDCODED
    static const float kd = 1.f; // HARDCODED
    // static const float ks = 1.f; // HARDCODED
    const float shininess = 0.5f;

    const glm::vec3 obj_A = glm::vec3(0.f, 0.5f, 0.f) * ka; // HARDCODED
    glm::vec3 obj_D = glm::vec3(1.f, 1.f, 1.f) * kd;        // HARDCODED
    // const glm::vec3 obj_S = glm::vec3(1.f, 1.f, 1.f) * ks; // HARDCODED

    const glm::vec3 obj_point = (glm::vec3(hit.point) - hit.sphere.center) / hit.sphere.radius / 2.f; // point in object space
    const glm::vec3 obj_normal = glm::normalize(obj_point);
    const float obj_side = hit.point[3];
    // const glm::mat3 transform = glm::inverse(glm::transpose(glm::mat3(hit.shape->ctm)));

    const uv uv_map = get_uv(obj_point);

    const bool enableBumpMap = false; // HARDCODED
    const glm::vec3 bump_normal =
        enableBumpMap && bump_map.width > 0
            ? get_bump_normal(bump_map, FilterType::Bilinear, uv_map, 0.01, obj_normal)
            : obj_normal;
    const glm::vec3 N = glm::normalize(/*transform **/ obj_normal);
    // const glm::vec3 V = glm::normalize(world_camera - glm::vec3(intersect));

    // Blend obj_D with textures

    const bool enableTextureMap = false;                   // HARDCODED
    const FilterType textureFilter = FilterType::Nearest; // HARDCODED


    if (enableTextureMap && texture.width > 0)
    {
        const float blend = 1.0; // HARDCODED
        const glm::vec3 texture_color = get_texture(texture, textureFilter, uv_map);
        obj_D = obj_D * (1.f - blend) + texture_color * blend;
    }

    glm::vec3 pixel_color(obj_A);

    for (const SceneLightData &light : lights)
    {
        // std::cout << light.pos[3] << " " << obj_side << std::endl;
        if (light.pos[3] != obj_side)
        {
            continue;
        }

        glm::vec3 phong(0.f);

        const glm::vec3 to_light = calc_light_vec(glm::vec3(hit.point), light);
        const glm::vec3 L_i = glm::normalize(to_light); // Vector from point to light
        const glm::vec3 R = glm::normalize(2.f * glm::dot(L_i, N) * N - L_i);

        const float D_intensity = std::max(glm::dot(N, L_i), 0.f);
        // const float S_intensity = std::pow(std::max(glm::dot(R, V), 0.f), shininess);

        const glm::vec3 D_color = glm::vec3(light.color) * obj_D;
        // const glm::vec3 S_color = glm::vec3(light.color) * obj_S;

        phong += D_color * D_intensity;
        // phong += S_color * S_intensity;

        if (glm::length(phong) < 1e-5f)
        {
            continue;
        }

        float light_intensity = 1.f;

        // Spotlight falloff
        if (light.type == LightType::LIGHT_SPOT)
        {
            const float angle = std::acos(glm::dot(-L_i, glm::vec3(light.dir)));
            const float inner = light.angle - light.penumbra;
            const float outer = light.angle;

            if (angle > outer)
            {
                light_intensity *= 0.f;
            }
            else if (angle > inner)
            {
                light_intensity *= 1.f - calc_falloff(inner, outer, angle);
            }
        }

        bool directional = light.type == LightType::LIGHT_DIRECTIONAL;

        // Attenuation effect
        // if (!directional && light_intensity > 1e-6f)
        // {
        //     const float dist = glm::distance(glm::vec3(obj_point), glm::vec3(light.pos));
        //     const float mult = 1.f / glm::dot(light.function, glm::vec3(1.f, dist, dist * dist));
        //     light_intensity *= std::min(1.f, mult);
        // }

        pixel_color += phong * light_intensity;
    }

    return glm::clamp(pixel_color, 0.f, 1.f);
}

glm::vec3 shadePixel(const Hit &hit, const ImageData &texture, const BumpMap &bump_map,
                     const std::vector<SceneLightData> &lights, PrimitiveType objectType, const glm::vec3 &surfaceNormal)
{
    static const float ka = 1.f; // HARDCODED
    static const float kd = 1.f; // HARDCODED
    // static const float ks = 1.f; // HARDCODED
    const float shininess = 0.5f;

    const glm::vec3 obj_A = glm::vec3(0.f, 0.5f, 0.f) * ka; // HARDCODED
    glm::vec3 obj_D = glm::vec3(1.f, 1.f, 1.f) * kd;        // HARDCODED
    // const glm::vec3 obj_S = glm::vec3(1.f, 1.f, 1.f) * ks; // HARDCODED

    const glm::vec3 obj_point = glm::vec3(hit.point) - hit.sphere.center; // point in object space
    const float obj_side = hit.point[3];
    // const glm::mat3 transform = glm::inverse(glm::transpose(glm::mat3(hit.shape->ctm)));

    uv uv_map;
    glm::vec3 obj_normal;

    if (objectType == PrimitiveType::Sphere)
    {
        const glm::vec3 normalized_point = glm::normalize(obj_point);
        uv_map = get_uv(normalized_point);

        const bool enableBumpMap = false; // HARDCODED
        obj_normal = enableBumpMap && bump_map.width > 0
                         ? get_bump_normal(bump_map, FilterType::Bilinear, uv_map, 0.01, normalized_point)
                         : normalized_point;
    }
    else
    {
        uv_map = get_cube_uv(obj_point, surfaceNormal);
        obj_normal = surfaceNormal;
    }

    const glm::vec3 N = glm::normalize(/*transform **/ obj_normal);
    // const glm::vec3 V = glm::normalize(world_camera - glm::vec3(intersect));

    // Blend obj_D with textures
    const bool enableTextureMap = true;                   // HARDCODED
    const FilterType textureFilter = FilterType::Bilinear; // HARDCODED

    if (enableTextureMap && texture.width > 0)
    {
        const float blend = 1.0; // HARDCODED
        const glm::vec3 texture_color = get_texture(texture, textureFilter, uv_map);
        obj_D = obj_D * (1.f - blend) + texture_color * blend;
    }

    glm::vec3 pixel_color(obj_A);

    for (const SceneLightData &light : lights)
    {
        // std::cout << light.pos[3] << " " << obj_side << std::endl;
        if (light.pos[3] != obj_side)
        {
            continue;
        }

        glm::vec3 phong(0.f);

        const glm::vec3 to_light = calc_light_vec(obj_point, light);
        const glm::vec3 L_i = glm::normalize(to_light); // Vector from point to light
        const glm::vec3 R = glm::normalize(2.f * glm::dot(L_i, N) * N - L_i);

        const float D_intensity = std::max(glm::dot(N, L_i), 0.f);
        // const float S_intensity = std::pow(std::max(glm::dot(R, V), 0.f), shininess);

        const glm::vec3 D_color = glm::vec3(light.color) * obj_D;
        // const glm::vec3 S_color = glm::vec3(light.color) * obj_S;

        phong += D_color * D_intensity;
        // phong += S_color * S_intensity;

        if (glm::length(phong) < 1e-5f)
        {
            continue;
        }

        float light_intensity = 1.f;

        // Spotlight falloff
        if (light.type == LightType::LIGHT_SPOT)
        {
            const float angle = std::acos(glm::dot(-L_i, glm::vec3(light.dir)));
            const float inner = light.angle - light.penumbra;
            const float outer = light.angle;

            if (angle > outer)
            {
                light_intensity *= 0.f;
            }
            else if (angle > inner)
            {
                light_intensity *= 1.f - calc_falloff(inner, outer, angle);
            }
        }

        bool directional = light.type == LightType::LIGHT_DIRECTIONAL;

        // Attenuation effect
        if (!directional && light_intensity > 1e-6f)
        {
            const float dist = glm::distance(glm::vec3(obj_point), glm::vec3(light.pos));
            const float mult = 1.f / glm::dot(light.function, glm::vec3(1.f, dist, dist * dist));
            light_intensity *= std::min(1.f, mult);
        }

        pixel_color += phong * light_intensity;
    }

    return glm::clamp(pixel_color, 0.f, 1.f);
}

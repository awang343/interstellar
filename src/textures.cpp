#include "textures.h"
#include "glm/fwd.hpp"
#include "utils/structs.h"
#include <iostream>

uv get_uv(const glm::vec3 &obj_point)
{
    const float theta = atan2(obj_point.z, obj_point.x);
    const float u = theta < 0 ? -(theta) / (2.f * M_PI) : 1.f - (theta) / (2.f * M_PI);
    const float v = asin(std::clamp(obj_point.y, -0.5f, 0.5f) / 0.5f) / M_PI + .5f;
    return uv{u, v};
}

uv get_cube_uv(const glm::vec3 &obj_point, const glm::vec3 &normal)
{
    float u, v;
    glm::vec3 absNormal = glm::abs(normal);

    if (absNormal.x > absNormal.y && absNormal.x > absNormal.z)
    {
        u = (obj_point.z / 0.5f + 1.0f) * 0.5f;
        v = (obj_point.y / 0.5f + 1.0f) * 0.5f;
    }
    else if (absNormal.y > absNormal.z)
    {
        u = (obj_point.x / 0.5f + 1.0f) * 0.5f;
        v = (obj_point.z / 0.5f + 1.0f) * 0.5f;
    }
    else
    {
        u = (obj_point.x / 0.5f + 1.0f) * 0.5f;
        v = (obj_point.y / 0.5f + 1.0f) * 0.5f;
    }

    return uv{u, v};
}

static std::pair<glm::vec3, glm::vec3> get_dpduv(const uv &uv_coords)
{
    const float r = 0.5f;

    const float theta = (1.f - uv_coords.u) * 2.f * M_PI;
    const float phi = (uv_coords.v - 0.5f) * M_PI;

    const float cosTheta = cos(theta);
    const float sinTheta = sin(theta);
    const float cosPhi = cos(phi);
    const float sinPhi = sin(phi);

    const float dtheta_du = -2.f * M_PI;
    glm::vec3 dpdtheta(-r * cosPhi * sinTheta, 0.f, r * cosPhi * cosTheta);
    glm::vec3 dpdu = dtheta_du * dpdtheta;

    const float dphi_dv = M_PI;
    glm::vec3 dpdphi(-r * sinPhi * cosTheta, r * cosPhi, -r * sinPhi * sinTheta);
    glm::vec3 dpdv = dphi_dv * dpdphi;

    return {dpdu, dpdv};
}

static std::pair<glm::vec3, glm::vec3> get_cube_dpduv(const uv &uv_coords, const glm::vec3 &normal)
{
    glm::vec3 dpdu, dpdv;
    float scale = 1.0f; // Assuming cube of size 1 unit

    glm::vec3 absNormal = glm::abs(normal);

    if (absNormal.x > absNormal.y && absNormal.x > absNormal.z)
    {
        dpdu = glm::vec3(0.0f, 0.0f, scale);
        dpdv = glm::vec3(0.0f, scale, 0.0f);
    }
    else if (absNormal.y > absNormal.z)
    {
        dpdu = glm::vec3(scale, 0.0f, 0.0f);
        dpdv = glm::vec3(0.0f, 0.0f, scale);
    }
    else
    {
        dpdu = glm::vec3(scale, 0.0f, 0.0f);
        dpdv = glm::vec3(0.0f, scale, 0.0f);
    }

    return {dpdu, dpdv};
}

static glm::vec2 bilinear_bump(const BumpMap &bmp, const uv &uv_map, const float m, const float n)
{
    const float w = bmp.width;
    const float h = bmp.height;

    const float x_left = uv_map.u * m * w - 0.5f;
    const float x_right = x_left + 1.f;
    const float y_top = (1.f - uv_map.v) * n * h - 0.5f;
    const float y_bottom = y_top + 1.f;

    const float c_left = std::fmod(std::floorf(x_left) + w, w);
    const float c_right = std::fmod(std::floorf(x_right) + w, w);
    const float r_top = std::fmod(std::floorf(y_top) + h, h);
    const float r_bottom = std::fmod(std::floorf(y_bottom) + h, h);

    const float alpha_x = x_left - std::floorf(x_left);
    const float alpha_y = y_top - std::floorf(y_top);

    const glm::vec2 tl = bmp.gradients[static_cast<int>(r_top * w + c_left)];
    const glm::vec2 tr = bmp.gradients[static_cast<int>(r_top * w + c_right)];
    const glm::vec2 bl = bmp.gradients[static_cast<int>(r_bottom * w + c_left)];
    const glm::vec2 br = bmp.gradients[static_cast<int>(r_bottom * w + c_right)];

    const glm::vec2 I_top = (1.f - alpha_x) * tl + alpha_x * tr;
    const glm::vec2 I_bottom = (1.f - alpha_x) * bl + alpha_x * br;

    glm::vec3 normal;
    return ((1.f - alpha_y) * I_top + alpha_y * I_bottom);
}

static glm::vec3 bilinear_texture(const ImageData &image, const uv &uv_map, const float m,
                                  const float n)
{
    const float w = image.width;
    const float h = image.height;
    const std::vector<RGBA> color_data = image.pixels;

    const float x_left = uv_map.u * m * w - 0.5f;
    const float x_right = x_left + 1.f;
    const float y_top = (1.f - uv_map.v) * n * h - 0.5f;
    const float y_bottom = y_top + 1.f;

    const float c_left = std::fmod(std::floorf(x_left) + w, w);
    const float c_right = std::fmod(std::floorf(x_right) + w, w);
    const float r_top = std::fmod(std::floorf(y_top) + h, h);
    const float r_bottom = std::fmod(std::floorf(y_bottom) + h, h);

    const float alpha_x = x_left - std::floorf(x_left);
    const float alpha_y = y_top - std::floorf(y_top);

    const RGBA tl = color_data[static_cast<int>(r_top * w + c_left)];
    const RGBA tr = color_data[static_cast<int>(r_top * w + c_right)];
    const RGBA bl = color_data[static_cast<int>(r_bottom * w + c_left)];
    const RGBA br = color_data[static_cast<int>(r_bottom * w + c_right)];

    const glm::vec3 I_top = (1.f - alpha_x) * glm::vec3(tl.r / 255.f, tl.g / 255.f, tl.b / 255.f) +
                            alpha_x * glm::vec3(tr.r / 255.f, tr.g / 255.f, tr.b / 255.f);
    const glm::vec3 I_bottom =
        (1.f - alpha_x) * glm::vec3(bl.r / 255.f, bl.g / 255.f, bl.b / 255.f) +
        alpha_x * glm::vec3(br.r / 255.f, br.g / 255.f, br.b / 255.f);

    return ((1.f - alpha_y) * I_top + alpha_y * I_bottom);
}

glm::vec3 get_texture(const ImageData &texture, const FilterType filter_type, const uv &uv_map)
{
    const float m = 1; // HARDCODED
    const float n = 1; // HARDCODED

    switch (filter_type)
    {
    case FilterType::Nearest:
    {
        const float w = texture.width;
        const float h = texture.height;
        const std::vector<RGBA> color_data = texture.pixels;

        const int c = std::fmod(std::floorf(uv_map.u * m * w), w);
        const int r = std::fmod(std::floorf((1 - uv_map.v) * n * h), h);

        const RGBA color = color_data[static_cast<int>(r * w + c)];
        return glm::vec3(color.r, color.g, color.b) / 255.f;
    }
    case FilterType::Bilinear:
    {
        return bilinear_texture(texture, uv_map, m, n);
    }
    case FilterType::Trilinear:
    {
        std::cerr << "Trilinear filtering not supported." << std::endl;
        exit(1);
        break;
    }
    }
}

glm::vec3 get_bump_normal(const BumpMap &bump_map, const FilterType filter_type, const uv &uv_map,
                          float bump_scale, const glm::vec3 &normal, const PrimitiveType objectType)
{
    const float w = bump_map.width;
    const float h = bump_map.height;
    const float m = 1; // HARDCODED
    const float n = 1; // HARDCODED

    const auto [dpdu, dpdv] =
        objectType == PrimitiveType::Sphere ? get_dpduv(uv_map) : get_cube_dpduv(uv_map, normal);

    glm::vec2 g;

    switch (filter_type)
    {
    case FilterType::Nearest:
    {
        const int c = std::fmod(std::floorf(uv_map.u * m * w) + w, w);
        const int r = std::fmod(std::floorf((1 - uv_map.v) * n * h) + h, h);

        g = bump_map.gradients[static_cast<int>(r * w + c)];
        break;
    }
    case FilterType::Bilinear:
    {
        g = bilinear_bump(bump_map, uv_map, m, n);
        break;
    }
    case FilterType::Trilinear:
    {
        std::cerr << "Trilinear filtering not supported." << std::endl;
        exit(1);
        break;
    }
    }

    const float dHdu = g.x * m * w;
    const float dHdv = g.y * n * h;

    glm::vec3 dNdu = bump_scale * dHdu * glm::cross(normal, dpdv);
    glm::vec3 dNdv = bump_scale * dHdv * glm::cross(dpdu, normal);

    glm::vec3 perturbed_normal = normal + dNdu + dNdv;
    return glm::normalize(perturbed_normal);
}

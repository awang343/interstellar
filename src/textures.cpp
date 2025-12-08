#include "textures.h"
#include "glm/fwd.hpp"
#include <iostream>

uv get_uv(const Hit &hit)
{
    const float theta = atan2(hit.point.z, hit.point.x);
    const float u = theta < 0 ? -(theta) / (2.f * M_PI) : 1.f - (theta) / (2.f * M_PI);
    const float v = asin(hit.point.y / 0.5f) / M_PI + .5f;
    return uv{u, v};
}

static std::pair<glm::vec3, glm::vec3> get_dpduv(const Hit &hit, const uv &uv_coords)
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

static glm::vec2 bilinear_bump(const BumpMap &bmp, const uv &uv_map, const float m, const float n)
{
    const float w = bmp.width;
    const float h = bmp.height;
    const glm::vec2 *bmp_data = bmp.gradients;

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

    const glm::vec2 tl = bmp_data[static_cast<int>(r_top * w + c_left)];
    const glm::vec2 tr = bmp_data[static_cast<int>(r_top * w + c_right)];
    const glm::vec2 bl = bmp_data[static_cast<int>(r_bottom * w + c_left)];
    const glm::vec2 br = bmp_data[static_cast<int>(r_bottom * w + c_right)];

    const glm::vec2 I_top = (1.f - alpha_x) * tl + alpha_x * tr;
    const glm::vec2 I_bottom = (1.f - alpha_x) * bl + alpha_x * br;

    return ((1.f - alpha_y) * I_top + alpha_y * I_bottom);
}

static glm::vec3 bilinear_texture(const ImageData &image, const uv &uv_map, const float m, const float n)
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
    const glm::vec3 I_bottom = (1.f - alpha_x) * glm::vec3(bl.r / 255.f, bl.g / 255.f, bl.b / 255.f) +
                               alpha_x * glm::vec3(br.r / 255.f, br.g / 255.f, br.b / 255.f);

    return ((1.f - alpha_y) * I_top + alpha_y * I_bottom);
}

glm::vec3 get_texture(const ImageData &texture, const Hit &hit, const FilterType filter_type, const uv &uv_map)
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

glm::vec3 get_bump_normal(const BumpMap &bump_map, const Hit &hit, const FilterType filter_type, const uv &uv_map,
                          float bump_scale)
{
    const float w = bump_map.width;
    const float h = bump_map.height;
    const float m = 1; // HARDCODED
    const float n = 1; // HARDCODED

    const auto [dpdu, dpdv] = get_dpduv(hit, uv_map);
    glm::vec2 g;

    switch (filter_type)
    {
    case FilterType::Nearest:
    {
        const glm::vec2 *gradients = bump_map.gradients;

        const int c = std::fmod(std::floorf(uv_map.u * m * w) + w, w);
        const int r = std::fmod(std::floorf((1 - uv_map.v) * n * h) + h, h);

        g = gradients[static_cast<int>(r * w + c)];
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

    glm::vec3 dNdu = bump_scale * dHdu * glm::cross(hit.normal, dpdv);
    glm::vec3 dNdv = bump_scale * dHdv * glm::cross(dpdu, hit.normal);

    glm::vec3 perturbed_normal = hit.normal + dNdu + dNdv;
    return glm::normalize(perturbed_normal);
}

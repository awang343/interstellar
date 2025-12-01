#include "textures.h"
#include "glm/fwd.hpp"
#include "raytracer.h"
#include "raytracescene.h"
#include "utils/ini_utils.h"

uv get_uv(const Hit &hit)
{
    switch (hit.shape->primitive.type)
    {
    case PrimitiveType::PRIMITIVE_CUBE:
    {
        uv uv_coords;
        if (fabs(hit.normal.x - 1) < 1e-6)
            uv_coords = {-hit.point.z, hit.point.y};
        else if (fabs(hit.normal.y - 1) < 1e-6)
            uv_coords = {hit.point.x, -hit.point.z};
        else if (fabs(hit.normal.z - 1) < 1e-6)
            uv_coords = {hit.point.x, hit.point.y};
        else if (fabs(hit.normal.x + 1) < 1e-6)
            uv_coords = {hit.point.z, hit.point.y};
        else if (fabs(hit.normal.y + 1) < 1e-6)
            uv_coords = {hit.point.x, hit.point.z};
        else if (fabs(hit.normal.z + 1) < 1e-6)
            uv_coords = {-hit.point.x, hit.point.y};

        uv_coords.u += 0.5;
        uv_coords.v += 0.5;

        return uv_coords;
    }

    case PrimitiveType::PRIMITIVE_SPHERE:
    {
        const float theta = atan2(hit.point.z, hit.point.x);
        const float u = theta < 0 ? -(theta) / (2.f * M_PI) : 1.f - (theta) / (2.f * M_PI);
        const float v = asin(hit.point.y / 0.5f) / M_PI + .5f;
        return uv{u, v};
    }

    case PrimitiveType::PRIMITIVE_CYLINDER:
    {
        float y = hit.point.y;
        float radius = 0.5f;

        if (fabs(y - 0.5f) < 1e-6f)
        { // top cap
            float u = hit.point.x / radius * 0.5f + 0.5f;
            float v = -hit.point.z / radius * 0.5f + 0.5f;
            return uv{u, v};
        }
        else if (fabs(y + 0.5f) < 1e-6f)
        { // bottom cap
            float u = hit.point.x / radius * 0.5f + 0.5f;
            float v = hit.point.z / radius * 0.5f + 0.5f;
            return uv{u, v};
        }
        else
        { // side
            float theta = atan2(-hit.point.z, -hit.point.x);
            float u = 0.5f - (theta) / (2.0f * M_PI);
            float v = y + 0.5f;
            return uv{u, v};
        }
    }

    case PrimitiveType::PRIMITIVE_CONE:
    {
        float y = hit.point.y;
        float radius = 0.5f;

        if (fabs(y + 0.5f) < 1e-6f)
        { // base cap
            float u = hit.point.x / radius * 0.5f + 0.5f;
            float v = hit.point.z / radius * 0.5f + 0.5f;
            return uv{u, v};
        }
        else
        { // side
            float theta = atan2(-hit.point.z, -hit.point.x);
            float u = 0.5f - (theta) / (2.0f * M_PI);
            float v = y + 0.5f; // apex at v=0, base at v=1
            return uv{u, v};
        }
    }

    default:
        return uv{0.f, 0.f};
    }
}

static std::pair<glm::vec3, glm::vec3> get_dpduv(const Hit &hit, const uv &uv_coords)
{
    switch (hit.shape->primitive.type)
    {
    case PrimitiveType::PRIMITIVE_CUBE:
    {
        if (fabs(hit.normal.x - 1) < 1e-6)
            return {-glm::vec3(0, 0, 1), glm::vec3(0, 1, 0)};
        else if (fabs(hit.normal.y - 1) < 1e-6)
            return {glm::vec3(1, 0, 0), -glm::vec3(0, 0, 1)};
        else if (fabs(hit.normal.z - 1) < 1e-6)
            return {glm::vec3(1, 0, 0), glm::vec3(0, 1, 0)};
        else if (fabs(hit.normal.x + 1) < 1e-6)
            return {glm::vec3(0, 0, 1), glm::vec3(0, 1, 0)};
        else if (fabs(hit.normal.y + 1) < 1e-6)
            return {glm::vec3(1, 0, 0), glm::vec3(0, 0, 1)};
        else if (fabs(hit.normal.z + 1) < 1e-6)
            return {-glm::vec3(1, 0, 0), glm::vec3(0, 1, 0)};
    }

    case PrimitiveType::PRIMITIVE_SPHERE:
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

    case PrimitiveType::PRIMITIVE_CYLINDER:
    {
        const float r = 0.5f;
        const float x = hit.point.x;
        const float y = hit.point.y;
        const float z = hit.point.z;

        if (fabs(y - 0.5f) < 1e-6f)
        {
            glm::vec3 dpdu(1.f, 0.f, 0.f);
            glm::vec3 dpdv(0.f, 0.f, -1.f);
            return {dpdu, dpdv};
        }
        else if (fabs(y + 0.5f) < 1e-6f)
        {
            glm::vec3 dpdu(1.f, 0.f, 0.f);
            glm::vec3 dpdv(0.f, 0.f, 1.f);
            return {dpdu, dpdv};
        }
        else
        {
            float theta = atan2(-z, -x);
            float dtheta_du = -2.f * M_PI;

            float cosT = cos(theta);
            float sinT = sin(theta);

            glm::vec3 dpdtheta(-r * sinT, 0.f, r * cosT);
            glm::vec3 dpdu = dtheta_du * dpdtheta;

            glm::vec3 dpdv(0.f, 1.f, 0.f);

            return {dpdu, dpdv};
        }
    }

    case PrimitiveType::PRIMITIVE_CONE:
    {
        const float r_base = 0.5f; // radius at base y = -0.5
        const float y = hit.point.y;

        if (fabs(y + 0.5f) < 1e-6f)
        {
            glm::vec3 dpdu(1.f, 0.f, 0.f);
            glm::vec3 dpdv(0.f, 0.f, 1.f);
            return {dpdu, dpdv};
        }
        else
        {
            float theta = atan2(-hit.point.z, -hit.point.x);
            const float r_at_y = r_base * (0.5f - y); // linear taper; verify sign with your orientation
            glm::vec3 dpdtheta(-r_at_y * sin(theta), 0.f, r_at_y * cos(theta));
            const float dtheta_du = -2.f * M_PI; // consistent with your get_uv
            glm::vec3 dpdu = dtheta_du * dpdtheta;

            const float dr_dy = -r_base;
            glm::vec3 dpdv(dr_dy * cos(theta), 1.f, dr_dy * sin(theta));

            return {dpdu, dpdv};
        }
    }

    default:
        return {glm::vec3(0, 0, 0), glm::vec3(0, 0, 0)}; // Placeholder
    }
}

static std::pair<glm::vec3, glm::vec3> get_duvdp(const Hit &hit)
{
    switch (hit.shape->primitive.type)
    {
    case PrimitiveType::PRIMITIVE_CUBE:
    {
        const float tol = 1e-5f;
        const glm::vec3 n = hit.normal;

        if (fabsf(n.x - 1.f) < tol) // +X face
            return {glm::vec3(0, 0, -1), glm::vec3(0, 1, 0)};
        else if (fabsf(n.x + 1.f) < tol) // -X face
            return {glm::vec3(0, 0, 1), glm::vec3(0, 1, 0)};
        else if (fabsf(n.y - 1.f) < tol) // +Y face
            return {glm::vec3(1, 0, 0), glm::vec3(0, 0, -1)};
        else if (fabsf(n.y + 1.f) < tol) // -Y face
            return {glm::vec3(1, 0, 0), glm::vec3(0, 0, 1)};
        else if (fabsf(n.z - 1.f) < tol) // +Z face
            return {glm::vec3(1, 0, 0), glm::vec3(0, 1, 0)};
        else // -Z face
            return {glm::vec3(-1, 0, 0), glm::vec3(0, 1, 0)};
    }

    case PrimitiveType::PRIMITIVE_SPHERE:
    {
        const float du_dx =
            (1.f / (2.f * M_PI)) * (hit.point.z / (hit.point.x * hit.point.x + hit.point.z * hit.point.z));
        const float du_dy = 0.f;
        const float du_dz =
            -(1.f / (2.f * M_PI)) * (hit.point.x / (hit.point.x * hit.point.x + hit.point.z * hit.point.z));

        const float dv_dx = 0.f;
        const float dv_dy = 1.f / (M_PI * std::sqrt(.5f * .5f - hit.point.y * hit.point.y));
        const float dv_dz = 0.f;

        return {glm::vec3(du_dx, du_dy, du_dz), glm::vec3(dv_dx, dv_dy, dv_dz)};
    }

    case PrimitiveType::PRIMITIVE_CYLINDER:
    {
        const float x = hit.point.x;
        const float y = hit.point.y;
        const float z = hit.point.z;
        const float r2 = x * x + z * z;

        if (fabs(y - 0.5f) < 1e-6f)
        { // top cap
            return {glm::vec3(1.f, 0.f, 0.f), glm::vec3(0.f, 0.f, -1.f)};
        }
        else if (fabs(y + 0.5f) < 1e-6f)
        { // bottom cap
            return {glm::vec3(1.f, 0.f, 0.f), glm::vec3(0.f, 0.f, 1.f)};
        }
        else
        { // side
            const float du_dx = -(z) / (2.f * M_PI * r2);
            const float du_dy = 0.f;
            const float du_dz = x / (2.f * M_PI * r2);

            return {glm::vec3(du_dx, du_dy, du_dz), glm::vec3(0, 1, 0)};
        }
    }

    case PrimitiveType::PRIMITIVE_CONE:
    {
        const float x = hit.point.x;
        const float y = hit.point.y;
        const float z = hit.point.z;
        const float r2 = x * x + z * z;

        // Side
        if (fabs(y + 0.5f) > 1e-6f)
        {
            const float du_dx = -(z) / (2.f * M_PI * r2);
            const float du_dy = 0.f;
            const float du_dz = x / (2.f * M_PI * r2);

            const float dv_dx = 0.f;
            const float dv_dy = 1.f;
            const float dv_dz = 0.f;

            return {glm::vec3(du_dx, du_dy, du_dz), glm::vec3(dv_dx, dv_dy, dv_dz)};
        }

        // Base cap (same as cylinder cap)
        glm::vec3 du = glm::vec3(1, 0, 0);
        glm::vec3 dv = glm::vec3(0, 0, 1);
        return {du, dv};
    }

    default:
        return {glm::vec3(0.f), glm::vec3(0.f)};
    }
}

static glm::vec2 bilinear_bump(std::shared_ptr<BumpMap> bmp, const float u, const float v, const float m, const float n)
{
    const float w = bmp->width;
    const float h = bmp->height;
    const glm::vec2 *bmp_data = bmp->gradients;

    const float x_left = u * m * w - 0.5f;
    const float x_right = x_left + 1.f;
    const float y_top = (1.f - v) * n * h - 0.5f;
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

static glm::vec3 bilinear_texture(std::shared_ptr<Image> image, const float u, const float v, const float m,
                                  const float n)
{
    const float w = image->width;
    const float h = image->height;
    const RGBA *color_data = image->data;

    const float x_left = u * m * w - 0.5f;
    const float x_right = x_left + 1.f;
    const float y_top = (1.f - v) * n * h - 0.5f;
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

float get_scale_factor(const RayTraceScene &scene, const Hit &hit, const glm::vec3 &norm, const glm::vec3 &d,
                       const uv &uv_map, const int samples)
{
    static const glm::mat3 inv_view = glm::inverse(scene.getCamera().getViewMatrix());
    const glm::vec3 norm_d = glm::normalize(d);

    const glm::vec3 dd_dx =
        inv_view * glm::vec3(scene.getCamera().getPlaneWidth() / scene.getWidth() / std::sqrt(samples), 0.f, 0.f);
    const glm::vec3 ddir_dx = (glm::dot(d, d) * dd_dx - glm::dot(d, dd_dx) * d) / std::pow(glm::dot(d, d), 1.5f);

    const glm::vec3 dd_dy =
        inv_view * glm::vec3(0.f, scene.getCamera().getPlaneHeight() / scene.getHeight() / std::sqrt(samples), 0.f);
    const glm::vec3 ddir_dy = (glm::dot(d, d) * dd_dy - glm::dot(d, dd_dy) * d) / std::pow(glm::dot(d, d), 1.5f);

    const float dtdx = -glm::dot(norm, hit.t * ddir_dx) / glm::dot(norm, norm_d);
    const float dtdy = -glm::dot(norm, hit.t * ddir_dy) / glm::dot(norm, norm_d);

    const glm::vec3 dpdx = hit.shape->inv_ctm * glm::vec4(hit.t * ddir_dx + dtdx * norm_d, 0.f);
    const glm::vec3 dpdy = hit.shape->inv_ctm * glm::vec4(hit.t * ddir_dy + dtdy * norm_d, 0.f);

    const auto [dudp, dvdp] = get_duvdp(hit);

    const float dsdu = (*hit.shape->texture_levels)[0]->width * hit.shape->primitive.material.textureMap.repeatU;
    const float dtdv = (*hit.shape->texture_levels)[0]->height * hit.shape->primitive.material.textureMap.repeatV;

    const float dSdx = glm::dot(dpdx, dudp) * dsdu;
    const float dTdx = glm::dot(dpdx, dvdp) * dtdv;
    const float dSdy = glm::dot(dpdy, dudp) * dsdu;
    const float dTdy = glm::dot(dpdy, dvdp) * dtdv;

    const float X = std::sqrtf(dSdx * dSdx + dTdx * dTdx);
    const float Y = std::sqrtf(dSdy * dSdy + dTdy * dTdy);

    return std::fmaxf(X, Y);
}

glm::vec3 get_texture(const Hit &hit, const FilterType filter_type, uv uv_map, float scale_factor)
{
    const float m = hit.shape->primitive.material.textureMap.repeatU;
    const float n = hit.shape->primitive.material.textureMap.repeatV;

    const float log_scale = std::log2f(scale_factor);
    const int mm_floor =
        std::clamp(static_cast<int>(log_scale), 0, static_cast<int>(std::size(*hit.shape->texture_levels)) - 1);
    const int mm_ceil = std::clamp(mm_floor + 1, 0, static_cast<int>(std::size(*hit.shape->texture_levels)) - 1);

    switch (filter_type)
    {
    case FilterType::Nearest:
    {
        const float w = (*hit.shape->texture_levels)[mm_ceil]->width;
        const float h = (*hit.shape->texture_levels)[mm_ceil]->height;
        const RGBA *color_data = (*hit.shape->texture_levels)[mm_ceil]->data;

        const int c = std::fmod(std::floorf(uv_map.u * m * w), w);
        const int r = std::fmod(std::floorf((1 - uv_map.v) * n * h), h);

        const RGBA color = color_data[static_cast<int>(r * w + c)];
        return glm::vec3(color.r, color.g, color.b) / 255.f;
    }
    case FilterType::Bilinear:
    {
        return bilinear_texture((*hit.shape->texture_levels)[mm_ceil], uv_map.u, uv_map.v, m, n);
    }
    case FilterType::Trilinear:
    {
        const float alpha = std::clamp(log_scale - mm_floor, 0.f, 1.f);
        return bilinear_texture((*hit.shape->texture_levels)[mm_floor], uv_map.u, uv_map.v, m, n) * (1.f - alpha) +
               bilinear_texture((*hit.shape->texture_levels)[mm_ceil], uv_map.u, uv_map.v, m, n) * alpha;
    }
    }
}

glm::vec3 get_bump_normal(const Hit &hit, const FilterType filter_type, uv uv_map, float bump_scale)
{
    const float w = hit.shape->bump_map->width;
    const float h = hit.shape->bump_map->height;
    const float m = hit.shape->primitive.material.bumpMap.repeatU;
    const float n = hit.shape->primitive.material.bumpMap.repeatV;

    const auto [dpdu, dpdv] = get_dpduv(hit, uv_map);
    glm::vec2 g;

    switch (filter_type)
    {
    case FilterType::Nearest:
    {
        const glm::vec2 *gradients = hit.shape->bump_map->gradients;

        const int c = std::fmod(std::floorf(uv_map.u * m * w) + w, w);
        const int r = std::fmod(std::floorf((1 - uv_map.v) * n * h) + h, h);

        g = gradients[static_cast<int>(r * w + c)];
        break;
    }
    case FilterType::Bilinear:
    {
        g = bilinear_bump(hit.shape->bump_map, uv_map.u, uv_map.v, m, n);
        break;
    }
    case FilterType::Trilinear:
    {
        std::cerr << "Trilinear filtering not supported for bump maps." << std::endl;
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

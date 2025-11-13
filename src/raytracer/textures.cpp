#include "textures.h"
#include "raytracer.h"
#include "raytracescene.h"

uv get_uv(const Hit &hit)
{
    switch (hit.shape->primitive.type)
    {
    case PrimitiveType::PRIMITIVE_CUBE:
    {
        float x, y;
        if (fabs(hit.normal.x - 1) < 1e-6)
        {
            x = -hit.point.z;
            y = hit.point.y;
        }
        else if (fabs(hit.normal.y - 1) < 1e-6)
        {
            x = hit.point.x;
            y = -hit.point.z;
        }
        else if (fabs(hit.normal.z - 1) < 1e-6)
        {
            x = hit.point.x;
            y = hit.point.y;
        }
        else if (fabs(hit.normal.x + 1) < 1e-6)
        {
            x = hit.point.z;
            y = hit.point.y;
        }
        else if (fabs(hit.normal.y + 1) < 1e-6)
        {
            x = hit.point.x;
            y = hit.point.z;
        }
        else if (fabs(hit.normal.z + 1) < 1e-6)
        {
            x = -hit.point.x;
            y = hit.point.y;
        }

        return uv{x + 0.5f, y + 0.5f};
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

static glm::vec3 bilinear(std::shared_ptr<Image> image, const float u, const float v, const float m, const float n)
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

glm::vec3 get_texture(const RenderShapeData *shape, const TextureFilterType filter_type, float u, float v,
                      float scale_factor)
{
    const float m = shape->primitive.material.textureMap.repeatU;
    const float n = shape->primitive.material.textureMap.repeatV;

    const float log_scale = std::log2f(scale_factor);
    const int mm_floor =
        std::clamp(static_cast<int>(log_scale), 0, static_cast<int>(std::size(*shape->texture_levels)) - 1);
    const int mm_ceil = std::clamp(mm_floor + 1, 0, static_cast<int>(std::size(*shape->texture_levels)) - 1);

    switch (filter_type)
    {
    case TextureFilterType::Nearest:
    {
        const float w = (*shape->texture_levels)[mm_ceil]->width;
        const float h = (*shape->texture_levels)[mm_ceil]->height;
        const RGBA *color_data = (*shape->texture_levels)[mm_ceil]->data;

        const int c = std::fmod(std::floorf(u * m * w), w);
        const int r = std::fmod(std::floorf((1 - v) * n * h), h);

        const RGBA color = color_data[static_cast<int>(r * w + c)];
        return glm::vec3(color.r, color.g, color.b) / 255.f;
    }
    case TextureFilterType::Bilinear:
    {
        return bilinear((*shape->texture_levels)[mm_ceil], u, v, m, n);
    }
    case TextureFilterType::Trilinear:
    {
        const float alpha = std::clamp(log_scale - mm_floor, 0.f, 1.f);
        return bilinear((*shape->texture_levels)[mm_floor], u, v, m, n) * (1.f - alpha) +
               bilinear((*shape->texture_levels)[mm_ceil], u, v, m, n) * alpha;
    }
    }
}

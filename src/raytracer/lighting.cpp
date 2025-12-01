#include "lighting.h"
#include "intersector.h"
#include "textures.h"

static glm::vec3 calc_light_vec(const RayTraceScene &scene, const glm::vec3 &point, const SceneLightData &light)
{
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

glm::vec3 RayTracer::shadePixel(const RayTraceScene &scene, const Hit &hit, const glm::vec3 &world_camera,
                                const glm::vec3 &world_dir, unsigned int &last_ray_id,
                                std::vector<unsigned int> &last_visited, int depth, float passthrough_scale_factor)
{
    static const float ka = scene.getGlobalData().ka;
    static const float kd = scene.getGlobalData().kd;
    static const float ks = scene.getGlobalData().ks;

    const glm::vec3 obj_A = hit.shape->primitive.material.cAmbient * ka;
    glm::vec3 obj_D = hit.shape->primitive.material.cDiffuse * kd;
    const glm::vec3 obj_S = hit.shape->primitive.material.cSpecular * ks;
    const glm::vec3 obj_refl = hit.shape->primitive.material.cReflective * ks;
    const float shininess = hit.shape->primitive.material.shininess;

    const glm::vec4 intersect = hit.shape->ctm * glm::vec4(hit.point, 1.f);
    const glm::mat3 transform = glm::inverse(glm::transpose(glm::mat3(hit.shape->ctm)));

    const uv uv_map = get_uv(hit);

    const glm::vec3 obj_normal = m_config.enableBumpMap && hit.shape->bump_map
                                     ? get_bump_normal(hit, m_config.bumpMapFilterType, uv_map, m_config.bumpScale)
                                     : hit.normal;
    const glm::vec3 N = glm::normalize(transform * obj_normal);
    const glm::vec3 V = glm::normalize(world_camera - glm::vec3(intersect));

    // Blend obj_D with textures
    float scale_factor = passthrough_scale_factor;
    if (m_config.enableTextureMap && hit.shape->texture_levels)
    {
        const float blend = hit.shape->primitive.material.blend;
        if (scale_factor < 0.f) // If no passthrough scale factor provided, compute it
        {
            const int spp = m_config.enableSuperSample ? m_config.samplesPerPixel : 1.f;
            scale_factor = m_config.enableMipMapping ? get_scale_factor(scene, hit, N, world_dir, uv_map, spp) : 1.f;
        }
        const glm::vec3 texture_color = get_texture(hit, m_config.textureFilterType, uv_map, scale_factor);
        obj_D = obj_D * (1.f - blend) + texture_color * blend;
    }

    glm::vec3 pixel_color(obj_A);

    for (const SceneLightData &light : scene.getLights())
    {
        glm::vec3 phong(0.f);

        const glm::vec3 to_light = calc_light_vec(scene, intersect, light);
        const glm::vec3 L_i = glm::normalize(to_light);
        const glm::vec3 R = glm::normalize(2.f * glm::dot(L_i, N) * N - L_i);

        const float D_intensity = std::max(glm::dot(N, L_i), 0.f);
        const float S_intensity = std::pow(std::max(glm::dot(R, V), 0.f), shininess);

        const glm::vec3 D_color = glm::vec3(light.color) * obj_D;
        const glm::vec3 S_color = glm::vec3(light.color) * obj_S;

        phong += D_color * D_intensity;
        phong += S_color * S_intensity;

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
            const float dist = glm::distance(glm::vec3(intersect), glm::vec3(light.pos));
            const float mult = 1.f / glm::dot(light.function, glm::vec3(1.f, dist, dist * dist));
            light_intensity *= std::min(1.f, mult);
        }

        // Shadow effect
        if (m_config.enableShadow && light_intensity > 1e-6f)
        {
            const glm::vec4 point = intersect + glm::vec4(L_i * 1e-3f, 0.f);
            const auto shadow_hit = traverseKDTree(scene, point, L_i, false, last_ray_id, last_visited);

            if (shadow_hit)
            {
                const glm::vec4 p = (*shadow_hit).shape->ctm * glm::vec4((*shadow_hit).point, 1.f);
                const float dist_to_obj = glm::length(p - intersect);
                const float dist_to_light = glm::length(light.pos - intersect);
                light_intensity *= directional || dist_to_obj < dist_to_light ? 0.f : 1.f;
            }
        }

        pixel_color += phong * light_intensity;
    }

    // Reflectivity
    if (m_config.enableReflection && glm::length(obj_refl) > 0.f && depth < m_config.maxRecursiveDepth)
    {
        const glm::vec3 refl_dir = glm::normalize(2.f * glm::dot(V, N) * N - V);
        const glm::vec4 point = intersect + glm::vec4(refl_dir * 1e-3f, 0.f);
        const auto refl_hit = traverseKDTree(scene, point, refl_dir, false, last_ray_id, last_visited);

        if (refl_hit)
        {
            pixel_color += obj_refl * shadePixel(scene, *refl_hit, point, refl_dir, last_ray_id, last_visited,
                                                 depth + 1, scale_factor);
        }
    }

    return glm::clamp(pixel_color, 0.f, 1.f);
}

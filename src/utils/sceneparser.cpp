#include "sceneparser.h"
#include "scenefilereader.h"

#include <glm/gtx/transform.hpp>
#include <iostream>
#include <memory>
#include <stack>

using ImageVector = std::vector<std::shared_ptr<Image>>;

inline int reflect(int x, int max)
{
    if (x < 0)
        return -x - 1;
    if (x >= max)
        return 2 * max - x - 1;
    return x;
}

std::shared_ptr<Image> downsample(std::shared_ptr<Image> base, const int factor)
{
    // factor must be a power of 2 for kernel generation to work
    std::vector<float> kernel(2 * factor);
    for (int i = 0; i < std::size(kernel); i++)
    {
        kernel[i] = (factor - fabs(i + 0.5f - factor)) / (factor * factor);
    }

    const int w = base->width;
    const int h = base->height;
    const RGBA *data = base->data;
    const int new_w = base->width / factor;
    const int new_h = base->height / factor;

    std::vector<RGBA> horz_filtered(new_w * h);
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < new_w; j++)
        {
            const float src_j = (j + .5f) * factor - .5f;
            glm::vec3 result(0.f);
            for (int o = 0; o < std::size(kernel); o++)
            {
                const int oj = reflect(static_cast<int>(src_j + o - factor + .5f), w - 1);

                const RGBA &px = data[i * w + oj];
                result.x += px.r / 255.f * kernel[o];
                result.y += px.g / 255.f * kernel[o];
                result.z += px.b / 255.f * kernel[o];
            }

            horz_filtered[i * new_w + j] =
                RGBA{static_cast<uint8_t>(result.x * 255.f), static_cast<uint8_t>(result.y * 255.f),
                     static_cast<uint8_t>(result.z * 255.f)};
        }
    }

    RGBA *filtered = new RGBA[new_w * new_h];
    for (int i = 0; i < new_h; i++)
    {
        for (int j = 0; j < new_w; j++)
        {
            const float src_i = (i + .5f) * factor - .5f;
            glm::vec3 result(0.f);
            for (int o = 0; o < std::size(kernel); o++)
            {
                const int oi = reflect(static_cast<int>(src_i + o - factor + .5f), h - 1);

                const RGBA &px = horz_filtered[oi * new_w + j];
                result.x += px.r / 255.f * kernel[o];
                result.y += px.g / 255.f * kernel[o];
                result.z += px.b / 255.f * kernel[o];
            }

            filtered[i * new_w + j] =
                RGBA{static_cast<uint8_t>(result.x * 255.f), static_cast<uint8_t>(result.y * 255.f),
                     static_cast<uint8_t>(result.z * 255.f)};
        }
    }

    return std::shared_ptr<Image>(new Image{filtered, new_w, new_h});
}

std::shared_ptr<ImageVector> generateMipmaps(std::string filename, std::string outputPath)
{
    auto vec = std::make_shared<ImageVector>();
    std::shared_ptr<Image> base(loadImageFromFile(filename));

    vec->push_back(base);

    int w = base->width;
    int h = base->height;
    int ds_factor = 2;

    while (w > 1 || h > 1)
    {
        auto new_image = downsample(base, ds_factor);
        w = new_image->width;
        h = new_image->height;

        std::string name = filename.substr(filename.find_last_of('/') + 1);
        std::string basename = name.substr(0, name.find_last_of('.'));
        std::string ext = name.substr(name.find_last_of('.'));
        std::string path = outputPath + basename + "_" + std::to_string(ds_factor) + ext;

        QImage(reinterpret_cast<const uchar *>(new_image->data), w, h,
               w * 4, // bytes per line
               QImage::Format_RGBA8888)
            .save(path.c_str());
        vec->push_back(new_image);

        ds_factor *= 2;
    }
    return vec;
}

float sample_height_wrap(const Image &img, int x, int y)
{
    x = (x % img.width + img.width) % img.width;
    y = (y % img.height + img.height) % img.height;

    return img.data[y * img.width + x].r / 255.0f; // grayscale in .r
}

std::shared_ptr<BumpMap> generateBumpMap(std::string filename)
{
    Image *bump_image = loadImageFromFile(filename);

    BumpMap bump_map;
    bump_map.width = bump_image->width;
    bump_map.height = bump_image->height;

    bump_map.gradients = new glm::vec2[bump_map.width * bump_map.height];

    for (int y = 0; y < bump_map.height; y++)
    {
        for (int x = 0; x < bump_map.width; x++)
        {

            float hL = sample_height_wrap(*bump_image, x - 1, y);
            float hR = sample_height_wrap(*bump_image, x + 1, y);
            float hD = sample_height_wrap(*bump_image, x, y + 1);
            float hU = sample_height_wrap(*bump_image, x, y - 1);

            float dX = (hR - hL) * 0.5f;
            float dY = (hD - hU) * 0.5f;

            bump_map.gradients[y * bump_map.width + x] = glm::vec2(dX, dY);
        }
    }

    return std::make_shared<BumpMap>(bump_map);
}

void scene_dfs(SceneNode *root, std::string texturepath, std::string mipspath, RenderData &renderData)
{
    std::stack<std::pair<SceneNode *, glm::mat4>> stack;
    std::unordered_map<std::string, std::shared_ptr<ImageVector>> textureCache;
    std::unordered_map<std::string, std::shared_ptr<BumpMap>> bumpMapCache;

    stack.push({root, glm::mat4(1.0f)}); // start with identity matrix

    while (!stack.empty())
    {
        SceneNode *node;
        glm::mat4 ctm;
        std::tie(node, ctm) = stack.top();
        stack.pop();

        // Apply node transformations
        for (SceneTransformation *transformation : node->transformations)
        {
            switch (transformation->type)
            {
            case TransformationType::TRANSFORMATION_TRANSLATE:
                ctm = ctm * glm::translate(transformation->translate);
                break;
            case TransformationType::TRANSFORMATION_SCALE:
                ctm = ctm * glm::scale(transformation->scale);
                break;
            case TransformationType::TRANSFORMATION_ROTATE:
                ctm = ctm * glm::rotate(transformation->angle, transformation->rotate);
                break;
            case TransformationType::TRANSFORMATION_MATRIX:
                ctm = ctm * transformation->matrix;
                break;
            default:
                std::cerr << "Unknown transformation type" << std::endl;
                break;
            }
        }

        // Process primitives
        for (ScenePrimitive *primitive : node->primitives)
        {
            const std::string &texture_name = primitive->material.textureMap.filename;
            const std::string &bumpmap_name = primitive->material.bumpMap.filename;

            std::shared_ptr<ImageVector> texture = nullptr;
            if (!texture_name.empty())
            {
                const std::string &path = texturepath + texture_name;

                // If already loaded, reuse it
                auto it = textureCache.find(texture_name);
                if (it != textureCache.end())
                {
                    texture = it->second;
                }
                else
                {
                    // Otherwise, load and store
                    texture = generateMipmaps(path, mipspath);
                    textureCache[texture_name] = texture;
                }
            }

            std::shared_ptr<BumpMap> bump_map = nullptr;
            if (!bumpmap_name.empty())
            {
                const std::string &path = texturepath + bumpmap_name;

                // If already loaded, reuse it
                auto it = bumpMapCache.find(bumpmap_name);
                if (it != bumpMapCache.end())
                {
                    bump_map = it->second;
                }
                else
                {
                    // Otherwise, load and store
                    bump_map = generateBumpMap(path);
                    bumpMapCache[bumpmap_name] = bump_map;
                }
            }

            renderData.shapes.push_back(
                RenderShapeData{std::size(renderData.shapes), *primitive, ctm, glm::inverse(ctm), texture, bump_map});
        }

        // Process lights
        for (SceneLight *light : node->lights)
        {
            glm::vec4 world_pos = ctm * glm::vec4(0.f, 0.f, 0.f, 1.f);
            glm::vec4 world_dir = glm::normalize(ctm * light->dir);

            SceneLightData sl_data{light->id, light->type,     light->color, light->function, world_pos,
                                   world_dir, light->penumbra, light->angle, light->width,    light->height};
            renderData.lights.push_back(sl_data);
        }

        // Push children onto stack
        // Pushing in reverse order preserves original traversal order
        for (auto it = node->children.rbegin(); it != node->children.rend(); ++it)
        {
            stack.push({*it, ctm});
        }
    }
}

bool SceneParser::parse(std::string scenepath, std::string texturepath, std::string mipspath, RenderData &renderData)
{
    ScenefileReader fileReader = ScenefileReader(scenepath);
    bool success = fileReader.readJSON();
    if (!success)
    {
        return false;
    }

    renderData.globalData = fileReader.getGlobalData();
    renderData.cameraData = fileReader.getCameraData();

    renderData.shapes.clear();
    renderData.lights.clear();
    scene_dfs(fileReader.getRootNode(), texturepath, mipspath, renderData);

    return true;
}

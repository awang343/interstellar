#pragma once

// a struct containing info of an image
#include "utils/rgba.h"
#include "utils/imagereader.h"
#include <iostream>
#include <qimage.h>
#include <vector>
#include <glm/glm.hpp>

struct ImageData {
    int width  = 0;
    int height = 0;
    std::vector<RGBA> pixels;  // size = width * height
};

struct BumpMap
{
    int width;
    int height;
    std::vector<glm::vec2> gradients;
};


// a struct containing info of a primitive
struct Sphere {
    ImageData textureFile; // should change to path to allow for hashing for multiple objects
    std::vector<std::vector<float>> points;
};

struct Cube {
    ImageData textureFile;
    std::vector<std::vector<float>> points;
    float side;
};


enum class PrimitiveType { Sphere, Cube };

struct Object {
    PrimitiveType type;

    ImageData textureFile;
    BumpMap bumpMapFile;
    std::vector<std::vector<float>> points;

    float radius = 0.0f;
    float side   = 0.0f;
};


// this function loads the two celestial sphere textures into image data
inline bool loadImageToStruct(const QString &path, ImageData &out)
{
    QImage img(path);
    if (img.isNull()) {
        std::cerr << "Failed to load image: " << path.toStdString() << "\n";
        return false;
    }

    QImage converted = img.convertToFormat(QImage::Format_RGBA8888);

    out.width  = converted.width();
    out.height = converted.height();
    out.pixels.resize(out.width * out.height);

    for (int y = 0; y < out.height; ++y) {
        const uchar *line = converted.constScanLine(y);
        for (int x = 0; x < out.width; ++x) {
            const uchar *src = line + 4 * x;
            RGBA &dst = out.pixels[y * out.width + x];
            dst.r = src[0];
            dst.g = src[1];
            dst.b = src[2];
            dst.a = src[3];
        }
    }
    return true;
}

float sample_height_wrap(const Image &img, int x, int y)
{
    x = (x % img.width + img.width) % img.width;
    y = (y % img.height + img.height) % img.height;

    return img.data[y * img.width + x].r / 255.0f; // grayscale in .r
}

// this function loads the two celestial sphere textures into image data
inline bool loadBumpMapToStruct(const QString &path, BumpMap &out)
{
    Image *bump_image = loadImageFromFile(path.toStdString());

    BumpMap bump_map;
    bump_map.width = bump_image->width;
    bump_map.height = bump_image->height;

    bump_map.gradients.resize(bump_map.width * bump_map.height);

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
    return true;
}

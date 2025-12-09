#pragma once

// a struct containing info of an image
#include "utils/rgba.h"
#include <iostream>
#include <qimage.h>
#include <vector>
struct ImageData {
    int width  = 0;
    int height = 0;
    std::vector<RGBA> pixels;  // size = width * height
};

// a struct containing info of a primitive
struct Sphere {
    ImageData textureFile; // should change to path to allow for hashing for multiple objects
    std::vector<std::vector<float>> points;
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

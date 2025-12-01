#pragma once

#include <QtCore>

enum class FilterType {
    Nearest = 0,
    Bilinear = 1,
    Trilinear = 2,
};

enum class SuperSamplerPattern {
    Grid = 0,
    Stratified = 1,
    Random = 2,
};

namespace IniUtils {
    FilterType filterTypeFromString(const QString& str);
    SuperSamplerPattern superSamplerPatternFromString(const QString& str);
} // namespace IniUtils

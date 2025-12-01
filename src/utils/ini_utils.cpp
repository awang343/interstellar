#include "ini_utils.h"

FilterType IniUtils::filterTypeFromString(const QString& str) {
    if (str == "nearest")
        return FilterType::Nearest;
    else if (str == "bilinear")
        return FilterType::Bilinear;
    else if (str == "trilinear")
        return FilterType::Trilinear;
    else
        throw std::runtime_error("Invalid filter type string.");
}

SuperSamplerPattern IniUtils::superSamplerPatternFromString(const QString& str) {
    if (str == "grid")
        return SuperSamplerPattern::Grid;
    else if (str == "stratified")
        return SuperSamplerPattern::Stratified;
    else if (str == "random")
        return SuperSamplerPattern::Random;
    else
        throw std::runtime_error("Invalid supersampler pattern string.");
}

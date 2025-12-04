#include "raytracescene.h"
#include "utils/sceneparser.h"

RayTraceScene::RayTraceScene(int width, int height, const RenderData &metaData)
    : w(width), h(height), data(metaData), camera(data.cameraData, width, height) {
    // Optional TODO: implement this. Store whatever you feel is necessary.
    w = width;
    h = height;
    data = metaData;
    camera = Camera(data.cameraData, width, height);
}

const int& RayTraceScene::width() const {
    return w;
}

const int& RayTraceScene::height() const {
    return h;
}

const SceneGlobalData& RayTraceScene::getGlobalData() const {
    return data.globalData;
}

const Camera& RayTraceScene::getCamera() const {
    return camera;
}

const std::vector<RenderShapeData>& RayTraceScene::getShapes() const {
    return data.shapes;
}

const std::vector<SceneLightData>& RayTraceScene::getLights() const {
    return data.lights;
}

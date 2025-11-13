#pragma once

#include <memory>

#include "camera/camera.h"
#include "utils/scenedata.h"
#include "utils/sceneparser.h"

#define PARALLEL_DEPTH 5

// A class representing a scene to be ray-traced

// Feel free to make your own design choices for RayTraceScene, the functions
// below are all optional / for your convenience. You can either implement and
// use these getters, or make your own design. If you decide to make your own
// design, feel free to delete these as TAs won't rely on them to grade your
// assignments.

struct AABB
{
    glm::vec3 min, max;
};

struct KDNode
{
    AABB bounds;
    int axis;
    float split;
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;
    std::vector<const RenderShapeData*> shapes;  // only for leaf nodes
    bool isLeaf() const { return left == nullptr && right == nullptr; }
};

class RayTraceScene
{
   public:
    RayTraceScene(int width, int height, const RenderData& metaData);

    // The getter of the width of the scene
    const int& getWidth() const;

    // The getter of the height of the scene
    const int& getHeight() const;

    // The getter of the global data of the scene
    const SceneGlobalData& getGlobalData() const;

    // The getter of the shared pointer to the camera instance of the scene
    const Camera& getCamera() const;

    const std::vector<SceneLightData>& getLights() const;
    const std::vector<RenderShapeData>& getShapes() const;
    const std::vector<glm::vec3>& getObjCams() const;

    void buildKDTreeSAH();
    const KDNode* getKDRoot() const;

   private:
    int width;
    int height;
    Camera camera;
    const SceneGlobalData& globalData;
    const std::vector<SceneLightData>& lights;
    const std::vector<RenderShapeData>& shapes;
    const std::vector<glm::vec3> obj_cams;
    std::unique_ptr<KDNode> kdTreeRoot;
};

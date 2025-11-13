#include "raytracescene.h"
#include "kdtree.h"

#include "utils/sceneparser.h"


Camera get_camera(int width, int height, const SceneCameraData &camera_data)
{
    glm::vec3 w = glm::normalize(glm::vec3(-camera_data.look));
    glm::vec3 v =
        glm::normalize(glm::vec3(camera_data.up) - glm::dot(glm::vec3(camera_data.up), w) * w);
    glm::vec3 u = glm::cross(v, w);

    glm::mat4 R(u[0], v[0], w[0], 0.f, u[1], v[1], w[1], 0.f, u[2], v[2], w[2], 0.f, 0.f, 0.f, 0.f,
                1.f);
    glm::mat4 T(1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, -camera_data.pos[0],
                -camera_data.pos[1], -camera_data.pos[2], 1.f);

    glm::mat4 view_matrix = R * T;

    float ratio = static_cast<float>(width) / height;

    return Camera(view_matrix, ratio, camera_data.heightAngle);
}

std::vector<glm::vec3> get_obj_cams(std::vector<RenderShapeData> shapes, Camera cam)
{
    const glm::mat4 view_matrix = cam.getViewMatrix();
    const glm::mat4 inverse_view = glm::inverse(view_matrix);
    const glm::vec4 world_camera(inverse_view * glm::vec4(0.f, 0.f, 0.f, 1.f));

    std::vector<glm::vec3> O_vec;
    O_vec.reserve(shapes.size());
    for (const auto &shape : shapes)
    {
        O_vec.push_back(shape.inv_ctm * world_camera);
    }

    return O_vec;
}

RayTraceScene::RayTraceScene(int width, int height, const RenderData &metaData)
    : width(width),
      height(height),
      globalData(metaData.globalData),
      lights(metaData.lights),
      shapes(metaData.shapes),
      camera(get_camera(width, height, metaData.cameraData)),
      obj_cams(get_obj_cams(shapes, camera))
{
}

const int &RayTraceScene::getWidth() const
{
    return width;
}

const int &RayTraceScene::getHeight() const
{
    return height;
}

const SceneGlobalData &RayTraceScene::getGlobalData() const
{
    return globalData;
}

const Camera &RayTraceScene::getCamera() const
{
    return camera;
}

const std::vector<SceneLightData> &RayTraceScene::getLights() const
{
    return lights;
}

const std::vector<RenderShapeData> &RayTraceScene::getShapes() const
{
    return shapes;
}

const std::vector<glm::vec3> &RayTraceScene::getObjCams() const
{
    return obj_cams;
}

const KDNode *RayTraceScene::getKDRoot() const
{
    return kdTreeRoot.get();
}

void RayTraceScene::buildKDTreeSAH()
{
    std::vector<ShapeWithBounds> shapesWithBounds;
    for (const auto &s : getShapes())
    {
        shapesWithBounds.push_back({&s, computeBounds(s)});
    }

    kdTreeRoot = recurseKDTreeSAH(shapesWithBounds);
}



#include "intersector.h"

#include <glm/gtx/component_wise.hpp>
#include <limits>
#include <stack>

// ------------------- Utilities -------------------

static inline std::optional<std::pair<float, float>> solveQuadratic(float a, float b, float c) {
    const float disc = b * b - 4 * a * c;
    if (disc < 0)
        return std::nullopt;
    const float s = std::sqrtf(disc);
    return std::make_pair((-b - s) / (2 * a), (-b + s) / (2 * a));
}

static bool intersectAABB(
    const AABB &box, const glm::vec3 &eye, const glm::vec3 &dir, const glm::vec3 &inv_dir,
    float &tmin, float &tmax
) {
    const glm::vec3 t1 = (box.min - glm::vec3(eye)) * inv_dir;
    const glm::vec3 t2 = (box.max - glm::vec3(eye)) * inv_dir;

    const glm::vec3 tminVec = glm::min(t1, t2);
    const glm::vec3 tmaxVec = glm::max(t1, t2);

    tmin = glm::compMax(tminVec);
    tmax = glm::compMin(tmaxVec);

    return tmax >= tmin && tmax >= 0.f;
}

// ------------------- Shape Intersections -------------------

static std::optional<Hit>
intersectCube(const glm::vec3 &O, const glm::vec3 &D, const RenderShapeData &shape) {
    float tnear = -std::numeric_limits<float>::infinity();
    float tfar = std::numeric_limits<float>::infinity();

    for (int i = 0; i < 3; i++) {
        if (fabsf(D[i]) > 1e-6f) {
            float t1 = (-0.5f - O[i]) / D[i];
            float t2 = (0.5f - O[i]) / D[i];
            tnear = std::max(tnear, std::min(t1, t2));
            tfar = std::min(tfar, std::max(t1, t2));
        } else if (O[i] < -0.5f || O[i] > 0.5f) {
            return std::nullopt; // parallel ray outside bounds
        }
    }

    if (tnear > tfar || tfar < 0)
        return std::nullopt;

    float t = (tnear >= 0) ? tnear : tfar;
    glm::vec3 p = O + t * D;

    glm::vec3 n(0.f);
    const float tol = 1e-5f;
    if (fabsf(p.x - 0.5f) < tol)
        n = {1, 0, 0};
    else if (fabsf(p.x + 0.5f) < tol)
        n = {-1, 0, 0};
    else if (fabsf(p.y - 0.5f) < tol)
        n = {0, 1, 0};
    else if (fabsf(p.y + 0.5f) < tol)
        n = {0, -1, 0};
    else if (fabsf(p.z - 0.5f) < tol)
        n = {0, 0, 1};
    else if (fabsf(p.z + 0.5f) < tol)
        n = {0, 0, -1};

    if (glm::dot(D, n) > 0.f)
        n = -1.f * n;

    return Hit{t, p, n, &shape};
}

static std::optional<Hit>
intersectSphere(const glm::vec3 &O, const glm::vec3 &D, const RenderShapeData &shape) {
    const float r = 0.5f;
    const float a = glm::dot(D, D);
    const float b = 2.f * glm::dot(O, D);
    const float c = glm::dot(O, O) - r * r;

    const auto roots = solveQuadratic(a, b, c);
    if (!roots)
        return std::nullopt;

    const auto [t1, t2] = *roots;
    const float t = (t1 >= 0.f) ? t1 : (t2 >= 0.f ? t2 : -1.f);
    if (t < 0.f)
        return std::nullopt;

    const glm::vec3 p = glm::clamp(O + t * D, -.5f, .5f);
    return Hit{t, p, p, &shape};
}

static std::optional<Hit>
intersectCylinder(const glm::vec3 &O, const glm::vec3 &D, const RenderShapeData &shape) {
    float r = 0.5f;
    float tClosest = std::numeric_limits<float>::infinity();
    glm::vec3 hitPoint, normal;

    // Side
    float a = D.x * D.x + D.z * D.z;
    if (fabs(a) > 1e-5f) {
        float b = 2.0f * (O.x * D.x + O.z * D.z);
        float c = O.x * O.x + O.z * O.z - r * r;
        auto roots = solveQuadratic(a, b, c);

        if (roots) {
            for (float tside : {roots->first, roots->second}) {
                if (tside < 0)
                    continue;
                glm::vec3 p = O + tside * D;
                if (p.y >= -0.5f && p.y <= 0.5f && tside < tClosest) {
                    tClosest = tside;
                    hitPoint = p;
                    normal = glm::vec3(p.x, 0.0f, p.z);
                }
            }
        }
    }

    // Caps
    if (D.y != 0.f) {
        for (float ycap : {-0.5f, 0.5f}) {
            float tcap = (ycap - O.y) / D.y;
            if (tcap < 0)
                continue;
            glm::vec3 p = O + tcap * D;
            if (p.x * p.x + p.z * p.z <= r * r && tcap < tClosest) {
                tClosest = tcap;
                hitPoint = p;
                normal = (ycap > 0) ? glm::vec3(0, 1, 0) : glm::vec3(0, -1, 0);
            }
        }
    }

    if (tClosest < std::numeric_limits<float>::infinity())
        return Hit{tClosest, hitPoint, normal, &shape};

    return std::nullopt;
}

static std::optional<Hit>
intersectCone(const glm::vec3 &O, const glm::vec3 &D, const RenderShapeData &shape) {
    float tClosest = std::numeric_limits<float>::infinity();
    glm::vec3 hitPoint, normal;

    float a = D.x * D.x + D.z * D.z - D.y * D.y / 4.f;
    if (fabs(a) > 1e-5f) {
        float b = 2.0f * (O.x * D.x + O.z * D.z - (O.y - 0.5f) * D.y / 4.f);
        float c = O.x * O.x + O.z * O.z - (O.y - 0.5f) * (O.y - 0.5f) / 4.f;

        auto roots = solveQuadratic(a, b, c);
        if (roots) {
            for (float tside : {roots->first, roots->second}) {
                if (tside < 0)
                    continue;
                glm::vec3 p = O + tside * D;
                if (p.y >= -0.5f && p.y <= 0.5f && tside < tClosest) {
                    tClosest = tside;
                    hitPoint = p;
                    normal = glm::vec3(p.x, -(p.y - 0.5f) / 4.f, p.z);
                }
            }
        }
    }

    // Base cap
    if (D.y != 0.0f) {
        float tCap = (-0.5f - O.y) / D.y;
        if (tCap >= 0) {
            glm::vec3 p = O + tCap * D;
            if (p.x * p.x + p.z * p.z <= 0.25f && tCap < tClosest) {
                tClosest = tCap;
                hitPoint = p;
                normal = {0, -1, 0};
            }
        }
    }

    if (tClosest < std::numeric_limits<float>::infinity())
        return Hit{tClosest, hitPoint, normal, &shape};

    return std::nullopt;
}

static std::optional<Hit>
intersectShape(const RenderShapeData &shape, const glm::vec3 &O, const glm::vec3 &dir) {
    const glm::vec3 D(glm::mat3(shape.inv_ctm) * dir);
    switch (shape.primitive.type) {
    case PrimitiveType::PRIMITIVE_CUBE:
        return intersectCube(O, D, shape);
    case PrimitiveType::PRIMITIVE_SPHERE:
        return intersectSphere(O, D, shape);
    case PrimitiveType::PRIMITIVE_CYLINDER:
        return intersectCylinder(O, D, shape);
    case PrimitiveType::PRIMITIVE_CONE:
        return intersectCone(O, D, shape);
    default:
        return std::nullopt;
    }
}

// ------------------- Main Functions -------------------

std::optional<Hit> checkIntersection(const RayTraceScene &scene, const glm::vec3 &dir) {
    std::optional<Hit> closest;
    float closestT = std::numeric_limits<float>::infinity();
    const std::vector<glm::vec3> &obj_cams = scene.getObjCams();
    const std::vector<RenderShapeData> &shapes = scene.getShapes();

    for (const auto &shape : shapes) {
        auto hit = intersectShape(shape, obj_cams[shape.id], dir);
        if (hit && hit->t < closestT) {
            closestT = hit->t;
            closest = hit;
        }
    }
    return closest;
}

std::optional<Hit> traverseKDTree(
    const RayTraceScene &scene, const glm::vec4 &start, const glm::vec3 &dir, bool cam_ray,
    unsigned int &ray_id, std::vector<unsigned int> &last_visited
) {
    static const std::vector<glm::vec3> &obj_cams = scene.getObjCams();

    ray_id += 1;
    std::optional<Hit> closestHit;
    float closestT = std::numeric_limits<float>::infinity();

    std::stack<const KDNode *> stack;
    stack.push(scene.getKDRoot());

    float tmin, tmax;
    const KDNode *node;

    const glm::vec3 inv_dir = 1.f / dir;

    while (!stack.empty()) {
        node = stack.top();
        stack.pop();

        if (!intersectAABB(node->bounds, start, dir, inv_dir, tmin, tmax))
            continue;
        if (tmin >= closestT)
            continue; // This node is worth looking at

        if (node->isLeaf()) {
            for (const auto *shape : node->shapes) {
                if (last_visited[shape->id] == ray_id)
                    continue;
                last_visited[shape->id] = ray_id;

                const glm::vec3 obj_start = cam_ray ? obj_cams[shape->id] : shape->inv_ctm * start;
                const auto hit = intersectShape(*shape, obj_start, dir);
                if (hit && hit->t < closestT) {
                    closestT = hit->t;
                    closestHit = hit;
                }
            }
            continue;
        }

        // Determine which child to visit first based on ray direction
        const KDNode *firstChild = node->left.get();
        const KDNode *secondChild = node->right.get();

        if (dir[node->axis] < 0.f)
            std::swap(firstChild, secondChild);

        const float dist_to_split = (node->split - start[node->axis]) / dir[node->axis];
        if (secondChild && dist_to_split < closestT)
            stack.push(secondChild);
        if (firstChild)
            stack.push(firstChild);
    }

    return closestHit;
}

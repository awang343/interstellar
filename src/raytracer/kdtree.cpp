#include "kdtree.h"
#include "raytracescene.h"

#include <algorithm>
#include <future>


static const inline float surfaceArea(const AABB &box)
{
    glm::vec3 d = box.max - box.min;
    return 2.0f * (d.x * d.y + d.y * d.z + d.z * d.x);
}

const AABB computeBounds(const RenderShapeData &shape)
{
    // We can assume that all the shapes are bounded by the 1x1x1 cube
    // centered at the origin in object space
    AABB box;
    box.min = glm::vec3(std::numeric_limits<float>::infinity());
    box.max = glm::vec3(-std::numeric_limits<float>::infinity());
    int signs[2] = {-1, 1};
    for (int x : signs)
        for (int y : signs)
            for (int z : signs)
            {
                glm::vec4 corner = shape.ctm * glm::vec4(0.5f * x, 0.5f * y, 0.5f * z, 1.0f);
                box.min = glm::min(box.min, glm::vec3(corner));
                box.max = glm::max(box.max, glm::vec3(corner));
            }
    return box;
}

std::unique_ptr<KDNode> recurseKDTreeSAH(const std::vector<ShapeWithBounds> &shapes, int depth)
{
    auto node = std::make_unique<KDNode>();

    // Compute bounds of all shapes in this node
    node->bounds.min = glm::vec3(std::numeric_limits<float>::infinity());
    node->bounds.max = glm::vec3(-std::numeric_limits<float>::infinity());
    for (auto &swb : shapes)
    {
        node->bounds.min = glm::min(node->bounds.min, swb.bounds.min);
        node->bounds.max = glm::max(node->bounds.max, swb.bounds.max);
    }

    // Leaf condition
    if (shapes.size() <= 10 || depth > 24)
    {
        node->shapes.reserve(shapes.size());
        for (auto &swb : shapes)
            node->shapes.push_back(swb.shape);
        return node;
    }

    float bestCost = std::numeric_limits<float>::infinity();
    int bestAxis = -1;
    float bestSplit = 0.0f;

    float parentArea = surfaceArea(node->bounds);

    struct Event
    {
        float pos;
        int index;
        bool isStart;
    };
    std::vector<Event> events;
    events.reserve(shapes.size() * 2);

    // Sweep SAH along each axis
    for (int axis = 0; axis < 3; axis++)
    {
        events.clear();
        for (int i = 0; i < (int) shapes.size(); i++)
        {
            events.push_back({shapes[i].bounds.min[axis], i, true});
            events.push_back({shapes[i].bounds.max[axis], i, false});
        }

        std::sort(events.begin(), events.end(), [](const Event &a, const Event &b)
                  { return a.pos < b.pos || (a.pos == b.pos && a.isStart && !b.isStart); });

        std::vector<AABB> prefix(events.size()), suffix(events.size());
        std::vector<int> prefixCount(events.size()), suffixCount(events.size());

        // Build prefix bounds/counts
        AABB curLeft{glm::vec3(INFINITY), glm::vec3(-INFINITY)};
        int nLeft = 0;
        for (size_t i = 0; i < events.size(); i++)
        {
            if (events[i].isStart)
            {
                nLeft++;
                curLeft.min = glm::min(curLeft.min, shapes[events[i].index].bounds.min);
                curLeft.max = glm::max(curLeft.max, shapes[events[i].index].bounds.max);
            }
            prefix[i] = curLeft;
            prefixCount[i] = nLeft;
        }

        // Build suffix bounds/counts
        AABB curRight{glm::vec3(INFINITY), glm::vec3(-INFINITY)};
        int nRight = 0;
        for (int i = (int) events.size() - 1; i >= 0; i--)
        {
            if (!events[i].isStart)
            {
                nRight++;
                curRight.min = glm::min(curRight.min, shapes[events[i].index].bounds.min);
                curRight.max = glm::max(curRight.max, shapes[events[i].index].bounds.max);
            }
            suffix[i] = curRight;
            suffixCount[i] = nRight;
        }

        // Evaluate candidate splits
        for (size_t i = 0; i + 1 < events.size(); i++)
        {
            int nL = prefixCount[i];
            int nR = suffixCount[i + 1];
            if (nL == 0 || nR == 0)
                continue;

            float cost = 1.0f + (surfaceArea(prefix[i]) / parentArea) * nL +
                         (surfaceArea(suffix[i + 1]) / parentArea) * nR;

            if (cost < bestCost)
            {
                bestCost = cost;
                bestAxis = axis;
                bestSplit = 0.5f * (events[i].pos + events[i + 1].pos);
            }
        }
    }

    // Partition shapes
    std::vector<ShapeWithBounds> leftShapes, rightShapes;
    for (auto &swb : shapes)
    {
        if (swb.bounds.min[bestAxis] <= bestSplit)
            leftShapes.push_back(swb);
        if (swb.bounds.max[bestAxis] >= bestSplit)
            rightShapes.push_back(swb);
    }

    // Avoid degenerate splits
    if (leftShapes.size() == shapes.size() || rightShapes.size() == shapes.size())
    {
        node->shapes.reserve(shapes.size());
        for (auto &swb : shapes)
            node->shapes.push_back(swb.shape);
        return node;
    }

    // Recurse
    node->axis = bestAxis;
    node->split = bestSplit;

    if (depth < PARALLEL_DEPTH)
    {
        auto leftFuture =
            std::async(std::launch::async, [&] { return recurseKDTreeSAH(leftShapes, depth + 1); });
        node->right = recurseKDTreeSAH(rightShapes, depth + 1);
        node->left = leftFuture.get();
    }
    else
    {
        node->left = recurseKDTreeSAH(leftShapes, depth + 1);
        node->right = recurseKDTreeSAH(rightShapes, depth + 1);
    }

    return node;
}


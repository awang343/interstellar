#include "raytracescene.h"

struct ShapeWithBounds
{
    const RenderShapeData *shape;
    AABB bounds;
};

const AABB computeBounds(const RenderShapeData &shape);

std::unique_ptr<KDNode> recurseKDTreeSAH(const std::vector<ShapeWithBounds> &shapes, int depth = 0);

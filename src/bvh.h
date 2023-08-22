#pragma once
#include "aabb.h"
struct BVHNode {
    AABB box;
    int left_node_id;
    int right_node_id;
    int primitive_id;
};

struct BBoxWithID {
    AABB box;
    int id;
};

int construct_bvh(const std::vector<BBoxWithID> &boxes,
                  std::vector<BVHNode> &node_pool);

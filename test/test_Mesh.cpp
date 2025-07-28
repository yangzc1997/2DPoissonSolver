#include <gtest/gtest.h>
#include "../src/Mesh.h"

namespace Poisson {

TEST(MeshTest, RectangleMesh) {
    std::vector<bool> edges = {true, false, true, false}; // AB和BC是边界
    Mesh mesh(2.0, 1.0, 2, 2, "4", edges);
    mesh.generate_mesh();
    
    // 验证节点数量
    EXPECT_EQ(mesh.nodes.size(), 9); // (2+1)*(2+1)=9
    
    // 验证边界标记
    // CD边界 (上边)[注意顶点处X方向有优先级,先排AD,BC，再排AB，CD]
    EXPECT_EQ(mesh.nodes[6].is_edges, 0); // (0,2)
    EXPECT_EQ(mesh.nodes[7].is_edges, 0); // (2,1)
    EXPECT_EQ(mesh.nodes[8].is_edges, 3); // (2,2)
    
    // BC边界 (右边)
    EXPECT_EQ(mesh.nodes[2].is_edges, 3); // (2,0)
    EXPECT_EQ(mesh.nodes[5].is_edges, 3); // (2,1)
    EXPECT_EQ(mesh.nodes[8].is_edges, 3); // (2,2) - 注意：角点只标记一个边界
    
    // 验证单元数量
    EXPECT_EQ(mesh.elements.size(), 4); // 2x2网格有4个四边形单元
}

TEST(MeshTest, TriangleMesh) {
    std::vector<bool> edges = {true, true, false, false}; // AB和AD是边界
    Mesh mesh(1.0, 1.0, 1, 1, "3", edges);
    mesh.generate_mesh();
    
    // 验证节点数量
    EXPECT_EQ(mesh.nodes.size(), 4); // (1+1)*(1+1)=4
    
    // 验证边界标记
    // CD边界 (上边)
    EXPECT_EQ(mesh.nodes[2].is_edges, 2); // (0,1)
    EXPECT_EQ(mesh.nodes[3].is_edges, 0); // (1,1)
    
    // AB边界 (下边)
    EXPECT_EQ(mesh.nodes[0].is_edges, 2); // (0,0)
    EXPECT_EQ(mesh.nodes[1].is_edges, 1); // (0,1) - 角点只标记一个边界
    
    // 验证单元数量
    EXPECT_EQ(mesh.elements.size(), 2); // 1个矩形分成2个三角形
}

TEST(MeshTest, ComplexBoundary) {
    std::vector<bool> edges = {true, true, true, true}; // 所有边界
    Mesh mesh(1.0, 1.0, 1, 1, "3", edges);
    mesh.generate_mesh();
    
    // 验证所有节点都是边界
    for (const auto& node : mesh.nodes) {
        EXPECT_NE(node.is_edges, 0);
    }
    
    // 验证角点标记
    EXPECT_EQ(mesh.nodes[0].is_edges, 2); // 左下角 - AD边界
    EXPECT_EQ(mesh.nodes[1].is_edges, 3); // 右下角 - BC边界
    EXPECT_EQ(mesh.nodes[2].is_edges, 2); // 左上角 - AD边界
    EXPECT_EQ(mesh.nodes[3].is_edges, 3); // 右上角 - BC边界
}

} // namespace Poisson
#include <gtest/gtest.h>
#include "Mesh.h"
#include <cmath>
#include <iostream>

namespace Poisson {

TEST(MeshTest, RectangleMesh) {
    // 创建 2x2 矩形网格
    Mesh mesh(2.0, 1.0, 2, 2, "4");
    
    // 验证节点数量
    EXPECT_EQ(mesh.nodes.size(), 9); // (2+1)*(2+1)=9
    
    // 验证节点坐标
    // 第0行
    EXPECT_DOUBLE_EQ(mesh.nodes[0]->x, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[0]->y, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[1]->x, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[1]->y, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[2]->x, 2.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[2]->y, 0.0);
    
    // 第1行
    EXPECT_DOUBLE_EQ(mesh.nodes[3]->x, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[3]->y, 0.5);
    EXPECT_DOUBLE_EQ(mesh.nodes[4]->x, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[4]->y, 0.5);
    EXPECT_DOUBLE_EQ(mesh.nodes[5]->x, 2.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[5]->y, 0.5);
    
    // 第2行
    EXPECT_DOUBLE_EQ(mesh.nodes[6]->x, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[6]->y, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[7]->x, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[7]->y, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[8]->x, 2.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[8]->y, 1.0);
    
    // 验证单元数量
    EXPECT_EQ(mesh.elements.size(), 4); // 2x2网格有4个四边形单元
    
    // 验证单元连接
    // 第一个单元 (左下角)
    EXPECT_EQ(mesh.elements[0].nodePtrs[0]->id, 0);
    EXPECT_EQ(mesh.elements[0].nodePtrs[1]->id, 1);
    EXPECT_EQ(mesh.elements[0].nodePtrs[2]->id, 4);
    EXPECT_EQ(mesh.elements[0].nodePtrs[3]->id, 3);
    
    // 第二个单元 (右下角)
    EXPECT_EQ(mesh.elements[1].nodePtrs[0]->id, 1);
    EXPECT_EQ(mesh.elements[1].nodePtrs[1]->id, 2);
    EXPECT_EQ(mesh.elements[1].nodePtrs[2]->id, 5);
    EXPECT_EQ(mesh.elements[1].nodePtrs[3]->id, 4);
    
    // 第三个单元 (左上角)
    EXPECT_EQ(mesh.elements[2].nodePtrs[0]->id, 3);
    EXPECT_EQ(mesh.elements[2].nodePtrs[1]->id, 4);
    EXPECT_EQ(mesh.elements[2].nodePtrs[2]->id, 7);
    EXPECT_EQ(mesh.elements[2].nodePtrs[3]->id, 6);
    
    // 第四个单元 (右上角)
    EXPECT_EQ(mesh.elements[3].nodePtrs[0]->id, 4);
    EXPECT_EQ(mesh.elements[3].nodePtrs[1]->id, 5);
    EXPECT_EQ(mesh.elements[3].nodePtrs[2]->id, 8);
    EXPECT_EQ(mesh.elements[3].nodePtrs[3]->id, 7);
}

TEST(MeshTest, TriangleMesh) {
    // 创建 1x1 三角形网格
    Mesh mesh(1.0, 1.0, 1, 1, "3");
    
    // 验证节点数量
    EXPECT_EQ(mesh.nodes.size(), 4); // (1+1)*(1+1)=4
    
    // 验证节点坐标
    EXPECT_DOUBLE_EQ(mesh.nodes[0]->x, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[0]->y, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[1]->x, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[1]->y, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[2]->x, 0.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[2]->y, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[3]->x, 1.0);
    EXPECT_DOUBLE_EQ(mesh.nodes[3]->y, 1.0);
    
    // 验证单元数量
    EXPECT_EQ(mesh.elements.size(), 2); // 1个矩形分成2个三角形
    
    // 验证单元连接
    // 第一个三角形 (左下角)
    EXPECT_EQ(mesh.elements[0].nodePtrs[0]->id, 0);
    EXPECT_EQ(mesh.elements[0].nodePtrs[1]->id, 1);
    EXPECT_EQ(mesh.elements[0].nodePtrs[2]->id, 2);
    
    // 第二个三角形 (右上角)
    EXPECT_EQ(mesh.elements[1].nodePtrs[0]->id, 1);
    EXPECT_EQ(mesh.elements[1].nodePtrs[1]->id, 3);
    EXPECT_EQ(mesh.elements[1].nodePtrs[2]->id, 2);
}

TEST(MeshTest, LargeGrid) {
    // 创建 10x10 网格
    const int Nx = 10, Ny = 10;
    Mesh mesh(1.0, 1.0, Nx, Ny, "4");
    
    // 验证节点数量
    EXPECT_EQ(mesh.nodes.size(), (Nx+1)*(Ny+1));
    
    // 验证单元数量
    EXPECT_EQ(mesh.elements.size(), Nx*Ny);
    
    // 验证节点坐标
    double dx = 1.0/Nx;
    double dy = 1.0/Ny;
    for (int j = 0; j <= Ny; j++) {
        for (int i = 0; i <= Nx; i++) {
            int id = j*(Nx+1) + i;
            EXPECT_DOUBLE_EQ(mesh.nodes[id]->x, i*dx);
            EXPECT_DOUBLE_EQ(mesh.nodes[id]->y, j*dy);
        }
    }
}

TEST(MeshTest, DefaultMeshType) {
    // 测试默认网格类型
    Mesh mesh(1.0, 1.0, 1, 1, "invalid");
    
    // 默认应为三角形网格
    EXPECT_EQ(mesh.elements.size(), 2);
}

TEST(MeshTest, MixedMeshTypes) {
    // 测试不同网格类型
    Mesh rectMesh(1.0, 1.0, 2, 2, "4");
    Mesh triMesh(1.0, 1.0, 2, 2, "3");
    
    // 矩形网格应有4个单元
    EXPECT_EQ(rectMesh.elements.size(), 4);
    
    // 三角形网格应有8个单元 (2x2网格分成8个三角形)
    EXPECT_EQ(triMesh.elements.size(), 8);
}

TEST(MeshTest, NodeIds) {
    // 测试节点ID是否正确分配
    Mesh mesh(1.0, 1.0, 2, 2, "4");
    
    for (size_t i = 0; i < mesh.nodes.size(); i++) {
        EXPECT_EQ(mesh.nodes[i]->id, static_cast<int>(i));
    }
}

TEST(MeshTest, ElementNodes) {
    // 测试单元节点指针是否正确
    Mesh mesh(1.0, 1.0, 1, 1, "4");
    
    // 检查第一个单元
    const auto& element = mesh.elements[0];
    EXPECT_EQ(element.nodePtrs[0], mesh.nodes[0]);
    EXPECT_EQ(element.nodePtrs[1], mesh.nodes[1]);
    EXPECT_EQ(element.nodePtrs[2], mesh.nodes[3]);
    EXPECT_EQ(element.nodePtrs[3], mesh.nodes[2]);
}

} // namespace Poisson
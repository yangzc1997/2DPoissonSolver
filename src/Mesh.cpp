// Mesh.cpp
#include "Mesh.h"
#include "Core_Export.h"
#include <iostream>

namespace Poisson{

// --------------------- 生成网格 ------------------ //
Mesh::Mesh(double lx_, double ly_, int Nx_, int Ny_, const std::string& mesh_type_)
    : lx(lx_), ly(ly_), Nx(Nx_), Ny(Ny_),mesh_type(mesh_type_)
{
    // === 1. 生成节点 ===
    const int total_nodes = (Nx + 1) * (Ny + 1);
    nodes.reserve(total_nodes);
    
    double dx = lx / Nx;
    double dy = ly / Ny;
    for (int j = 0; j <= Ny; j++) {
        for (int i = 0; i <= Nx; i++) {
            int id = j * (Nx + 1) + i;
            nodes.emplace_back(std::make_shared<Node>(i * dx, j * dy, id));
        }
    }
    
    // === 2. 生成单元 ===
        // 根据网格类型生成单元
    if (mesh_type == "4" || mesh_type == "rectangle") {
        elements.reserve(Nx * Ny);  //预分配节点内存
        generate_rectangle_elements();
    } else if (mesh_type == "3" || mesh_type == "triangle") {
        elements.reserve(2 * Nx * Ny);  //预分配节点内存
        generate_triangle_elements();
    } else {
        std::cout << "未定义的网格类型: " << mesh_type << std::endl;
        std::cout << "默认使用三角形单元" << std::endl;
        elements.reserve(2 * Nx * Ny);
        generate_triangle_elements();  //补充其他网格类型
    }
}


/// @brief 生成四边形网格节点连接信息
void Mesh::generate_rectangle_elements() {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            int n1 = i * (Nx + 1) + j;
            int n2 = n1 + 1;
            int n3 = n1 + (Nx + 1);
            int n4 = n3 + 1;

            // 创建四边形单元
            elements.emplace_back(std::vector<std::shared_ptr<Node>>{
                nodes[n1], nodes[n2], nodes[n4], nodes[n3]});
        }
    }
}

/// @brief 生成三角形网格节点连接信息
void Mesh::generate_triangle_elements() {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            int n1 = i * (Nx + 1) + j;
            int n2 = n1 + 1;
            int n3 = n1 + (Nx + 1);
            int n4 = n3 + 1;

            // 一个矩形可以划分为两个三角形
            elements.emplace_back(std::vector<std::shared_ptr<Node>>{
                nodes[n1], nodes[n2], nodes[n3]
            });
            elements.emplace_back(std::vector<std::shared_ptr<Node>>{
                nodes[n2], nodes[n4], nodes[n3]
            });
        }
    }
}


/// @brief 打印网格信息
void Mesh::printInfo() const {
    std::cout << "网格信息:\n";
    std::cout << "  类型: " << (mesh_type == "rectangle" || mesh_type == "4" ? "四边形" : "三角形") << "\n";
    std::cout << "  网格密度: Nx=" << Nx << "  ;  Ny= " << Ny << "\n";
    std::cout << "  节点数: " << nodes.size() << "\n";
    std::cout << "  单元数: " << elements.size() << "\n";
    
    std::cout << "  【节点信息】: " << "\n";
    for (const auto& node : nodes) {
        std::cout << "NodePosition: (" << node->x << ", " << node->y << "); Node Id: " << node->id << std::endl;
    }
}

} //namespace Poisson

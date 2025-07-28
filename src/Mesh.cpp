// 网格生成
#include "Mesh.h"
#include "Core_Export.h"
#include <iostream>

namespace Poisson{

/// @brief 获取网格参数方便生成网格和连接信息
/// @param lx x方向区域长度
/// @param ly y方向区域长度
/// @param Nx x方向划分份数
/// @param Ny x方向划分份数
/// @param mesh_type 网格类型（三角形["3" or "triangle"]或四边形["4" or "rectangle"]）
/// @param edge_ABCD bool型vector，记录了AB,AD,BC,CD 边是否属于固定边界条件
Mesh::Mesh(double lx, double ly, int Nx, int Ny, const std::string& mesh_type, const std::vector<bool>&  edge_ABCD) 
    : lx(lx), ly(ly), Nx(Nx), Ny(Ny), mesh_type(mesh_type), edge_ABCD(edge_ABCD) {}


// --------------------- 生成网格 ------------------ //
/// @brief 生成网格
void Mesh::generate_mesh() {
    double dx = lx / Nx;
    double dy = ly / Ny;

    // 清空节点和单元
    nodes.clear();
    elements.clear();

    // 预分配内存
    nodes.reserve((Nx + 1) * (Ny + 1));

    // 生成节点
    for (int i = 0; i <= Ny; i++) {
        for (int j = 0; j <= Nx; j++) {
            int is_edge = 0;  // 0表示非边界或者自由边界；1-4分别表示AB, AD, BC, CD边

            // 标记边界条件
            if  (j == 0 && edge_ABCD[1]) {
                is_edge = 2;  // AD边界
            } else if (j == Nx && edge_ABCD[2]) {
                is_edge = 3;  // BC边界
            }  else if (i == 0 && edge_ABCD[0]) {
                is_edge = 1;  // AB边界
            } else if (i == Ny && edge_ABCD[3]) {
                is_edge = 4;  // CD边界
            }

            // 添加节点
            nodes.push_back(Node{j * dx, i * dy, is_edge});  
        }
    }

    // 根据网格类型生成单元
    if (mesh_type == "4") {
        elements.reserve(Nx * Ny);  //预分配节点内存
        generate_rectangle_elements();
    } else if (mesh_type == "3") {
        elements.reserve(2 * Nx * Ny);  //预分配节点内存
        generate_triangle_elements();
    } else {
        std::cout << "未定义的网格类型: " << mesh_type << std::endl;
        std::cout << "默认使用三角形单元" << std::endl;
        elements.reserve(2 * Nx * Ny);
        generate_rectangle_elements();  //补充其他网格类型
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
            elements.push_back({n1, n2, n4, n3});
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
            elements.push_back({n1, n2, n3});
            elements.push_back({n2, n4, n3});
        }
    }
}

/// @brief 生成复合网格节点连接信息
void Mesh::generate_complex_elements() {
    // TODO: 生成复合网格节点连接信息
}


/// @brief 打印网格信息
void Mesh::print_mesh() const {
    std::cout << "Nodes: " << nodes.size() << std::endl;
    for (const auto& node : nodes) {
        std::cout << "Node: (" << node.x << ", " << node.y << "), Edge type: " << node.is_edges << std::endl;
    }

    std::cout << "Mesh Type: " << mesh_type << std::endl;
    std::cout << "Elements: " << elements.size() << std::endl;
    int n_element_nodes = 0 ;
    for (const auto& elem : elements) {
        n_element_nodes++;
        std::cout << " number:" << n_element_nodes;
        std::cout << "    Element_nodes: ";
        for (const auto& node : elem) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Total Nodes: " << nodes.size() << std::endl;
    std::cout << "Mesh Type: " << mesh_type << std::endl;
}

} //namespace Poisson

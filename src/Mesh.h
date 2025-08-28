// Mesh.h
#ifndef MESH_H
#define MESH_H

#include "Core_Export.h"
#include <vector>


namespace Poisson{

/// @brief 单个网格节点（每个节点坐标）
struct Node {
    double x; ///< x坐标
    double y; ///< y坐标
    int id;    ///< 节点全局ID

    Node(double x_, double y_, int id_)
        : x(x_), y(y_), id(id_) {}
};

/// @brief 单个网格单元
struct Element {
   std::vector<Node*> nodePtrs; ///< 每个单元中节点的指针数组:长度为3或4
    
    Element(std::vector<Node*> ptrs)
        : nodePtrs(ptrs) {}
    
    /// @brief 获取单元节点数
    int get_num_nodes() const {
        return nodePtrs.size();
    }
};

/// @brief 网格类
struct Mesh {
public:
    Mesh()= default;

    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;
    Mesh(Mesh&&) = default;
    Mesh& operator=(Mesh&&) = default;
    ~Mesh() = default;

    std::vector<Node> nodes;         ///< 网格节点
    std::vector<Element> elements;    ///< 网格单元
    
    /// @brief 获取节点数量
    size_t get_num_nodes() const { return nodes.size(); }
    
    /// @brief 获取单元数量
    size_t get_num_elements() const { return elements.size(); }

    /// @brief 输出网格信息
    void print_mesh_information(bool nodesInf = false) const;

};

} // namespace Poisson

#endif // MESH_H

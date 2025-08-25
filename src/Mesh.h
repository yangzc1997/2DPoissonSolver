// Mesh.h
#ifndef MESH_H
#define MESH_H

#include "Core_Export.h"
#include <vector>
#include <string>
#include <memory>

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
   std::vector<std::shared_ptr<Node>> nodePtrs; ///< 每个单元中节点的指针数组:长度为3或4
    
    Element(std::vector<std::shared_ptr<Node>> ptrs)
        : nodePtrs(ptrs) {}
    
    /// @brief 获取单元节点数
    int numNodes() const {
        return nodePtrs.size();
    }

    /// @brief 获取单元中节点的指针 
    const std::vector<std::shared_ptr<Node>>& getNodes() const { return nodePtrs; }
};

/// @brief 网格类
struct POISSONCORE_API Mesh {
public:
    /// @brief 构造函数，直接生成网格
    /// @param lx 区域长度
    /// @param ly 区域宽度
    /// @param Nx x方向单元数
    /// @param Ny y方向单元数
    /// @param type 网格类型 ("triangle" 或 "rectangle")
    Mesh(double lx_, double ly_, int Nx_, int Ny_, const std::string& mesh_type_);
    
    // 禁用拷贝
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;
    ~Mesh() = default;

    std::vector<std::shared_ptr<Node>> nodes;         ///< 网格节点
    std::vector<Element> elements;    ///< 网格单元
    
    /// @brief 获取节点数量
    size_t num_nodes() const { return nodes.size(); }
    
    /// @brief 获取单元数量
    size_t num_elements() const { return elements.size(); }

    /// @brief 输出网格信息
    void printInfo() const;

private:
    double lx, ly;   ///< 网格的长度和宽度
    int Nx, Ny;      ///< 网格划分的数量
    std::string mesh_type="3"; ///< 网格类型: "rectangle=4" 或 "triangle=3"

    /// @brief 生成三角形和四边形单元信息
    void generate_triangle_elements();
    void generate_rectangle_elements();
};

}
#endif // MESH_H

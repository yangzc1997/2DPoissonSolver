// Mesh.h
#ifndef MESH_H
#define MESH_H

#include "Core_Export.h"
#include <vector>
#include <string>

namespace POISSONCORE_API Poisson{

/// @brief 节点结构体（每个节点坐标和是否为边界）
struct POISSONCORE_API Node {
    double x, y;  // 节点坐标
    int is_edges;  // 0表示非边界或自由边界，1-4分别表示AB,AD,BC,CD边界

    // 添加带参数的构造函数
    Node(double x_val, double y_val, int edge_type)
        : x(x_val), y(y_val), is_edges(edge_type) {}
};

/// @brief 构建有限元网格
class POISSONCORE_API Mesh {
public:
    //---------------------- 公有成员变量 ------------------//
    std::vector<Node> nodes;                            ///< 网格节点
    std::vector<std::vector<int>> elements;             ///< 最终单元（方便后续调用)

    /// @brief 构造函数,用于读取网格参数
    Mesh(double lx, double ly, int Nx, int Ny, const std::string& mesh_type, const std::vector<bool>&  edge_ABCD);
    
    /// @brief 根据网格类型生成网格
    void generate_mesh();

    /// @brief 打印网格信息
    void print_mesh() const;

private:
    // --------------------- 私有成员变量 ------------------- //
    double lx, ly;   ///< 网格的长度和宽度
    int Nx, Ny;      ///< 网格划分的数量
    std::string mesh_type; ///< 网格类型: "rectangle=4" 或 "triangle=3"
    std::vector<bool>  edge_ABCD; ///< 边界是否属于固定边界

    // -------------- 私有成员函数 ------------- //
    /// @brief 生成四边形网格
    void generate_rectangle_elements();

    /// @brief 生成三角形网格
    void generate_triangle_elements();

    /// @brief 生成复合网络
    void generate_complex_elements();

};


} //namespace Poisson

#endif // MESH_H

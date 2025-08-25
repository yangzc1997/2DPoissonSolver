// ReadInputData.h
// 读取输入文件
#ifndef READINPUTDATA_H
#define READINPUTDATA_H

#include "Core_Export.h"
#include <vector>
#include <stdexcept>
#include <cctype> // 用于isspace()
#include <array>
#include <iostream>
#include <filesystem>

namespace Poisson{

/// @brief 表示单条边界信息的结构体
struct SingleEdge {
    std::string value = "";  // 边界函数的表达式
    std::array<double, 2> range = {0.0, 1.0}; // 边界作用范围【min,max】

    SingleEdge() = default;

    SingleEdge(const std::string& val, const std::array<double, 2>& rng)
        : value(val), range(rng){}
};

/// @brief 长方形区域中各种边界条件集合的结构体
struct BoundaryConditionData {
    SingleEdge AB;  // 底边 (x轴方向)
    SingleEdge AD;  // 左边 (y轴方向)
    SingleEdge BC;  // 右边 (y轴方向)
    SingleEdge CD;  // 顶边 (x轴方向)

    BoundaryConditionData() = default;

    BoundaryConditionData(const std::string& ab, const std::array<double, 2>& ab_range,
                         const std::string& ad, const std::array<double, 2>& ad_range,
                         const std::string& bc, const std::array<double, 2>& bc_range,
                         const std::string& cd, const std::array<double, 2>& cd_range)
        : AB(ab, ab_range), 
          AD(ad, ad_range), 
          BC(bc, bc_range), 
          CD(cd, cd_range) {}
};


/// @brief 处理输入文件，读取计算参数
struct ReadInputData {

    // 网格参数
    double lx = 1.0;                   ///< 区域长度和
    double ly = 1.0;                   ///< 区域宽度
    int Nx = 10;                       ///< 网格划分数量(x方向)
    int Ny = 10;                       ///< 网格划分数量(y方向)
    std::string mesh_type = "triangle";///< 网格类型

    // 边界条件
    BoundaryConditionData bc;             ///< 边界条件

    // 函数定义
    std::string initial_guess = "1.0"; ///< 初始猜测
    std::string source = "0.0";         ///< 源函数
    std::string source_derivatives = "0.0"; ///< 源函数的导数

    // 求解器设置
    std::string solver_type = "sparse"; ///< "sparse" 或 "dense"
    double rel_tol = 1.0e-7;           ///< 牛顿法相对容差
    double abs_tol = 1.0e-7;           ///< 牛顿法绝对容差
    int max_iter = 100;                ///< 最大迭代次数
    
    // 输出设置
    std::string output_file = "./SolutionResults.vtk"; ///< 计算结果输出路径和文件名

    ReadInputData(const std::filesystem::path& jsonPath);
    ~ReadInputData() = default;

    /// @brief 输出读取到的计算参数
    void printParameters() const;
};

} // namespace Poisson

#endif  //READINPUTDATA_H
// ReadInputData.h
// 读取输入文件
#ifndef READINPUTDATA_H
#define READINPUTDATA_H

#include "Core_Export.h"
#include "BoundaryData.h"
#include <filesystem>
#include <string>
#include <fstream>

namespace Poisson{

/// @brief 处理输入文件，读取计算参数
struct ReadInputData {

    // 网格参数
    double lx = 1.0;                   ///< 区域长度和
    double ly = 1.0;                   ///< 区域宽度
    int Nx = 10;                       ///< 网格划分数量(x方向)
    int Ny = 10;                       ///< 网格划分数量(y方向)
    std::string mesh_type = "triangle";///< 网格类型

    // 边界条件
    BoundaryConditionInfo bc;             ///< 边界条件

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
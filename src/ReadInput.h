// 读取输入文件
#ifndef READINPUT_H
#define READINPUT_H

#include "Core_Export.h"
#include <vector>
#include <stdexcept>
#include <cctype> // 用于isspace()

namespace Poisson{

/// @brief 处理输入文件，读取计算参数
class POISSONCORE_API ReadInput {
public:
    // 成员变量用于存储输入参数
    double lx = 1.0;                   ///< 区域长度和
    double ly = 1.0;                   ///< 区域宽度
    int Nx = 10;                        ///< 网格划分数量(x方向)
    int Ny = 10;                        ///< 网格划分数量(y方向)
    std::string mesh_type = "triangle";     ///< 网格类型
    std::string uAB = "", uCD = "", uAD = "", uBC = ""; ///< 边界条件(空字符串表示自由边界)
    std::vector<bool> edge_ABCD;   ///< 分别表示AB, AD, BC, CD是否属于边界（固定边界）
    std::string initial_guess = "1.0";              ///< 初始猜测
    std::string source = "0.0";               ///< 源函数
    std::string source_derivatives = "0.0";   ///< 源函数的导数
    std::string solver_type = "sparse" ; // "sparse" 或 "dense"
    double rel_tol = 1.0e-7, abs_tol = 1.0e-7;  ///< 牛顿法容差
    int max_iter = 100;         ///< 最大迭代次数
    std::string output_file = "./SolutionResults.vtk"; ///< 计算结果输出路径和文件名
 
    // 构造函数
    ReadInput(const std::string& filename);

    /// @brief 输出读取到的计算参数
    void print_parameters();
};

} // namespace Poisson


#endif  //READINPUT_H
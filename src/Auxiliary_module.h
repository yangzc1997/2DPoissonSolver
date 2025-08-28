#ifndef AUXILIARY_MODULE_H
#define AUXILIARY_MODULE_H

#include <string>
#include <functional>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "Mesh.h"
#include "BoundaryData.h"
#include "PoissonSolver_FiniteElementData.h"

namespace Poisson {

namespace Auxiliary {

// 网格生成模块
namespace MeshGenerator {
    Mesh generate_mesh(double lx, double ly, int Nx, int Ny, const std::string& mesh_type);
}

// 表达式解析模块
namespace ExpressionParser {
    using fuxy = std::function<double(double, double, double)>;
    using fu = std::function<double(double)>;
    using fxy = std::function<double(double, double)>;

    fuxy parse_expression_fuxy(const std::string& expr_str);
    fu parse_expression_fu(const std::string& expr_str);
    fxy parse_expression_fxy(const std::string& expr_str);
}

// 求解参数预处理模块
namespace SolvingParameterPreprocessing {
    using FiniteElementDataSet = std::vector<FiniteElementData>;
    using vec_t = Eigen::VectorXd;
    using fuxy = std::function<double(double, double, double)>;
    using fxy = std::function<double(double, double)>;

    // 利用网格重新组织需要计算的有限单元
    FiniteElementDataSet getFiniteElementDataSet(const class Mesh& mesh);
    
    // 获得D边界NodeID
    std::vector<int> get_dirichlet_nodeIDs(
        const Mesh& mesh, 
        const BoundaryConditionInfo& bc,
        double lx, double ly);
        
    // 获得泊松方程的初始解
    vec_t get_initial_value_of_u(
        const Mesh& mesh, 
        const BoundaryConditionInfo& bc,
        const std::string& initial_guess_expr,
        double lx, double ly);
        
    // 解析（字符串型）源函数和及其导数为函数类型
    fuxy create_source_function(const std::string& expr_str);
    fuxy create_source_derivative(const std::string& expr_str);
}

} // namespace Auxiliary

} // namespace Poisson

#endif // AUXILIARY_MODULE_H
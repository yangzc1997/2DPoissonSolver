// Auxiliary_SolvingParameterPreprocessing.h

#ifndef AUXILIARY_SOLVING_PARAMETER_PREPROCESSING_H
#define AUXILIARY_SOLVING_PARAMETER_PREPROCESSING_H

#include "Mesh.h"
#include "BoundaryData.h"
#include "PoissonSolver_FiniteElementData.h"
#include "UsingAlias.h"
#include <vector>
#include <string>

namespace Poisson {

namespace Auxiliary {

// 求解参数预处理模块
namespace SolvingParameterPreprocessing {
    using FiniteElementDataSet = std::vector<FiniteElementData>;

    // 利用网格重新组织需要计算的有限单元
    FiniteElementDataSet generate_FiniteElementDataSet(const class Mesh& mesh);
    
    // 获得D边界NodeID
    std::vector<int> generate_dirichlet_nodeIDs(
        const FiniteElementDataSet& feData,  
        const BoundaryConditionInfo& bc,
        double lx, double ly);
        
    // 获得泊松方程的初始解
    vec_t generate_initial_value_of_u(
        const FiniteElementDataSet& feData, 
        const BoundaryConditionInfo& bc,
        const std::string& initial_guess_expr,
        double lx, double ly);
        
    // 解析（字符串型）源函数和及其导数为函数类型
    fuxy create_source_function(const std::string& expr_str);
    fuxy create_source_derivative(const std::string& expr_str);
} // namespace SolvingParameterPreprocessing

} // namespace Auxiliary

} // namespace Poisson

#endif // AUXILIARY_SOLVING_PARAMETER_PREPROCESSING_H
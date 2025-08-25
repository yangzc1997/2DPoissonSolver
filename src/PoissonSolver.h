// PoissonSolver.h
#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "Core_Export.h"
#include "ElementCalculator.h"
#include "ExpressionEvaluator.h"
#include "Mesh.h"
#include "ReadInputData.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <memory>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace Poisson {

    using smat_t = Eigen::SparseMatrix<double>;
    using vec_t = Eigen::VectorXd;
    using mat_t = Eigen::MatrixXd;

/// @class PoissonSolverBase
/// @brief 泊松求解器
class POISSONCORE_API PoissonSolver {
public:
    PoissonSolver(const ReadInputData& read_input, const Mesh& mesh,const vec_t& u0, 
                  const std::vector<int>& DirichletNodeIDs);
    ~PoissonSolver() = default;

    // 牛顿迭代求解
    bool solveByNewtonMethod();

    // 输出相关的函数
    void output_results() const;
    void print_results(int max_display = 20) const;

private:
    const ReadInputData& read_input; ///< 输入参数
    const Mesh& mesh;           ///< 网格数据
    vec_t u;                    ///< 解向量
    vec_t f;                    ///< 载荷向量
    smat_t K;                   ///< 刚度矩阵
    std::vector<int> dirichlet_nodes; ///< Dirichlet边界节点集合
    
    // 源函数及其导数
    ExpressionEvaluator func_source;
    ExpressionEvaluator func_source_derivatives;
    std::unique_ptr<ElementCalculator> element_calculator;
    
    int integration_order = 7; ///< 积分阶数

    // 全局向量和矩阵的计算与组装
    void calAndAssembleGlobalSystem();

    // 单元向量和矩阵组装
    void applyBoundaryCondition();
    
    // 输出迭代信息
    void Iter_print(int iter, double step_abs_tol, double step_rel_tol, const std::string& header) const;
};

} // namespace Poisson

#endif // POISSON_SOLVER_H
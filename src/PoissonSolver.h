// PoissonSolver.h
#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "FiniteElementCalculator.h"
#include "Mesh.h"
#include "PoissonSolver_FiniteElementData.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <memory>
#include <functional>
#include <iomanip>
#include "UsingAlias.h"
#include "Auxiliary_Integrate2DLegendre.h"

namespace Poisson {

    using FiniteElementDataSet = std::vector<FiniteElementData>;
    using GaussPoints = std::vector<Auxiliary::Integrate2DLegendre::GaussPoint>;

/// @class PoissonSolverBase
/// @brief 泊松求解器
class PoissonSolver {
public:
    PoissonSolver(
            const FiniteElementDataSet& feDataSet_, const vec_t& u_, 
            const std::vector<int>& dirichlet_nodes_,
            const fuxy& source_func, const fuxy& source_deriv_func,
            int max_iter_, double rel_tol_,
            double abs_tol_, const std::string& output_file
        );

    PoissonSolver(const PoissonSolver&) = delete;
    PoissonSolver& operator=(const PoissonSolver&) = delete;
    ~PoissonSolver() = default;

    // 牛顿迭代求解
    bool solveByNewtonMethod();

    // 输出相关的函数
    void output_results() const;
    void print_results(int max_display = 20) const;

private:
    const FiniteElementDataSet& feDataSet;
    vec_t u;  
    const std::vector<int> dirichlet_nodes; ///< Dirichlet边界节点集合,
    fuxy sourceFunc;
    fuxy sourceDerivFunc;
    int   max_iter;
    double rel_tol;
    double abs_tol;
    const std::string output_file;

    vec_t f;                    ///< 载荷向量
    smat_t K;                   ///< 刚度矩阵

    std::unique_ptr<FiniteElementCalculator> element_calculator = nullptr;

    // 全局向量和矩阵的计算与组装
    void calAndAssembleGlobalSystem(const GaussPoints& triangleGaussPoints, 
    const GaussPoints& rectangleGaussPoints);

    // 单元向量和矩阵组装
    void applyBoundaryCondition();
    
    // 输出迭代信息
    void Iter_print(int iter, double step_abs_tol, double step_rel_tol, const std::string& header) const;
};

} // namespace Poisson

#endif // POISSON_SOLVER_H
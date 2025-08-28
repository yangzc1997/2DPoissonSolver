// PoissonSolver.h
#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "Core_Export.h"
#include "FiniteElementCalculator.h"
#include "Mesh.h"
#include "PoissonSolver_FiniteElementData.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <memory>
#include <functional>
#include <iomanip>

namespace Poisson {

    using namespace Eigen;
    using smat_t = Eigen::SparseMatrix<double>;
    using vec_t = Eigen::VectorXd;
    using mat_t = Eigen::MatrixXd;
    using NodeCoords = std::vector<Eigen::Vector2d>;
    using fuxy = std::function<double(double u_val, double x_val, double y_val)>;

/// @class PoissonSolverBase
/// @brief 泊松求解器
class POISSONCORE_API PoissonSolver {
public:
    PoissonSolver(
            const std::string& mesh_type_,
            const Mesh& mesh_, const vec_t& u_, 
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
    const std::string mesh_type;
    const Mesh& mesh;
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

    std::unique_ptr<FiniteElementCalculator> element_calculator;
    
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
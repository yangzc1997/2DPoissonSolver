// PoissonSolverBase.h
#ifndef POISSON_SOLVER_BASE_H
#define POISSON_SOLVER_BASE_H

#include "Core_Export.h"
#include "ElementCalculator.h"
#include "ExpressionEvaluator.h"
#include "Mesh.h"
#include "ReadInput.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <memory>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace Poisson {

/// @class PoissonSolverBase
/// @brief 泊松求解器基类，提供公共接口和实现
class POISSONCORE_API PoissonSolverBase {
public:
    PoissonSolverBase(const ReadInput& read_input, const Mesh& mesh);
    virtual ~PoissonSolverBase() = default;
    
    // 公有接口
    virtual void assemble() = 0;
    virtual void solveNewton() = 0;
    virtual void initialize_u();
    
    // 公共实现（不依赖于矩阵类型）
    void output_results();
    void print_results();

protected:
    // 公共成员
    const ReadInput& read_input;
    const Mesh& mesh;
    
    // 边界条件和源项
    ExpressionEvaluator func_uAB;
    ExpressionEvaluator func_uCD;
    ExpressionEvaluator func_uAD;
    ExpressionEvaluator func_uBC;
    ExpressionEvaluator func_guess;
    ExpressionEvaluator func_source;
    ExpressionEvaluator func_source_derivatives;

    std::unique_ptr<ElementCalculator> element_calculator;
    Eigen::VectorXd u;           // 解向量
    Eigen::VectorXd f;           // 载荷向量
    
    // 公共辅助函数
    std::pair<std::vector<double>, std::vector<double>> get_element_coords(const std::vector<int>& element) const;
    std::pair<double, double> transform_coordinates(const std::vector<double>& x_coords, const std::vector<double>& y_coords, 
                                                    double p1, double p2, const std::vector<int>& element) const;
    void Iter_print(const int& iter, const double& step_abs_tol, const double& step_rel_tol) const;
    
    // 局部计算函数
    virtual double compute_localK(int i, int j,
                                 const std::vector<double>& x_coords,
                                 const std::vector<double>& y_coords,
                                 double p1, double p2,
                                 const std::vector<int>& element) const;
                                 
    virtual double compute_localF(int i,
                          const std::vector<double>& x_coords,
                          const std::vector<double>& y_coords,
                          double p1, double p2,
                          const std::vector<int>& element) const;
};

} // namespace Poisson

#endif // POISSON_SOLVER_BASE_H
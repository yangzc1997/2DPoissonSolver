// PoissonSolverDense.h
#ifndef POISSON_SOLVER_DENSE_H
#define POISSON_SOLVER_DENSE_H

#include "Core_Export.h"
#include "PoissonSolverBase.h"
//#include <eigen3/Eigen/Dense>

namespace Poisson {

/// @class PoissonSolverDense
/// @brief 稠密矩阵求解器
class POISSONCORE_API PoissonSolverDense : public PoissonSolverBase {
public:
    PoissonSolverDense(const ReadInput& read_input, const Mesh& mesh);
    
    void assemble() override;
    void solveNewton() override;

private:
    Eigen::MatrixXd K; // 刚度矩阵
    
};

} // namespace Poisson

#endif // POISSON_SOLVER_DENSE_H
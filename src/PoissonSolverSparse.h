// PossionSolverSparse.h
#ifndef POISSON_SOLVER_SPARSE_H
#define POISSON_SOLVER_SPARSE_H

#include "Core_Export.h"
#include "PoissonSolverBase.h"
// #include <eigen3/Eigen/Sparse>

namespace Poisson {

/// @class PoissonSolverSparse
/// @brief 稀疏矩阵求解器
class POISSONCORE_API PoissonSolverSparse : public PoissonSolverBase {
public:
    PoissonSolverSparse(const ReadInput& read_input, const Mesh& mesh);

    void assemble() override;
    void solveNewton() override;

private:
    Eigen::SparseMatrix<double> K; // 刚度矩阵
    
};

} // namespace Poisson

#endif // POISSON_SOLVER_SPARSE_H
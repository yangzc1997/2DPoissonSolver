// PoissonSolver_FiniteElementData.h

#ifndef POISSON_SOLVER_FINITE_ELEMENT_DATA_H
#define POISSON_SOLVER_FINITE_ELEMENT_DATA_H

#include <vector>
#include <eigen3/Eigen/Dense>

namespace Poisson {

// 这里面放的是有限单元的节点坐标及其编号
struct FiniteElementData {
    std::vector<Eigen::Vector2d> NodeCoords;
    std::vector<int> NodeIndexs;
};

} // namespace Poisson

#endif // POISSON_SOLVER_FINITE_ELEMENT_DATA_H


// PoissonSolver_FiniteElementData.h

#ifndef POISSON_SOLVER_FINITE_ELEMENT_DATA_H
#define POISSON_SOLVER_FINITE_ELEMENT_DATA_H

#include <vector>
#include <eigen3/Eigen/Dense>

namespace Poisson {

// 这里面放的是有限单元的节点坐标及其编号
struct FiniteElementData {
    std::vector<Eigen::Vector2d> NodeCoords; // N X 2 的数组[[x1,y1],[x2,y2],...]
    std::vector<int> DofIndexs;  // 自由度编号
};

} // namespace Poisson

#endif // POISSON_SOLVER_FINITE_ELEMENT_DATA_H


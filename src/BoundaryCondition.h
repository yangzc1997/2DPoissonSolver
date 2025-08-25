// BoundaryCondition.h

#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "ReadInputData.h"
#include "Core_Export.h"
#include "Mesh.h"
#include <vector>
#include <memory>
#include <eigen3/Eigen/Dense>

namespace Poisson {

/// @class BoundaryCondition
/// @brief 边界条件处理类
class POISSONCORE_API BoundaryCondition {
public:
    /// @brief 构造函数
    /// @param read_input 输入参数
    /// @param mesh 网格数据
    BoundaryCondition(const ReadInputData& read_input, const Mesh& mesh);

    ~BoundaryCondition() = default;
    
    /// @brief 获取初始解向量
    Eigen::VectorXd getInitialSolution() const { return initialSolution; }
    
    /// @brief 获取第一类（D）边界的NodeID
    const std::vector<int> getDirichletNodeID() const { return DirichletNodeID; }

private:
    const ReadInputData& readInput;
    const Mesh& mesh;
    Eigen::VectorXd initialSolution; // 初始解向量
    std::vector<int> DirichletNodeID;  // 第一类边界条件的NodeID

    /// @brief 处理边界条件并给出初始解
    void processBoundaryConditions();
};

} // namespace Poisson


#endif // BOUNDARY_CONDITION_H

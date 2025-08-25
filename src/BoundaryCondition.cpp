// BoundaryCondition.cpp
#include "BoundaryCondition.h"
#include "ExpressionEvaluator.h"
#include <cmath>
#include <iostream>

namespace Poisson {

BoundaryCondition::BoundaryCondition(const ReadInputData& read_input, const Mesh& mesh)
    : readInput(read_input), mesh(mesh) 
{
    processBoundaryConditions();
}

void BoundaryCondition::processBoundaryConditions() 
{
    // 将字符串解析为函数
    ExpressionEvaluator func_uAB(readInput.bc.AB.value);
    ExpressionEvaluator func_uBC(readInput.bc.BC.value);
    ExpressionEvaluator func_uAD(readInput.bc.AD.value);
    ExpressionEvaluator func_uCD(readInput.bc.CD.value);
    ExpressionEvaluator func_guess(readInput.initial_guess);

    // 标记每个Node的边界条件，并给定u的初始值
    const int num_nodes = mesh.nodes.size();
    constexpr double tolerance = 1e-8;

    initialSolution.resize(num_nodes);
    DirichletNodeID.clear();

    // 辅助函数：检查点是否在边界上
    auto isOnBoundary = [&](const Node& node, double fixedValue, bool isXAxis, 
                            const std::array<double, 2> range, 
                            const ExpressionEvaluator& func) -> bool 
    {
        // 检查坐标是否在固定值附近
        double coord = isXAxis ? node.x : node.y;
        if (std::fabs(coord - fixedValue) > tolerance) {
            return false;
        }
        
        // 检查是否在范围内
        double varyingCoord = isXAxis ? node.y : node.x;
        if (varyingCoord < range[0] || varyingCoord > range[1]) {
            return false;
        }
        
        return true;
    };

    // 边界定义数组
    struct BoundaryInfo {
        double fixedValue;                 // 固定坐标值
        bool isXAxis;                      // 是否是x轴固定
        const std::array<double, 2> range; // 作用范围
        ExpressionEvaluator& func; // 边界函数
        bool isDefined;           // 是否定义
    };
    
    std::array<BoundaryInfo, 4> boundaries = {
        BoundaryInfo{0.0, true, readInput.bc.AD.range, func_uAD, !readInput.bc.AD.value.empty()}, // 左边 (AD)
        BoundaryInfo{readInput.lx, true, readInput.bc.BC.range, func_uBC, !readInput.bc.BC.value.empty()},  // 右边 (BC)
        BoundaryInfo{0.0, false, readInput.bc.AB.range, func_uAB, !readInput.bc.AB.value.empty()}, // 底边 (AB)
        BoundaryInfo{readInput.ly, false, readInput.bc.CD.range, func_uCD, !readInput.bc.CD.value.empty()}  // 顶边 (CD)
    };

    for (int i = 0; i < num_nodes; i++) {
        const auto& node = *(mesh.nodes[i]);
        bool isBoundaryNode = false;
        
        // 检查所有边界
        for (const auto& boundary : boundaries) {
            if (!boundary.isDefined) continue;
            
            if (isOnBoundary(node, boundary.fixedValue, boundary.isXAxis, 
                             boundary.range, boundary.func)) 
            {
                DirichletNodeID.push_back(node.id);
                initialSolution(i) = boundary.func.evaluate_xy(node.x, node.y);
                isBoundaryNode = true;
                break; // 找到边界后退出循环
            }
        }
        
        // 如果不是边界节点
        if (!isBoundaryNode) {
            initialSolution(i) = func_guess.evaluate_xy(node.x, node.y);
        }
    }
}

} // namespace Poisson

// Auxiliary_SolvingParameterPreprocessing.cpp
#include "Auxiliary_SolvingParameterPreprocessing.h"
#include "Auxiliary_ExpressionParser.h"
#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <set>


namespace Poisson {

namespace Auxiliary {

namespace SolvingParameterPreprocessing {

    // 利用网格重新组织需要计算的有限单元
    FiniteElementDataSet generate_FiniteElementDataSet(const Mesh& mesh)
    {
        FiniteElementDataSet dataSet;
        
        // 预处理每个单元的数据
        for (auto element : mesh.elements) {
            FiniteElementData elemData;
            
            // 提取节点坐标和索引
            const int num_nodes = element.get_num_nodes();
            elemData.NodeCoords.reserve(num_nodes);
            elemData.DofIndexs.reserve(num_nodes);
            
            for (int i = 0; i < num_nodes; ++i) {
                const Node* node = element.nodePtrs[i];
                elemData.NodeCoords.emplace_back(node->x, node->y);
                elemData.DofIndexs.push_back(node->id);
            }

           dataSet.push_back(elemData);
        }
        
        return dataSet;
    
    }

    // 获取属于D边界的nodeID (自由度ID)
    std::vector<int> generate_dirichlet_nodeIDs(
        const FiniteElementDataSet& feData, 
        const BoundaryConditionInfo& bc,
        double lx, double ly) 
    {
        std::vector<int> dirichlet_nodeIds;
        constexpr double tolerance = 1e-8;
        
        // 边界定义数组
        struct BoundaryDef {
            double fixed_value;                 // 固定坐标值
            bool is_x_axis;                     // 是否是x轴固定
            std::array<double, 2> range;        // 作用范围
            bool is_defined;                    // 是否定义
        };
        
        std::array<BoundaryDef, 4> boundaries = {
            BoundaryDef{0.0, true, bc.AD.range, !bc.AD.value.empty()},   // 左边 (AD)
            BoundaryDef{lx, true, bc.BC.range, !bc.BC.value.empty()},    // 右边 (BC)
            BoundaryDef{0.0, false, bc.AB.range, !bc.AB.value.empty()},  // 底边 (AB)
            BoundaryDef{ly, false, bc.CD.range, !bc.CD.value.empty()}   // 顶边 (CD)
        };
        
        // 收集所有节点信息
        std::unordered_map<int, Eigen::Vector2d> allNodes;
        for (const auto& elemData : feData) {
            for (size_t i = 0; i < elemData.DofIndexs.size(); ++i) {
                int nodeId = elemData.DofIndexs[i];
                const auto& coord = elemData.NodeCoords[i];
                allNodes[nodeId] = coord;
            }
        }
        
        // 遍历所有节点
        for (const auto& [nodeId, coord] : allNodes) {
            double x = coord.x();
            double y = coord.y();
            
            for (const auto& boundary : boundaries) {
                if (!boundary.is_defined) continue;
                
                // 检查坐标是否在固定值附近
                double coord_val = boundary.is_x_axis ? x : y;
                if (std::fabs(coord_val - boundary.fixed_value) > tolerance) {
                    continue;
                }
                
                // 检查是否在范围内
                double varying_coord = boundary.is_x_axis ? y : x;
                if (varying_coord < boundary.range[0] || varying_coord > boundary.range[1]) {
                    continue;
                }
                
                // 找到边界节点
                dirichlet_nodeIds.push_back(nodeId);
                break; // 找到边界后退出循环
            }
        }
        
        return dirichlet_nodeIds;
    }
        

    // 获取解向量的初始值
    vec_t generate_initial_value_of_u(
        const FiniteElementDataSet& feData, 
        const BoundaryConditionInfo& bc,
        const std::string& initial_guess_expr,
        double lx, double ly)
    {
        // 创建表达式求值器
        auto initial_guess_func = ExpressionParser::parse_expression_fxy(initial_guess_expr);
        auto ab_func = ExpressionParser::parse_expression_fxy(bc.AB.value);
        auto bc_func = ExpressionParser::parse_expression_fxy(bc.BC.value);
        auto cd_func = ExpressionParser::parse_expression_fxy(bc.CD.value);
        auto ad_func = ExpressionParser::parse_expression_fxy(bc.AD.value);
        
        // 收集所有节点信息
        std::unordered_map<int, Eigen::Vector2d> allNodes;
        for (const auto& elemData : feData) {
            for (size_t i = 0; i < elemData.DofIndexs.size(); ++i) {
                int nodeId = elemData.DofIndexs[i];
                const auto& coord = elemData.NodeCoords[i];
                allNodes[nodeId] = coord;
            }
        }
        
        // 确定节点总数
        int maxNodeId = 0;
        for (const auto& [nodeId, _] : allNodes) {
            if (nodeId > maxNodeId) maxNodeId = nodeId;
        }
        const int numNodes = maxNodeId + 1;
        
        // 创建初始解向量
        vec_t u0 = vec_t::Zero(numNodes);
        constexpr double tolerance = 1e-8;
        
        // 边界定义数组
        struct BoundaryDef {
            double fixed_value;                 // 固定坐标值
            bool is_x_axis;                     // 是否是x轴固定
            std::array<double, 2> range;        // 作用范围
            fxy func;                           // 边界函数
            bool is_defined;                    // 是否定义
        };
        
        std::array<BoundaryDef, 4> boundaries = {
            BoundaryDef{0.0, true, bc.AD.range, ad_func, !bc.AD.value.empty()},   // 左边 (AD)
            BoundaryDef{lx, true, bc.BC.range, bc_func, !bc.BC.value.empty()},    // 右边 (BC)
            BoundaryDef{0.0, false, bc.AB.range, ab_func, !bc.AB.value.empty()},  // 底边 (AB)
            BoundaryDef{ly, false, bc.CD.range, cd_func, !bc.CD.value.empty()}   // 顶边 (CD)
        };
        
        // 遍历所有节点
        for (const auto& [nodeId, coord] : allNodes) {
            double x = coord.x();
            double y = coord.y();
            bool is_boundary_node = false;
            
            for (const auto& boundary : boundaries) {
                if (!boundary.is_defined) continue;
                
                // 检查坐标是否在固定值附近
                double coord_val = boundary.is_x_axis ? x : y;
                if (std::fabs(coord_val - boundary.fixed_value) > tolerance) {
                    continue;
                }
                
                // 检查是否在范围内
                double varying_coord = boundary.is_x_axis ? y : x;
                if (varying_coord < boundary.range[0] || varying_coord > boundary.range[1]) {
                    continue;
                }
                
                // 边界节点：使用边界条件
                u0(nodeId) = boundary.func(x, y);
                is_boundary_node = true;
                break;
            }
            
            // 内部节点：使用初始猜测
            if (!is_boundary_node) {
                u0(nodeId) = initial_guess_func(x, y);
            }
        }
        
        return u0;
    }


    // 生成源函数的表达式
    fuxy create_source_function(const std::string& expr_str) {
        return ExpressionParser::parse_expression_fuxy(expr_str);
    }

    // 生成源函数导数的表达式
    fuxy create_source_derivative(const std::string& expr_str) {
        return ExpressionParser::parse_expression_fuxy(expr_str);
    }

}  // SolvingParameterPreprocessing

} // namespace Auxiliary

} // namespace Poisson

// Auxiliary_module.cpp
#include "Auxiliary_module.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <array>
#include <memory>
#include "../lib/exprtk.hpp" 

namespace Poisson {

namespace Auxiliary {


namespace MeshGenerator{

    Mesh generate_mesh(double lx, double ly, int Nx, int Ny, const std::string& mesh_type)
    {
        Mesh mesh;

        // 生成节点
        const int total_nodes = (Nx + 1) * (Ny + 1);
        mesh.nodes.reserve(total_nodes);
        
        double dx = lx / Nx;
        double dy = ly / Ny;
        for (int j = 0; j <= Ny; j++) {
            for (int i = 0; i <= Nx; i++) {
                int id = j * (Nx + 1) + i;
                mesh.nodes.emplace_back(i * dx, j * dy, id);
            }
        }
        
        // 生成单元
        if (mesh_type == "4" || mesh_type == "rectangle") {
            // 生成四边形单元
            mesh.elements.reserve(Nx * Ny);  // 预分配单元内存
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int n1 = j * (Nx + 1) + i;
                    int n2 = n1 + 1;
                    int n3 = n1 + (Nx + 1);
                    int n4 = n3 + 1;

                    mesh.elements.emplace_back(std::vector<Node*>{
                        &mesh.nodes[n1], &mesh.nodes[n2], &mesh.nodes[n4], &mesh.nodes[n3]
                    });
                }
            }
        } else {
            // 生成三角形单元
            mesh.elements.reserve(2 * Nx * Ny);
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int n1 = j * (Nx + 1) + i;
                    int n2 = n1 + 1;
                    int n3 = n1 + (Nx + 1);
                    int n4 = n3 + 1;

                    // 一个矩形可以划分为两个三角形
                    mesh.elements.emplace_back(std::vector<Node*>{
                        &mesh.nodes[n1], &mesh.nodes[n2], &mesh.nodes[n3]
                    });
                    mesh.elements.emplace_back(std::vector<Node*>{
                        &mesh.nodes[n2], &mesh.nodes[n4], &mesh.nodes[n3]
                    });
                }
            }
        }

        return mesh;
    }

}  // namespace MeshGenerator


namespace ExpressionParser{
    using fuxy = std::function<double(double, double, double)>;
    using fu = std::function<double(double)>;
    using fxy = std::function<double(double, double)>;

    fuxy parse_expression_fuxy(const std::string& expr_str){
        typedef double T;
        exprtk::symbol_table<T> symbol_table;
        exprtk::expression<T> expression;
        exprtk::parser<T> parser;
        
        // 需要单独给u,x,y变量加入到堆之中
        struct ExpressionVars {
            T u = 0.0;
            T x = 0.0;
            T y = 0.0;
        };
        // 这里需要使用c++智能指针自动管理内存，因为我自己很难确定什么时候该释放掉这里的内存
        auto vars = std::make_shared<ExpressionVars>();  

        // 添加变量绑定
        symbol_table.add_variable("u", vars->u);
        symbol_table.add_variable("x", vars->x);
        symbol_table.add_variable("y", vars->y);
        symbol_table.add_constants(); // 添加常量如 pi, e

        // 注册符号表并编译表达式
        expression.register_symbol_table(symbol_table);
        if (!parser.compile(expr_str, expression)) {
            if (expr_str != ""){
                std::cerr << "警告: 表达式: { " << expr_str << " } 解析失败\n";
                std::cerr << "提示: 只能使用 u, x, y 作为变量，支持标准数学函数\n";
            }
            return [](double, double, double) { 
                return std::numeric_limits<double>::quiet_NaN(); 
            };
        }
        
        // 返回可执行的函数
        return [expression, vars](double u_val, double x_val, double y_val) mutable {
            vars->u = u_val;
            vars->x = x_val;
            vars->y = y_val;
            return expression.value();
        };
    }

    fu parse_expression_fu(const std::string& expr_str) {
        auto func = parse_expression_fuxy(expr_str);
        return [func](double u) { return func(u, 0.0, 0.0); };
    }

    fxy parse_expression_fxy(const std::string& expr_str) {
        auto func = parse_expression_fuxy(expr_str);
        return [func](double x, double y) { return func(0.0, x, y); };
    }
}  // namespace ExpressionParser


namespace SolvingParameterPreprocessing {

    // 利用网格重新组织需要计算的有限单元
    FiniteElementDataSet getFiniteElementDataSet(const class Mesh& mesh)
    {
        FiniteElementDataSet dataSet;
        
        // 预处理每个单元的数据
        for (auto element : mesh.elements) {
            FiniteElementData elemData;
            
            // 提取节点坐标和索引
            const int num_nodes = element.get_num_nodes();
            elemData.NodeCoords.reserve(num_nodes);
            elemData.NodeIndexs.reserve(num_nodes);
            
            for (int i = 0; i < num_nodes; ++i) {
                const Node* node = element.nodePtrs[i];
                elemData.NodeCoords.emplace_back(node->x, node->y);
                elemData.NodeIndexs.push_back(node->id);
            }
        }
        
        return dataSet;
    
    }

    // 获取属于D边界的nodeID
    std::vector<int> get_dirichlet_nodeIDs(
        const Mesh& mesh, 
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
        
        // 遍历所有节点
        for (const auto& node : mesh.nodes) {
            for (const auto& boundary : boundaries) {
                if (!boundary.is_defined) continue;
                
                // 检查坐标是否在固定值附近
                double coord = boundary.is_x_axis ? node.x : node.y;
                if (std::fabs(coord - boundary.fixed_value) > tolerance) {
                    continue;
                }
                
                // 检查是否在范围内
                double varying_coord = boundary.is_x_axis ? node.y : node.x;
                if (varying_coord < boundary.range[0] || varying_coord > boundary.range[1]) {
                    continue;
                }
                
                // 找到边界节点
                dirichlet_nodeIds.push_back(node.id);
                break; // 找到边界后退出循环
            }
        }
        
        return dirichlet_nodeIds;
    }
        

    // 获取解向量的初始值
    vec_t get_initial_value_of_u(
        const Mesh& mesh, 
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
        
        // 创建初始解向量
        vec_t u0(mesh.nodes.size());
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
        for (const auto& node : mesh.nodes) {
            bool is_boundary_node = false;
            
            for (const auto& boundary : boundaries) {
                if (!boundary.is_defined) continue;
                
                // 检查坐标是否在固定值附近
                double coord = boundary.is_x_axis ? node.x : node.y;
                if (std::fabs(coord - boundary.fixed_value) > tolerance) {
                    continue;
                }
                
                // 检查是否在范围内
                double varying_coord = boundary.is_x_axis ? node.y : node.x;
                if (varying_coord < boundary.range[0] || varying_coord > boundary.range[1]) {
                    continue;
                }
                
                // 边界节点：使用边界条件
                u0(node.id) = boundary.func(node.x, node.y);
                is_boundary_node = true;
                break;
            }
            
            // 内部节点：使用初始猜测
            if (!is_boundary_node) {
                u0(node.id) = initial_guess_func(node.x, node.y);
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

// PoissonSolverBase.cpp

#include "PoissonSolverBase.h"
#include "Core_Export.h"
#include <iostream>
#include <cmath>

namespace Poisson {

/// @brief 使用有限元方法求解二维非线性泊松方程的求解器（稠密矩阵）
/// @param read_input 输入文件读取类实例
/// @param mesh 网格生成类实例
PoissonSolverBase::PoissonSolverBase(const ReadInput& read_input, const Mesh& mesh)
    : read_input(read_input),
      mesh(mesh),
      func_uAB(read_input.uAB),
      func_uCD(read_input.uCD),
      func_uAD(read_input.uAD),
      func_uBC(read_input.uBC),
      func_guess(read_input.initial_guess),
      func_source(read_input.source),
      func_source_derivatives(read_input.source_derivatives),
      element_calculator(nullptr)
{
    // 根据网格类型选择计算器
    if (read_input.mesh_type == "4") {
        element_calculator = std::make_unique<RectangleCalculator>();
    } else {
        element_calculator = std::make_unique<TriangleCalculator>();
    }
    
    // 初始化解向量和载荷向量
    const int num_nodes = mesh.nodes.size();
    u.resize(num_nodes);
    f.resize(num_nodes);
    u.setZero();
    f.setZero();
}

void PoissonSolverBase::initialize_u() 
{
    for (int i = 0; i < mesh.nodes.size(); i++) {
        switch (mesh.nodes[i].is_edges) {
        case 1: u(i) = func_uAB.evaluate_xy(mesh.nodes[i].x, mesh.nodes[i].y); break;
        case 2: u(i) = func_uAD.evaluate_xy(mesh.nodes[i].x, mesh.nodes[i].y); break;
        case 3: u(i) = func_uBC.evaluate_xy(mesh.nodes[i].x, mesh.nodes[i].y); break;
        case 4: u(i) = func_uCD.evaluate_xy(mesh.nodes[i].x, mesh.nodes[i].y); break;
        default: u(i) = func_guess.evaluate_xy(mesh.nodes[i].x, mesh.nodes[i].y); break;
        }
    }
}

// 组装单元坐标
std::pair<std::vector<double>, std::vector<double>> 
PoissonSolverBase::get_element_coords(const std::vector<int>& element) const
{
    std::vector<double> x_coords, y_coords;
    for (int node_idx : element) {
        x_coords.push_back(mesh.nodes[node_idx].x);
        y_coords.push_back(mesh.nodes[node_idx].y);
    }
    return {x_coords, y_coords};
}


// 坐标变换辅助函数
std::pair<double, double> PoissonSolverBase::transform_coordinates(
    const std::vector<double>& x_coords,
    const std::vector<double>& y_coords,
    double p1, double p2,
    const std::vector<int>& element) const
{
    if (element.size() == 4) {
        return {
            ((x_coords[1]-x_coords[0])*p1 + (x_coords[1]+x_coords[0])) / 2.0,
            ((y_coords[2]-y_coords[0])*p2 + (y_coords[2]+y_coords[0])) / 2.0
        };
    } else {
        return {
            x_coords[0] + (x_coords[1]-x_coords[0])*p1 + (x_coords[2]-x_coords[0])*p2,
            y_coords[0] + (y_coords[1]-y_coords[0])*p1 + (y_coords[2]-y_coords[0])*p2
        };
    }
}



// 计算局部刚度矩阵
double PoissonSolverBase::compute_localK(int i, int j,
                        const std::vector<double>& x_coords,
                        const std::vector<double>& y_coords,
                        double p1, double p2,
                        const std::vector<int>& element) const
{
    // 计算位置坐标
    auto [x, y] = transform_coordinates(x_coords, y_coords, p1, p2, element);

    // 计算当前解值
    const double u_val = element_calculator->u_value(u, element, p1, p2);

    // 获取形函数和其梯度
    double N_i = element_calculator->shape_function(i, p1, p2);
    double N_j = element_calculator->shape_function(j, p1, p2);
     
    Eigen::Vector2d dN_i = element_calculator->shape_gradient(i, x_coords, y_coords, p1, p2);
    Eigen::Vector2d dN_j = element_calculator->shape_gradient(j, x_coords, y_coords, p1, p2);
    
    // 计算刚度矩阵元
    return -1.0 * (dN_i.dot(dN_j) + N_i * N_j * func_source_derivatives.evaluate_uxy(u_val, x, y));
}


// 计算局部载荷向量
double PoissonSolverBase::compute_localF(int i,
                          const std::vector<double>& x_coords,
                          const std::vector<double>& y_coords,
                          double p1, double p2,
                          const std::vector<int>& element) const 
{
    // 计算位置坐标
    auto [x, y] = transform_coordinates(x_coords, y_coords, p1, p2, element);
    
    // 计算当前解值和其梯度
    const double u_val = element_calculator->u_value(u, element, p1, p2);
    const Eigen::Vector2d grad_u = element_calculator->u_gradient(u, element, x_coords, y_coords, p1, p2);
    
    // 获取形函数和其梯度
    const double N_i = element_calculator->shape_function(i, p1, p2);
    const Eigen::Vector2d dN_i = element_calculator->shape_gradient(i, x_coords, y_coords, p1, p2);
    
    // 计算载荷向量元
    return N_i * func_source.evaluate_uxy(u_val, x, y) + grad_u.dot(dN_i);
}



// 输出迭代信息
void PoissonSolverBase::Iter_print(const int& iter, const double& step_abs_tol, const double& step_rel_tol) const
{
    // 创建表格头部
    if (iter == 0) {
        const std::string header(60, '-');
        std::cout << "" << header << "\n"
                  << std::left << std::setw(8) << "Iter" 
                  << std::setw(16) << "Abs. Error"
                  << std::setw(16) << "Rel. Error";
        std::cout << "\n" << header << std::endl;
    }
    
    // 输出迭代信息
    std::cout << std::left << std::setw(8) << iter 
              << std::setw(16) << std::scientific << std::setprecision(5) << step_abs_tol;

    if(iter > 0) {
        std::cout << std::setw(16) << step_rel_tol;  
    } else {
        std::cout << std::setw(16) << "-";
    }
    std::cout << std::endl;
}


// 输出计算结果到文件中
void PoissonSolverBase::output_results()
{
    const std::string& filename = read_input.output_file; 
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    std::cout << "\n计算结果储存文件为: " << filename << std::endl;

    // 1. 写入VTK文件头
    file << "# vtk DataFile Version 3.0\n";
    file << "Poisson Solver Output\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // 2. 写入点数据
    file << "POINTS " << mesh.nodes.size() << " double\n";
    for (const auto& point : mesh.nodes) {
        file << point.x << " " << point.y << " 0.0" << "\n";
    }

    // 3. 写入单元数据
    const auto& cells = mesh.elements;
    size_t total_cell_size = 0;
    for (const auto& cell : cells) {
        total_cell_size += cell.size() + 1;
    }
    
    file << "CELLS " << cells.size() << " " << total_cell_size << "\n";
    for (const auto& cell : cells) {
        file << cell.size();
        for (int node : cell) {
            file << " " << node;
        }
        file << "\n";
    }
    
    // 4. 写入单元类型
    file << "CELL_TYPES " << cells.size() << "\n";
    for (size_t i = 0; i < cells.size(); i++) {
        if (cells[i].size() == 4) {
            file << "9\n";
        } else {
            file << "5\n";
        }
    }

    // 5. 写入点数据（解u）
    file << "POINT_DATA " << u.size() << "\n";
    file << "SCALARS u double 1\n";
    file << "LOOKUP_TABLE default\n";
    file << std::scientific << std::setprecision(8);
    for (int i = 0; i < u.size(); i++) {
        file << u(i) << "\n";
    }
    
    file.close();
}


// 打印结果到文件中
void PoissonSolverBase::print_results()
{
    const int max_display = 20;
    const int step = (u.size() > max_display) ? u.size() / max_display : 1;
    
    std::cout << std::endl;
    std::cout << "================ 计算结果摘要 ================\n";
    std::cout << std::setw(10) << "节点ID" 
              << std::setw(12) << "X坐标" 
              << std::setw(12) << "Y坐标" 
              << std::setw(15) << "U(X,Y)" 
              << std::setw(15) << "边界状态\n";
    
    for (int i = 0; i < u.size(); i += step) {
        std::cout << std::setw(10) << i
                  << std::setw(12) << std::fixed << std::setprecision(4) << mesh.nodes[i].x
                  << std::setw(12) << mesh.nodes[i].y
                  << std::setw(15) << std::scientific << std::setprecision(6) << u(i)
                  << std::setw(15) << mesh.nodes[i].is_edges << std::endl;
    }
    
    if (step > 1) {
        std::cout << "... 已省略 " << (u.size() - max_display) << " 个节点结果 ...\n";
    }
    
    double max_u = u.maxCoeff();
    double min_u = u.minCoeff();
    std::cout << "\n解范围: [" << min_u << ", " << max_u << "]\n";
    std::cout << "边界点数量: " 
              << std::count_if(mesh.nodes.begin(), mesh.nodes.end(), 
                    [](const Node& node) { return node.is_edges != 0; })
              << " / " << mesh.nodes.size() << std::endl;
    std::cout << "============================================\n";
}

} // namespace Poisson
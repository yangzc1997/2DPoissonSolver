// PoissonSolver.cpp

#include "PoissonSolver.h"
#include "Core_Export.h"
#include <iostream>
#include <cmath>
#include <unordered_set>

namespace Poisson {

    using namespace Eigen;
    using smat_t = Eigen::SparseMatrix<double>;
    using vec_t = Eigen::VectorXd;
    using mat_t = Eigen::MatrixXd;
    using NodeCoords = std::vector<Eigen::Vector2d>;

/// @brief 使用有限元方法求解二维非线性泊松方程的求解器
/// @param read_input 输入文件读取类实例
/// @param mesh 网格生成类实例
PoissonSolver::PoissonSolver(const ReadInputData& read_input, const Mesh& mesh, const vec_t& u0, const std::vector<int>& DirichletNodeIDs)
    : read_input(read_input),
      mesh(mesh),
      u(u0),
      dirichlet_nodes(DirichletNodeIDs),
      func_source(read_input.source),
      func_source_derivatives(read_input.source_derivatives)
{
    // 根据网格类型选择计算器
    if (read_input.mesh_type == "4") {
        element_calculator = std::make_unique<RectangleCalculator>();
    } else {
        element_calculator = std::make_unique<TriangleCalculator>();
    }
    element_calculator->setIntegrationOrder(integration_order); // 设置单元积分阶数
    
    // 初始化载荷向量
    const int num_nodes = mesh.nodes.size();
    f.setZero(num_nodes);

    // 初始化刚度矩阵
    K.resize(num_nodes, num_nodes);
    // K.reserve(VectorXi::Constant(num_nodes, 8)); // 三角形网格点最多有8个近邻
    std::vector<Triplet<double>> triplets;
    triplets.reserve(num_nodes*8);
    for (const auto& element : mesh.elements) {
        const int n_elemNode = element.numNodes();
        for (int i = 0; i < n_elemNode; i++){
            const int i_idx = element.nodePtrs[i]->id;
            for (int j = 0; j < n_elemNode; j++){
                const int j_idx = element.nodePtrs[j]->id;
                triplets.emplace_back(i_idx, j_idx, 1.0);
            }
        }
    }
    K.setFromTriplets(triplets.begin(), triplets.end());
}

// 牛顿迭代求解
bool PoissonSolver::solveByNewtonMethod() 
{
    std::cout << "\n牛顿迭代求解中..." << std::endl;
    int iter = 0;
    vec_t delta_u = vec_t::Zero(u.size());
    const std::string header(60, '-');
    
    // 牛顿迭代过程
    while (iter <= read_input.max_iter) {
        calAndAssembleGlobalSystem();
        
        // 求增量
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(K);
        delta_u = solver.solve(f);
        u += delta_u;
        
        // 误差计算
        double step_abs_tol = f.norm();
        double step_rel_tol = delta_u.norm() / (read_input.rel_tol*1e-2 + u.norm());
        
        Iter_print(iter, step_abs_tol, step_rel_tol, header);
        
        if (step_rel_tol < read_input.rel_tol && step_abs_tol < read_input.abs_tol) {
            std::cout << "\n注意：牛顿迭代已收敛！迭代次数: " << iter << std::endl;
            std::cout << header << std::endl;
            return true;
        }
        iter++;
    }

    if (iter > read_input.max_iter) {
        std::cout << "\n注意：未收敛！！！达到设定最大迭代次数(" << read_input.max_iter << ")" << std::endl;  
        std::cout << header << std::endl;
    }

    return false;
}


// 全局矩阵和向量的计算与组装
void PoissonSolver::calAndAssembleGlobalSystem() 
{
    K.setZero();
    f.setZero();

    // 创建函数对象(方便解耦单元类和字符串转函数类)
    func_uxy sourceFunc = [&](double u_val, double x_val, double y_val) {
        return func_source.evaluate_uxy(u_val, x_val, y_val);
    };
    
    func_uxy sourceDerivFunc = [&](double u_val, double x_val, double y_val) {
        return func_source_derivatives.evaluate_uxy(u_val, x_val, y_val);
    };
    
    // 遍历所有单元
    for (const auto& element : mesh.elements) {
        const int num_nodes = element.numNodes();
        
        // 获取单元节点坐标和局部解向量
        NodeCoords coords;
        vec_t local_u(num_nodes);
        for (int i = 0; i < num_nodes; ++i) {
            const int node_idx = element.nodePtrs[i]->id;
            coords.push_back(Vector2d(element.nodePtrs[i]->x, element.nodePtrs[i]->y));
            local_u(i) = u(node_idx);
        }
        
        // 计算单元矩阵和向量
        auto [K_local, F_local] = element_calculator->computeElementMatrixAndVector(
            coords, sourceFunc, sourceDerivFunc, local_u
        );
        // mat_t K_local = element_calculator->computeStiffnessMatrix(coords, sourceDerivFunc, local_u);
        // vec_t F_local = element_calculator->computeLoadVector(coords, sourceFunc, local_u);
   
        // 组装到全局系统
        for (int i_local = 0;  i_local < num_nodes; ++i_local) {
            const int global_i = element.nodePtrs[i_local]->id;            
            f(global_i) += F_local(i_local);
            for (int j_local = 0; j_local < num_nodes; ++j_local) {
                const int global_j = element.nodePtrs[j_local]->id;
                K.coeffRef(global_i, global_j) += K_local(i_local, j_local);  
            }
        }
    }

    // 最后统一处理D边界
    applyBoundaryCondition();
}


// 施加边界条件
void PoissonSolver::applyBoundaryCondition(){
    // 标记需要修改的行同时设置载荷向量为0
    vec_t diag_mod = vec_t::Zero(K.rows());
    for (int index : dirichlet_nodes) {
        diag_mod(index) = 1.0;
        f(index) = 0.0;
    }
    
    // 应用边界条件
    for (int i = 0; i < K.outerSize(); ++i) {
        for (smat_t::InnerIterator it(K, i); it; ++it) {
            const int row = it.row();
            const int col = it.col();
            
            // 如果是边界节点所在的行或列
            if (diag_mod(row) > 0 || diag_mod(col) > 0) {
                if (row == col) {
                    // 对角元素：如果是边界节点设为1，否则不变
                    it.valueRef() = diag_mod(row) > 0 ? 1.0 : it.value();
                } else {
                    // 非对角元素：如果行或列是边界节点设为0
                    it.valueRef() = (diag_mod(row) > 0 || diag_mod(col) > 0) ? 0.0 : it.value();
                }
            }
        }
    }
}


// 输出迭代信息
void PoissonSolver::Iter_print(int iter, double step_abs_tol, double step_rel_tol,const std::string& header) const
{
    // 创建表格头部
    if (iter == 0) {
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
void PoissonSolver::output_results() const
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
    for (const auto& node_ptr : mesh.nodes) {
        const auto& node = *node_ptr; // 解引用智能指针
        file << node.x << " " << node.y << " 0.0" << "\n";
    }

    // 3. 写入单元数据
    const auto& cells = mesh.elements;
    size_t total_cell_size = 0;
    for (const auto& cell : cells) {
        total_cell_size += cell.numNodes() + 1;
    }
    
    file << "CELLS " << cells.size() << " " << total_cell_size << "\n";
    for (const auto& cell : cells) {
        file << cell.numNodes();
        for (int i = 0; i < cell.numNodes(); i++) {
            const int node_id = cell.nodePtrs[i]->id; // 使用节点ID
            file << " " << node_id;
        }
        file << "\n";
    }
    
    // 4. 写入单元类型
    file << "CELL_TYPES " << cells.size() << "\n";
    for (size_t i = 0; i < cells.size(); i++) {
        if (cells[i].numNodes() == 4) {
            file << "9\n"; // VTK_QUAD
        } else {
            file << "5\n"; // VTK_TRIANGLE
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
void PoissonSolver::print_results(int max_display) const 
{
    int step = (u.size() > max_display) ? u.size() / max_display : 1;
    
    std::cout << std::endl;
    std::cout << "================ 计算结果摘要 ================\n";
    std::cout << std::setw(10) << "节点ID" 
              << std::setw(12) << "X坐标" 
              << std::setw(12) << "Y坐标" 
              << std::setw(15) << "U(X,Y)" 
              << std::setw(15) << "边界状态\n";

    std::unordered_set<int> dirichlet_nodes_set(dirichlet_nodes.begin(), dirichlet_nodes.end()); // 方便后续查找

    for (int i = 0; i < u.size(); i += step) {
        const auto& node = *(mesh.nodes[i]);
        bool is_boundary = (dirichlet_nodes_set.find(i) != dirichlet_nodes_set.end());
        
        std::cout << std::setw(10) << i
                  << std::setw(12) << std::fixed << std::setprecision(4) << node.x
                  << std::setw(12) << node.y
                  << std::setw(15) << std::scientific << std::setprecision(6) << u(i)
                  << std::setw(15) << (is_boundary ? "Dirichlet" : "Free") << std::endl;
    }
    
    if (step > 1) {
        std::cout << "... 已省略 " << (u.size() - max_display) << " 个节点结果 ...\n";
    }
    
    double max_u = u.maxCoeff();
    double min_u = u.minCoeff();
    std::cout << "\n解范围: [" << min_u << ", " << max_u << "]\n";
    
    // 统计边界节点数量
    std::cout << "边界点数量: " << dirichlet_nodes.size()
              << " / " << mesh.nodes.size() << std::endl;
    std::cout << "============================================\n";
}

} // namespace Poisson
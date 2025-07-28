// PoissonSolverSparse.cpp

#include "PoissonSolverSparse.h"
#include "Core_Export.h"
#include <iostream>
#include <vector>

namespace Poisson {

PoissonSolverSparse::PoissonSolverSparse(const ReadInput& read_input, const Mesh& mesh)
    : PoissonSolverBase(read_input, mesh)
{
    K.resize(mesh.nodes.size(), mesh.nodes.size());
    K.reserve(Eigen::VectorXi::Constant(mesh.nodes.size(), 10));
}


void PoissonSolverSparse::assemble() 
{
    K.setZero();
    f.setZero();

    // 遍历所有单元
    for (const auto& element : mesh.elements) {
        auto coords = get_element_coords(element);
        auto& x_coords = coords.first;
        auto& y_coords = coords.second;
        
        // 遍历单元节点
        for (int i = 0; i < element.size(); i++) {
            const int ii = element[i];
            
            // 边界节点处理
            if (mesh.nodes[ii].is_edges != 0) {
                K.coeffRef(ii, ii) = 1.0;
                continue;
            }

            // 计算刚度矩阵部分
            for (int j = i; j < element.size(); j++) {
                const int jj = element[j];
                if (mesh.nodes[jj].is_edges == 0) {               
                    // 创建被积函数
                    auto funcK = [this, i, j, &element, &x_coords, &y_coords](double p1, double p2) {
                        return compute_localK(i, j, x_coords, y_coords, p1, p2, element);
                    };

                    // 积分计算局部刚度
                    const double local_K = element_calculator->integrate(x_coords, y_coords, funcK);

                    if (i == j) {
                        K.coeffRef(ii, ii) += local_K;
                    } else {
                        K.coeffRef(ii, jj) += local_K;
                        K.coeffRef(jj, ii) += local_K;
                    }
                }
            }
            
            // 计算载荷向量部分
            auto funcF = [this, i, &element, &x_coords, &y_coords](double p1, double p2) {
                return compute_localF(i, x_coords, y_coords, p1, p2, element);
            };
            
            f(ii) += element_calculator->integrate(x_coords, y_coords, funcF);
        }
    }
}


// 牛顿迭代求解
void PoissonSolverSparse::solveNewton() 
{
    std::cout << "\n牛顿迭代求解中... (稀疏矩阵)" << std::endl;
    int iter = 0;
    Eigen::VectorXd delta_u(u.size());
    delta_u.setZero();
    
    // 牛顿迭代过程
    while (iter <= read_input.max_iter) {
        assemble();
        
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(K);
        solver.factorize(K);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("矩阵分解失败");
        }
        delta_u = solver.solve(f);
        u += delta_u;
        
        // 误差计算
        double step_abs_tol = f.norm();
        double step_rel_tol = delta_u.norm() / (read_input.rel_tol*1e-2 + u.norm());
        
        Iter_print(iter, step_abs_tol, step_rel_tol);
        
        if (step_rel_tol < read_input.rel_tol && step_abs_tol < read_input.abs_tol) {
            std::cout << "\n注意：牛顿迭代已收敛！迭代次数: " << iter << std::endl;
            break;
        }
        iter++;
    }

    const std::string header(60, '-');
    if (iter > read_input.max_iter) {
        std::cout << "\n注意：未收敛！！！达到设定最大迭代次数(" << read_input.max_iter << ")" << std::endl;  
    }
    std::cout << header << std::endl;
}

} // namespace Poisson

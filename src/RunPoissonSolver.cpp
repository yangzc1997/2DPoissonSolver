// RunPoissonSolver.cpp
#include "Core_Export.h"
#include "Mesh.h"
#include "RunPoissonSolver.h"
#include "Auxiliary_MeshGenerator.h"
#include "Auxiliary_SolvingParameterPreprocessing.h"
#include "PoissonSolver_FiniteElementData.h"
#include "PoissonSolver.h"
#include <fstream>
#include <iostream>
#include <eigen3/Eigen/Dense>

namespace Poisson{

RunPoissonSolver::RunPoissonSolver(const std::filesystem::path& jsonPath)
    : m_jsonPath(jsonPath){}

bool RunPoissonSolver::readFromJson(){
    try{
        std::cout << "配置文件为：" << m_jsonPath << std::endl;
        m_readInputData = std::make_unique<ReadInputData>(m_jsonPath);
        m_readInputData->printParameters();
        return true;

    } catch(const std::exception& e){
        std::cerr << " 读取配置文件出现问题: " << e.what() << "\n";
        std::cout << " 提示: 请使用'-h'查看配置文件格式\n";
        return false;
    }
}


bool RunPoissonSolver::simulate() {
    try {
        m_timer.start(); // 开始计时

        // 网格生成
        std::cout << "\n网格生成中..." << std::endl;
        Mesh m_mesh = Auxiliary::MeshGenerator::generate_mesh(
                m_readInputData->lx, m_readInputData->ly, 
                m_readInputData->Nx, m_readInputData->Ny, 
                m_readInputData->mesh_type
            );
        m_mesh.print_mesh_information();

        // 牛顿迭代求解前的预处理
        std::cout << "\n求解预处理中..." << std::endl;

        // 从网格中生成有限元
        auto FiniteElements = Auxiliary::SolvingParameterPreprocessing::generate_FiniteElementDataSet(m_mesh);

        // 获取初始解
        auto u0 = Auxiliary::SolvingParameterPreprocessing::generate_initial_value_of_u(
            FiniteElements, m_readInputData->bc, m_readInputData->initial_guess,
            m_readInputData->lx, m_readInputData->ly
        );
        
        // 生成D边界有限元节点编号（自由度）
        auto dirichletNodes = Auxiliary::SolvingParameterPreprocessing::generate_dirichlet_nodeIDs(
            FiniteElements, m_readInputData->bc,
            m_readInputData->lx, m_readInputData->ly
        );

        // 源函数及其导数
        auto source_func = Auxiliary::SolvingParameterPreprocessing::create_source_function(m_readInputData->source);
        auto source_deriv_func = Auxiliary::SolvingParameterPreprocessing::create_source_derivative(m_readInputData->source_derivatives);

        // 创建求解器
        PoissonSolver m_solver(
            FiniteElements, u0, 
            dirichletNodes,
            source_func, source_deriv_func,
            m_readInputData->max_iter, m_readInputData->rel_tol,
            m_readInputData->abs_tol, m_readInputData->output_file
        );

        // 求解
        bool success = m_solver.solveByNewtonMethod();

        /// 统计计算时长
        m_timer.pause();
        m_timer.print_time();

        if (success) {
            // 输出结果
            m_solver.output_results();
            // m_solver.print_results();
            return true;
        }
        return false;

    } catch (const std::exception& e) {
        std::cerr << "模拟过程中出错: " << e.what() << std::endl;
        return false;
    }
}

} // namespace Poisson

// RunPoissonSolver.cpp
#include "Core_Export.h"
#include "RunPoissonSolver.h"
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

        // 创建网格
        std::cout << "\n开始网格生成..." << std::endl;
        m_mesh = std::make_unique<Mesh>(
                m_readInputData->lx, m_readInputData->ly, 
                m_readInputData->Nx, m_readInputData->Ny, 
                m_readInputData->mesh_type
            );
        // m_mesh->printInfo();

        // 处理边界条件
        std::cout << "\n处理边界条件..." << std::endl;
        m_boundaryCondition = std::make_unique<BoundaryCondition>(
            *m_readInputData, *m_mesh
        );

        // 获取初始解和边界节点
        const Eigen::VectorXd initialSolution = m_boundaryCondition->getInitialSolution();
        const std::vector<int> dirichletNodes = m_boundaryCondition->getDirichletNodeID();
    
        // 创建求解器
        m_solver = std::make_unique<PoissonSolver>(
            *m_readInputData, *m_mesh, initialSolution, dirichletNodes
        );

        // 求解
        bool success = m_solver->solveByNewtonMethod();

        /// 统计计算时长
        m_timer.pause();
        m_timer.print_time();

        if (success) {
            // 输出结果
            m_solver->output_results();
            // m_solver->print_results();
            return true;
        }
        return false;

    } catch (const std::exception& e) {
        std::cerr << "模拟过程中出错: " << e.what() << std::endl;
        return false;
    }
}

} // namespace Poisson

// main.cpp
#include "Core_Export.h"
#include "Help.h"
#include "ReadInput.h"
#include "Timer.h"
#include "Mesh.h"
#include "PoissonSolverSparse.h"
#include "PoissonSolverDense.h"

using namespace Poisson;

/// @brief  有限元方法求解二维泊松方程c++实现
/// @param argv [-i input.json]
/// @return 0表示成功
int main(int argc, char* argv[]) {
    // 输出开头信息
    std::string outstr(80,'=');
    std::cout << outstr << std::endl;
    std::cout << "欢迎使用有限单元法求解非线性泊松方程程序!" << std::endl;

    //====================== Step0 处理初始化信息 =======================//
    // 调用 check_help 函数来检查是否需要显示帮助信息                          
    std::string infilename = "input.json";
    if (argc > 1){
        infilename = check_help(argc, argv);
    } else { 
        std::cout << "使用默认输入文件: input.json" << std::endl; 
    } 

    // 使用计时器
    Timer timer;

    //=========================== step1 读取输入文件信息 =========================//
    // 读取输入文件
    std::cout << "配置文件读取中..." << std::endl;
    // 创建一个ReadInput对象并传入输入文件
    ReadInput read_input(infilename);
    read_input.print_parameters();  // 显示计算参数

    // 开始计时
    timer.start();

    //=========================== 网格剖分 =========================//
    // 网格和单元划分
    std::cout << "网格生成中..." << std::endl;
    Mesh meshes(read_input.lx, read_input.ly, read_input.Nx, read_input.Ny, read_input.mesh_type, read_input.edge_ABCD);
    meshes.generate_mesh();
    // meshes.print_mesh();
   
    //=============== 根据配置选择求解器 ==============//
    std::unique_ptr<PoissonSolverBase> poissonsolver;
    if (read_input.solver_type == "sparse") {
        std::cout << "使用稀疏矩阵求解器" << std::endl;
        poissonsolver = std::make_unique<PoissonSolverSparse>(read_input, meshes);
    } else {
        std::cout << "使用稠密矩阵求解器" << std::endl;
        poissonsolver = std::make_unique<PoissonSolverDense>(read_input, meshes);
    }

    //========== 求解过程 =============//
    poissonsolver->initialize_u();
    poissonsolver->solveNewton();

    //========== 输出计算结果 =============//
    poissonsolver->output_results();   // 输出到文件
    // poissonsolver->print_results();  // 输出到屏幕
    
    /// 统计计算时长
    timer.pause();
    timer.print_time();

    std::cout << outstr << std::endl;

    return 0;
}

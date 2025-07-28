#include "ReadInput.h"
#include "../lib/json.hpp"
#include "Core_Export.h"
#include <iostream>
#include <sstream>
#include <fstream>

using json = nlohmann::json;

namespace Poisson{

/// @brief 读取输入文件内容
/// @param filename 输入文件名
ReadInput::ReadInput(const std::string& filename) {
    try{
        // 打开输入文件
        std::cout << "配置文件为：" << filename << std::endl;
        std::ifstream infile(filename);
        if (!infile) {
            throw std::runtime_error("无法打开输入文件!");
        }

        json  config;
        try{
            infile >> config;
        } catch (const json::parse_error& e){
            throw std::runtime_error("JSON解析错误: " + std::string(e.what()));
        }

        // === 1. 解析网格参数 ===
        if (!config.contains("mesh") || !config["mesh"].is_object()) {
            throw std::runtime_error("缺失或无效的'mesh'配置块");
        }
        lx = config["mesh"].value("lx", 1.0);
        ly = config["mesh"].value("ly", 1.0);
        Nx = config["mesh"].value("Nx", 10);
        Ny = config["mesh"].value("Ny", 10);
        mesh_type = config["mesh"].value("type", "triangle");

        if (mesh_type == "triangle" || mesh_type == "3" || mesh_type == "Triangle"){
            mesh_type =  "3";
        } else if (mesh_type == "rectangle" || mesh_type == "4" || mesh_type == "Rectangle") {
            mesh_type =  "4";
        } else {
            mesh_type =  "3";  // 将来可能添加混合单元类型
        }

        // 网格参数验证
        if (lx <= 0 || ly <= 0) throw std::runtime_error("lx和ly必须大于0");
        if (Nx < 2 || Ny < 2) throw std::runtime_error("Nx和Ny必须至少为2");

        // === 2. 解析边界条件 ===
        if (config.contains("boundary_conditions") && config["boundary_conditions"].is_object()) {
                const auto& bc = config["boundary_conditions"];
                uAB = bc.value("AB", "");
                uAD = bc.value("AD", "");
                uBC = bc.value("BC", "");
                uCD = bc.value("CD", "");
            } else {
                throw std::runtime_error("Error: 无边界条件设置！！！"); 
            }
        // 设置边界激活标志
        edge_ABCD = {!uAB.empty(), !uAD.empty(), !uBC.empty(), !uCD.empty()};

        // === 3. 解析函数 ===
        if (!config.contains("functions") || !config["functions"].is_object()) {
            throw std::runtime_error("缺失或无效的'functions'配置块");
        }
        const auto& functions = config["functions"];
        initial_guess = functions.value("initial_guess", "1");  // 默认值=1
        source = functions.value("source", "0");                // 默认无源项
        source_derivatives = functions.value("source_derivatives", "0");

        // === 4. 解析求解器设置 ===
        // 解析求解器类型
        if (config.contains("solver_setting") && config["solver_setting"].is_object()) {
            const auto& solver = config["solver_setting"];
            solver_type = solver.value("solver_type", "sparse"); // 默认为稀疏求解器
        } else {
            solver_type = "dense";
        }
        if (solver_type == "sparse" || solver_type == "Sparse"){
            solver_type == "sparse";
        } else{
            solver_type = "dense";
        }
        // 容差设置
        if (config.contains("solver_setting") && config["solver_setting"].is_object()) {
            const auto& solver = config["solver_setting"];
            rel_tol = solver.value("rel_tol", rel_tol);
            abs_tol = solver.value("abs_tol", abs_tol);
            max_iter = solver.value("max_iter", max_iter);
        } else {
            std::cout << "   Warning：未设置求解容差，使用默认值！" <<  std::endl;
        }
        // 验证求解器参数
        if (rel_tol <= 0 || abs_tol <= 0) throw std::runtime_error("容差必须大于0");
        if (max_iter <= 0) throw std::runtime_error("最大迭代次数必须大于0");

        // === 5. 结果输出文件设置 ===
        if (config.contains("output_file") || config["output_file"].is_string()) {
            output_file = config["output_file"].get<std::string>();
        } else {
            output_file = "./SolutionResults.vtk";
        }

        infile.close();

    } 
    catch(const std::exception& e){
        std::cerr << " 读取输入配置出现问题: " << e.what() << "\n";
        std::cout << " 提示: 请使用'-h'查看配置文件格式\n";
    }
}


/// @brief 移除json文件中的注释



/// @brief 输出input文件的参数设置
void ReadInput::print_parameters() {
    auto safePrint = [](const std::string& title, auto value) {
        try {
            std::cout << "  " << title << ": " << value << std::endl;
        } catch (...) {
            std::cout << "  " << title << ": [INVALID]" << std::endl;
        }
    };

    std::cout << "----------------计算参数为 (请检查输入配置是否正确)----------------------\n";
    std::cout << "[网格参数]" << std::endl;
    safePrint("lx        ", lx);
    safePrint("ly        ", ly);
    safePrint("Nx        ", Nx);
    safePrint("Ny        ", Ny);
    safePrint("类型      ", (mesh_type == "3") ? "triangle":"rectangle");

    std::cout << "\n[边界条件]" << std::endl;
    safePrint("AB (上边) ", uAB.empty() ? "自由边界" : uAB);
    safePrint("BC (右边) ", uBC.empty() ? "自由边界" : uBC);
    safePrint("CD (下边) ", uCD.empty() ? "自由边界" : uCD);
    safePrint("AD (左边) ", uAD.empty() ? "自由边界" : uAD);

    std::cout << "\n[函数定义]" << std::endl;
    safePrint("初值表达式", initial_guess);
    safePrint("源项函数  ", source);
    safePrint("源项导数  ", source_derivatives);

    std::cout << "\n[求解器设置]" << std::endl;
    safePrint("求解方式  ", solver_type);
    safePrint("相对容差  ", rel_tol);
    safePrint("绝对容差  ", abs_tol);
    safePrint("最大迭代  ", max_iter);

    std::cout << "\n[输出设置]" << std::endl;
    safePrint("输出结果  ", output_file);
    std::cout << "------------------------------------------------------------------\n" << std::endl;
}

} //namespace Poisson
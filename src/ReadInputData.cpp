// ReadInputData.cpp
#include "ReadInputData.h"
#include "../lib/json.hpp"
#include "Core_Export.h"
#include <stdexcept>
#include <iostream>
#include <array>

using json = nlohmann::json;

namespace Poisson{

/// @brief 读取输入文件内容
/// @param filename 输入文件名
ReadInputData::ReadInputData(const std::filesystem::path& jsonPath) {
    try{
        // 打开输入文件
        std::ifstream infile(jsonPath);
        if (!infile) {
            throw std::runtime_error("无法打开输入文件!");
        }

        json  config;
        try{
            infile >> config;
        } catch (const json::parse_error& e){
            throw std::runtime_error("JSON解析错误: " + std::string(e.what()));
        }

        // === 解析网格参数 ===
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
        // 使用默认范围作为后备值
        const std::array<double, 2> default_x_range = {0.0, lx};
        const std::array<double, 2> default_y_range = {0.0, ly};

        if (config.contains("boundary_conditions") && config["boundary_conditions"].is_object()) {
            const auto& boundary_obj = config["boundary_conditions"];
            
            // 辅助函数：解析单一边界
            auto parse_edge = [this](const json& obj, const std::string& key, 
                                    const std::array<double, 2>& default_range,
                                    bool is_x_direction) 
                -> std::pair<std::string, std::array<double, 2>> 
            {
                std::array<double, 2> valid_range = default_range;
                double max_val = is_x_direction ? this->lx : this->ly;
                
                if (obj.contains(key)) {
                    if (obj[key].is_string()) {
                        return {obj[key].get<std::string>(), valid_range};
                    } else if (obj[key].is_object()) {
                        std::string val = obj[key].value("value", "");
                        std::array<double, 2> range = valid_range; // 初始化为默认值
                        
                        // 尝试直接读取为 array
                        if (obj[key].contains("range") && obj[key]["range"].is_array()) {
                            try {
                                // 直接尝试解析为 array<double, 2>
                                range = obj[key]["range"].get<std::array<double, 2>>();
                            } catch (const json::exception& e) {
                                std::cerr << "警告: " << key << "边界范围格式错误【正确格式[min,max]】，使用默认范围" << std::endl;
                            }
                        }
                        
                        // 确保范围值在有效区间内
                        // 调整最小值
                        if (range[0] < 0.0) {
                            std::cerr << "警告: " << key << "边界范围最小值小于0，自动调整为0" << std::endl;
                            range[0] = 0.0;
                        } else if (range[0] > max_val) {
                            std::cerr << "警告: " << key << "边界范围最小值大于" << max_val 
                                    << "，自动调整为" << max_val << std::endl;
                            range[0] = max_val;
                        }
                        
                        // 调整最大值
                        if (range[1] < 0.0) {
                            std::cerr << "警告: " << key << "边界范围最大值小于0，自动调整为0" << std::endl;
                            range[1] = 0.0;
                        } else if (range[1] > max_val) {
                            std::cerr << "警告: " << key << "边界范围最大值大于" << max_val 
                                    << "，自动调整为" << max_val << std::endl;
                            range[1] = max_val;
                        }
                        
                        // 确保最小值小于最大值
                        if (range[0] > range[1]) {
                            std::cerr << "警告: " << key << "边界范围无效（最小值大于最大值），使用默认范围" << std::endl;
                            range = valid_range;
                        }
                        
                        return {val, range};
                    }
                }
                return {"", valid_range}; // 默认为自由边界
            };

            // 解析各条边界
            auto [ab_val, ab_range] = parse_edge(boundary_obj, "AB", default_x_range, true); // x方向
            auto [ad_val, ad_range] = parse_edge(boundary_obj, "AD", default_y_range, false); // y方向
            auto [bc_val, bc_range] = parse_edge(boundary_obj, "BC", default_y_range, false); // y方向
            auto [cd_val, cd_range] = parse_edge(boundary_obj, "CD", default_x_range, true); // x方向
            
            // 创建边界条件对象
            bc = BoundaryConditionInfo(ab_val, ab_range, ad_val, ad_range, bc_val, bc_range, cd_val, cd_range);
        
        } else {
            bc = BoundaryConditionInfo("", {0.0, lx}, "", {0.0, ly}, "", {0.0, ly}, "", {0.0, lx});
            std::cout << "  警告：未配置边界条件，所有边界均为自由边界！" << std::endl;
        }

        // === 解析求解函数设置 ===
        if (!config.contains("functions") || !config["functions"].is_object()) {
            throw std::runtime_error("缺失或无效的'functions'配置块");
        }
        const auto& functions = config["functions"];
        initial_guess = functions.value("initial_guess", "1");  // 默认值=1
        source = functions.value("source", "0");                // 默认无源项
        source_derivatives = functions.value("source_derivatives", "0");

        // === 解析求解器设置 ===
        // 解析求解器类型
        if (config.contains("solver_setting") && config["solver_setting"].is_object()) {
            const auto& solver = config["solver_setting"];
            solver_type = solver.value("solver_type", "sparse"); // 默认为稀疏求解器
            if (solver_type == "sparse" || solver_type == "Sparse"){
                solver_type = "sparse";
            } else{
                solver_type = "dense";
            }
            
            // 容差设置
            rel_tol = solver.value("rel_tol", rel_tol);
            abs_tol = solver.value("abs_tol", abs_tol);
            max_iter = solver.value("max_iter", max_iter);
        } else {
            std::cout << "  警告：未设置求解器配置，使用默认值！" <<  std::endl;
        }
        // 验证求解器参数
        if (rel_tol <= 0 || abs_tol <= 0) throw std::runtime_error("容差必须大于0");
        if (max_iter < 0) throw std::runtime_error("最大迭代次数必须大于等于0");

        // === 结果输出文件设置 ===
        if (config.contains("output_file") && config["output_file"].is_string()) {
            output_file = config["output_file"].get<std::string>();
        } else {
            output_file = "./SolutionResults.vtk";
        }

        infile.close();
    } 
    catch(const std::exception& e){
        std::cerr << " 读取输入配置出现问题: " << e.what() << "\n";
        std::cout << " 提示: 请使用'-h'查看配置文件格式\n";
        throw;
    }
}

/// @brief 输出input文件的参数设置
void ReadInputData::printParameters() const {
    auto safePrint = [](const std::string& title, auto value) {
        try {
            std::cout << "  " << title << ": " << value << std::endl;
        } catch (...) {
            std::cout << "  " << title << ": [INVALID]" << std::endl;
        }
    };

    auto printRange = [](const std::array<double, 2>& range) -> std::string {
        if (range.size() != 2) return "[无效范围]";
        return "[" + std::to_string(range[0]) + ", " + std::to_string(range[1]) + "]";
    };

    std::cout << "----------------计算参数为 (请检查输入配置是否正确)----------------------\n";
    std::cout << "[网格参数]" << std::endl;
    safePrint("lx        ", lx);
    safePrint("ly        ", ly);
    safePrint("Nx        ", Nx);
    safePrint("Ny        ", Ny);
    safePrint("类型      ", (mesh_type == "3") ? "triangle":"rectangle");

    std::cout << "\n[边界条件]" << std::endl;
    safePrint("AB (下边) ", bc.AB.value.empty() ? "自由边界" : bc.AB.value + " " + printRange(bc.AB.range));
    safePrint("BC (右边) ", bc.BC.value.empty() ? "自由边界" : bc.BC.value + " " + printRange(bc.BC.range));
    safePrint("CD (上边) ", bc.CD.value.empty() ? "自由边界" : bc.CD.value + " " + printRange(bc.CD.range));
    safePrint("AD (左边) ", bc.AD.value.empty() ? "自由边界" : bc.AD.value + " " + printRange(bc.AD.range));

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
    std::cout << "------------------------------------------------------------------" << std::endl;
}

} //namespace Poisson
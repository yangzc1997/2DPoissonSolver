#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include "../lib/json.hpp"
#include "../src/ReadInputData.h" // 更新头文件引用

using json = nlohmann::json;
namespace fs = std::filesystem;

namespace Poisson {

TEST(ReadInputTest, JsonLibrary) {
    const fs::path filename = "temp_input.json";
    
    json j;
    j["mesh"]["lx"] = 1.0;
    j["mesh"]["ly"] = 2.0;
    j["mesh"]["Nx"] = 10;
    j["mesh"]["Ny"] = 20;
    j["mesh"]["type"] = "triangle";
    
    j["boundary_conditions"]["AB"] = {{"value", "x+y"}, {"range", {0.0, 1.0}}};
    j["boundary_conditions"]["CD"] = {{"value", "x*y"}, {"range", {0.5, 2.0}}};
    
    j["functions"]["initial_guess"] = "1.0";
    j["functions"]["source"] = "sin(x)*cos(y)";
    j["functions"]["source_derivatives"] = "cos(x)*cos(y)";
    
    j["solver_setting"]["solver_type"] = "sparse";
    j["solver_setting"]["rel_tol"] = 1e-8;
    j["solver_setting"]["abs_tol"] = 1e-10;
    j["solver_setting"]["max_iter"] = 100;
    
    j["output_file"] = "results.vtk";
    
    // 写入文件
    std::ofstream file(filename);
    file << j.dump(4);
    file.close();
    
    // 测试
    ReadInputData input(filename);

    // 验证网格参数
    EXPECT_DOUBLE_EQ(input.lx, 1.0);
    EXPECT_DOUBLE_EQ(input.ly, 2.0);
    EXPECT_EQ(input.Nx, 10);
    EXPECT_EQ(input.Ny, 20);
    EXPECT_EQ(input.mesh_type, "3"); // 三角形网格应为 "3"
    
    // 验证边界条件
    EXPECT_EQ(input.bc.AB.value, "x+y");
    EXPECT_EQ(input.bc.CD.value, "x*y");
    
    // 验证边界范围
    EXPECT_DOUBLE_EQ(input.bc.AB.range[0], 0.0);
    EXPECT_DOUBLE_EQ(input.bc.AB.range[1], 1.0);
    EXPECT_DOUBLE_EQ(input.bc.CD.range[0], 0.5);
    EXPECT_DOUBLE_EQ(input.bc.CD.range[1], 1.0);
    
    // 验证函数
    EXPECT_EQ(input.initial_guess, "1.0");
    EXPECT_EQ(input.source, "sin(x)*cos(y)");
    EXPECT_EQ(input.source_derivatives, "cos(x)*cos(y)");
    
    // 验证求解器设置
    EXPECT_EQ(input.solver_type, "sparse");
    EXPECT_DOUBLE_EQ(input.rel_tol, 1e-8);
    EXPECT_DOUBLE_EQ(input.abs_tol, 1e-10);
    EXPECT_EQ(input.max_iter, 100);
    
    // 验证输出文件
    EXPECT_EQ(input.output_file, "results.vtk");
    
    // 删除临时文件
    fs::remove(filename);
}

TEST(ReadInputTest, DefaultValues) {
    const fs::path filename = "temp_default.json";
    
    json j;
    j["mesh"]["lx"] = 1.0;
    j["mesh"]["ly"] = 1.0;
    j["mesh"]["Nx"] = 10;
    j["mesh"]["Ny"] = 10;
    
    // 添加必要的配置块
    j["functions"] = json::object();
    j["functions"]["initial_guess"] = "1.0";
    j["functions"]["source"] = "0.0";
    j["functions"]["source_derivatives"] = "0.0";

    // 写入文件
    std::ofstream file(filename);
    file << j.dump(4);
    file.close();
    
    // 测试
    ReadInputData input(filename);
    
    // 验证默认值
    EXPECT_EQ(input.mesh_type, "3"); // 默认三角形网格
    EXPECT_EQ(input.bc.AB.value, ""); // 默认无边界条件
    EXPECT_EQ(input.initial_guess, "1.0"); // 默认初始猜测
    EXPECT_EQ(input.source, "0.0"); // 默认无源项
    EXPECT_EQ(input.source_derivatives, "0.0"); // 默认无源项导数
    EXPECT_EQ(input.solver_type, "sparse"); // 默认稀疏求解器
    EXPECT_DOUBLE_EQ(input.rel_tol, 1.0e-7); // 默认相对容差
    EXPECT_DOUBLE_EQ(input.abs_tol, 1.0e-7); // 默认绝对容差
    EXPECT_EQ(input.max_iter, 100); // 默认最大迭代次数
    EXPECT_EQ(input.output_file, "./SolutionResults.vtk"); // 默认输出文件
    
    // 删除临时文件
    fs::remove(filename);
}

TEST(ReadInputTest, InvalidInput) {
    // 测试无效输入（文件不存在）
    EXPECT_THROW(ReadInputData("nonexistent.json"), std::runtime_error);
    
    // 创建无效网格文件
    const fs::path filename = "temp_invalid.json";
    
    json j;
    j["mesh"]["lx"] = -1.0; // 无效lx
    j["mesh"]["ly"] = 1.0;
    j["mesh"]["Nx"] = 1; // 无效Nx（至少为2）
    j["mesh"]["Ny"] = 10;
    
    // 写入文件
    std::ofstream file(filename);
    file << j.dump(4);
    file.close();
    
    // 测试
    EXPECT_THROW(ReadInputData("temp_invalid.json"), std::runtime_error);
    
    // 删除临时文件
    fs::remove(filename);
}

TEST(ReadInputTest, BoundaryRangeAdjustment) {
    const fs::path filename = "temp_range.json";
    
    json j;
    j["mesh"]["lx"] = 1.0;
    j["mesh"]["ly"] = 1.0;
    j["mesh"]["Nx"] = 10;
    j["mesh"]["Ny"] = 10;
    
    // 创建无效范围边界条件
    j["boundary_conditions"]["AB"] = {{"value", "x"}, {"range", {-1.0, 2.0}}}; // 超出范围
    j["boundary_conditions"]["BC"] = {{"value", "y"}, {"range", {0.5, 0.2}}}; // 最小值大于最大值
    
    // 添加必要的配置块
    j["functions"] = json::object();
    j["functions"]["initial_guess"] = "1.0";
    j["functions"]["source"] = "0.0";
    j["functions"]["source_derivatives"] = "0.0";
    
    // 写入文件
    std::ofstream file(filename);
    file << j.dump(4);
    file.close();
    
    // 测试
    ReadInputData input(filename);
    
    // 验证范围调整
    EXPECT_DOUBLE_EQ(input.bc.AB.range[0], 0.0); // 最小值调整为0
    EXPECT_DOUBLE_EQ(input.bc.AB.range[1], 1.0); // 最大值调整为1
    
    EXPECT_DOUBLE_EQ(input.bc.BC.range[0], 0.0); // 最小值调整为0
    EXPECT_DOUBLE_EQ(input.bc.BC.range[1], 1.0); // 最大值调整为1
    
    // 删除临时文件
    fs::remove(filename);
}

TEST(ReadInputTest, MeshTypeConversion) {
    const fs::path filename = "temp_mesh_type.json";
    
    json j;
    j["mesh"]["lx"] = 1.0;
    j["mesh"]["ly"] = 1.0;
    j["mesh"]["Nx"] = 10;
    j["mesh"]["Ny"] = 10;
    
    // 添加必要的配置块
    j["functions"] = json::object();
    j["functions"]["initial_guess"] = "1.0";
    j["functions"]["source"] = "0.0";
    j["functions"]["source_derivatives"] = "0.0";
    
    // 测试不同网格类型表示
    j["mesh"]["type"] = "rectangle";
    std::ofstream file1(filename);
    file1 << j.dump(4);
    file1.close();
    ReadInputData input1(filename);
    EXPECT_EQ(input1.mesh_type, "4");
    
    j["mesh"]["type"] = "4";
    std::ofstream file2(filename);
    file2 << j.dump(4);
    file2.close();
    ReadInputData input2(filename);
    EXPECT_EQ(input2.mesh_type, "4");
    
    j["mesh"]["type"] = "triangle";
    std::ofstream file3(filename);
    file3 << j.dump(4);
    file3.close();
    ReadInputData input3(filename);
    EXPECT_EQ(input3.mesh_type, "3");
    
    j["mesh"]["type"] = "3";
    std::ofstream file4(filename);
    file4 << j.dump(4);
    file4.close();
    ReadInputData input4(filename);
    EXPECT_EQ(input4.mesh_type, "3");
    
    j["mesh"]["type"] = "invalid";
    std::ofstream file5(filename);
    file5 << j.dump(4);
    file5.close();
    ReadInputData input5(filename);
    EXPECT_EQ(input5.mesh_type, "3"); // 无效类型默认为三角形
    
    // 删除临时文件
    fs::remove(filename);
}

} // namespace Poisson
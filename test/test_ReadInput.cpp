#include <gtest/gtest.h>
#include <fstream>
#include "../lib/json.hpp"
#include "../src/ReadInput.h"
using json = nlohmann::json;

namespace Poisson {

TEST(ReadInputTest, JsonLibrary) {
    const std::string filename = "temp_input.json";
    
    json j;
    j["mesh"]["lx"] = 1.0;
    j["mesh"]["ly"] = 2.0;
    j["mesh"]["Nx"] = 10;
    j["mesh"]["Ny"] = 20;
    j["mesh"]["type"] = "triangle";
    
    j["boundary_conditions"]["AB"] = "x+y";
    j["boundary_conditions"]["CD"] = "x*y";
    
    j["functions"]["initial_guess"] = "1.0";
    j["functions"]["source"] = "sin(x)*cos(y)";
    
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
    ReadInput input(filename);

    // 验证网格参数
    EXPECT_DOUBLE_EQ(input.lx, 1.0);
    EXPECT_DOUBLE_EQ(input.ly, 2.0);
    EXPECT_EQ(input.Nx, 10);
    EXPECT_EQ(input.Ny, 20);
    EXPECT_EQ(input.mesh_type, "3");
    
    // 验证边界条件
    EXPECT_EQ(input.uAB, "x+y");
    EXPECT_EQ(input.uCD, "x*y");
    EXPECT_TRUE(input.edge_ABCD[0]); // AB边界激活
    EXPECT_FALSE(input.edge_ABCD[1]); // AD边界未激活
    EXPECT_FALSE(input.edge_ABCD[2]); // BC边界未激活
    EXPECT_TRUE(input.edge_ABCD[3]); // CD边界激活
    
    // 验证函数
    EXPECT_EQ(input.initial_guess, "1.0");
    EXPECT_EQ(input.source, "sin(x)*cos(y)");
    
    // 验证求解器设置
    EXPECT_EQ(input.solver_type, "sparse");
    EXPECT_DOUBLE_EQ(input.rel_tol, 1e-8);
    EXPECT_DOUBLE_EQ(input.abs_tol, 1e-10);
    EXPECT_EQ(input.max_iter, 100);
    
    // 验证输出文件
    EXPECT_EQ(input.output_file, "results.vtk");
    
    // 删除临时文件
    std::remove(filename.c_str());
}

} // namespace Poisson
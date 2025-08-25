// main.cpp
#include "Core_Export.h"
#include "Help.h"
#include "RunPoissonSolver.h"
#include <iostream>
#include <string>
#include <filesystem>

using namespace Poisson;

int main(int argc, char* argv[]) {
    try {
        // 输出开头信息
        std::string outstr(80,'=');
        std::cout << outstr << std::endl;
        std::cout << "欢迎使用有限单元法求解非线性泊松方程程序!" << std::endl;

        // 处理初始化信息                        
        std::filesystem::path jsonPath = "input.json";
        if (argc > 1){
            jsonPath = check_help(argc, argv);
        } else {
            std::cout << "使用默认配置文件: input.json" << std::endl; 
        }

        // 创建求解任务
        RunPoissonSolver runTask(jsonPath);
        
        // 读取配置文件
        try {
            bool readJsonSucceed = runTask.readFromJson();

            if (!readJsonSucceed){
                std::cout << "ERROR: Loading .json file failed!" << std::endl;
                std::cout << outstr << std::endl;
                return 1; // Json读取失败
            }
        } catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            std::cout << outstr << std::endl;
            return 1; // Json读取失败
        }

        // 执行模拟 
        bool simulationSucceed = runTask.simulate();

        if (!simulationSucceed){
            std::cerr << "错误: 模拟失败!" << std::endl;
            std::cout << outstr << std::endl;
            return 2; // 求解失败
        }
    
        std::cout << "\n模拟已成功完成!" << std::endl;
        std::cout << outstr << std::endl;
        return 0; // 模拟成功
    
    } catch (const std::exception& e) {
        std::cerr << "\n 程序发生致命错误: " << e.what() << std::endl;
        return -1;  // 整体发生错误
    }
}

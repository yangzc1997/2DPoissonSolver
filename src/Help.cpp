#include "Help.h"
#include <sys/stat.h> // 用于目录操作
#include "../lib/json.hpp"
#include "Core_Export.h"

using json = nlohmann::json;

namespace Poisson {

/// @brief 打印帮助信息
void print_help() {
    std::cout << "------------------------------------------------\n";
    std::cout << "用法：\n";
    std::cout << "  program.exe [-h] [-i <input_file>] [-s [<sample_file>]]\n\n";
    std::cout << "选项：\n";
    std::cout << "  -h, --help          显示此帮助信息\n";
    std::cout << "  -i <input_file>     指定输入配置文件路径(JSON格式)\n";
    std::cout << "  -s [<sample_file>]  生成JSON配置样例文件\n";
    std::cout << "\n输入文件必须是有效的JSON格式，包含以下结构：\n";
    
    // 打印简化的JSON结构参考
    std::cout << R"({
  "mesh": {
    "lx": <正浮点数>,       // x方向长度
    "ly": <正浮点数>,       // y方向长度
    "Nx": <整数(≥2)>,     // x方向单元数
    "Ny": <整数(≥2)>,     // y方向单元数
    "type": <字符串>       // 网格类型("3"-三角形/"4"-四边形)
  },
  "boundary_conditions": {
    "AB": <表达式或数字>,  // 上边界(可选)
    "BC": <表达式或数字>,  // 右边界(可选)
    "CD": <表达式或数字>,  // 下边界(可选)
    "AD": <表达式或数字>   // 左边界(可选)
  },
  "functions": {
    "initial_guess": <表达式>,  // 初始值函数
    "source": <表达式>,         // 源项函数
    "source_derivatives": <表达式> // 源项导数
  },
  "solver_setting": {
    "solver_type": <求解器矩阵类型> // sparse或dense
    "rel_tol": <正浮点数>,      // 相对容差
    "abs_tol": <正浮点数>,      // 绝对容差
    "max_iter": <正整数>       // 最大迭代次数
  },
  "output_file": "<文件路径>"  // 输出VTK文件路径
}
)" << std::endl;

    std::cout << "------------------------------------------------\n";
}

/// @brief 创建目录（跨平台方法）
bool create_directory(const std::string& path) {
    // 尝试创建目录（支持Windows和POSIX系统）
    #ifdef _WIN32
        return _mkdir(path.c_str()) == 0;
    #else
        return mkdir(path.c_str(), 0755) == 0;
    #endif
}

/// @brief 确保目录存在
void ensure_directory_exists(const std::string& path) {
    // 查找最后一个路径分隔符
    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) {
        // 没有目录路径，直接在当前目录创建文件
        return;
    }
    
    std::string dir_path = path.substr(0, pos);
    if (!dir_path.empty()) {
        // 尝试创建目录（如果不存在）
        create_directory(dir_path);
    }
}

/// @brief 生成JSON格式的配置样例文件
void output_input_sample(const std::string& filename) {
    // 创建JSON样例对象
    json sample;
    
    // 网格设置
    sample["mesh"]["lx"] = 2.0;
    sample["mesh"]["ly"] = 1.0;
    sample["mesh"]["Nx"] = 20;
    sample["mesh"]["Ny"] = 10;
    sample["mesh"]["type"] = "4";  // 四边形网格
    
    // 边界条件
    sample["boundary_conditions"]["AD"] = "1.0";  // 左边值
    sample["boundary_conditions"]["AB"] = "1.0";  // 下边值
    sample["boundary_conditions"]["BC"] = "1.0";  // 右边值
    sample["boundary_conditions"]["CD"] = "1.0";  // 上边值
    
    // 函数定义
    sample["functions"]["initial_guess"] = "exp(-x-y)";
    sample["functions"]["source"] = "exp(-u)";
    sample["functions"]["source_derivatives"] = "-exp(-u)";
    
    // 求解器设置
    sample["solver_setting"]["rel_tol"] = 1e-5;
    sample["solver_setting"]["abs_tol"] = 1e-5;
    sample["solver_setting"]["max_iter"] = 100;
    
    // 输出设置
    sample["output_file"] = "./results.vtk";
    
    try {
        // 确保目录存在
        ensure_directory_exists(filename);
        
        // 写入文件
        std::ofstream outfile(filename);
        if (!outfile) {
            throw std::runtime_error("无法创建文件: " + filename);
        }
        
        // 缩进4个空格
        outfile << sample.dump(4);
        outfile.close();
        
        std::cout << "  JSON配置样例已生成: " << filename << "\n";
        std::cout << "  请编辑此文件并按需修改参数后使用\n";
    } catch (const std::exception& e) {
        std::cerr << "  生成配置样例失败: " << e.what() << std::endl;
    }
}

/// @brief 检查命令行参数
std::string check_help(int argc, char* argv[]) {
    std::string inputFile;
    bool generateSample = false;
    std::string sampleFile = "input_sample.json";

    // 遍历命令行参数
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_help();
            std::string outstr(80,'=');
            std::cout << outstr << std::endl;
            exit(0);
        }
        else if (arg == "-s" || arg == "--sample") {
            generateSample = true;
            // 检查下一个参数是否是文件名
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                sampleFile = argv[i + 1];
                i++;  // 跳过文件名参数
            }
        }
        else if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) {
                inputFile = argv[i + 1];
                i++;  // 跳过文件名参数
            } else {
                std::cerr << "错误: -i 参数后缺少文件名\n";
                print_help();
                exit(1);
            }
        }
        else if (arg[0] == '-') {
            std::cerr << "错误: 未知参数 " << arg << "\n";
            print_help();
            exit(1);
        }
        else {
            // 未指定参数的输入文件
            if (inputFile.empty()) {
                inputFile = arg;
            } else {
                std::cerr << "错误: 检测到多余参数 " << arg << "\n";
                print_help();
                exit(1);
            }
        }
    }

    // 处理生成样例请求
    if (generateSample) {
        output_input_sample(sampleFile);
        std::string outstr(80,'=');
        std::cout << outstr << std::endl;
        exit(0);
    }

    // 设置默认输入文件
    if (inputFile.empty()) {
        inputFile = "input.json";
        std::cout << "  未指定输入文件，默认使用: " << inputFile << "\n";
    }

    return inputFile;
}

} // namespace Poisson
#include "Help.h"
#include "../lib/json.hpp"
#include "Core_Export.h"

using json = nlohmann::json;

namespace Poisson {

/// @brief 打印帮助信息
void print_help() {
    std::cout << "------------------------------------------------\n";
    std::cout << "Usage:\n";
    std::cout << "  program.exe [-h] [-i <input_file>] [-s [<sample_file>]]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help          Show this help message\n";
    std::cout << "  -i <input_file>     Specify input configuration file (JSON format)\n";
    std::cout << "  -s [<sample_file>]  Generate JSON configuration sample file\n";
    std::cout << "\nInput file must be valid JSON with the following structure:\n";
    
    // 打印简化的JSON结构参考（使用英文）
    std::cout << R"({
  "mesh": {
    "lx": <positive float>,   // Domain length in x-direction
    "ly": <positive float>,   // Domain length in y-direction
    "Nx": <integer (>=2)>,   // Number of elements in x-direction
    "Ny": <integer (>=2)>,   // Number of elements in y-direction
    "type": <string>         // Mesh type ("3" for triangle ; "4" for rectangle)
  },
  "boundary_conditions": {
    "AB": <string or object>, // Bottom edge (x-direction)
    "AD": <string or object>, // Left edge (y-direction)
    "BC": <string or object>, // Right edge (y-direction)
    "CD": <string or object>  // Top edge (x-direction)
  },
  "functions": {
    "initial_guess": <expression>,  // Initial guess function
    "source": <expression>,         // Source function
    "source_derivatives": <expression> // Source derivatives
  },
  "solver_setting": {
    "solver_type": <string>,       // "sparse" or "dense"
    "rel_tol": <positive float>,   // Relative tolerance
    "abs_tol": <positive float>,   // Absolute tolerance
    "max_iter": <positive integer> // Maximum iterations
  },
  "output_file": "<file path>"    // Output VTK file path
}

Boundary conditions can be:
- Simple string expression:
    "AB": " sin(x) "
- Object form <with range>:
    "AB": {
        "value": " sin(x) ",
        "range": [0.5, 1.0]  // Boundary range
    }
}
)" << std::endl;

    std::cout << "------------------------------------------------\n";
}

/// @brief 生成JSON格式的配置样例文件
void output_input_sample(const std::filesystem::path& filename) {
    // 创建JSON样例对象
    json jsSample;
    
    // 网格设置
    jsSample["mesh"]["lx"] = 2.0;
    jsSample["mesh"]["ly"] = 1.0;
    jsSample["mesh"]["Nx"] = 20;
    jsSample["mesh"]["Ny"] = 10;
    jsSample["mesh"]["type"] = "4";  // Quadrilateral mesh
    
    // 边界条件 - 使用对象形式展示范围
    jsSample["boundary_conditions"]["AB"] = {
        {"value", "0.0"},
        {"range", {0.0, 2.0}}
    };
    jsSample["boundary_conditions"]["AD"] = {
        {"value", "x*y"},
        {"range", {0.0, 1.0}}
    };
    jsSample["boundary_conditions"]["BC"] = {
        {"value", "exp(-y)"},
        {"range", {0.0, 0.8}}
    };
    jsSample["boundary_conditions"]["CD"] = {
        {"value", "sin(x)"},
        {"range", {0.5, 2.0}}
    };
    
    // 函数定义
    jsSample["functions"]["initial_guess"] = "exp(-x-y)";
    jsSample["functions"]["source"] = "exp(-u)";
    jsSample["functions"]["source_derivatives"] = "-exp(-u)";
    
    // 求解器设置
    jsSample["solver_setting"]["solver_type"] = "sparse";
    jsSample["solver_setting"]["rel_tol"] = 1e-5;
    jsSample["solver_setting"]["abs_tol"] = 1e-5;
    jsSample["solver_setting"]["max_iter"] = 100;
    
    // 输出设置
    jsSample["output_file"] = "./results.vtk";
    
    try {
        // 创建目录
        std::filesystem::path dir_path(filename);
        if (!std::filesystem::exists(dir_path.parent_path())) {
            // 如果父目录不存在，则创建它
            std::filesystem::create_directories(dir_path.parent_path());
        }
        
        // 写入文件
        std::ofstream outfile(filename);
        if (!outfile) {
            throw std::runtime_error("Cannot create file: " + filename.string());
        }
        
        // 缩进4个空格
        outfile << jsSample.dump(4);
        outfile.close();
        
        std::cout << "  JSON sample configuration generated: " << filename << "\n";
        std::cout << "  Please edit this file and modify parameters as needed\n";
    } catch (const std::exception& e) {
        std::cerr << "  Failed to generate sample configuration: " << e.what() << std::endl;
    }
}

/// @brief 检查命令行参数
std::filesystem::path check_help(int argc, char* argv[]) {
    std::filesystem::path inputFile;
    bool generateSample = false;
    std::filesystem::path sampleFile = "input_sample.json";

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
                std::cerr << "Error: -i requires a filename\n";
                print_help();
                exit(1);
            }
        }
        else if (arg[0] == '-') {
            std::cerr << "Error: Unknown option " << arg << "\n";
            print_help();
            exit(1);
        }
        else {
            // 未指定参数的输入文件
            if (inputFile.empty()) {
                inputFile = arg;
            } else {
                std::cerr << "Error: Extra argument " << arg << "\n";
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
        std::cout << "  No input file specified, using default: " << inputFile << "\n";
    }

    return inputFile;
}

} // namespace Poisson
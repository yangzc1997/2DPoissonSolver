#ifndef HELP_H
#define HELP_H

#include "Core_Export.h"
#include <iostream>
#include <fstream>
#include <string>

namespace Poisson{

/// @brief 打印帮助信息
POISSONCORE_API void print_help();

/// @brief 检查程序交互参数，确认是否打印帮助信息，并获取输入文件名
/// @param argc main函数系统交互参数
/// @param argv main函数系统交互参数
/// @return 输入文件名
POISSONCORE_API std::string check_help(int argc, char* argv[]);

/// @brief 输出输入文件样例
POISSONCORE_API void output_input_sample(const std::string& filename);

/// @brief 创建目录（跨平台方法）
POISSONCORE_API bool create_directory(const std::string& path);

} // namespace Poisson

#endif
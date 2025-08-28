#ifndef HELP_H
#define HELP_H

#include "Core_Export.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

namespace Poisson{

/// @brief 打印帮助信息
POISSONCORE_API void print_help();

/// @brief 检查程序交互参数，确认是否打印帮助信息，并获取输入文件名
/// @param argc main函数系统交互参数
/// @param argv main函数系统交互参数
/// @return 输入文件名
POISSONCORE_API std::filesystem::path check_help(int argc, char* argv[]);

/// @brief 输出输入文件样例
POISSONCORE_API void output_input_sample(const std::filesystem::path& filename);

} // namespace Poisson

#endif
// ExpressionEvaluator.h
#ifndef EXPRESSION_EVALUATOR_H
#define EXPRESSION_EVALUATOR_H

#include "Core_Export.h"
#include <functional>
#include <string>
#include <cmath>
#include <limits>
#include "../lib/exprtk.hpp"  // 引入 exprtk 库


namespace Poisson{

/// @brief 将字符串表达式转化为可直接使用的函数表达式对象，用于转化输入文件中的源函数及其导数、以及边界条件
class POISSONCORE_API ExpressionEvaluator {
public:
    /// @brief 构造函数，传入表达式字符串和函数类型标识
    /// @param expr_str 字符串表达式
    /// @param expression_type 表达式类型(u/xy/uxy表达式)
    ExpressionEvaluator(const std::string& expr_str);

    ~ExpressionEvaluator() = default;

    /// @brief 计算表达式值
    double evaluate(double u = 0.0, double x = 0.0, double y = 0.0) const;
    
    // 简化接口
    double evaluate_u(double u) const { return evaluate(u); }
    double evaluate_xy(double x, double y) const { return evaluate(0.0, x, y); }
    double evaluate_uxy(double u, double x, double y) const { return evaluate(u, x, y); }

private:
    // 表达式求值函数
    std::function<double(double, double, double)> parse_expression(const std::string& expr_str);

    std::function<double(double, double, double)> expression_func;  // 关于 u, x, y 的表达式
    
    // 变量存储
    double u = 0.0, x = 0.0, y = 0.0;
};

} //namespace Poisson

#endif // EXPRESSION_EVALUATOR_H

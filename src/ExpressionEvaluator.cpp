// 根据字符串表达式返回函数表达式
#include "ExpressionEvaluator.h"
#include "Core_Export.h"
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace Poisson{

/// @brief 将字符串表达式构建为解析表达式
/// @param expr_str 字符串表达式
/// @param expression_type 表达式类型(u/xy/uxy表达式)
ExpressionEvaluator::ExpressionEvaluator(const std::string& expr_str) {
    expression_func = parse_expression(expr_str);
}

/// @brief 计算关于 u,x,y 的表达式值
double ExpressionEvaluator::evaluate(double u_val, double x_val, double y_val) const {
    return expression_func(u_val, x_val, y_val);
}

/// @brief 使用 exprtk 库解析关于 u 的字符串表达式
std::function<double(double, double, double)> ExpressionEvaluator::parse_expression(const std::string& expr_str) {
    typedef double T;
    exprtk::symbol_table<T> symbol_table;
    exprtk::expression<T> expression;
    exprtk::parser<T> parser;
    
    // 添加变量绑定
    symbol_table.add_variable("u", u);
    symbol_table.add_variable("x", x);
    symbol_table.add_variable("y", y);
    
    // 注册符号表并编译表达式
    expression.register_symbol_table(symbol_table);
    if (!parser.compile(expr_str, expression)) {
        if (!expr_str.empty()){
            std::cerr << "警告: 请确认源函数/边界条件/初始解函数是否正确:['" << expr_str << "]'\n";
            std::cerr << "提示: 只能使用 u, x, y 作为变量，支持标准数学函数\n";
        }
        // return [] (double u_val, double x_val, double y_val) { return 0.0; };  // 返回一个值为0的默认函数
       return [] (double u_val, double x_val, double y_val) {return std::numeric_limits<double>::quiet_NaN(); }; // 返回 NAN
    }
    
        // 返回可执行的函数
    return [this, expression](double u_val, double x_val, double y_val) mutable -> double {
        this->u = u_val;
        this->x = x_val;
        this->y = y_val;
        return expression.value();
    };
}

} //namespace Poisson
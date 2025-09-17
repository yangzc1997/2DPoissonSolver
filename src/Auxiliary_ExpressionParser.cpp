// Auxiliary_ExpressionParser.cpp
#include "Auxiliary_ExpressionParser.h"
#include <iostream>
#include <limits>
#include <stdexcept>
#include <functional>
#include <memory>
#include "../lib/exprtk.hpp"

namespace Poisson {

namespace Auxiliary {

namespace ExpressionParser{

    fuxy parse_expression_fuxy(const std::string& expr_str){
        typedef double T;
        exprtk::symbol_table<T> symbol_table;
        exprtk::expression<T> expression;
        exprtk::parser<T> parser;
        
        // 需要单独给u,x,y变量加入到堆之中
        struct ExpressionVars {
            T u = 0.0;
            T x = 0.0;
            T y = 0.0;
        };
        // 这里需要使用c++智能指针自动管理内存，因为我自己很难确定什么时候该释放掉这里的内存
        auto vars = std::make_shared<ExpressionVars>();  

        // 添加变量绑定
        symbol_table.add_variable("u", vars->u);
        symbol_table.add_variable("x", vars->x);
        symbol_table.add_variable("y", vars->y);

        // 添加常量
        symbol_table.add_constant("pi", M_PI);
        symbol_table.add_constant("e", M_E);
        symbol_table.add_constants(); // 添加常量如 pi, e

        // 注册符号表并编译表达式
        expression.register_symbol_table(symbol_table);
        if (!parser.compile(expr_str, expression)) {
            if (expr_str != ""){
                std::cerr << "警告: 表达式: { " << expr_str << " } 解析失败\n";
                std::cerr << "提示: 只能使用 u, x, y 作为变量，支持标准数学函数\n";
            }
            return [](double, double, double) { 
                return std::numeric_limits<double>::quiet_NaN(); 
            };
        }
        
        // 返回可执行的函数
        return [expression, vars](double u_val, double x_val, double y_val) mutable {
            vars->u = u_val;
            vars->x = x_val;
            vars->y = y_val;
            return expression.value();
        };
    }

    fu parse_expression_fu(const std::string& expr_str) {
        auto func = parse_expression_fuxy(expr_str);
        return [func](double u) { return func(u, 0.0, 0.0); };
    }

    fxy parse_expression_fxy(const std::string& expr_str) {
        auto func = parse_expression_fuxy(expr_str);
        return [func](double x, double y) { return func(0.0, x, y); };
    }
    
}  // namespace ExpressionParser

} // namespace Auxiliary

} // namespace Poisson

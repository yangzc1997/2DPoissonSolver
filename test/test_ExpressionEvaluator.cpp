#include <gtest/gtest.h>
#include "../src/Auxiliary_ExpressionParser.h"
#include <cmath> 
#include <iostream>

namespace Poisson {

using namespace Auxiliary::ExpressionParser;

TEST(ExpressionParserTest, BasicExpression) {
    auto func = parse_expression_fxy("x + y");
    EXPECT_DOUBLE_EQ(func(1.0, 2.0), 3.0);
    EXPECT_DOUBLE_EQ(func(3.0, 4.0), 7.0);
}

TEST(ExpressionParserTest, FunctionExpression) {
    auto func = parse_expression_fxy("sin(x) + cos(y)");
    EXPECT_NEAR(func(0.0, 0.0), 1.0, 1e-9);
    EXPECT_NEAR(func(M_PI/2, 0.0), 2.0, 1e-9);
}

TEST(ExpressionParserTest, UxyExpression) {
    auto func = parse_expression_fuxy("u + x + y");
    EXPECT_DOUBLE_EQ(func(1.0, 2.0, 3.0), 6.0);
    EXPECT_DOUBLE_EQ(func(4.0, 5.0, 6.0), 15.0);
}

TEST(ExpressionParserTest, InvalidExpression) {
    // 无效表达式应该返回 NaN
    auto func = parse_expression_fxy("invalid!expression");
    double result = func(1.0, 2.0);
    EXPECT_TRUE(std::isnan(result));
}

TEST(ExpressionParserTest, BoundaryPiecewiseFunction) {
    // 定义分段函数：当 x < 0.5 时返回 sin(x)，否则返回 cos(y)
    std::string piecewise_expr = "if (x < 0.5, sin(x), cos(y))";
    auto func = parse_expression_fxy(piecewise_expr);
    
    // 测试 x < 0.5 的情况
    EXPECT_NEAR(func(0.3, 0.0), std::sin(0.3), 1e-9);
    EXPECT_NEAR(func(0.4, 0.5), std::sin(0.4), 1e-9);
    
    // 测试 x >= 0.5 的情况
    EXPECT_NEAR(func(0.6, 0.0), std::cos(0.0), 1e-9);
    EXPECT_NEAR(func(0.7, M_PI/2), std::cos(M_PI/2), 1e-9);
    
    // 测试边界值 x = 0.5
    EXPECT_NEAR(func(0.5, 0.0), std::cos(0.0), 1e-9);
}

TEST(ExpressionParserTest, ComplexBoundaryCondition) {
    // 更复杂的分段函数：不同区域不同表达式
    std::string complex_expr = 
        "if ( (x < 0.5) and (y < 0.5), sin(x)*cos(y), "
        "if ( (x < 0.5) and (y >= 0.5), x*y, "
        "if ( (x >= 0.5) and (y < 0.5), x+y, "
        "exp(-(x*x + y*y)) )))";
    
    auto func = parse_expression_fxy(complex_expr);
    
    // 区域1: x < 0.5 && y < 0.5
    EXPECT_NEAR(func(0.3, 0.4), std::sin(0.3)*std::cos(0.4), 1e-9);
    
    // 区域2: x < 0.5 && y >= 0.5
    EXPECT_NEAR(func(0.3, 0.6), 0.3 * 0.6, 1e-9);
    
    // 区域3: x >= 0.5 && y < 0.5
    EXPECT_NEAR(func(0.6, 0.4), 0.6+0.4, 1e-9);
    
    // 区域4: x >= 0.5 && y >= 0.5
    double x = 0.7, y = 0.8;
    EXPECT_NEAR(func(x, y), std::exp(-(x*x + y*y)), 1e-9);
    
    // 边界情况
    EXPECT_NEAR(func(0.5, 0.5), std::exp(-(0.5 * 0.5 + 0.5 * 0.5)), 1e-9);
}

TEST(ExpressionParserTest, ConstantExpression) {
    auto func = parse_expression_fxy("pi + e");
    double expected = M_PI + M_E;
    EXPECT_DOUBLE_EQ(func(0.0, 0.0), expected);
}

TEST(ExpressionParserTest, VariableDependence) {
    auto func_xy = parse_expression_fxy("x * y");
    auto func_uxy = parse_expression_fuxy("u * x * y");
    
    // fxy 应该忽略 u 值
    EXPECT_DOUBLE_EQ(func_xy(2.0, 3.0), 6.0);
    EXPECT_DOUBLE_EQ(func_xy(4.0, 5.0), 20.0);
    
    // fuxy 应该使用所有变量
    EXPECT_DOUBLE_EQ(func_uxy(2.0, 3.0, 4.0), 24.0);
    EXPECT_DOUBLE_EQ(func_uxy(1.0, 2.0, 3.0), 6.0);
}

TEST(ExpressionParserTest, EmptyExpression) {
    auto func = parse_expression_fxy("");
    double result = func(1.0, 2.0);
    EXPECT_TRUE(std::isnan(result));
}

TEST(ExpressionParserTest, AdvancedMathFunctions) {
    auto func = parse_expression_fxy("log(x) + sqrt(y)");
    EXPECT_NEAR(func(1.0, 1.0), 0.0 + 1.0, 1e-9);
    EXPECT_NEAR(func(M_E, 4.0), 1.0 + 2.0, 1e-9);
}

} // namespace Poisson
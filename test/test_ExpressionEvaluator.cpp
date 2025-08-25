#include <gtest/gtest.h>
#include <cmath> 
#include "../src/ExpressionEvaluator.h"

namespace Poisson {

TEST(ExpressionEvaluatorTest, BasicExpression) {
    ExpressionEvaluator eval("x + y");
    EXPECT_DOUBLE_EQ(eval.evaluate_xy(1.0, 2.0), 3.0);
    EXPECT_DOUBLE_EQ(eval.evaluate_xy(3.0, 4.0), 7.0);
}

TEST(ExpressionEvaluatorTest, FunctionExpression) {
    ExpressionEvaluator eval("sin(x) + cos(y)");
    EXPECT_NEAR(eval.evaluate_xy(0.0, 0.0), 1.0, 1e-9);
    EXPECT_NEAR(eval.evaluate_xy(M_PI/2, 0.0), 2.0, 1e-9);
}

TEST(ExpressionEvaluatorTest, UxyExpression) {
    ExpressionEvaluator eval("u + x + y");
    EXPECT_DOUBLE_EQ(eval.evaluate_uxy(1.0, 2.0, 3.0), 6.0);
    EXPECT_DOUBLE_EQ(eval.evaluate_uxy(4.0, 5.0, 6.0), 15.0);
}

TEST(ExpressionEvaluatorTest, InvalidExpression) {
    // 无效表达式应该返回 NaN
    ExpressionEvaluator eval("invalid!expression");
    double result = eval.evaluate_xy(1.0, 2.0);
    EXPECT_TRUE(std::isnan(result));
}

TEST(ExpressionEvaluatorTest, BoundaryPiecewiseFunction) {
    // 定义分段函数：当 x < 0.5 时返回 sin(x)，否则返回 cos(y)
    std::string piecewise_expr = "if (x < 0.5, sin(x), cos(y))";
    ExpressionEvaluator eval(piecewise_expr);
    
    // 测试 x < 0.5 的情况
    EXPECT_NEAR(eval.evaluate_xy(0.3, 0.0), std::sin(0.3), 1e-9);
    EXPECT_NEAR(eval.evaluate_xy(0.4, 0.5), std::sin(0.4), 1e-9);
    
    // 测试 x >= 0.5 的情况
    EXPECT_NEAR(eval.evaluate_xy(0.6, 0.0), std::cos(0.0), 1e-9);
    EXPECT_NEAR(eval.evaluate_xy(0.7, M_PI/2), std::cos(M_PI/2), 1e-9);
    
    // 测试边界值 x = 0.5
    EXPECT_NEAR(eval.evaluate_xy(0.5, 0.0), std::cos(0.0), 1e-9);
}

TEST(ExpressionEvaluatorTest, ComplexBoundaryCondition) {
    // 更复杂的分段函数：不同区域不同表达式
    std::string complex_expr = 
        "if ( (x < 0.5) and (y < 0.5), sin(x)*cos(y), "
        "if ( (x < 0.5) and (y >= 0.5), x*y, "
        "if ( (x >= 0.5) and (y < 0.5), x+y, "
        "exp(-(x*x + y*y)) )))";
    
    ExpressionEvaluator eval(complex_expr);
    
    // 区域1: x < 0.5 && y < 0.5
    EXPECT_NEAR(eval.evaluate_xy(0.3, 0.4), std::sin(0.3)*std::cos(0.4), 1e-9);
    
    // 区域2: x < 0.5 && y >= 0.5
    EXPECT_NEAR(eval.evaluate_xy(0.3, 0.6), 0.3 * 0.6, 1e-9);
    
    // 区域3: x >= 0.5 && y < 0.5
    EXPECT_NEAR(eval.evaluate_xy(0.6, 0.4), 0.6+0.4, 1e-9);
    
    // 区域4: x >= 0.5 && y >= 0.5
    double x = 0.7, y = 0.8;
    EXPECT_NEAR(eval.evaluate_xy(x, y), std::exp(-(x*x + y*y)), 1e-9);
    
    // 边界情况
    EXPECT_NEAR(eval.evaluate_xy(0.5, 0.5), std::exp(-(0.5 * 0.5 + 0.5 * 0.5)), 1e-9);
}

} // namespace Poisson
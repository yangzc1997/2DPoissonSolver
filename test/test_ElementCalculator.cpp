#include <gtest/gtest.h>
#include "ElementCalculator.h"
#include <cmath>
#include <iostream>

namespace Poisson {

// 测试三角形单元的高斯积分点
TEST(GaussPointsTest, TrianglePoints) {
    TriangleCalculator calculator;
    
    // 测试不同阶数的积分点
    for (int order = 1; order <= 10; order++) {
        calculator.setIntegrationOrder(order);
        auto points = calculator.gaussPoints();
        
        // 验证积分点数量
        int expected_points = 0;
        switch (order) {
            case 1: expected_points = 1; break;
            case 2: expected_points = 3; break;
            case 3: expected_points = 4; break;
            case 4: expected_points = 6; break;
            case 5: expected_points = 7; break;
            case 6: expected_points = 12; break;
            case 7: expected_points = 13; break;
            case 8: expected_points = 16; break;
            case 9: expected_points = 19; break;
            case 10: expected_points = 25; break;
            default: expected_points = 25; break;
        }
        EXPECT_EQ(points.size(), expected_points) << "Order " << order << " has wrong number of points";
        
        // 验证权重和（应为0.5）
        double sum_weights = 0.0;
        for (const auto& gp : points) {
            sum_weights += gp.weight;
        }
        EXPECT_NEAR(sum_weights, 0.5, 1e-10) << "Order " << order << " has incorrect weight sum";
        
        // 验证点坐标在三角形内
        for (const auto& gp : points) {
            EXPECT_GE(gp.s, 0.0) << "Order " << order << " point s < 0";
            EXPECT_GE(gp.t, 0.0) << "Order " << order << " point t < 0";
            EXPECT_LE(gp.s + gp.t, 1.0 + 1e-10) << "Order " << order << " point s+t > 1";
        }
    }
}

// 测试四边形单元的高斯积分点
TEST(GaussPointsTest, RectanglePoints) {
    RectangleCalculator calculator;
    
    // 测试不同阶数的积分点
    for (int order = 1; order <= 10; order++) {
        calculator.setIntegrationOrder(order);
        auto points = calculator.gaussPoints();
        
        // 验证积分点数量
        int expected_points = order * order;
        EXPECT_EQ(points.size(), expected_points) << "Order " << order << " has wrong number of points";
        
        // 验证权重和（应为4.0）
        double sum_weights = 0.0;
        for (const auto& gp : points) {
            sum_weights += gp.weight;
        }
        EXPECT_NEAR(sum_weights, 4.0, 1e-10) << "Order " << order << " has incorrect weight sum";
        
        // 验证点坐标在[-1,1]范围内
        for (const auto& gp : points) {
            EXPECT_GE(gp.s, -1.0) << "Order " << order << " point s < -1";
            EXPECT_LE(gp.s, 1.0) << "Order " << order << " point s > 1";
            EXPECT_GE(gp.t, -1.0) << "Order " << order << " point t < -1";
            EXPECT_LE(gp.t, 1.0) << "Order " << order << " point t > 1";
        }
    }
}

// 测试积分精度 - 三角形单元
TEST(GaussPointsTest, TriangleIntegrationAccuracy) {
    TriangleCalculator calculator;
    
    // 测试多项式积分
    auto poly2d = [](double s, double t) { return s*s + t*t; };
    double exact_integral = 1.0/6.0; // ∫∫(s² + t²) ds dt over triangle = 1/6
    
    for (int order = 1; order <= 10; order++) {
        calculator.setIntegrationOrder(order);
        auto points = calculator.gaussPoints();
        
        double integral = 0.0;
        for (const auto& gp : points) {
            integral += gp.weight * poly2d(gp.s, gp.t);
        }
        
        // 根据阶数设置不同的精度要求
        double tolerance = 1e-6;
        if (order >= 3) tolerance = 1e-8;
        
        EXPECT_NEAR(integral, exact_integral, tolerance) 
            << "Order " << order << " failed to integrate polynomial";
    }
}

// 测试积分精度 - 四边形单元
TEST(GaussPointsTest, RectangleIntegrationAccuracy) {
    RectangleCalculator calculator;
    
    // 测试多项式积分
    auto poly2d = [](double s, double t) { return s*s*s + t*t*t; };
    double exact_integral = 0.0; // ∫∫(s³ + t³) ds dt over [-1,1]^2 = 0
    
    for (int order = 1; order <= 10; order++) {
        calculator.setIntegrationOrder(order);
        auto points = calculator.gaussPoints();
        
        double integral = 0.0;
        for (const auto& gp : points) {
            integral += gp.weight * poly2d(gp.s, gp.t);
        }
        
        // 根据阶数设置不同的精度要求
        double tolerance = 1e-6;
        if (order >= 3) tolerance = 1e-10;
        if (order >= 5) tolerance = 1e-12;
        
        EXPECT_NEAR(integral, exact_integral, tolerance) 
            << "Order " << order << " failed to integrate polynomial";
    }
}

// 测试常数函数积分 - 三角形单元
TEST(GaussPointsTest, TriangleConstantIntegration) {
    TriangleCalculator calculator;
    
    // 常数函数积分应等于面积0.5
    auto constant_func = [](double s, double t) { return 1.0; };
    double exact_integral = 0.5;
    
    for (int order = 1; order <= 10; order++) {
        calculator.setIntegrationOrder(order);
        auto points = calculator.gaussPoints();
        
        double integral = 0.0;
        for (const auto& gp : points) {
            integral += gp.weight * constant_func(gp.s, gp.t);
        }
        
        EXPECT_NEAR(integral, exact_integral, 1e-12) 
            << "Order " << order << " failed to integrate constant function";
    }
}

// 测试常数函数积分 - 四边形单元
TEST(GaussPointsTest, RectangleConstantIntegration) {
    RectangleCalculator calculator;
    
    // 常数函数积分应等于面积4.0
    auto constant_func = [](double s, double t) { return 1.0; };
    double exact_integral = 4.0;
    
    for (int order = 1; order <= 10; order++) {
        calculator.setIntegrationOrder(order);
        auto points = calculator.gaussPoints();
        
        double integral = 0.0;
        for (const auto& gp : points) {
            integral += gp.weight * constant_func(gp.s, gp.t);
        }
        
        EXPECT_NEAR(integral, exact_integral, 1e-12) 
            << "Order " << order << " failed to integrate constant function";
    }
}

// 测试积分点唯一性
TEST(GaussPointsTest, UniquePoints) {
    TriangleCalculator tri_calculator;
    RectangleCalculator rect_calculator;
    
    for (int order = 1; order <= 10; order++) {
        tri_calculator.setIntegrationOrder(order);
        auto tri_points = tri_calculator.gaussPoints();
        
        rect_calculator.setIntegrationOrder(order);
        auto rect_points = rect_calculator.gaussPoints();
        
        // 检查三角形积分点是否唯一
        for (size_t i = 0; i < tri_points.size(); i++) {
            for (size_t j = i + 1; j < tri_points.size(); j++) {
                double dist = std::hypot(tri_points[i].s - tri_points[j].s, 
                                        tri_points[i].t - tri_points[j].t);
                EXPECT_GT(dist, 1e-8) << "Order " << order << " has duplicate triangle points";
            }
        }
        
        // 检查四边形积分点是否唯一
        for (size_t i = 0; i < rect_points.size(); i++) {
            for (size_t j = i + 1; j < rect_points.size(); j++) {
                double dist = std::hypot(rect_points[i].s - rect_points[j].s, 
                                        rect_points[i].t - rect_points[j].t);
                EXPECT_GT(dist, 1e-8) << "Order " << order << " has duplicate rectangle points";
            }
        }
    }
}

} // namespace Poisson
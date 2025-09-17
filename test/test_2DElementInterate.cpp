#include <gtest/gtest.h>
#include "../src/Auxiliary_Integrate2DLegendre.h"
#include <cmath>
#include <iostream>

namespace Poisson {
namespace Auxiliary {
namespace Integrate2DLegendre {

// 测试三角形单元的高斯积分点
TEST(Integrate2DLegendreTest, TrianglePoints) {
    // 测试不同阶数的积分点
    for (int order = 1; order <= 10; order++) {
        auto points = generateTriangleGaussPoints(order);
        
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
TEST(Integrate2DLegendreTest, RectanglePoints) {
    // 测试不同阶数的积分点
    for (int order = 1; order <= 10; order++) {
        auto points = generateRectangleGaussPoints(order);
        
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
TEST(Integrate2DLegendreTest, TriangleIntegrationAccuracy) {
    // 测试多项式积分
    auto poly2d = [](double s, double t) { return s*s + t*t; };
    double exact_integral = 1.0/6.0; // ∫∫(s² + t²) ds dt over triangle = 1/6
    
    for (int order = 3; order <= 10; order++) {
        auto points = generateTriangleGaussPoints(order);
        
        double integral = 0.0;
        for (const auto& gp : points) {
            integral += gp.weight * poly2d(gp.s, gp.t);
        }
        
        // 根据阶数设置不同的精度要求
        double tolerance = 1e-6;
        
        EXPECT_NEAR(integral, exact_integral, tolerance) 
            << "Order " << order << " failed to integrate polynomial";
    }
}

// 测试积分精度 - 四边形单元
TEST(Integrate2DLegendreTest, RectangleIntegrationAccuracy) {
    // 测试多项式积分
    auto poly2d = [](double s, double t) { return s*s*s + t*t*t; };
    double exact_integral = 0.0; // ∫∫(s³ + t³) ds dt over [-1,1]^2 = 0
    
    for (int order = 1; order <= 10; order++) {
        auto points = generateRectangleGaussPoints(order);
        
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
TEST(Integrate2DLegendreTest, TriangleConstantIntegration) {
    // 常数函数积分应等于面积0.5
    auto constant_func = [](double s, double t) { return 1.0; };
    double exact_integral = 0.5;
    
    for (int order = 1; order <= 10; order++) {
        auto points = generateTriangleGaussPoints(order);
        
        double integral = 0.0;
        for (const auto& gp : points) {
            integral += gp.weight * constant_func(gp.s, gp.t);
        }
        
        EXPECT_NEAR(integral, exact_integral, 1e-12) 
            << "Order " << order << " failed to integrate constant function";
    }
}

// 测试常数函数积分 - 四边形单元
TEST(Integrate2DLegendreTest, RectangleConstantIntegration) {
    // 常数函数积分应等于面积4.0
    auto constant_func = [](double s, double t) { return 1.0; };
    double exact_integral = 4.0;
    
    for (int order = 1; order <= 10; order++) {
        auto points = generateRectangleGaussPoints(order);
        
        double integral = 0.0;
        for (const auto& gp : points) {
            integral += gp.weight * constant_func(gp.s, gp.t);
        }
        
        EXPECT_NEAR(integral, exact_integral, 1e-12) 
            << "Order " << order << " failed to integrate constant function";
    }
}

// 测试积分点唯一性
TEST(Integrate2DLegendreTest, UniquePoints) {
    for (int order = 1; order <= 10; order++) {
        auto tri_points = generateTriangleGaussPoints(order);
        auto rect_points = generateRectangleGaussPoints(order);
        
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

// 测试默认阶数
TEST(Integrate2DLegendreTest, DefaultOrder) {
    // 测试三角形默认阶数
    auto tri_points_default = generateTriangleGaussPoints();
    auto tri_points_7 = generateTriangleGaussPoints(7);
    EXPECT_EQ(tri_points_default.size(), tri_points_7.size());
    
    // 测试四边形默认阶数
    auto rect_points_default = generateRectangleGaussPoints();
    auto rect_points_7 = generateRectangleGaussPoints(7);
    EXPECT_EQ(rect_points_default.size(), rect_points_7.size());
}

// 测试高阶积分点
TEST(Integrate2DLegendreTest, HighOrderIntegration) {
    // 测试高阶多项式积分 - 三角形
    auto poly2d_tri = [](double s, double t) { return s*s*s*s + t*t*t*t; };
    double exact_integral_tri = 1.0/15.0; // ∫∫(s⁴ + t⁴) ds dt = 1/15
    
    auto tri_points = generateTriangleGaussPoints(10);
    double integral_tri = 0.0;
    for (const auto& gp : tri_points) {
        integral_tri += gp.weight * poly2d_tri(gp.s, gp.t);
    }
    EXPECT_NEAR(integral_tri, exact_integral_tri, 1e-8);
    
    // 测试高阶多项式积分 - 四边形
    auto poly2d_rect = [](double s, double t) { return s*s*s*s + t*t*t*t; };
    double exact_integral_rect = 8.0/5.0; // ∫∫(s⁴ + t⁴) ds dt over [-1,1]^2 = 32/5
    
    auto rect_points = generateRectangleGaussPoints(10);
    double integral_rect = 0.0;
    for (const auto& gp : rect_points) {
        integral_rect += gp.weight * poly2d_rect(gp.s, gp.t);
    }
    EXPECT_NEAR(integral_rect, exact_integral_rect, 1e-8);
}

} // namespace Integrate2DLegendre
} // namespace Auxiliary
} // namespace Poisson
// ElementCalculator.cpp
#include "ElementCalculator.h"
#include "Core_Export.h"
#include <boost/math/quadrature/gauss.hpp>
#include <iostream>

namespace Poisson {

using namespace Eigen;
using namespace boost::math::quadrature;


// ==================== 单元计算器实现 ====================== //
// 空间某点处的解值u(x,y)
double ElementCalculator::u_value(const VectorXd& global_u, 
                                   const std::vector<int>& element, 
                                   double x, double y) const 
{
    double result = 0.0;
    for (int i = 0; i < element.size(); i++) {
        result += global_u(element[i]) * shape_function(i, x, y);
    }
    return result;
}


// u在某点的导数u'(x,y)
Vector2d ElementCalculator::u_gradient(const VectorXd& global_u, 
                                       const std::vector<int>& element,
                                       const std::vector<double>& x_coords, 
                                       const std::vector<double>& y_coords,
                                       double x, double y) const 
{
    Vector2d grad(0, 0);
    for (int i = 0; i < element.size(); i++) {
        grad += global_u(element[i]) * shape_gradient(i, x_coords, y_coords, x, y);
    }
    return grad;
}


// ====================== 三角形计算器实现 ====================== //
// 三角形单元雅可比行列式
double TriangleCalculator::jacobian(const std::vector<double>& x_coords, 
                                    const std::vector<double>& y_coords) const 
{
    return (x_coords[1]-x_coords[0])*(y_coords[2]-y_coords[0]) - 
           (x_coords[2]-x_coords[0])*(y_coords[1]-y_coords[0]);
}


// 三角形单元的形函数
double TriangleCalculator::shape_function(int i, double x, double y) const 
{
    static const std::array<std::function<double(double, double)>, 3> shape_funcs = {
        [](double x, double y) { return 1 - x - y; },  // N1 implementation
        [](double x, double y) { return x;}, // N2 implementation
        [](double x, double y) { return y;} // N3 implementation
    };
    return shape_funcs[i](x, y);
}

// 三角形单元的形函数梯度
Vector2d TriangleCalculator::shape_gradient(int i, const std::vector<double>& x_coords, 
    const std::vector<double>& y_coords, double s, double t) const 
{
    double A2 = jacobian(x_coords, y_coords);
    static const std::array<std::function<Vector2d(const std::vector<double>, const std::vector<double>, double, double, const double)>, 3> gradients = {
        [](const std::vector<double>& x_coords, const std::vector<double>& y_coords, double s, double t, const double& A2) {
            return Vector2d((y_coords[1] - y_coords[2])/A2, (x_coords[2] - x_coords[1])/A2);
        }, // N1 gradient
        [](const std::vector<double>& x_coords, const std::vector<double>& y_coords, double s, double t, const double& A2) {
            return Vector2d((y_coords[2] - y_coords[0])/A2, (x_coords[0] - x_coords[2])/A2);
        }, // N2 gradient
        [](const std::vector<double>& x_coords, const std::vector<double>& y_coords, double s, double t, const double& A2) {
            return Vector2d( (y_coords[0] - y_coords[1])/A2, (x_coords[1] - x_coords[0])/A2);
        } // N3 gradient
    };
    return gradients[i](x_coords, y_coords, s, t, A2);
}


double TriangleCalculator::integrate(const std::vector<double>& x_coords,
                                    const std::vector<double>& y_coords,
                                    const std::function<double(double, double)>& func) const 
{
    const double jac = jacobian(x_coords, y_coords);
    // 使用通用积分模板
    return integrate2D<7, 7>( // 使用7点高斯积分
        0.0, 1.0, // 外层积分范围 (s)
        [](double s) { return 0.0; }, // 内层积分下限 (s)
        [](double s) { return 1.0 - s; }, // 内层积分上限 (s)
        [jac, &func](double s, double t) {
            return func(s, t) * jac; // 应用雅可比变换
        },
        false // 不检查边界
    );
}

// ====================== 四边形计算器实现 ====================== //
// 四边形单元雅可比行列式
double RectangleCalculator::jacobian(const std::vector<double>& x_coords, 
                                     const std::vector<double>& y_coords) const 
{
    return (x_coords[1]-x_coords[0])*(y_coords[2]-y_coords[0])/4.0;
}


// 四边形形函数
double RectangleCalculator::shape_function(int i, double s, double t) const 
{
    static const std::array<std::function<double(double, double)>, 4> shape_funcs = {
        [](double s, double t) { return (1.0-s)*(1.0-t)/4.0; }, // N1
        [](double s, double t) { return (1.0+s)*(1.0-t)/4.0; }, // N2
        [](double s, double t) { return (1.0+s)*(1.0+t)/4.0; }, // N3
        [](double s, double t) { return (1.0-s)*(1.0+t)/4.0; }  // N4
    };
    return shape_funcs[i](s, t);
}

// 四边形单元梯度
Vector2d RectangleCalculator::shape_gradient(int i, const std::vector<double>& x_coords,
    const std::vector<double>& y_coords, double s, double t) const 
{
    static const std::array<std::function<Vector2d(const std::vector<double>, const std::vector<double>, double, double)>, 4> gradients = {
        [](const std::vector<double>& x_coords, const std::vector<double>& y_coords, double s, double t) {
            return Vector2d(-(1.0-t)/(2.0*(x_coords[1]-x_coords[0])), -(1.0-s)/(2.0*(y_coords[2]-y_coords[0])));
        }, // N1 gradient
        [](const std::vector<double>& x_coords, const std::vector<double>& y_coords, double s, double t) {
            return Vector2d((1.0-t)/(2.0*(x_coords[1]-x_coords[0])), -(1.0+s)/(2.0*(y_coords[2]-y_coords[0])));
        }, // N2 gradient
        [](const std::vector<double>& x_coords, const std::vector<double>& y_coords, double s, double t) {
            return Vector2d((1.0+t)/(2.0*(x_coords[1]-x_coords[0])), (1.0+s)/(2.0*(y_coords[2]-y_coords[0])));
        }, // N3 gradient
        [](const std::vector<double>& x_coords, const std::vector<double>& y_coords, double s, double t) {
            return Vector2d( -(1.0+t)/(2.0*(x_coords[1]-x_coords[0])), (1.0-s)/(2.0*(y_coords[2]-y_coords[0])));
        } // N4 gradient
    };
    return gradients[i](x_coords, y_coords, s, t);
}

// 四边形单元积分实现
double RectangleCalculator::integrate(const std::vector<double>& x_coords,
                                    const std::vector<double>& y_coords,
                                    const std::function<double(double, double)>& func) const 
{
    const double jac = jacobian(x_coords, y_coords);
    
    // 使用通用积分模板
    return integrate2D<7, 7>( // 使用5点高斯积分
        -1.0, 1.0, // 外层积分范围 (t)
        [](double t) { return -1.0; }, // 内层积分下限 (s)
        [](double t) { return 1.0; }, // 内层积分上限 (s)
        [jac, &func](double s, double t) {
            return func(s, t) * jac; // 应用雅可比变换
        }, false // 不检查边界
    );
}

} // namespace Poisson
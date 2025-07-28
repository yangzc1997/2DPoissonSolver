// ElementCalculator.h
#ifndef ELEMENT_CALCULATOR_H
#define ELEMENT_CALCULATOR_H

#include "Core_Export.h"
#include <vector>
#include <functional>
#include <eigen3/Eigen/Dense>
#include <boost/math/quadrature/gauss.hpp>

namespace Poisson {

/// @class ElementCalculator
/// @brief 单元计算器基类
class POISSONCORE_API ElementCalculator {
public:
    virtual ~ElementCalculator() = default;
    
    // 纯虚函数接口
    virtual double jacobian(const std::vector<double>& x_coords, 
                          const std::vector<double>& y_coords) const = 0;
    
    virtual double u_value(const Eigen::VectorXd& global_u, 
                          const std::vector<int>& element, 
                          double param1, double param2) const;
                          
    virtual Eigen::Vector2d u_gradient(const Eigen::VectorXd& global_u, 
                                     const std::vector<int>& element,
                                     const std::vector<double>& x_coords, 
                                     const std::vector<double>& y_coords,
                                     double param1, double param2) const;
    
    virtual double shape_function(int i, double param1, double param2) const = 0;
                                
    virtual Eigen::Vector2d shape_gradient(int i, 
                                         const std::vector<double>& x_coords, 
                                         const std::vector<double>& y_coords,
                                         double param1, double param2) const = 0;
    
    // 积分方法
    virtual double integrate(const std::vector<double>& x_coords,
                            const std::vector<double>& y_coords,
                            const std::function<double(double, double)>& func) const = 0;

    
protected:
    //提供一个二维积分器
    template<unsigned N1 = 20, unsigned N2 = 20>
    double integrate2D(
        double x0,
        double x1,
        const std::function<double(double)>& y0,
        const std::function<double(double)>& y1,
        const std::function<double(double, double)>& func,
        bool check_bounds = false
    ) const
    {
        using OuterRule = boost::math::quadrature::gauss<double, N1>;
        using InnerRule = boost::math::quadrature::gauss<double, N2>;

        auto inner_integral = [&](double x) {
            const double a = y0(x);
            const double b = y1(x);
            if (check_bounds && a >= b)
                return 0.0;

            return InnerRule::integrate(
                [x, &func](double y) {return func(x, y);},a,b
            );
        };

        return OuterRule::integrate(inner_integral, x0, x1);
    }
    
};


// ====================== 三角形计算器 ====================== //
class POISSONCORE_API TriangleCalculator : public ElementCalculator {
public:
    // 实现所有虚函数
    double jacobian(const std::vector<double>& x_coords, 
                   const std::vector<double>& y_coords) const override;
    
    double shape_function(int i, double x, double y) const override;
                         
    Eigen::Vector2d shape_gradient(int i, 
                                  const std::vector<double>& x_coords, 
                                  const std::vector<double>& y_coords,
                                  double x, double y) const override;
    
    double integrate(const std::vector<double>& x_coords,
                    const std::vector<double>& y_coords,
                    const std::function<double(double, double)>& func) const override;
};


// ====================== 四边形计算器 ====================== //
class POISSONCORE_API RectangleCalculator : public ElementCalculator {
public:
    double jacobian(const std::vector<double>& x_coords, 
                    const std::vector<double>& y_coords) const override;
    
    double shape_function(int i, double s, double t) const override;
                                
    Eigen::Vector2d shape_gradient(int i, 
                                 const std::vector<double>& x_coords, 
                                 const std::vector<double>& y_coords,
                                 double s, double t) const override;
    
    double integrate(const std::vector<double>& x_coords,
                    const std::vector<double>& y_coords,
                    const std::function<double(double, double)>& func) const override;
};

} // namespace Poisson

#endif // ELEMENT_CALCULATOR_H

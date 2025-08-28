// FiniteElementCalculator.h
#ifndef FINITE_ELEMENT_CALCULATOR_H
#define FINITE_ELEMENT_CALCULATOR_H

#include "Core_Export.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <functional>

namespace Poisson {

using namespace Eigen;
using vec_t = Eigen::VectorXd;
using vec_t2 = Eigen::Vector2d;
using mat_t = Eigen::MatrixXd;
using mat_t2 = Eigen::Matrix<double, Eigen::Dynamic, 2>;
using NodeCoords = std::vector<Eigen::Vector2d>;
using fuxy = std::function<double(double u_val, double x_val, double y_val)>; // 源函数及其导数

// 高斯积分点结构体
struct GaussPoint {
    double s;
    double t;
    double weight;
};

/// @class FiniteElementCalculator
/// @brief 单元计算器基类
class POISSONCORE_API FiniteElementCalculator {
public:
    FiniteElementCalculator(int defaultOrder) : integrationOrder(defaultOrder) {}
    
    virtual ~FiniteElementCalculator() = default;
    
    /// @brief 计算形函数向量
    virtual vec_t shapeFunction(double s, double t) const = 0;
    
    /// @brief 计算形函数梯度矩阵
    virtual mat_t2 shapeFunctionGradient(double s, double t) const = 0;
    
    /// @brief 计算实空间坐标
    virtual vec_t2 getPhysicalCoordinates(const NodeCoords& coords,double s, double t) const = 0;
    
    /// @brief 计算坐标变换雅可比矩阵
    virtual mat_t2 jacobian(const NodeCoords& coords) const = 0;
    
    /// @brief 设置积分阶数
    virtual void setIntegrationOrder(int order) = 0;

    /// @brief 获取高斯积分点
    std::vector<GaussPoint> gaussPoints() const;

    /// @brief 计算单元刚度矩阵
    mat_t computeStiffnessMatrix( const NodeCoords& coords, const fuxy& source_derivativesFunc, const vec_t& local_u) const;
    
    /// @brief 计算单元载荷向量
    vec_t computeLoadVector(const NodeCoords& coords, const fuxy& sourceFunc, const vec_t& local_u) const;

    // @brief 同时计算单元载荷向量和刚度矩阵
    std::pair<mat_t, vec_t> computeElementMatrixAndVector(
        const NodeCoords& coords,
        const fuxy& sourceFunc,
        const fuxy& source_derivativesFunc,
        const vec_t& local_u) const;

protected:
    int integrationOrder;
    std::vector<GaussPoint> gaussPointsCache;

    // 生成积分点
    virtual void generateGaussPoints2D(int order) = 0;
};

// ====================== 三角形计算器 ====================== //
class POISSONCORE_API TriangleCalculator : public FiniteElementCalculator {
public:
    TriangleCalculator();

    vec_t shapeFunction(double s, double t) const override;

    mat_t2 shapeFunctionGradient(double s, double t) const override;

    vec_t2 getPhysicalCoordinates(const NodeCoords& coords, double s, double t) const override;

    mat_t2 jacobian(const NodeCoords& coords) const override;

    void setIntegrationOrder(int order) override;

private:
    void generateGaussPoints2D(int order) override;
};

// ====================== 四边形计算器 ====================== //
class POISSONCORE_API RectangleCalculator : public FiniteElementCalculator {
public:
    RectangleCalculator(); 

    vec_t shapeFunction(double s, double t) const override;
    
    mat_t2 shapeFunctionGradient(double s, double t) const override;

    vec_t2 getPhysicalCoordinates(const NodeCoords& coords, double s, double t) const override;
    
    mat_t2 jacobian(const NodeCoords& coords) const override;

    void setIntegrationOrder(int order) override;
    
private:
    void generateGaussPoints2D(int order) override;
};

} // namespace Poisson

#endif // FINITE_ELEMENT_CALCULATOR_H
// FiniteElementCalculator.h
#ifndef FINITE_ELEMENT_CALCULATOR_H
#define FINITE_ELEMENT_CALCULATOR_H

#include <eigen3/Eigen/Dense>
#include <vector>
#include <functional>
#include "UsingAlias.h"
#include "Auxiliary_Integrate2DLegendre.h"


namespace Poisson {
    using GaussPoint = Auxiliary::Integrate2DLegendre::GaussPoint;

/// @class FiniteElementCalculator
/// @brief 单元计算器基类
class FiniteElementCalculator {
public:
    FiniteElementCalculator() = default;
    
    virtual ~FiniteElementCalculator() = default;
    
    /// @brief 计算形函数向量
    virtual vec_t shapeFunction(double s, double t) const = 0;
    
    /// @brief 计算形函数梯度矩阵
    virtual mat_t2 shapeFunctionGradient(double s, double t) const = 0;
    
    /// @brief 计算实空间坐标
    virtual vec_t2 getPhysicalCoordinates(const NodeCoords& coords,double s, double t) const = 0;
    
    /// @brief 计算坐标变换雅可比矩阵
    virtual mat_t2 jacobian(const NodeCoords& coords) const = 0;

    /// @brief 计算单元刚度矩阵
    mat_t computeStiffnessMatrix( const NodeCoords& coords, const fuxy& source_derivativesFunc, const vec_t& local_u, const std::vector<GaussPoint>& gaussPoints) const;
    
    /// @brief 计算单元载荷向量
    vec_t computeLoadVector(const NodeCoords& coords, const fuxy& sourceFunc, const vec_t& local_u, const std::vector<GaussPoint>& gaussPoints) const;

    // @brief 同时计算单元载荷向量和刚度矩阵
    std::pair<mat_t, vec_t> computeElementMatrixAndVector(
        const NodeCoords& coords,
        const fuxy& sourceFunc,
        const fuxy& source_derivativesFunc,
        const vec_t& local_u,
        const std::vector<GaussPoint>& gaussPoints
        ) const;
};

// ====================== 三角形计算器 ====================== //
class TriangleCalculator : public FiniteElementCalculator {
public:
    TriangleCalculator() = default;

    vec_t shapeFunction(double s, double t) const override;

    mat_t2 shapeFunctionGradient(double s, double t) const override;

    vec_t2 getPhysicalCoordinates(const NodeCoords& coords, double s, double t) const override;

    mat_t2 jacobian(const NodeCoords& coords) const override;

};

// ====================== 四边形计算器 ====================== //
class RectangleCalculator : public FiniteElementCalculator {
public:
    RectangleCalculator() = default; 

    vec_t shapeFunction(double s, double t) const override;
    
    mat_t2 shapeFunctionGradient(double s, double t) const override;

    vec_t2 getPhysicalCoordinates(const NodeCoords& coords, double s, double t) const override;
    
    mat_t2 jacobian(const NodeCoords& coords) const override;
};

} // namespace Poisson

#endif // FINITE_ELEMENT_CALCULATOR_H
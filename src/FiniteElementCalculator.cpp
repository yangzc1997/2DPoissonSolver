// FiniteElementCalculator.cpp

#include "FiniteElementCalculator.h"
#include <cmath>

namespace Poisson {
    using namespace Eigen;

// ================== 刚度矩阵和载荷向量的统一计算 ================ //
std::pair<mat_t, vec_t> FiniteElementCalculator::computeElementMatrixAndVector(
    const NodeCoords& coords,
    const fuxy& sourceFunc,
    const fuxy& source_derivativesFunc,
    const vec_t& local_u,
    const std::vector<GaussPoint>& gaussPoints
    ) const
{
    const int num_nodes = coords.size();
    mat_t K_local = mat_t::Zero(num_nodes, num_nodes);
    vec_t F_local = vec_t::Zero(num_nodes);
    
    // 计算公共量
    const mat_t2 Jac = jacobian(coords);
    const mat_t2 Jac_inverse = Jac.inverse();
    const double abs_jac = Jac.determinant();

    // 遍历所有高斯点
    for (const auto& gp : gaussPoints) {
        const double s = gp.s;
        const double t = gp.t;
        const double w = gp.weight;
    
        // 计算物理坐标
        vec_t2 xy = getPhysicalCoordinates(coords, s, t);
        double x = xy.x(), y = xy.y();

        // 计算形函数及其梯度
        vec_t N = shapeFunction(s, t);
        mat_t2 N_grad = shapeFunctionGradient(s, t) * Jac_inverse;
        
        // 计算当前解值和解的梯度
        double u_val = N.dot(local_u);
        vec_t2 grad_u = N_grad.transpose() * local_u;

        // 计算刚度矩阵所需
        double dFdu = source_derivativesFunc(u_val, x, y);
        
        // 计算载荷向量所需
        double F = sourceFunc(u_val, x, y);

        // 计算刚度矩阵贡献
        mat_t grad_dot = N_grad * N_grad.transpose();
        mat_t N_Nt = N * N.transpose();
        mat_t K_contrib = -grad_dot - dFdu * N_Nt;
        
        // 计算载荷向量贡献
        vec_t dN_du = N_grad * grad_u;
        vec_t F_N = F * N;
        vec_t F_contrib = dN_du + F_N;
        
        double w_absJac = w * abs_jac;

        K_local += w_absJac * K_contrib;
        F_local += w_absJac * F_contrib;
    }
    
    return {K_local, F_local};
}

// 计算单元刚度矩阵
mat_t FiniteElementCalculator::computeStiffnessMatrix(
    const NodeCoords& coords,
    const fuxy& source_derivativesFunc,
    const vec_t& local_u,
    const std::vector<GaussPoint>& gaussPoints
    ) const
{
    const int num_nodes = coords.size();
    mat_t K_local = mat_t::Zero(num_nodes, num_nodes);
    const mat_t2 Jac = jacobian(coords);
    const mat_t2 Jac_inverse = Jac.inverse();
    const double abs_jac = Jac.determinant();
    
    for (const auto& gp : gaussPoints) {
        const double s = gp.s;
        const double t = gp.t;
        const double w = gp.weight;
        
        // 计算物理坐标
        vec_t2 xy = getPhysicalCoordinates(coords, s, t);
        double x = xy.x(), y = xy.y();
        
        // 计算形函数及其梯度
        vec_t N = shapeFunction(s, t);
        mat_t2 N_grad = shapeFunctionGradient(s, t) * Jac_inverse;
        
        // 计算当前解值
        double u_val = N.dot(local_u);
        
        // 计算源函数导数
        double dFdu = source_derivativesFunc(u_val, x, y);
        
        mat_t dN_dN = N_grad * N_grad.transpose();
        mat_t dFdu_NN = dFdu * (N * N.transpose());
        mat_t K_contrib = -dN_dN - dFdu_NN;

        // 计算刚度矩阵
        K_local += w * abs_jac * K_contrib;
    }
    
    return K_local;
}

// 单元载荷向量
vec_t FiniteElementCalculator::computeLoadVector(
    const NodeCoords& coords,
    const fuxy& sourceFunc,
    const vec_t& local_u,
    const std::vector<GaussPoint>& gaussPoints
    ) const 
{
    const int num_nodes = coords.size();
    vec_t F_local = vec_t::Zero(num_nodes);
    const mat_t2 Jac = jacobian(coords);
    const mat_t2 Jac_inverse = Jac.inverse();;
    const double abs_jac = Jac.determinant();
    
    for (const auto& gp : gaussPoints) {
        const double s = gp.s;
        const double t = gp.t;
        const double w = gp.weight;
        
        // 计算物理坐标
        vec_t2 xy = getPhysicalCoordinates(coords, s, t);
        double x = xy.x(), y = xy.y();
        
        // 计算形函数及其梯度
        vec_t N = shapeFunction(s, t);
        mat_t2 N_grad = shapeFunctionGradient(s, t) * Jac_inverse;
        
        // 计算当前解值
        double u_val = N.dot(local_u);
        
        // 计算解梯度
        vec_t2 grad_u = N_grad.transpose() * local_u;
        
        // 计算源函数
        double F = sourceFunc(u_val, x, y);
        
        vec_t dN_du = N_grad * grad_u;
        vec_t F_N = F * N;
        vec_t F_contrib = dN_du + F_N;

        // 计算载荷向量
        F_local += w * abs_jac * F_contrib;
    }
    
    return F_local;
}


// ====================== 三角形计算器实现 ====================== //
vec_t TriangleCalculator::shapeFunction(double s, double t) const {
    vec_t N(3);
    N << (1 - s - t),
         s,
         t;
    return N;
}

mat_t2 TriangleCalculator::shapeFunctionGradient(double s, double t) const {
    mat_t2 N_grad(3, 2);
    N_grad << -1, -1,
               1,  0,
               0,  1;
    return N_grad;
}

vec_t2 TriangleCalculator::getPhysicalCoordinates(const NodeCoords& coords,
    double s, double t) const 
{
    vec_t2 xy;
    xy << coords[0].x() + s * (coords[1].x() - coords[0].x()) + t * (coords[2].x() - coords[0].x()),
          coords[0].y() + s * (coords[1].y() - coords[0].y()) + t * (coords[2].y() - coords[0].y());
    return xy;
}

mat_t2 TriangleCalculator::jacobian(const NodeCoords& coords) const 
{
    mat_t2 J(2, 2);
    J << coords[1].x() - coords[0].x(), coords[2].x() - coords[0].x(),
         coords[1].y() - coords[0].y(), coords[2].y() - coords[0].y();

    return J;
}


// ====================== 四边形计算器实现 ====================== //
vec_t RectangleCalculator::shapeFunction(double s, double t) const {
    vec_t N(4);
    N << 0.25 * (1 - s) * (1 - t),
         0.25 * (1 + s) * (1 - t),
         0.25 * (1 + s) * (1 + t),
         0.25 * (1 - s) * (1 + t);
    return N;
}

mat_t2 RectangleCalculator::shapeFunctionGradient(double s, double t) const {
    mat_t2 N_grad(4, 2);
    N_grad << -0.25*(1-t), -0.25*(1-s),
               0.25*(1-t), -0.25*(1+s),
               0.25*(1+t),  0.25*(1+s),
              -0.25*(1+t),  0.25*(1-s);
    return N_grad;
}


vec_t2 RectangleCalculator::getPhysicalCoordinates(
    const NodeCoords& coords,
    double s, double t) const 
{
    vec_t2 xy;
    xy << 0.5 * ((coords[1].x() - coords[0].x()) * s + (coords[1].x() + coords[0].x())),
          0.5 * ((coords[3].y() - coords[0].y()) * t + (coords[3].y() + coords[0].y()));
    return xy;
}


mat_t2 RectangleCalculator::jacobian(const NodeCoords& coords) const 
{
    mat_t2 J(2, 2);
    J << 0.5 * (coords[1].x() - coords[0].x()), 0.0,
         0.0, 0.5 * (coords[3].y() - coords[0].y());
    return J;
}

} // namespace Poisson
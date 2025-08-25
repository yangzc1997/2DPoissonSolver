// ElementCalculator.cpp

#include "ElementCalculator.h"
#include "Core_Export.h"
#include <cmath>

namespace Poisson {

// ================== 刚度矩阵和载荷向量的统一计算 ================ //
std::pair<mat_t, vec_t> ElementCalculator::computeElementMatrixAndVector(
    const NodeCoords& coords,
    const func_uxy& sourceFunc,
    const func_uxy& source_derivativesFunc,
    const vec_t& local_u) const
{
    const int num_nodes = coords.size();
    mat_t K_local = mat_t::Zero(num_nodes, num_nodes);
    vec_t F_local = vec_t::Zero(num_nodes);
    
    // 计算公共量
    const mat_t2 Jac = jacobian(coords);
    const mat_t2 Jac_inverse = Jac.inverse();
    const double abs_jac = Jac.determinant();

    // 遍历所有高斯点
    for (const auto& gp : gaussPointsCache) {
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
mat_t ElementCalculator::computeStiffnessMatrix(
    const NodeCoords& coords,
    const func_uxy& source_derivativesFunc,
    const vec_t& local_u) const
{
    const int num_nodes = coords.size();
    mat_t K_local = mat_t::Zero(num_nodes, num_nodes);
    const mat_t2 Jac = jacobian(coords);
    const mat_t2 Jac_inverse = Jac.inverse();
    const double abs_jac = Jac.determinant();
    
    for (const auto& gp : gaussPointsCache) {
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
vec_t ElementCalculator::computeLoadVector(
    const NodeCoords& coords,
    const func_uxy& sourceFunc,
    const vec_t& local_u) const 
{
    const int num_nodes = coords.size();
    vec_t F_local = vec_t::Zero(num_nodes);
    const mat_t2 Jac = jacobian(coords);
    const mat_t2 Jac_inverse = Jac.inverse();;
    const double abs_jac = Jac.determinant();
    
    for (const auto& gp : gaussPointsCache) {
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

// 获取高斯积分点
std::vector<GaussPoint> ElementCalculator::gaussPoints() const {
    return gaussPointsCache;
}


// ====================== 三角形计算器实现 ====================== //
TriangleCalculator::TriangleCalculator(): ElementCalculator(5) { 
    setIntegrationOrder(integrationOrder);
}

void TriangleCalculator::setIntegrationOrder(int order) {
    integrationOrder = order;
    generateGaussPoints2D(order);
}


// 注意，这里的权重是归一化的，实际计算需要再乘以0.5 
void TriangleCalculator::generateGaussPoints2D(int order) {
    gaussPointsCache.clear();
    
    switch (order) {
    case 1: // 1点积分 (阶数1)
        gaussPointsCache.push_back({1.0/3.0, 1.0/3.0, 1.0/2.0});
        break;
        
    case 2: // 3点积分 (阶数2)
        gaussPointsCache.push_back({1.0/6.0, 1.0/6.0, 1.0/3.0/2.0});
        gaussPointsCache.push_back({2.0/3.0, 1.0/6.0, 1.0/3.0/2.0});
        gaussPointsCache.push_back({1.0/6.0, 2.0/3.0, 1.0/3.0/2.0});
        break;
        
    case 3: // 4点积分 (阶数3)
        gaussPointsCache.push_back({1.0/3.0, 1.0/3.0, -27.0/48.0/2.0});
        gaussPointsCache.push_back({0.2, 0.2, 25.0/48.0/2.0});
        gaussPointsCache.push_back({0.2, 0.6, 25.0/48.0/2.0});
        gaussPointsCache.push_back({0.6, 0.2, 25.0/48.0/2.0});
        break;
        
    case 4: // 6点积分 (阶数4)
        gaussPointsCache.push_back({0.44594849091597, 0.44594849091597, 0.22338158967801/2.0});
        gaussPointsCache.push_back({0.44594849091597, 0.10810301816807, 0.22338158967801/2.0});
        gaussPointsCache.push_back({0.10810301816807, 0.44594849091597, 0.22338158967801/2.0});
        gaussPointsCache.push_back({0.09157621350977, 0.09157621350977, 0.10995174365532/2.0});
        gaussPointsCache.push_back({0.09157621350977, 0.81684757298046, 0.10995174365532/2.0});
        gaussPointsCache.push_back({0.81684757298046, 0.09157621350977, 0.10995174365532/2.0});
        break;
        
    case 5: // 7点积分 (阶数5)
        gaussPointsCache.push_back({1.0/3.0, 1.0/3.0, 0.225/2.0});
        gaussPointsCache.push_back({0.47014206410511, 0.47014206410511, 0.13239415278851/2.0});
        gaussPointsCache.push_back({0.47014206410511, 0.05971587178977, 0.13239415278851/2.0});
        gaussPointsCache.push_back({0.05971587178977, 0.47014206410511, 0.13239415278851/2.0});
        gaussPointsCache.push_back({0.10128650732346, 0.10128650732346, 0.12593918054483/2.0});
        gaussPointsCache.push_back({0.10128650732346, 0.79742698535309, 0.12593918054483/2.0});
        gaussPointsCache.push_back({0.79742698535309, 0.10128650732346, 0.12593918054483/2.0});
        break;
        
    case 6: // 12点积分 (阶数6)
        gaussPointsCache.push_back({0.24928674517091, 0.24928674517091, 0.11678627572638/2.0});
        gaussPointsCache.push_back({0.24928674517091, 0.50142650965818, 0.11678627572638/2.0});
        gaussPointsCache.push_back({0.50142650965818, 0.24928674517091, 0.11678627572638/2.0});
        gaussPointsCache.push_back({0.06308901449150, 0.06308901449150, 0.05084490637021/2.0});
        gaussPointsCache.push_back({0.06308901449150, 0.87382197101700, 0.05084490637021/2.0});
        gaussPointsCache.push_back({0.87382197101700, 0.06308901449150, 0.05084490637021/2.0});
        gaussPointsCache.push_back({0.31035245103378, 0.63650249912140, 0.08285107561837/2.0});
        gaussPointsCache.push_back({0.63650249912140, 0.05314504984482, 0.08285107561837/2.0});
        gaussPointsCache.push_back({0.05314504984482, 0.31035245103378, 0.08285107561837/2.0});
        gaussPointsCache.push_back({0.63650249912140, 0.31035245103378, 0.08285107561837/2.0});
        gaussPointsCache.push_back({0.31035245103378, 0.05314504984482, 0.08285107561837/2.0});
        gaussPointsCache.push_back({0.05314504984482, 0.63650249912140, 0.08285107561837/2.0});
        break;
        
    case 7: // 13点积分 (阶数7)
        gaussPointsCache.push_back({1.0/3.0, 1.0/3.0, -0.14957004446767/2.0});
        gaussPointsCache.push_back({0.26034596607904, 0.26034596607904, 0.17561525743321/2.0});
        gaussPointsCache.push_back({0.26034596607904, 0.47930806784192, 0.17561525743321/2.0});
        gaussPointsCache.push_back({0.47930806784192, 0.26034596607904, 0.17561525743321/2.0});
        gaussPointsCache.push_back({0.06513010290222, 0.06513010290222, 0.05334723560884/2.0});
        gaussPointsCache.push_back({0.06513010290222, 0.86973979419556, 0.05334723560884/2.0});
        gaussPointsCache.push_back({0.86973979419556, 0.06513010290222, 0.05334723560884/2.0});
        gaussPointsCache.push_back({0.31286549600487, 0.63844418856981, 0.07711376089026/2.0});
        gaussPointsCache.push_back({0.63844418856981, 0.04869031542532, 0.07711376089026/2.0});
        gaussPointsCache.push_back({0.04869031542532, 0.31286549600487, 0.07711376089026/2.0});
        gaussPointsCache.push_back({0.63844418856981, 0.31286549600487, 0.07711376089026/2.0});
        gaussPointsCache.push_back({0.31286549600487, 0.04869031542532, 0.07711376089026/2.0});
        gaussPointsCache.push_back({0.04869031542532, 0.63844418856981, 0.07711376089026/2.0});
        break;
        
    case 8: // 16点积分 (阶数8)
        gaussPointsCache.push_back({1.0/3.0, 1.0/3.0, 0.14431560767779/2.0});
        gaussPointsCache.push_back({0.17056930775176, 0.17056930775176, 0.10321737053472/2.0});
        gaussPointsCache.push_back({0.17056930775176, 0.65886138449648, 0.10321737053472/2.0});
        gaussPointsCache.push_back({0.65886138449648, 0.17056930775176, 0.10321737053472/2.0});
        gaussPointsCache.push_back({0.05054722831703, 0.05054722831703, 0.03245849762320/2.0});
        gaussPointsCache.push_back({0.05054722831703, 0.89890554336594, 0.03245849762320/2.0});
        gaussPointsCache.push_back({0.89890554336594, 0.05054722831703, 0.03245849762320/2.0});
        gaussPointsCache.push_back({0.45929258829272, 0.45929258829272, 0.09509163426728/2.0});
        gaussPointsCache.push_back({0.45929258829272, 0.08141482341456, 0.09509163426728/2.0});
        gaussPointsCache.push_back({0.08141482341456, 0.45929258829272, 0.09509163426728/2.0});
        gaussPointsCache.push_back({0.26311282963464, 0.72849239295540, 0.02723031417443/2.0});
        gaussPointsCache.push_back({0.72849239295540, 0.00839477740996, 0.02723031417443/2.0});
        gaussPointsCache.push_back({0.00839477740996, 0.26311282963464, 0.02723031417443/2.0});
        gaussPointsCache.push_back({0.72849239295540, 0.26311282963464, 0.02723031417443/2.0});
        gaussPointsCache.push_back({0.26311282963464, 0.00839477740996, 0.02723031417443/2.0});
        gaussPointsCache.push_back({0.00839477740996, 0.72849239295540, 0.02723031417443/2.0});
        break;
        
    case 9: // 19点积分 (阶数9)
        gaussPointsCache.push_back({1.0/3.0, 1.0/3.0, 0.09713579628280/2.0});
        gaussPointsCache.push_back({0.48968251919874, 0.48968251919874, 0.03133470022714/2.0});
        gaussPointsCache.push_back({0.48968251919874, 0.02884473323268, 0.03133470022714/2.0});
        gaussPointsCache.push_back({0.02884473323268, 0.48968251919874, 0.03133470022714/2.0});
        gaussPointsCache.push_back({0.43708959149294, 0.43708959149294, 0.07782754100474/2.0});
        gaussPointsCache.push_back({0.43708959149294, 0.12582081701412, 0.07782754100474/2.0});
        gaussPointsCache.push_back({0.12582081701412, 0.43708959149294, 0.07782754100474/2.0});
        gaussPointsCache.push_back({0.18820353561903, 0.18820353561903, 0.07964773892721/2.0});
        gaussPointsCache.push_back({0.18820353561903, 0.62359292876194, 0.07964773892721/2.0});
        gaussPointsCache.push_back({0.62359292876194, 0.18820353561903, 0.07964773892721/2.0});
        gaussPointsCache.push_back({0.04472951339445, 0.04472951339445, 0.02557767565870/2.0});
        gaussPointsCache.push_back({0.04472951339445, 0.91054097321110, 0.02557767565870/2.0});
        gaussPointsCache.push_back({0.91054097321110, 0.04472951339445, 0.02557767565870/2.0});
        gaussPointsCache.push_back({0.22196298916077, 0.74119859878478, 0.04328353937729/2.0});
        gaussPointsCache.push_back({0.74119859878478, 0.03683841205445, 0.04328353937729/2.0});
        gaussPointsCache.push_back({0.03683841205445, 0.22196298916077, 0.04328353937729/2.0});
        gaussPointsCache.push_back({0.74119859878478, 0.22196298916077, 0.04328353937729/2.0});
        gaussPointsCache.push_back({0.22196298916077, 0.03683841205445, 0.04328353937729/2.0});
        gaussPointsCache.push_back({0.03683841205445, 0.74119859878478, 0.04328353937729/2.0});
        break;
        
    case 10: // 25点积分 (阶数10)
        gaussPointsCache.push_back({1.0/3.0, 1.0/3.0, 0.09081799038275/2.0});
        gaussPointsCache.push_back({0.48557763338366, 0.48557763338366, 0.03672595775647/2.0});
        gaussPointsCache.push_back({0.48557763338366, 0.02884473323268, 0.03672595775647/2.0});
        gaussPointsCache.push_back({0.02884473323268, 0.48557763338366, 0.03672595775647/2.0});
        gaussPointsCache.push_back({0.10948157548504, 0.10948157548504, 0.04532105943553/2.0});
        gaussPointsCache.push_back({0.10948157548504, 0.78103684902992, 0.04532105943553/2.0});
        gaussPointsCache.push_back({0.78103684902992, 0.10948157548504, 0.04532105943553/2.0});
        gaussPointsCache.push_back({0.30793983876412, 0.55035294182100, 0.07275791684542/2.0});
        gaussPointsCache.push_back({0.55035294182100, 0.14170721941488, 0.07275791684542/2.0});
        gaussPointsCache.push_back({0.14170721941488, 0.30793983876412, 0.07275791684542/2.0});
        gaussPointsCache.push_back({0.55035294182100, 0.30793983876412, 0.07275791684542/2.0});
        gaussPointsCache.push_back({0.30793983876412, 0.14170721941488, 0.07275791684542/2.0});
        gaussPointsCache.push_back({0.14170721941488, 0.55035294182100, 0.07275791684542/2.0});
        gaussPointsCache.push_back({0.24667256063990, 0.72832390459741, 0.02832724253106/2.0});
        gaussPointsCache.push_back({0.72832390459741, 0.02500353476269, 0.02832724253106/2.0});
        gaussPointsCache.push_back({0.02500353476269, 0.24667256063990, 0.02832724253106/2.0});
        gaussPointsCache.push_back({0.72832390459741, 0.24667256063990, 0.02832724253106/2.0});
        gaussPointsCache.push_back({0.24667256063990, 0.02500353476269, 0.02832724253106/2.0});
        gaussPointsCache.push_back({0.02500353476269, 0.72832390459741, 0.02832724253106/2.0});
        gaussPointsCache.push_back({0.06665406347960, 0.92365593358750, 0.00942166696373/2.0});
        gaussPointsCache.push_back({0.92365593358750, 0.00969000293290, 0.00942166696373/2.0});
        gaussPointsCache.push_back({0.00969000293290, 0.06665406347960, 0.00942166696373/2.0});
        gaussPointsCache.push_back({0.92365593358750, 0.06665406347960, 0.00942166696373/2.0});
        gaussPointsCache.push_back({0.06665406347960, 0.00969000293290, 0.00942166696373/2.0});
        gaussPointsCache.push_back({0.00969000293290, 0.92365593358750, 0.00942166696373/2.0});
        break;
            
        default:
            // 默认使用阶数10
            generateGaussPoints2D(10);
            break;
    }
}

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
RectangleCalculator::RectangleCalculator(): ElementCalculator(7) {
    setIntegrationOrder(integrationOrder);
}

void RectangleCalculator::setIntegrationOrder(int order) {
    integrationOrder = order;
    generateGaussPoints2D(order);
}

// Golub–Welsch 算法生成[-1,1]区间的积分点和权重
void RectangleCalculator::generateGaussPoints2D(int n)  {
        // 1D 节点和权重
    std::vector<double> x(n), w(n);

    // 构造 Jacobi 三对角矩阵
    vec_t beta(n-1);
    for (int k = 1; k < n; ++k)
        beta(k-1) = k / std::sqrt(4.0*k*k - 1.0);

    mat_t J = mat_t::Zero(n,n);
    for (int i = 0; i < n-1; ++i) {
        J(i,i+1) = beta(i);
        J(i+1,i) = beta(i);
    }

    Eigen::SelfAdjointEigenSolver<mat_t> es(J);
    const auto& evals = es.eigenvalues();
    const auto& evecs = es.eigenvectors();

    for (int i = 0; i < n; ++i) {
        x[i] = evals(i);
        w[i] = 2.0 * evecs(0,i) * evecs(0,i);
    }

    gaussPointsCache.clear();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            gaussPointsCache.push_back({x[i], x[j], w[i] * w[j]});
        }
    }
}

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
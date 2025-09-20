// Auxiliary_Integrate2DLegendre.h

#ifndef AUXILIARY_INTEGRATE_2D_LEGENDRE_H
#define AUXILIARY_INTEGRATE_2D_LEGENDRE_H

#include <vector>
#include <eigen3/Eigen/Dense>

namespace Poisson {
namespace Auxiliary {

/// @namespace Integrate2DLegendre
/// @brief 二维高斯-勒让德积分点生成器
namespace Integrate2DLegendre {

    // 高斯积分点结构体
    struct GaussPoint {
        double s;
        double t;
        double weight;
    };

    /// @brief 生成三角形单元的高斯积分点
    std::vector<GaussPoint> generateTriangleGaussPoints(int order=7);
    
    /// @brief 生成四边形单元的高斯积分点
    std::vector<GaussPoint> generateRectangleGaussPoints(int order=7);

} // namespace Integrate2DLegendre

} // namespace Auxiliary

} // namespace Poisson

#endif // AUXILIARY_INTEGRATE_2D_LEGENDRE_H
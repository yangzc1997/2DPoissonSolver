// BoundaryData.h
#ifndef BOUNDARY_DATA_H
#define BOUNDARY_DATA_H

#include <array>
#include <string>

namespace Poisson {

/// @brief 表示单条边界信息的结构体
struct SingleEdgeInfo {
    std::string value = "";  // 边界函数的表达式
    std::array<double, 2> range = {0.0, 1.0}; // 边界作用范围【min,max】

    SingleEdgeInfo() = default;

    SingleEdgeInfo(const std::string& val, const std::array<double, 2>& rng)
        : value(val), range(rng){}
};

/// @brief 长方形区域中各种边界条件集合的结构体
struct BoundaryConditionInfo {
    SingleEdgeInfo AB;  // 底边 (x轴方向)
    SingleEdgeInfo AD;  // 左边 (y轴方向)
    SingleEdgeInfo BC;  // 右边 (y轴方向)
    SingleEdgeInfo CD;  // 顶边 (x轴方向)

    BoundaryConditionInfo() = default;

    BoundaryConditionInfo(const std::string& ab, const std::array<double, 2>& ab_range,
                         const std::string& ad, const std::array<double, 2>& ad_range,
                         const std::string& bc, const std::array<double, 2>& bc_range,
                         const std::string& cd, const std::array<double, 2>& cd_range)
        : AB(ab, ab_range), 
          AD(ad, ad_range), 
          BC(bc, bc_range), 
          CD(cd, cd_range) {}
};

} // namespace Poisson

#endif // BOUNDARY_DATA_H
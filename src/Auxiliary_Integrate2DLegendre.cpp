// Auxiliary_Integrate2DLegendre.h

#include "Auxiliary_Integrate2DLegendre.h"
#include <eigen3/Eigen/Dense>
#include <vector>

namespace Poisson {

namespace Auxiliary {

namespace Integrate2DLegendre {

    // 生成三角形单元的高斯积分点
    std::vector<GaussPoint> generateTriangleGaussPoints(int order) {
        
        std::vector<GaussPoint> gaussPoints;
        
        switch (order) {
        case 1: // 1点积分 (阶数1)
            gaussPoints.push_back({1.0/3.0, 1.0/3.0, 1.0/2.0});
            break;
        
        case 2: // 3点积分 (阶数2)
            gaussPoints.push_back({1.0/6.0, 1.0/6.0, 1.0/3.0/2.0});
            gaussPoints.push_back({2.0/3.0, 1.0/6.0, 1.0/3.0/2.0});
            gaussPoints.push_back({1.0/6.0, 2.0/3.0, 1.0/3.0/2.0});
            break;
            
        case 3: // 4点积分 (阶数3)
            gaussPoints.push_back({1.0/3.0, 1.0/3.0, -27.0/48.0/2.0});
            gaussPoints.push_back({0.2, 0.2, 25.0/48.0/2.0});
            gaussPoints.push_back({0.2, 0.6, 25.0/48.0/2.0});
            gaussPoints.push_back({0.6, 0.2, 25.0/48.0/2.0});
            break;
            
        case 4: // 6点积分 (阶数4)
            gaussPoints.push_back({0.44594849091597, 0.44594849091597, 0.22338158967801/2.0});
            gaussPoints.push_back({0.44594849091597, 0.10810301816807, 0.22338158967801/2.0});
            gaussPoints.push_back({0.10810301816807, 0.44594849091597, 0.22338158967801/2.0});
            gaussPoints.push_back({0.09157621350977, 0.09157621350977, 0.10995174365532/2.0});
            gaussPoints.push_back({0.09157621350977, 0.81684757298046, 0.10995174365532/2.0});
            gaussPoints.push_back({0.81684757298046, 0.09157621350977, 0.10995174365532/2.0});
            break;
            
        case 5: // 7点积分 (阶数5)
            gaussPoints.push_back({1.0/3.0, 1.0/3.0, 0.225/2.0});
            gaussPoints.push_back({0.47014206410511, 0.47014206410511, 0.13239415278851/2.0});
            gaussPoints.push_back({0.47014206410511, 0.05971587178977, 0.13239415278851/2.0});
            gaussPoints.push_back({0.05971587178977, 0.47014206410511, 0.13239415278851/2.0});
            gaussPoints.push_back({0.10128650732346, 0.10128650732346, 0.12593918054483/2.0});
            gaussPoints.push_back({0.10128650732346, 0.79742698535309, 0.12593918054483/2.0});
            gaussPoints.push_back({0.79742698535309, 0.10128650732346, 0.12593918054483/2.0});
            break;
            
        case 6: // 12点积分 (阶数6)
            gaussPoints.push_back({0.24928674517091, 0.24928674517091, 0.11678627572638/2.0});
            gaussPoints.push_back({0.24928674517091, 0.50142650965818, 0.11678627572638/2.0});
            gaussPoints.push_back({0.50142650965818, 0.24928674517091, 0.11678627572638/2.0});
            gaussPoints.push_back({0.06308901449150, 0.06308901449150, 0.05084490637021/2.0});
            gaussPoints.push_back({0.06308901449150, 0.87382197101700, 0.05084490637021/2.0});
            gaussPoints.push_back({0.87382197101700, 0.06308901449150, 0.05084490637021/2.0});
            gaussPoints.push_back({0.31035245103378, 0.63650249912140, 0.08285107561837/2.0});
            gaussPoints.push_back({0.63650249912140, 0.05314504984482, 0.08285107561837/2.0});
            gaussPoints.push_back({0.05314504984482, 0.31035245103378, 0.08285107561837/2.0});
            gaussPoints.push_back({0.63650249912140, 0.31035245103378, 0.08285107561837/2.0});
            gaussPoints.push_back({0.31035245103378, 0.05314504984482, 0.08285107561837/2.0});
            gaussPoints.push_back({0.05314504984482, 0.63650249912140, 0.08285107561837/2.0});
            break;
            
        case 7: // 13点积分 (阶数7)
            gaussPoints.push_back({1.0/3.0, 1.0/3.0, -0.14957004446767/2.0});
            gaussPoints.push_back({0.26034596607904, 0.26034596607904, 0.17561525743321/2.0});
            gaussPoints.push_back({0.26034596607904, 0.47930806784192, 0.17561525743321/2.0});
            gaussPoints.push_back({0.47930806784192, 0.26034596607904, 0.17561525743321/2.0});
            gaussPoints.push_back({0.06513010290222, 0.06513010290222, 0.05334723560884/2.0});
            gaussPoints.push_back({0.06513010290222, 0.86973979419556, 0.05334723560884/2.0});
            gaussPoints.push_back({0.86973979419556, 0.06513010290222, 0.05334723560884/2.0});
            gaussPoints.push_back({0.31286549600487, 0.63844418856981, 0.07711376089026/2.0});
            gaussPoints.push_back({0.63844418856981, 0.04869031542532, 0.07711376089026/2.0});
            gaussPoints.push_back({0.04869031542532, 0.31286549600487, 0.07711376089026/2.0});
            gaussPoints.push_back({0.63844418856981, 0.31286549600487, 0.07711376089026/2.0});
            gaussPoints.push_back({0.31286549600487, 0.04869031542532, 0.07711376089026/2.0});
            gaussPoints.push_back({0.04869031542532, 0.63844418856981, 0.07711376089026/2.0});
            break;
            
        case 8: // 16点积分 (阶数8)
            gaussPoints.push_back({1.0/3.0, 1.0/3.0, 0.14431560767779/2.0});
            gaussPoints.push_back({0.17056930775176, 0.17056930775176, 0.10321737053472/2.0});
            gaussPoints.push_back({0.17056930775176, 0.65886138449648, 0.10321737053472/2.0});
            gaussPoints.push_back({0.65886138449648, 0.17056930775176, 0.10321737053472/2.0});
            gaussPoints.push_back({0.05054722831703, 0.05054722831703, 0.03245849762320/2.0});
            gaussPoints.push_back({0.05054722831703, 0.89890554336594, 0.03245849762320/2.0});
            gaussPoints.push_back({0.89890554336594, 0.05054722831703, 0.03245849762320/2.0});
            gaussPoints.push_back({0.45929258829272, 0.45929258829272, 0.09509163426728/2.0});
            gaussPoints.push_back({0.45929258829272, 0.08141482341456, 0.09509163426728/2.0});
            gaussPoints.push_back({0.08141482341456, 0.45929258829272, 0.09509163426728/2.0});
            gaussPoints.push_back({0.26311282963464, 0.72849239295540, 0.02723031417443/2.0});
            gaussPoints.push_back({0.72849239295540, 0.00839477740996, 0.02723031417443/2.0});
            gaussPoints.push_back({0.00839477740996, 0.26311282963464, 0.02723031417443/2.0});
            gaussPoints.push_back({0.72849239295540, 0.26311282963464, 0.02723031417443/2.0});
            gaussPoints.push_back({0.26311282963464, 0.00839477740996, 0.02723031417443/2.0});
            gaussPoints.push_back({0.00839477740996, 0.72849239295540, 0.02723031417443/2.0});
            break;
            
        case 9: // 19点积分 (阶数9)
            gaussPoints.push_back({1.0/3.0, 1.0/3.0, 0.09713579628280/2.0});
            gaussPoints.push_back({0.48968251919874, 0.48968251919874, 0.03133470022714/2.0});
            gaussPoints.push_back({0.48968251919874, 0.02884473323268, 0.03133470022714/2.0});
            gaussPoints.push_back({0.02884473323268, 0.48968251919874, 0.03133470022714/2.0});
            gaussPoints.push_back({0.43708959149294, 0.43708959149294, 0.07782754100474/2.0});
            gaussPoints.push_back({0.43708959149294, 0.12582081701412, 0.07782754100474/2.0});
            gaussPoints.push_back({0.12582081701412, 0.43708959149294, 0.07782754100474/2.0});
            gaussPoints.push_back({0.18820353561903, 0.18820353561903, 0.07964773892721/2.0});
            gaussPoints.push_back({0.18820353561903, 0.62359292876194, 0.07964773892721/2.0});
            gaussPoints.push_back({0.62359292876194, 0.18820353561903, 0.07964773892721/2.0});
            gaussPoints.push_back({0.04472951339445, 0.04472951339445, 0.02557767565870/2.0});
            gaussPoints.push_back({0.04472951339445, 0.91054097321110, 0.02557767565870/2.0});
            gaussPoints.push_back({0.91054097321110, 0.04472951339445, 0.02557767565870/2.0});
            gaussPoints.push_back({0.22196298916077, 0.74119859878478, 0.04328353937729/2.0});
            gaussPoints.push_back({0.74119859878478, 0.03683841205445, 0.04328353937729/2.0});
            gaussPoints.push_back({0.03683841205445, 0.22196298916077, 0.04328353937729/2.0});
            gaussPoints.push_back({0.74119859878478, 0.22196298916077, 0.04328353937729/2.0});
            gaussPoints.push_back({0.22196298916077, 0.03683841205445, 0.04328353937729/2.0});
            gaussPoints.push_back({0.03683841205445, 0.74119859878478, 0.04328353937729/2.0});
            break;
            
        case 10: // 25点积分 (阶数10)
            gaussPoints.push_back({1.0/3.0, 1.0/3.0, 0.09081799038275/2.0});
            gaussPoints.push_back({0.48557763338366, 0.48557763338366, 0.03672595775647/2.0});
            gaussPoints.push_back({0.48557763338366, 0.02884473323268, 0.03672595775647/2.0});
            gaussPoints.push_back({0.02884473323268, 0.48557763338366, 0.03672595775647/2.0});
            gaussPoints.push_back({0.10948157548504, 0.10948157548504, 0.04532105943553/2.0});
            gaussPoints.push_back({0.10948157548504, 0.78103684902992, 0.04532105943553/2.0});
            gaussPoints.push_back({0.78103684902992, 0.10948157548504, 0.04532105943553/2.0});
            gaussPoints.push_back({0.30793983876412, 0.55035294182100, 0.07275791684542/2.0});
            gaussPoints.push_back({0.55035294182100, 0.14170721941488, 0.07275791684542/2.0});
            gaussPoints.push_back({0.14170721941488, 0.30793983876412, 0.07275791684542/2.0});
            gaussPoints.push_back({0.55035294182100, 0.30793983876412, 0.07275791684542/2.0});
            gaussPoints.push_back({0.30793983876412, 0.14170721941488, 0.07275791684542/2.0});
            gaussPoints.push_back({0.14170721941488, 0.55035294182100, 0.07275791684542/2.0});
            gaussPoints.push_back({0.24667256063990, 0.72832390459741, 0.02832724253106/2.0});
            gaussPoints.push_back({0.72832390459741, 0.02500353476269, 0.02832724253106/2.0});
            gaussPoints.push_back({0.02500353476269, 0.24667256063990, 0.02832724253106/2.0});
            gaussPoints.push_back({0.72832390459741, 0.24667256063990, 0.02832724253106/2.0});
            gaussPoints.push_back({0.24667256063990, 0.02500353476269, 0.02832724253106/2.0});
            gaussPoints.push_back({0.02500353476269, 0.72832390459741, 0.02832724253106/2.0});
            gaussPoints.push_back({0.06665406347960, 0.92365593358750, 0.00942166696373/2.0});
            gaussPoints.push_back({0.92365593358750, 0.00969000293290, 0.00942166696373/2.0});
            gaussPoints.push_back({0.00969000293290, 0.06665406347960, 0.00942166696373/2.0});
            gaussPoints.push_back({0.92365593358750, 0.06665406347960, 0.00942166696373/2.0});
            gaussPoints.push_back({0.06665406347960, 0.00969000293290, 0.00942166696373/2.0});
            gaussPoints.push_back({0.00969000293290, 0.92365593358750, 0.00942166696373/2.0});
            break;
                
        default:
            // 默认使用阶数10
            return generateTriangleGaussPoints(10);
        }
        
        return gaussPoints;
    }

    // 生成四边形单元的高斯积分点(Golub–Welsch 算法)
    std::vector<GaussPoint> generateRectangleGaussPoints(int order) {
        std::vector<GaussPoint> gaussPoints;

        std::vector<double> x(order), w(order);

        // 构造 Jacobi 三对角矩阵
        Eigen::VectorXd beta(order-1);
        for (int k = 1; k < order; ++k){
            beta(k-1) = k / std::sqrt(4.0*k*k - 1.0);
        }

        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(order,order);
        for (int i = 0; i < order-1; ++i) {
            J(i,i+1) = beta(i);
            J(i+1,i) = beta(i);
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
        const auto& evals = es.eigenvalues();
        const auto& evecs = es.eigenvectors();

        for (int i = 0; i < order; ++i) {
            x[i] = evals(i);
            w[i] = 2.0 * evecs(0,i) * evecs(0,i);
        }

        for (int i = 0; i < order; ++i) {
            for (int j = 0; j < order; ++j) {
                gaussPoints.push_back({x[i], x[j], w[i] * w[j]});
            }
        }

        return gaussPoints;
    }

} // namespace Integrate2DLegendre

} // namespace Auxiliary

} // namespace Poisson
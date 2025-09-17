// UsingAlias.h
#ifndef USING_ALIAS_H
#define USING_ALIAS_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <functional>

namespace Poisson {
    using fuxy = std::function<double(double u_val, double x_val, double y_val)>;
    using fu = std::function<double(double u_val)>;
    using fxy = std::function<double(double x_val, double y_val)>;

    using vec_t = Eigen::VectorXd;
    using vec_t2 = Eigen::Vector2d;
    using mat_t = Eigen::MatrixXd;
    using mat_t2 = Eigen::Matrix<double, Eigen::Dynamic, 2>;
    using smat_t = Eigen::SparseMatrix<double>;

    using NodeCoords = std::vector<Eigen::Vector2d>;
}

#endif  //USING_ALIAS_H


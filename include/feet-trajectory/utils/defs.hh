#pragma once
#include <Eigen/Core>

namespace feettrajectory
{
typedef Eigen::Ref<Eigen::MatrixXd> RefMat;
typedef Eigen::Ref<const Eigen::MatrixXd>& ConstRefMat;

template <typename T>
using RefVecTmpl = Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, 1> >;
typedef RefVecTmpl<double> RefVec;
typedef RefVecTmpl<int> RefVecInt;

template <typename T>
using ConstRefVecTmpl = Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1> >&;
typedef ConstRefVecTmpl<double> ConstRefVec;
typedef ConstRefVecTmpl<int> ConstRefVecInt;

typedef Eigen::Ref<Eigen::Vector3d> RefVec3d;
typedef const Eigen::Ref<Eigen::Vector3d>& ConstRefVec3d;

typedef Eigen::Ref<Eigen::Vector4d> RefVec4d;
typedef const Eigen::Ref<Eigen::Vector4d>& ConstRefVec4d;

typedef Eigen::VectorXd::Index Index;
typedef Eigen::Transpositions<Eigen::Dynamic, Eigen::Dynamic, int>
    Transpositions_t;

namespace fmt
{
static const Eigen::IOFormat custom(3, 0, ", ", "\n", "[", "]");
static const Eigen::IOFormat matlab(14, 0, ", ", ";\n", "[", "]", "[", "];");
static const Eigen::IOFormat matlabVector(14, 0, ", ", ";", "", "", "[", "];");
static const Eigen::IOFormat commaInit(14, 0, ", ", ", ", "", "", " << ", ";");
}  // end of namespace defs

} /* feettrajectory */


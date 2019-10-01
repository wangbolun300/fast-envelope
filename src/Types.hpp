#pragma once


#include <Eigen/Dense>

namespace fastEnvelope {
#ifdef FAST_ENVELOPE_USE_FLOAT
    typedef float Scalar;
#define SCALAR_ZERO 1e-6
#define SCALAR_ZERO_2 1e-12
#define SCALAR_ZERO_3 1e-18
#else
    typedef double Scalar;
#define SCALAR_ZERO 1e-8
#define SCALAR_ZERO_2 1e-16
#define SCALAR_ZERO_3 1e-24
#endif

    typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;

    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef Eigen::Matrix<Scalar, 2, 1> Vector2;


    typedef Eigen::Matrix<int, 4, 1> Vector4i;
    typedef Eigen::Matrix<int, 3, 1> Vector3i;
    typedef Eigen::Matrix<int, 2, 1> Vector2i;

    typedef Eigen::Matrix < Vector3, Eigen::Dynamic, 1, 0, 12, 1> HalfPlanes;
    //HalfPlanes hp(4, 1); hp(11, 1); no: hp(13, 1);
    //hp.size() 4, 11
}

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


	static const bool OUT_PRISM = 1;
	static const bool IN_PRISM = 0;

	static const int FE_CUT_COPLANAR = 4;
	static const int FE_CUT_EMPTY = -1;
	static const int FE_CUT_FACE = 3;

	static const int NOT_DEGENERATED = 0;
	static const int NERLY_DEGENERATED = 1;
	static const int DEGENERATED_SEGMENT = 2;
	static const int DEGENERATED_POINT = 3;

    typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;

    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef Eigen::Matrix<Scalar, 2, 1> Vector2;


    typedef Eigen::Matrix<int, 4, 1> Vector4i;
    typedef Eigen::Matrix<int, 3, 1> Vector3i;
    typedef Eigen::Matrix<int, 2, 1> Vector2i;

}

// Uncomment the following to distinguish between UNCERTAIN_POSITION and UNCERTAIN_CONFIGURATION
//#define USE_FILTER_FOR_DEGENERATE_CONFIGS

enum Filtered_Orientation {
	POSITIVE = 1,					// Intersection point is OVER the reference plane
	NEGATIVE = -1,					// Intersection point is UNDER the reference plane
	UNCERTAIN_POSITION = 0,			// Precision is not sufficient (intersection point too close to coplanarity)
#ifdef USE_FILTER_FOR_DEGENERATE_CONFIGS
	UNCERTAIN_CONFIGURATION = -2	// Precision is not sufficient (point configuration too close to degeneracy)
#else
	UNCERTAIN_CONFIGURATION = 0	// Precision is not sufficient (point configuration too close to degeneracy)
#endif
};

//
// orient3D_LPI - Indirect exact 3D orientation predicate with floating point filter.
//
// Input: eight points p, q, r, s, t, a, b, c
// where:
// <p, q> define a straight line L
// <r, s, t> define a plane P1
// <a, b, c> define a plane P2
//
// Output:
// Let i be the exact intersection point between L and P1.
// orient3D_LPI returns POSITIVE if the tetrahedron (i, a, b, c) has positive volume.
// Returns NEGATIVE if it has negative volume.
// Returns UNCERTAIN_CONFIGURATION if point configuration is degenerate (L and P1 are 
// parallel, <r, s, t> are collinear, ...) or nearly so (FP precision not enough to guarantee)
// Returns UNCERTAIN_POSITION if one of the following conditions holds:
// - i is exactly on P2
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

int orient3D_LPI_filtered(
	double px, double py, double pz, double qx, double qy, double qz,
	double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz,
	double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz);

//
// orient3D_TPI - Indirect exact 3D orientation predicate with floating point filter.
//
// Input: twelve points v1, v2, v3, w1, w2, w3, u1, u2, u3, q1, q2, q3
// where:
// the three vi define a plane PV
// the three wi define a plane PW
// the three ui define a plane PU
// the three qi define a plane PQ
//
// Output:
// Let i be the exact intersection point of the three planes PV, PW, PU
// orient3D_TPI returns POSITIVE if the tetrahedron (i, q1, q2, q3) has positive volume.
// Returns NEGATIVE if it has negative volume.
// Returns UNCERTAIN_CONFIGURATION if point configuration is degenerate (i is undefined) 
// or nearly so (FP precision not enough to guarantee).
// Returns UNCERTAIN_POSITION if one of the following conditions holds:
// - i is exactly on PQ
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

int orient3D_TPI_filtered(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z);


//////////////////////////////////////////////////////////////////////////////////////////
//
// Split versions - use these functions if the same intersection point 'i' must be checked
// against several reference planes.
//
//
// Example for LPI (checks the orientation of 'i' wrt <a1,b1,c1> and <a2,b2,c2>):
//
//  /* The following coordinates implicitly define the intersection point 'i' */
//	double px, py, pz;
//	double qx, qy, qz;
//	double rx, ry, rz;
//	double sx, sy, sz;
//	double tx, ty, tz;
//
//  /* This is the first reference plane */
//	double a1x, a1y, a1z;
//	double b1x, b1y, b1z;
//	double c1x, c1y, c1z;
//
//  /* This is the second reference plane */
//	double a2x, a2y, a2z;
//	double b2x, b2y, b2z;
//	double c2x, c2y, c2z;
//
//  ... initialize the above variables somehow ...
//
//  /* These are support variables that must be declared and will be initialized by orient3D_LPI_prefilter() */
//	double a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
//
//  int orientation_1, orientation_2;
//	if (orient3D_LPI_prefilter(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5))
//  {
//   orientation_1 = orient3D_LPI_postfilter(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5, px, py, pz, a1x, a1y, a1z, b1x, b1y, b1z, c1x, c1y, c1z);
//   orientation_2 = orient3D_LPI_postfilter(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5, px, py, pz, a2x, a2y, a2z, b2x, b2y, b2z, c2x, c2y, c2z);
//  }
//	else orientation_1 = orientation_2 = Filtered_Orientation::UNCERTAIN_CONFIGURATION;
//
//
//////////////////////////////////////////////////////////////////////////////////////////

//
// orient3D_LPI_prefilter
//
// Input: five points p, q, r, s, t
// where:
// <p, q> define a straight line L
// <r, s, t> define a plane P1
//
// Output:
// ten double values a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5
// These values must be passed as input to subsequent calls of orient3D_LPI_postfilter()
// Returns TRUE on success. FALSE if either overflow occurs or precision is insufficient.

bool orient3D_LPI_prefilter(
	const double& px, const double& py, const double& pz,
	const double& qx, const double& qy, const double& qz,
	const double& rx, const double& ry, const double& rz,
	const double& sx, const double& sy, const double& sz,
	const double& tx, const double& ty, const double& tz,
	double& a11, double& a12, double& a13, double& d, double& fa11, double& fa12, double& fa13, double& max1, double& max2, double& max5);

//
// orient3D_LPI_postfilter
//
// Input: ten double values a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5 and four points p, a, b, c
// where:
// the ten double values are those initialized by orient3D_LPI_prefilter()
// p is the first point defining the line L in orient3D_LPI_prefilter()
// <a, b, c> is the reference plane against which 'i' is checked for orientation
//
// Output:
// Returns POSITIVE if the tetrahedron (i, a, b, c) has positive volume.
// Returns NEGATIVE if it has negative volume.
// Returns UNCERTAIN_POSITION if one of the following conditions holds:
// - i is exactly on the plane by <a, b, c>
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

int orient3D_LPI_postfilter(
	const double& a11, const double& a12, const double& a13, const double& d,
	const double& fa11, const double& fa12, const double& fa13,
	const double& max1, const double& max2, const double& max5,
	const double& px, const double& py, const double& pz,
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz,
	const double& cx, const double& cy, const double& cz);


//
// orient3D_TPI_prefilter
//
// Input: nine points v1, v2, v3, w1, w2, w3, u1, u2, u3
// where:
// the three vi define a plane PV
// the three wi define a plane PW
// the three ui define a plane PU
//
// Output:
// eleven double values d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7
// These values must be passed as input to subsequent calls of orient3D_TPI_postfilter()
// Returns TRUE on success. FALSE if either overflow occurs or precision is insufficient.

bool orient3D_TPI_prefilter(
	const double& v1x, const double& v1y, const double& v1z, const double& v2x, const double& v2y, const double& v2z, const double& v3x, const double& v3y, const double& v3z,
	const double& w1x, const double& w1y, const double& w1z, const double& w2x, const double& w2y, const double& w2z, const double& w3x, const double& w3y, const double& w3z,
	const double& u1x, const double& u1y, const double& u1z, const double& u2x, const double& u2y, const double& u2z, const double& u3x, const double& u3y, const double& u3z,
	double& d, double& n1, double& n2, double& n3, double& max1, double& max2, double& max3, double& max4, double& max5, double& max6, double& max7);


//
// orient3D_TPI_postfilter
//
// Input: eleven double values d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7
//        and three points q1, q2, q3
// where:
// the eleven double values are those initialized by orient3D_TPI_prefilter()
// the three points define a reference plane against which 'i' is checked for orientation.
//
// Output:
// Returns POSITIVE if the tetrahedron (i, q1, q2, q3) has positive volume.
// Returns NEGATIVE if it has negative volume.
// Returns UNCERTAIN_POSITION if one of the following conditions holds:
// - i is exactly on the plane by <q1, q2, q3>
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

int orient3D_TPI_postfilter(
	const double& d, const double& n1, const double& n2, const double& n3,
	const double& max1, const double& max2, const double& max3, const double& max4, const double& max5, const double& max6, const double& max7,
	const double& q1x, const double& q1y, const double& q1z, const double& q2x, const double& q2y, const double& q2z, const double& q3x, const double& q3y, const double& q3z);

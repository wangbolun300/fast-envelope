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

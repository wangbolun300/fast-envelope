/****************************************************************************
* Indirect exact 3D orientation predicates with floating point filter       *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2019: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU Lesser General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or (at  *
* your option) any later version.                                           *
*                                                                           *
* This program is distributed in the hope that it will be useful, but       *
* WITHOUT ANY WARRANTY; without even the implied warranty of                *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  *
* General Public License for more details.                                  *
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  *
* along with this program.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/

/****************************************************************************

This code implements indirect 3D orientation predicates.

A 3D orientation predicate takes four points (i, a, b, c) and
returns the sign of the volume of the tetrahedron formed by i,a,b and c.
Traditionally, each of the four points is explicitly represented by its 
three coordinates.

If a,b,c are not aligned they define an oriented plane.
In this case, the sign of orient3d(i,a,b,c) states whether 'i' is above
this plane (>0), on the plane (=0), or below the plane (<0).

In this code, the point 'i' is not represented by its three cordinates.
Rather, it is implicitly represented as the intersection of other geometric
entities.

Currently, 'i' can be represented as a:
1) line-plane intersection (LPI)
2) three planes intersection (TPI)

In turn, lines and planes are represented by pairs and triples of explicit 
points respectively. Here, 'explicit' means that the point is represented
by its three coordinates.

Before using any other function from this code, initFPU() must be called
to properly set the FPU configuration. This guarantees that the exact
versions of the functions produce the correct result on IEEE-754 compliant
machines.

Then, the easiest way to use this code is by calling either orient3D_LPI()
or orient3D_TPI().

Three versions are provided for each of the above two functions 
(* = LPI / TPI).

1) orient3D_*_Filtered
Fast. Might return an UNCERTAIN answer in near-degenerate cases.

2) orient3D_*_Exact
Slow. Employs floating point expansions to ensure an exact answer in any case

3) orient3D_*
Fast and always exact. Uses orient3D_Filtered for most of the cases, while
switching to the Exact version only when necessary.


When an implicit point 'i' must be checked against more than one reference 
plane, avoid repeating part of the calculations by using the split versions
of the predicates:

orient3D_*_prefilter
orient3D_*_postfilter
orient3D_*_pre_exact
orient3D_*_post_exact

where the "orient3D_*_pre*" function generates common support variables to
be used by the "orient3D_*_post*" predicate to do the actual check for all 
the reference planes (see detailed comments below for usage examples).

****************************************************************************/

// Uncomment the following to use multi-stage filters instead of almost-static filter only
//#define USE_MULTISTAGE_FILTERS

//
// initFPU - Make sure to call this function before using the other functions in your code.
// This appropriately sets the FPU mode to make sure that the error calculations are correct.
// It is sufficient to call this function once at the beginning of your program.
//

void initFPU();

#ifdef USE_MULTISTAGE_FILTERS
// Interval_number
class interval_number
{
	typedef union error_approx_type_t
	{
		double d;
		uint64_t u;

		inline error_approx_type_t() {}
		inline error_approx_type_t(double a) : d(a) {}
		inline uint64_t is_negative() const { return u >> 63; }
	} casted_double;

	double low, high;

public:
	inline interval_number() {}
	inline interval_number(double a) : low(a), high(a) {}
	inline interval_number(double a, double b) : low(a), high(b) {}
	inline interval_number(const interval_number& b) : low(b.low), high(b.high) {}

	inline bool signIsReliable() const { return (low>0 || high <0 || (low == 0 && high == 0)); }
	inline int sign() const { return (low > 0) ? (1) : ((low < 0) ? (-1) : (0)); }

	inline interval_number& operator=(const interval_number& b) { low = b.low; high = b.high; return *this; }
	inline interval_number operator+(const interval_number& b) const { return interval_number(-((-low) - b.low), high + b.high); }
	inline interval_number operator-(const interval_number& b) const { return interval_number(-(b.high - low), high - b.low); }
	interval_number operator*(const interval_number& b) const;
};
#endif

enum Filtered_Orientation {
	POSITIVE = 1,					// Intersection point is OVER the reference plane
	NEGATIVE = -1,					// Intersection point is UNDER the reference plane
	UNCERTAIN = 0					// Precision is not sufficient
};

//
// orient3D_LPI_filtered
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
// Returns UNCERTAIN if point configuration is degenerate (L and P1 are 
// parallel, <r, s, t> are collinear, ...) or nearly so (FP precision not enough to guarantee)
// or if one of the following conditions holds:
// - i is exactly on P2
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

int orient3D_LPI_filtered(
	double px, double py, double pz, double qx, double qy, double qz,
	double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz,
	double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz);

//
// orient3D_TPI_filtered
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
// Returns UNCERTAIN if point configuration is degenerate (i is undefined) 
// or nearly so (FP precision not enough to guarantee).
// or if one of the following conditions holds:
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
//	LPI_filtered_suppvars s;
//
//  int orientation_1, orientation_2;
//	if (orient3D_LPI_prefilter(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, s))
//  {
//   orientation_1 = orient3D_LPI_postfilter(s, px, py, pz, a1x, a1y, a1z, b1x, b1y, b1z, c1x, c1y, c1z);
//   orientation_2 = orient3D_LPI_postfilter(s, px, py, pz, a2x, a2y, a2z, b2x, b2y, b2z, c2x, c2y, c2z);
//  }
//	else orientation_1 = orientation_2 = Filtered_Orientation::UNCERTAIN;
//
//
//////////////////////////////////////////////////////////////////////////////////////////

class LPI_filtered_suppvars
{
public:
	double a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
#ifdef USE_MULTISTAGE_FILTERS
	interval_number ia11, ia12, ia13, id;
	double qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz;
#endif

	LPI_filtered_suppvars() {}
};

//
// orient3D_LPI_prefilter
//
// Input: five points p, q, r, s, t
// where:
// <p, q> define a straight line L
// <r, s, t> define a plane P1
//
// Output:
// Inits the support object s.
// This so-initialized object must be passed as input to subsequent calls of orient3D_LPI_postfilter()
// Returns TRUE on success. FALSE if either overflow occurs or precision is insufficient.

bool orient3D_LPI_prefilter(
	const double& px, const double& py, const double& pz,
	const double& qx, const double& qy, const double& qz,
	const double& rx, const double& ry, const double& rz,
	const double& sx, const double& sy, const double& sz,
	const double& tx, const double& ty, const double& tz,
	LPI_filtered_suppvars& s);

//
// orient3D_LPI_postfilter
//
// Input: a support object s and four points p, a, b, c
// where:
// s was previously initialized by orient3D_LPI_prefilter()
// p is the first point defining the line L in orient3D_LPI_prefilter()
// <a, b, c> is the reference plane against which 'i' is checked for orientation
//
// Output:
// Returns POSITIVE if the tetrahedron (i, a, b, c) has positive volume.
// Returns NEGATIVE if it has negative volume.
// Returns UNCERTAIN if one of the following conditions holds:
// - i is exactly on the plane by <a, b, c>
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

int orient3D_LPI_postfilter(
	const LPI_filtered_suppvars& s,
	const double& px, const double& py, const double& pz,
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz,
	const double& cx, const double& cy, const double& cz);

class TPI_filtered_suppvars
{
public:
	double d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7;
#ifdef USE_MULTISTAGE_FILTERS
	interval_number id, in1, in2, in3;
	double v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
		w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
		u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z;
#endif

	TPI_filtered_suppvars() {}
};

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
// Inits the support object s.
// This so-initialized object must be passed as input to subsequent calls of orient3D_TPI_postfilter()
// Returns TRUE on success. FALSE if either overflow occurs or precision is insufficient.

bool orient3D_TPI_prefilter(
	const double& v1x, const double& v1y, const double& v1z, const double& v2x, const double& v2y, const double& v2z, const double& v3x, const double& v3y, const double& v3z,
	const double& w1x, const double& w1y, const double& w1z, const double& w2x, const double& w2y, const double& w2z, const double& w3x, const double& w3y, const double& w3z,
	const double& u1x, const double& u1y, const double& u1z, const double& u2x, const double& u2y, const double& u2z, const double& u3x, const double& u3y, const double& u3z,
	TPI_filtered_suppvars& s);


//
// orient3D_TPI_postfilter
//
// Input: a support object s and three points q1, q2, q3
// where:
// s was previously initialized by orient3D_TPI_prefilter()
// the three points define a reference plane against which 'i' is checked for orientation.
//
// Output:
// Returns POSITIVE if the tetrahedron (i, q1, q2, q3) has positive volume.
// Returns NEGATIVE if it has negative volume.
// Returns UNCERTAIN if one of the following conditions holds:
// - i is exactly on the plane by <q1, q2, q3>
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

int orient3D_TPI_postfilter(
	const TPI_filtered_suppvars& s,
	const double& q1x, const double& q1y, const double& q1z, const double& q2x, const double& q2y, 
	const double& q2z, const double& q3x, const double& q3y, const double& q3z);



//
// orient3D_LPI_exact
//
// Input: eight points p, q, r, s, t, a, b, c
// where:
// <p, q> define a straight line L
// <r, s, t> define a plane P1
// <a, b, c> define a plane P2
//
// Output:
// Let i be the exact intersection point between L and P1.
// orient3D_LPI_exact returns 1 if the tetrahedron T = (i, a, b, c) has positive volume,
// -1 if T has negative volume.
// Zero is returned if one of the following conditions holds:
// - degenerate input (L and P1 are parallel, <r, s, t> are collinear, ...)
// - i is exactly on P2 (i.e. T has zero volume).

int orient3D_LPI_exact(
	double px, double py, double pz, double qx, double qy, double qz,
	double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz,
	double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz);

//
// orient3D_TPI_exact
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
// orient3D_TPI_exact returns 1 if the tetrahedron T = (i, q1, q2, q3) has positive volume.
// Returns -1 if it has negative volume.
// Zero is returned if one of the following conditions holds:
// - degenerate input (i is undefined)
// - i is exactly on PQ (i.e. T has zero volume).

int orient3D_TPI_exact(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z);


//////////////////////////////////////////////////////////////////////////////////////////
//
// Orient3D_LPI_exact (split version) - use these functions if the same intersection point 
// 'i' must be checked against several reference planes.
//
//
// Example (checks the orientation of 'i' wrt <a1,b1,c1> and <a2,b2,c2>):
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
//  /* The following object contains support variables */
//  LPI_exact_suppvars s;
//
//  int orientation_1, orientation_2;
//	if (orient3D_LPI_pre_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, s))
//  {
//   orientation_1 = orient3D_LPI_post_exact(s, px, py, pz, 
//											 a1x, a1y, a1z, b1x, b1y, b1z, c1x, c1y, c1z);
//   orientation_2 = orient3D_LPI_post_exact(s, px, py, pz, 
//											 a2x, a2y, a2z, b2x, b2y, b2z, c2x, c2y, c2z);
//  }
//	else orientation_1 = orientation_2 = 0; // 'i' is undefined
//
//
//////////////////////////////////////////////////////////////////////////////////////////

class LPI_exact_suppvars
{
public:
	double a11[2], a12[2], a13[2], d[192], n[192];
	int dl, nl;

	LPI_exact_suppvars() {}
};

bool orient3D_LPI_pre_exact(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	LPI_exact_suppvars& s);

int orient3D_LPI_post_exact(
	LPI_exact_suppvars& s,
	double px, double py, double pz,
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz);


//////////////////////////////////////////////////////////////////////////////////////////
//
// Orient3D_TPI_exact (split version) - use these functions if the same intersection point 
// 'i' must be checked against several reference planes.
//
//
// Example (checks the orientation of 'i' wrt <a1,b1,c1> and <a2,b2,c2>):
//
//  /* The following coordinates implicitly define the intersection point 'i' */
//  double v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z; // First intersecting plane
//  double w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z; // Second intersecting plane
//  double u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z; // Third intersecting plane
//
//  /* This is the first reference plane */
//	double a1x, a1y, a1z, b1x, b1y, b1z, c1x, c1y, c1z;
//
//  /* This is the second reference plane */
//	double a2x, a2y, a2z, b2x, b2y, b2z, c2x, c2y, c2z;
//
//  ... initialize the above variables somehow ...
//
//  /* The following object contains support variables */
//  TPI_exact_suppvars s;
//
//  int orientation_1, orientation_2;
//	if (orient3D_TPI_pre_exact(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, 
//							   w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z, 
//							   u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, 
//							   s))
//  {
//   orientation_1 = orient3D_TPI_post_exact(s, 
//											 a1x, a1y, a1z, b1x, b1y, b1z, c1x, c1y, c1z);
//   orientation_2 = orient3D_TPI_post_exact(s, 
//											 a2x, a2y, a2z, b2x, b2y, b2z, c2x, c2y, c2z);
//  }
//	else orientation_1 = orientation_2 = 0; // 'i' is undefined
//  
//  /* PLEASE NOTE  */
//  /* Re-using 's' for another orient3D_TPI_pre_exact() call may cause memory leaks. */
//  /* Create another one, or dynamically allocate 's' and destroy it by 'new/delete' */
//  /* for each call.                                                                 */
//
//
//////////////////////////////////////////////////////////////////////////////////////////

class TPI_exact_suppvars
{
public:
	double dp[256], n1p[256], n2p[256], n3p[256];
	double *d, *n1, *n2, *n3;
	int dl, n1l, n2l, n3l;

	TPI_exact_suppvars();
	~TPI_exact_suppvars();
};

bool orient3D_TPI_pre_exact(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	TPI_exact_suppvars& s
	);

int orient3D_TPI_post_exact(
	TPI_exact_suppvars& s,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z
	);




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
// orient3D_LPI returns 1 if the tetrahedron T = (i, a, b, c) has positive volume,
// -1 if T has negative volume.
// Zero is returned if one of the following conditions holds:
// - degenerate input (L and P1 are parallel, <r, s, t> are collinear, ...)
// - i is exactly on P2 (i.e. T has zero volume).

int orient3D_LPI(
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
// orient3D_TPI returns 1 if the tetrahedron T = (i, q1, q2, q3) has positive volume.
// Returns -1 if it has negative volume.
// Zero is returned if one of the following conditions holds:
// - degenerate input (i is undefined)
// - i is exactly on PQ (i.e. T has zero volume).

int orient3D_TPI(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z);


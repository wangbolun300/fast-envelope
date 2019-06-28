#include <math.h>
#include <fenv.h>
#include <stdio.h>
#include <FastEnvelope/ip_filtered.h>
#include<arbitraryprecision/fprecision.h>

// Uncomment the following to activate a second filter on the determinants 'd' in both the following functions.
// If uncommented, the two functions may return -2 if 'd' is too close to zero, meaning that
// the intersection point is close to infinity.
#define USE_FILTER_FOR_DEGENERATE_CONFIGS
namespace fastEnvelope

{
	// Indirect 3D orientation predicate with floating point filter.
	//
	// Input: eight points p, q, r, s, t, a, b, c
	// where:
	// <p, q> define a straight line L
	// <r, s, t> define a plane P1
	// <a, b, c> define a plane P2
	//
	// Output:
	// Let i be the exact intersection point between L and P1.
	// orient3D_LPI returns 1 if the tetrahedron (i, a, b, c) has positive volume.
	// Returns -1 if it has negative volume.
	// Zero is returned if one of the following conditions holds:
	// - degenerate input (L and P1 are parallel, <r, s, t> are collinear, ...)
	// - i is exactly on P2
	// - the floating point precision is insufficient to guarantee an exact answer
	// - the input coordinates cause an under/overflow during the computation.

	int ip_filtered::orient3D_LPI_filtered(
		double px, double py, double pz,
		double qx, double qy, double qz,
		double rx, double ry, double rz,
		double sx, double sy, double sz,
		double tx, double ty, double tz,
		double ax, double ay, double az,
		double bx, double by, double bz,
		double cx, double cy, double cz)
	{
		::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

		double a11, a12, a13, a21, a22, a23, a31, a32, a33, d21, d31, d22, d32, d23, d33;
		double px_rx, py_ry, pz_rz, px_cx, py_cy, pz_cz;
		double a2233, a2133, a2132;
		double d, n;
		double d11, d12, d13;
		double d2233, d2332, d2133, d2331, d2132, d2231;
		double det;

		a11 = (px - qx);
		a12 = (py - qy);
		a13 = (pz - qz);
		a21 = (sx - rx);
		a22 = (sy - ry);
		a23 = (sz - rz);
		a31 = (tx - rx);
		a32 = (ty - ry);
		a33 = (tz - rz);
		a2233 = ((a22 * a33) - (a23 * a32));
		a2133 = ((a21 * a33) - (a23 * a31));
		a2132 = ((a21 * a32) - (a22 * a31));
		d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));

		// The almost static filter for 'd' might be moved here

		px_rx = px - rx;
		py_ry = py - ry;
		pz_rz = pz - rz;
		n = ((((py_ry)* a2133) - ((px_rx)* a2233)) - ((pz_rz)* a2132));

		px_cx = px - cx;
		py_cy = py - cy;
		pz_cz = pz - cz;

		d11 = (d * px_cx) + (a11 * n);
		d21 = (ax - cx);
		d31 = (bx - cx);
		d12 = (d * py_cy) + (a12 * n);
		d22 = (ay - cy);
		d32 = (by - cy);
		d13 = (d * pz_cz) + (a13 * n);
		d23 = (az - cz);
		d33 = (bz - cz);

		d2233 = d22 * d33;
		d2332 = d23 * d32;
		d2133 = d21 * d33;
		d2331 = d23 * d31;
		d2132 = d21 * d32;
		d2231 = d22 * d31;

		det = d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231);

		if (::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID)) return 0; // Fast reject in case of under/overflow

		// An additional static filter here may be advantageous... possibly initialized based on input coordinates

		// Almost static filter
		double fa11 = fabs(a11);
		double fa21 = fabs(a21);
		double fa31 = fabs(a31);
		double fa12 = fabs(a12);
		double fa22 = fabs(a22);
		double fa32 = fabs(a32);
		double fa13 = fabs(a13);
		double fa23 = fabs(a23);
		double fa33 = fabs(a33);

		double max1, max2, max5;

		max1 = fa23;
		if (max1 < fa13) max1 = fa13;
		if (max1 < fa33) max1 = fa33;

		max2 = fa12;
		if (max2 < fa22) max2 = fa22;
		if (max2 < fa32) max2 = fa32;

		max5 = fa11;
		if (max5 < fa21) max5 = fa21;
		if (max5 < fa31) max5 = fa31;

#ifdef USE_FILTER_FOR_DEGENERATE_CONFIGS
		double deps = 5.1107127829973299e-015 * max1 * max2 * max5;
		if (d <= deps && d >= -deps) return -2;
#endif

		double fpxrx = fabs(px_rx);
		double fpyry = fabs(py_ry);
		double fpzrz = fabs(pz_rz);
		double fd11 = fabs(d11);
		double fd21 = fabs(d21);
		double fd31 = fabs(d31);
		double fd12 = fabs(d12);
		double fd22 = fabs(d22);
		double fd32 = fabs(d32);
		double fd13 = fabs(d13);
		double fd23 = fabs(d23);
		double fd33 = fabs(d33);
		double fpxcx = fabs(px_cx);
		double fpycy = fabs(py_cy);
		double fpzcz = fabs(pz_cz);

		if (max1 < fpzrz) max1 = fpzrz;
		if (max2 < fpyry) max2 = fpyry;
		if (max5 < fpxrx) max5 = fpxrx;

		double max3, max4, max6;

		max3 = fa12;
		if (max3 < fd32) max3 = fd32;
		if (max3 < fd22) max3 = fd22;
		if (max3 < fpycy) max3 = fpycy;

		max4 = fa13;
		if (max4 < fpzcz) max4 = fpzcz;
		if (max4 < fd33) max4 = fd33;
		if (max4 < fd23) max4 = fd23;

		max6 = fa11;
		if (max6 < fd21) max6 = fd21;
		if (max6 < fd31) max6 = fd31;
		if (max6 < fpxcx) max6 = fpxcx;

		double eps = 1.3865993466947057e-013 * max2 * max1 * max5 * max6 * max3 * max4;

		if ((det > eps)) return (d > 0) ? (1) : (-1);
		if ((det < -eps)) return (d > 0) ? (-1) : (1);
		return 0;
	}


	// Indirect 3D orientation predicate with floating point filter.
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
	// orient3D_TPI returns 1 if the tetrahedron (i, q1, q2, q3) has positive volume.
	// Returns -1 if it has negative volume.
	// Zero is returned if one of the following conditions holds:
	// - degenerate input (i is undefined)
	// - i is exactly on PQ
	// - the floating point precision is insufficient to guarantee an exact answer
	// - the input coordinates cause an under/overflow during the computation.

	int ip_filtered::orient3D_TPI_filtered(
		double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
		double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
		double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
		double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z
	)
	{
		::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

		v3x -= v2x;
		v3y -= v2y;
		v3z -= v2z;
		v2x -= v1x;
		v2y -= v1y;
		v2z -= v1z;
		double nvx = v2y * v3z - v2z * v3y;
		double nvy = v3x * v2z - v3z * v2x;
		double nvz = v2x * v3y - v2y * v3x;

		w3x -= w2x;
		w3y -= w2y;
		w3z -= w2z;
		w2x -= w1x;
		w2y -= w1y;
		w2z -= w1z;
		double nwx = w2y * w3z - w2z * w3y;
		double nwy = w3x * w2z - w3z * w2x;
		double nwz = w2x * w3y - w2y * w3x;

		u3x -= u2x;
		u3y -= u2y;
		u3z -= u2z;
		u2x -= u1x;
		u2y -= u1y;
		u2z -= u1z;
		double nux = u2y * u3z - u2z * u3y;
		double nuy = u3x * u2z - u3z * u2x;
		double nuz = u2x * u3y - u2y * u3x;

		double nwyuz = nwy * nuz - nwz * nuy;
		double nwxuz = nwx * nuz - nwz * nux;
		double nwxuy = nwx * nuy - nwy * nux;

		double nvyuz = nvy * nuz - nvz * nuy;
		double nvxuz = nvx * nuz - nvz * nux;
		double nvxuy = nvx * nuy - nvy * nux;

		double nvywz = nvy * nwz - nvz * nwy;
		double nvxwz = nvx * nwz - nvz * nwx;
		double nvxwy = nvx * nwy - nvy * nwx;

		double d = nvx * nwyuz - nvy * nwxuz + nvz * nwxuy;

		double p1 = nvx * v1x + nvy * v1y + nvz * v1z;
		double p2 = nwx * w1x + nwy * w1y + nwz * w1z;
		double p3 = nux * u1x + nuy * u1y + nuz * u1z;

		double n1 = p1 * nwyuz - p2 * nvyuz + p3 * nvywz;
		double n2 = p2 * nvxuz - p3 * nvxwz - p1 * nwxuz;
		double n3 = p3 * nvxwy - p2 * nvxuy + p1 * nwxuy;

		double dq3x = d * q3x;
		double dq3y = d * q3y;
		double dq3z = d * q3z;

		double a11 = n1 - dq3x;
		double a12 = n2 - dq3y;
		double a13 = n3 - dq3z;
		double a21 = q1x - q3x;
		double a22 = q1y - q3y;
		double a23 = q1z - q3z;
		double a31 = q2x - q3x;
		double a32 = q2y - q3y;
		double a33 = q2z - q3z;

		double det = a11 * (a22*a33 - a23 * a32) - a12 * (a21*a33 - a23 * a31) + a13 * (a21*a32 - a22 * a31);

		if (::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID)) return 0; // Fast reject in case of under/overflow

		// Almost static filter

		double fv2x = fabs(v2x);
		double fv2y = fabs(v2y);
		double fv2z = fabs(v2z);
		double fv3x = fabs(v3x);
		double fv3y = fabs(v3y);
		double fv3z = fabs(v3z);

		double fw2x = fabs(w2x);
		double fw2y = fabs(w2y);
		double fw2z = fabs(w2z);
		double fw3x = fabs(w3x);
		double fw3y = fabs(w3y);
		double fw3z = fabs(w3z);

		double fu2x = fabs(u2x);
		double fu3x = fabs(u3x);
		double fu2y = fabs(u2y);
		double fu2z = fabs(u2z);
		double fu3y = fabs(u3y);
		double fu3z = fabs(u3z);

		double max2, max4, max5, max7;

		max4 = fv2y;
		if (max4 < fv3y) max4 = fv3y;
		if (max4 < fw3y) max4 = fw3y;
		if (max4 < fw2y) max4 = fw2y;
		max2 = fv3x;
		if (max2 < fv2x) max2 = fv2x;
		if (max2 < fw2x) max2 = fw2x;
		if (max2 < fw3x) max2 = fw3x;
		max5 = fv2z;
		if (max5 < fv3z) max5 = fv3z;
		if (max5 < fw3z) max5 = fw3z;
		if (max5 < fw2z) max5 = fw2z;
		max7 = fu2x;
		if (max7 < fu3x) max7 = fu3x;
		if (max7 < fw2x) max7 = fw2x;
		if (max7 < fw3x) max7 = fw3x;

#ifdef USE_FILTER_FOR_DEGENERATE_CONFIGS
		double max9 = fu2y;
		if (max9 < fu3y) max9 = fu3y;
		if (max9 < fw2y) max9 = fw2y;
		if (max9 < fw3y) max9 = fw3y;
		double max10 = fu2y;
		if (max10 < fu3z) max10 = fu3z;
		if (max10 < fw2z) max10 = fw2z;
		if (max10 < fw3z) max10 = fw3z;

		double deps = 8.8881169117764924e-014 * (((((max4 * max5) * max2) * max10) * max7) * max9);
		if (d <= deps && d >= -deps) return -2;
#endif

		double fa21 = fabs(a21);
		double fa22 = fabs(a22);
		double fa23 = fabs(a23);
		double fa31 = fabs(a31);
		double fa32 = fabs(a32);
		double fa33 = fabs(a33);

		double max1, max3, max6, max8;

		if (max4 < fu2y) max4 = fu2y;
		if (max4 < fu3y) max4 = fu3y;
		if (max2 < fu2x) max2 = fu2x;
		if (max2 < fu3x) max2 = fu3x;
		if (max5 < fu2z) max5 = fu2z;
		if (max5 < fu3z) max5 = fu3z;
		if (max7 < fa21) max7 = fa21;
		if (max7 < fa31) max7 = fa31;

		max1 = max4;
		if (max1 < max2) max1 = max2;
		max3 = max5;
		if (max3 < max4) max3 = max4;
		max6 = fu2x;
		if (max6 < fu3x) max6 = fu3x;
		if (max6 < fu2z) max6 = fu2z;
		if (max6 < fw3y) max6 = fw3y;
		if (max6 < fw2x) max6 = fw2x;
		if (max6 < fw3z) max6 = fw3z;
		if (max6 < fw2y) max6 = fw2y;
		if (max6 < fw2z) max6 = fw2z;
		if (max6 < fu2y) max6 = fu2y;
		if (max6 < fu3z) max6 = fu3z;
		if (max6 < fu3y) max6 = fu3y;
		if (max6 < fw3x) max6 = fw3x;
		if (max6 < fa22) max6 = fa22;
		if (max6 < fa32) max6 = fa32;
		max8 = fa22;
		if (max8 < fa23) max8 = fa23;
		if (max8 < fa33) max8 = fa33;
		if (max8 < fa32) max8 = fa32;

		double eps = 3.4025182954957945e-012 * (((((((max1 * max3) * max2) * max5) * max7) * max4) * max6) * max8);

		if ((det > eps)) return (d > 0) ? (1) : (-1);
		if ((det < -eps)) return (d > 0) ? (-1) : (1);
		return 0;
	}



	//using multiple precision

//#include<arbitraryprecision/intervalprecision.h>
	int ip_filtered::orient3D_LPI_filtered_multiprecision(
		double px, double py, double pz,
		double qx, double qy, double qz,
		double rx, double ry, double rz,
		double sx, double sy, double sz,
		double tx, double ty, double tz,
		double ax, double ay, double az,
		double bx, double by, double bz,
		double cx, double cy, double cz)
	{
		::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

		double a11, a12, a13, a21, a22, a23, a31, a32, a33, d21, d31, d22, d32, d23, d33;
		double px_rx, py_ry, pz_rz, px_cx, py_cy, pz_cz;
		double a2233, a2133, a2132;
		double d, n;
		double d11, d12, d13;
		double d2233, d2332, d2133, d2331, d2132, d2231;
		double det;

		a11 = (px - qx);
		a12 = (py - qy);
		a13 = (pz - qz);
		a21 = (sx - rx);
		a22 = (sy - ry);
		a23 = (sz - rz);
		a31 = (tx - rx);
		a32 = (ty - ry);
		a33 = (tz - rz);
		a2233 = ((a22 * a33) - (a23 * a32));
		a2133 = ((a21 * a33) - (a23 * a31));
		a2132 = ((a21 * a32) - (a22 * a31));
		d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));

		// The almost static filter for 'd' might be moved here

		px_rx = px - rx;
		py_ry = py - ry;
		pz_rz = pz - rz;
		n = ((((py_ry)* a2133) - ((px_rx)* a2233)) - ((pz_rz)* a2132));

		px_cx = px - cx;
		py_cy = py - cy;
		pz_cz = pz - cz;

		d11 = (d * px_cx) + (a11 * n);
		d21 = (ax - cx);
		d31 = (bx - cx);
		d12 = (d * py_cy) + (a12 * n);
		d22 = (ay - cy);
		d32 = (by - cy);
		d13 = (d * pz_cz) + (a13 * n);
		d23 = (az - cz);
		d33 = (bz - cz);

		d2233 = d22 * d33;
		d2332 = d23 * d32;
		d2133 = d21 * d33;
		d2331 = d23 * d31;
		d2132 = d21 * d32;
		d2231 = d22 * d31;

		det = d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231);

		if (::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID)) return 0; // Fast reject in case of under/overflow

		// An additional static filter here may be advantageous... possibly initialized based on input coordinates

		// Almost static filter
		double fa11 = fabs(a11);
		double fa21 = fabs(a21);
		double fa31 = fabs(a31);
		double fa12 = fabs(a12);
		double fa22 = fabs(a22);
		double fa32 = fabs(a32);
		double fa13 = fabs(a13);
		double fa23 = fabs(a23);
		double fa33 = fabs(a33);

		double max1, max2, max5;

		max1 = fa23;
		if (max1 < fa13) max1 = fa13;
		if (max1 < fa33) max1 = fa33;

		max2 = fa12;
		if (max2 < fa22) max2 = fa22;
		if (max2 < fa32) max2 = fa32;

		max5 = fa11;
		if (max5 < fa21) max5 = fa21;
		if (max5 < fa31) max5 = fa31;

#ifdef USE_FILTER_FOR_DEGENERATE_CONFIGS
		double deps = 5.1107127829973299e-015 * max1 * max2 * max5;
		if (d <= deps && d >= -deps) return -2;
#endif

		double fpxrx = fabs(px_rx);
		double fpyry = fabs(py_ry);
		double fpzrz = fabs(pz_rz);
		double fd11 = fabs(d11);
		double fd21 = fabs(d21);
		double fd31 = fabs(d31);
		double fd12 = fabs(d12);
		double fd22 = fabs(d22);
		double fd32 = fabs(d32);
		double fd13 = fabs(d13);
		double fd23 = fabs(d23);
		double fd33 = fabs(d33);
		double fpxcx = fabs(px_cx);
		double fpycy = fabs(py_cy);
		double fpzcz = fabs(pz_cz);

		if (max1 < fpzrz) max1 = fpzrz;
		if (max2 < fpyry) max2 = fpyry;
		if (max5 < fpxrx) max5 = fpxrx;

		double max3, max4, max6;

		max3 = fa12;
		if (max3 < fd32) max3 = fd32;
		if (max3 < fd22) max3 = fd22;
		if (max3 < fpycy) max3 = fpycy;

		max4 = fa13;
		if (max4 < fpzcz) max4 = fpzcz;
		if (max4 < fd33) max4 = fd33;
		if (max4 < fd23) max4 = fd23;

		max6 = fa11;
		if (max6 < fd21) max6 = fd21;
		if (max6 < fd31) max6 = fd31;
		if (max6 < fpxcx) max6 = fpxcx;

		double eps = 1.3865993466947057e-013 * max2 * max1 * max5 * max6 * max3 * max4;

		if ((det > eps)) return (d > 0) ? (1) : (-1);
		if ((det < -eps)) return (d > 0) ? (-1) : (1);
		return 0;
	}



}
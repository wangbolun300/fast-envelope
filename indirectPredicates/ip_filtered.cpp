
///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// To compile on MSVC: use /fp:strict
// On GNU GCC: use -frounding-math
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <fenv.h>
#include "ip_filtered.h"

#pragma intrinsic(fabs)

#ifndef WIN32

#include <fenv.h>

#if __APPLE__
#define _FPU_SETCW(cw) // nothing https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/float.3.html
#else
#include <fpu_control.h>
#endif

#endif

#ifdef USE_MULTISTAGE_FILTERS

/////////////////// Dynamic filters with interval arithmetic /////////////

bool orient3D_LPI_pre_dfilter(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	LPI_filtered_suppvars &svs)
{
	setFPUModeToRoundUP();

	svs.ia11 = (px - qx);
	svs.ia12 = (py - qy);
	svs.ia13 = (pz - qz);
	interval_number a21(sx - rx);
	interval_number a22(sy - ry);
	interval_number a23(sz - rz);
	interval_number a31(tx - rx);
	interval_number a32(ty - ry);
	interval_number a33(tz - rz);
	interval_number a2233((a22 * a33) - (a23 * a32));
	interval_number a2133((a21 * a33) - (a23 * a31));
	interval_number a2132((a21 * a32) - (a22 * a31));
	svs.id = (((svs.ia11 * a2233) - (svs.ia12 * a2133)) + (svs.ia13 * a2132));

	if (!svs.id.signIsReliable())
	{
		setFPUModeToRoundNEAR();
		return false;
	}

	interval_number px_rx(px - rx);
	interval_number py_ry(py - ry);
	interval_number pz_rz(pz - rz);
	interval_number n(((((py_ry)*a2133) - ((px_rx)*a2233)) - ((pz_rz)*a2132)));

	svs.ia11 = (svs.ia11 * n);
	svs.ia12 = (svs.ia12 * n);
	svs.ia13 = (svs.ia13 * n);

	setFPUModeToRoundNEAR();

	if (!isfinite(svs.ia11) || !isfinite(svs.ia12) || !isfinite(svs.ia13))
		return false;

	svs.d = NAN; // This makes it possible to know which filter succeeded

	return true;
}

int orient3D_LPI_post_dfilter(const LPI_filtered_suppvars &svs,
							  double px, double py, double pz,
							  double ax, double ay, double az,
							  double bx, double by, double bz,
							  double cx, double cy, double cz)
{
	setFPUModeToRoundUP();

	interval_number px_cx(px - cx);
	interval_number py_cy(py - cy);
	interval_number pz_cz(pz - cz);

	interval_number d11((svs.id * px_cx) + (svs.ia11));
	interval_number d21((ax - cx));
	interval_number d31((bx - cx));
	interval_number d12((svs.id * py_cy) + (svs.ia12));
	interval_number d22((ay - cy));
	interval_number d32((by - cy));
	interval_number d13((svs.id * pz_cz) + (svs.ia13));
	interval_number d23((az - cz));
	interval_number d33((bz - cz));

	interval_number d2233(d22 * d33);
	interval_number d2332(d23 * d32);
	interval_number d2133(d21 * d33);
	interval_number d2331(d23 * d31);
	interval_number d2132(d21 * d32);
	interval_number d2231(d22 * d31);

	interval_number det(d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231));
	setFPUModeToRoundNEAR();

	if (!det.signIsReliable())
		return Filtered_Orientation::UNCERTAIN;

	return det.sign() * svs.id.sign();
}

bool orient3D_TPI_pre_dfilter(
	double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z,
	double ow1x, double ow1y, double ow1z, double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z,
	double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z, double ou3x, double ou3y, double ou3z,
	TPI_filtered_suppvars &svs)
{
	setFPUModeToRoundUP();

	interval_number v3x(ov3x - ov2x);
	interval_number v3y(ov3y - ov2y);
	interval_number v3z(ov3z - ov2z);
	interval_number v2x(ov2x - ov1x);
	interval_number v2y(ov2y - ov1y);
	interval_number v2z(ov2z - ov1z);
	interval_number w3x(ow3x - ow2x);
	interval_number w3y(ow3y - ow2y);
	interval_number w3z(ow3z - ow2z);
	interval_number w2x(ow2x - ow1x);
	interval_number w2y(ow2y - ow1y);
	interval_number w2z(ow2z - ow1z);
	interval_number u3x(ou3x - ou2x);
	interval_number u3y(ou3y - ou2y);
	interval_number u3z(ou3z - ou2z);
	interval_number u2x(ou2x - ou1x);
	interval_number u2y(ou2y - ou1y);
	interval_number u2z(ou2z - ou1z);

	interval_number nvx(v2y * v3z - v2z * v3y);
	interval_number nvy(v3x * v2z - v3z * v2x);
	interval_number nvz(v2x * v3y - v2y * v3x);

	interval_number nwx(w2y * w3z - w2z * w3y);
	interval_number nwy(w3x * w2z - w3z * w2x);
	interval_number nwz(w2x * w3y - w2y * w3x);

	interval_number nux(u2y * u3z - u2z * u3y);
	interval_number nuy(u3x * u2z - u3z * u2x);
	interval_number nuz(u2x * u3y - u2y * u3x);

	interval_number nwyuz(nwy * nuz - nwz * nuy);
	interval_number nwxuz(nwx * nuz - nwz * nux);
	interval_number nwxuy(nwx * nuy - nwy * nux);

	svs.id = (nvx * nwyuz - nvy * nwxuz + nvz * nwxuy);

	if (!svs.id.signIsReliable())
	{
		setFPUModeToRoundNEAR();
		return false;
	}

	interval_number nvyuz(nvy * nuz - nvz * nuy);
	interval_number nvxuz(nvx * nuz - nvz * nux);
	interval_number nvxuy(nvx * nuy - nvy * nux);

	interval_number nvywz(nvy * nwz - nvz * nwy);
	interval_number nvxwz(nvx * nwz - nvz * nwx);
	interval_number nvxwy(nvx * nwy - nvy * nwx);

	interval_number p1(nvx * ov1x + nvy * ov1y + nvz * ov1z);
	interval_number p2(nwx * ow1x + nwy * ow1y + nwz * ow1z);
	interval_number p3(nux * ou1x + nuy * ou1y + nuz * ou1z);

	svs.in1 = (p1 * nwyuz - p2 * nvyuz + p3 * nvywz);
	svs.in2 = (p2 * nvxuz - p3 * nvxwz - p1 * nwxuz);
	svs.in3 = (p3 * nvxwy - p2 * nvxuy + p1 * nwxuy);

	setFPUModeToRoundNEAR();
	if (!isfinite(svs.in1) || !isfinite(svs.in2) || !isfinite(svs.in3))
		return false;

	svs.d = NAN; // This makes it possible to know which filter succeeded

	return true;
}

int orient3D_TPI_post_dfilter(
	const TPI_filtered_suppvars &svs,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z)
{
	setFPUModeToRoundUP();

	interval_number dq3x(svs.id * q3x);
	interval_number dq3y(svs.id * q3y);
	interval_number dq3z(svs.id * q3z);

	interval_number a11(svs.in1 - dq3x);
	interval_number a12(svs.in2 - dq3y);
	interval_number a13(svs.in3 - dq3z);
	interval_number a21(q1x - q3x);
	interval_number a22(q1y - q3y);
	interval_number a23(q1z - q3z);
	interval_number a31(q2x - q3x);
	interval_number a32(q2y - q3y);
	interval_number a33(q2z - q3z);

	interval_number det(a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31));
	setFPUModeToRoundNEAR();

	if (!det.signIsReliable())
		return Filtered_Orientation::UNCERTAIN;

	return det.sign() * svs.id.sign();
}

int orient3D_LPI_dfiltered(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz)
{
	LPI_filtered_suppvars svs;

	if (!orient3D_LPI_pre_dfilter(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, svs))
		return Filtered_Orientation::UNCERTAIN;
	return orient3D_LPI_post_dfilter(svs, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
}

int orient3D_TPI_dfiltered(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z)
{
	TPI_filtered_suppvars svs;

	if (!orient3D_TPI_pre_dfilter(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z, u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, svs))
		return Filtered_Orientation::UNCERTAIN;
	return orient3D_TPI_post_dfilter(svs, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}
#endif

/////////////////// Multi-stage filters (almost-static + dynamic) ///////////////
// Almost-static only if USE_MULTISTAGE_FILTERS is not defined

bool orient3D_LPI_prefilter(
	const double &px, const double &py, const double &pz,
	const double &qx, const double &qy, const double &qz,
	const double &rx, const double &ry, const double &rz,
	const double &sx, const double &sy, const double &sz,
	const double &tx, const double &ty, const double &tz,
	LPI_filtered_suppvars &svs)
{
	double a21, a22, a23, a31, a32, a33;

	double &a11 = svs.a11;
	double &a12 = svs.a12;
	double &a13 = svs.a13;
	double &d = svs.d;
	double &fa11 = svs.fa11;
	double &fa12 = svs.fa12;
	double &fa13 = svs.fa13;
	double &max1 = svs.max1;
	double &max2 = svs.max2;
	double &max5 = svs.max5;

	a11 = (px - qx);
	a12 = (py - qy);
	a13 = (pz - qz);
	a21 = (sx - rx);
	a22 = (sy - ry);
	a23 = (sz - rz);
	a31 = (tx - rx);
	a32 = (ty - ry);
	a33 = (tz - rz);
	double a2233 = ((a22 * a33) - (a23 * a32));
	double a2133 = ((a21 * a33) - (a23 * a31));
	double a2132 = ((a21 * a32) - (a22 * a31));
	d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));

	fa11 = fabs(a11);
	double fa21 = fabs(a21);
	double fa31 = fabs(a31);
	fa12 = fabs(a12);
	double fa22 = fabs(a22);
	double fa32 = fabs(a32);
	fa13 = fabs(a13);
	double fa23 = fabs(a23);
	double fa33 = fabs(a33);

	max1 = fa23;
	if (max1 < fa13)
		max1 = fa13;
	if (max1 < fa33)
		max1 = fa33;

	max2 = fa12;
	if (max2 < fa22)
		max2 = fa22;
	if (max2 < fa32)
		max2 = fa32;

	max5 = fa11;
	if (max5 < fa21)
		max5 = fa21;
	if (max5 < fa31)
		max5 = fa31;

	double deps = 5.1107127829973299e-015 * max1 * max2 * max5;
	if (!isfinite(d) || (d <= deps && d >= -deps))
#ifdef USE_MULTISTAGE_FILTERS
		return orient3D_LPI_pre_dfilter(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, svs);
#else
		return false;
#endif

	double px_rx = px - rx;
	double py_ry = py - ry;
	double pz_rz = pz - rz;

	double n = ((((py_ry)*a2133) - ((px_rx)*a2233)) - ((pz_rz)*a2132));

	a11 *= n;
	a12 *= n;
	a13 *= n;

	if (!isfinite(a11) || !isfinite(a12) || !isfinite(a13))
#ifdef USE_MULTISTAGE_FILTERS
		return orient3D_LPI_pre_dfilter(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, svs);
#else
		return false;
#endif

	double fpxrx = fabs(px_rx);
	double fpyry = fabs(py_ry);
	double fpzrz = fabs(pz_rz);

	if (max1 < fpzrz)
		max1 = fpzrz;
	if (max2 < fpyry)
		max2 = fpyry;
	if (max5 < fpxrx)
		max5 = fpxrx;

	//#ifdef USE_MULTISTAGE_FILTERS
	//	svs.qx = qx;
	//	svs.qy = qy;
	//	svs.qz = qz;
	//	svs.rx = rx;
	//	svs.ry = ry;
	//	svs.rz = rz;
	//	svs.sx = sx;
	//	svs.sy = sy;
	//	svs.sz = sz;
	//	svs.tx = tx;
	//	svs.ty = ty;
	//	svs.tz = tz;
	//#endif

	return true;
}

int orient3D_LPI_postfilter(const LPI_filtered_suppvars &svs,
							const double &px, const double &py, const double &pz,
							const double &ax, const double &ay, const double &az,
							const double &bx, const double &by, const double &bz,
							const double &cx, const double &cy, const double &cz)
{
#ifdef USE_MULTISTAGE_FILTERS
	if (svs.d == NAN)
		return orient3D_LPI_post_dfilter(svs, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
#endif
	double px_cx, py_cy, pz_cz;
	double d11, d12, d13, d21, d31, d22, d32, d23, d33;
	double d2233, d2332, d2133, d2331, d2132, d2231;
	double det;

	const double &a11 = svs.a11, &a12 = svs.a12, &a13 = svs.a13;
	const double &d = svs.d;
	const double &fa11 = svs.fa11, &fa12 = svs.fa12, &fa13 = svs.fa13;
	const double &max1 = svs.max1, &max2 = svs.max2, &max5 = svs.max5;

	px_cx = px - cx;
	py_cy = py - cy;
	pz_cz = pz - cz;

	d11 = (d * px_cx) + (a11);
	d21 = (ax - cx);
	d31 = (bx - cx);
	d12 = (d * py_cy) + (a12);
	d22 = (ay - cy);
	d32 = (by - cy);
	d13 = (d * pz_cz) + (a13);
	d23 = (az - cz);
	d33 = (bz - cz);

	d2233 = d22 * d33;
	d2332 = d23 * d32;
	d2133 = d21 * d33;
	d2331 = d23 * d31;
	d2132 = d21 * d32;
	d2231 = d22 * d31;

	det = d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231);

	if (!isfinite(det))
		//#ifdef USE_MULTISTAGE_FILTERS
		//		return orient3D_LPI_dfiltered(px, py, pz, svs.qx, svs.qy, svs.qz, svs.rx, svs.ry, svs.rz,
		//		        svs.sx, svs.sy, svs.sz, svs.tx, svs.ty, svs.tz,	ax, ay, az,	bx, by, bz,	cx, cy, cz);
		//#else
		return Filtered_Orientation::UNCERTAIN; // Fast reject in case of under/overflow
												//#endif

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

	double max3, max4, max6;

	max3 = fa12;
	if (max3 < fd32)
		max3 = fd32;
	if (max3 < fd22)
		max3 = fd22;
	if (max3 < fpycy)
		max3 = fpycy;

	max4 = fa13;
	if (max4 < fpzcz)
		max4 = fpzcz;
	if (max4 < fd33)
		max4 = fd33;
	if (max4 < fd23)
		max4 = fd23;

	max6 = fa11;
	if (max6 < fd21)
		max6 = fd21;
	if (max6 < fd31)
		max6 = fd31;
	if (max6 < fpxcx)
		max6 = fpxcx;

	double eps = 1.3865993466947057e-013 * max1 * max2 * max5 * max6 * max3 * max4;

	if ((det > eps))
		return (d > 0) ? (Filtered_Orientation::POSITIVE) : (Filtered_Orientation::NEGATIVE);
	if ((det < -eps))
		return (d > 0) ? (Filtered_Orientation::NEGATIVE) : (Filtered_Orientation::POSITIVE);

	//#ifdef USE_MULTISTAGE_FILTERS
	//	return orient3D_LPI_dfiltered(px, py, pz, svs.qx, svs.qy, svs.qz, svs.rx, svs.ry, svs.rz,
	//		svs.sx, svs.sy, svs.sz, svs.tx, svs.ty, svs.tz, ax, ay, az, bx, by, bz, cx, cy, cz);
	//#else
	return Filtered_Orientation::UNCERTAIN; // Fast reject in case of under/overflow
	//#endif
}

bool orient3D_TPI_prefilter(
	const double &ov1x, const double &ov1y, const double &ov1z, const double &ov2x, const double &ov2y, const double &ov2z, const double &ov3x, const double &ov3y, const double &ov3z,
	const double &ow1x, const double &ow1y, const double &ow1z, const double &ow2x, const double &ow2y, const double &ow2z, const double &ow3x, const double &ow3y, const double &ow3z,
	const double &ou1x, const double &ou1y, const double &ou1z, const double &ou2x, const double &ou2y, const double &ou2z, const double &ou3x, const double &ou3y, const double &ou3z,
	TPI_filtered_suppvars &svs)
{
	double &d = svs.d, &n1 = svs.n1, &n2 = svs.n2, &n3 = svs.n3;
	double &max1 = svs.max1, &max2 = svs.max2, &max3 = svs.max3, &max4 = svs.max4, &max5 = svs.max5, &max6 = svs.max6, &max7 = svs.max7;

	double v3x = ov3x - ov2x;
	double v3y = ov3y - ov2y;
	double v3z = ov3z - ov2z;
	double v2x = ov2x - ov1x;
	double v2y = ov2y - ov1y;
	double v2z = ov2z - ov1z;
	double w3x = ow3x - ow2x;
	double w3y = ow3y - ow2y;
	double w3z = ow3z - ow2z;
	double w2x = ow2x - ow1x;
	double w2y = ow2y - ow1y;
	double w2z = ow2z - ow1z;
	double u3x = ou3x - ou2x;
	double u3y = ou3y - ou2y;
	double u3z = ou3z - ou2z;
	double u2x = ou2x - ou1x;
	double u2y = ou2y - ou1y;
	double u2z = ou2z - ou1z;

	double nvx = v2y * v3z - v2z * v3y;
	double nvy = v3x * v2z - v3z * v2x;
	double nvz = v2x * v3y - v2y * v3x;

	double nwx = w2y * w3z - w2z * w3y;
	double nwy = w3x * w2z - w3z * w2x;
	double nwz = w2x * w3y - w2y * w3x;

	double nux = u2y * u3z - u2z * u3y;
	double nuy = u3x * u2z - u3z * u2x;
	double nuz = u2x * u3y - u2y * u3x;

	double nwyuz = nwy * nuz - nwz * nuy;
	double nwxuz = nwx * nuz - nwz * nux;
	double nwxuy = nwx * nuy - nwy * nux;

	d = nvx * nwyuz - nvy * nwxuz + nvz * nwxuy;

	// Almost static filter for d
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

	max4 = fv2y;
	if (max4 < fv3y)
		max4 = fv3y;
	if (max4 < fw3y)
		max4 = fw3y;
	if (max4 < fw2y)
		max4 = fw2y;
	max2 = fv3x;
	if (max2 < fv2x)
		max2 = fv2x;
	if (max2 < fw2x)
		max2 = fw2x;
	if (max2 < fw3x)
		max2 = fw3x;
	max5 = fv2z;
	if (max5 < fv3z)
		max5 = fv3z;
	if (max5 < fw3z)
		max5 = fw3z;
	if (max5 < fw2z)
		max5 = fw2z;
	max7 = fu2x;
	if (max7 < fu3x)
		max7 = fu3x;
	if (max7 < fw2x)
		max7 = fw2x;
	if (max7 < fw3x)
		max7 = fw3x;

	double max9 = fu2y;
	if (max9 < fu3y)
		max9 = fu3y;
	if (max9 < fw2y)
		max9 = fw2y;
	if (max9 < fw3y)
		max9 = fw3y;
	double max10 = fu2z;
	if (max10 < fu3z)
		max10 = fu3z;
	if (max10 < fw2z)
		max10 = fw2z;
	if (max10 < fw3z)
		max10 = fw3z;

	double deps = 8.8881169117764924e-014 * (((((max4 * max5) * max2) * max10) * max7) * max9);
	if (!isfinite(d) || (d <= deps && d >= -deps))
#ifdef USE_MULTISTAGE_FILTERS
		return orient3D_TPI_pre_dfilter(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z,
										ow1x, ow1y, ow1z, ow2x, ow2y, ow2z, ow3x, ow3y, ow3z, ou1x, ou1y, ou1z, ou2x, ou2y, ou2z,
										ou3x, ou3y, ou3z, svs);
#else
		return false;
#endif

	double nvyuz = nvy * nuz - nvz * nuy;
	double nvxuz = nvx * nuz - nvz * nux;
	double nvxuy = nvx * nuy - nvy * nux;

	double nvywz = nvy * nwz - nvz * nwy;
	double nvxwz = nvx * nwz - nvz * nwx;
	double nvxwy = nvx * nwy - nvy * nwx;

	double p1 = nvx * ov1x + nvy * ov1y + nvz * ov1z;
	double p2 = nwx * ow1x + nwy * ow1y + nwz * ow1z;
	double p3 = nux * ou1x + nuy * ou1y + nuz * ou1z;

	n1 = p1 * nwyuz - p2 * nvyuz + p3 * nvywz;
	n2 = p2 * nvxuz - p3 * nvxwz - p1 * nwxuz;
	n3 = p3 * nvxwy - p2 * nvxuy + p1 * nwxuy;

	if (!isfinite(n1) || !isfinite(n2) || !isfinite(n3))
#ifdef USE_MULTISTAGE_FILTERS
		return orient3D_TPI_pre_dfilter(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z,
										ow1x, ow1y, ow1z, ow2x, ow2y, ow2z, ow3x, ow3y, ow3z, ou1x, ou1y, ou1z, ou2x, ou2y, ou2z,
										ou3x, ou3y, ou3z, svs);
#else
		return false;
#endif

	if (max4 < fu2y)
		max4 = fu2y;
	if (max4 < fu3y)
		max4 = fu3y;
	if (max2 < fu2x)
		max2 = fu2x;
	if (max2 < fu3x)
		max2 = fu3x;
	if (max5 < fu2z)
		max5 = fu2z;
	if (max5 < fu3z)
		max5 = fu3z;

	max1 = max4;
	if (max1 < max2)
		max1 = max2;
	max3 = max5;
	if (max3 < max4)
		max3 = max4;
	max6 = fu2x;
	if (max6 < fu3x)
		max6 = fu3x;
	if (max6 < fu2z)
		max6 = fu2z;
	if (max6 < fw3y)
		max6 = fw3y;
	if (max6 < fw2x)
		max6 = fw2x;
	if (max6 < fw3z)
		max6 = fw3z;
	if (max6 < fw2y)
		max6 = fw2y;
	if (max6 < fw2z)
		max6 = fw2z;
	if (max6 < fu2y)
		max6 = fu2y;
	if (max6 < fu3z)
		max6 = fu3z;
	if (max6 < fu3y)
		max6 = fu3y;
	if (max6 < fw3x)
		max6 = fw3x;

	//#ifdef USE_MULTISTAGE_FILTERS
	//	svs.v1x = ov1x; svs.v1y = ov1y; svs.v1z = ov1z; svs.v2x = ov2x; svs.v2y = ov2y; svs.v2z = ov2z; svs.v3x = ov3x; svs.v3y = ov3y; svs.v3z = ov3z;
	//	svs.w1x = ow1x; svs.w1y = ow1y; svs.w1z = ow1z; svs.w2x = ow2x; svs.w2y = ow2y; svs.w2z = ow2z; svs.w3x = ow3x; svs.w3y = ow3y; svs.w3z = ow3z;
	//	svs.u1x = ou1x; svs.u1y = ou1y; svs.u1z = ou1z; svs.u2x = ou2x; svs.u2y = ou2y; svs.u2z = ou2z; svs.u3x = ou3x; svs.u3y = ou3y; svs.u3z = ou3z;
	//#endif
	return true;
}

int orient3D_TPI_postfilter(
	const TPI_filtered_suppvars &svs,
	const double &q1x, const double &q1y, const double &q1z, const double &q2x, const double &q2y, const double &q2z, const double &q3x, const double &q3y, const double &q3z)
{
#ifdef USE_MULTISTAGE_FILTERS
	if (svs.d == NAN)
		return orient3D_TPI_post_dfilter(svs, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
#endif

	const double &d = svs.d, &n1 = svs.n1, &n2 = svs.n2, &n3 = svs.n3;
	const double &max1 = svs.max1, &max2 = svs.max2, &max3 = svs.max3, &max4 = svs.max4, &max5 = svs.max5, &max6 = svs.max6, &max7 = svs.max7;

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

	double det = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
	bool infinite = true;
	infinite = isfinite(det);
	if (!infinite)
		//#ifdef USE_MULTISTAGE_FILTERS
		//	return orient3D_TPI_dfiltered(svs.v1x, svs.v1y, svs.v1z, svs.v2x, svs.v2y, svs.v2z, svs.v3x, svs.v3y, svs.v3z,
		//		svs.w1x, svs.w1y, svs.w1z, svs.w2x, svs.w2y, svs.w2z, svs.w3x, svs.w3y, svs.w3z, svs.u1x, svs.u1y, svs.u1z, svs.u2x, svs.u2y, svs.u2z,
		//		svs.u3x, svs.u3y, svs.u3z, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
		//#else
		return Filtered_Orientation::UNCERTAIN;
	//#endif

	double fa21 = fabs(a21);
	double fa22 = fabs(a22);
	double fa23 = fabs(a23);
	double fa31 = fabs(a31);
	double fa32 = fabs(a32);
	double fa33 = fabs(a33);

	double nmax7 = max7;
	if (nmax7 < fa21)
		nmax7 = fa21;
	if (nmax7 < fa31)
		nmax7 = fa31;

	double nmax6 = max6;
	if (nmax6 < fa22)
		nmax6 = fa22;
	if (nmax6 < fa32)
		nmax6 = fa32;

	double max8 = fa22;
	if (max8 < fa23)
		max8 = fa23;
	if (max8 < fa33)
		max8 = fa33;
	if (max8 < fa32)
		max8 = fa32;

	double eps = 3.4025182954957945e-012 * (((((((max1 * max3) * max2) * max5) * nmax7) * max4) * nmax6) * max8);

	if ((det > eps))
		return (d > 0) ? (Filtered_Orientation::POSITIVE) : (Filtered_Orientation::NEGATIVE);
	if ((det < -eps))
		return (d > 0) ? (Filtered_Orientation::NEGATIVE) : (Filtered_Orientation::POSITIVE);
	//#ifdef USE_MULTISTAGE_FILTERS
	//	return orient3D_TPI_dfiltered(svs.v1x, svs.v1y, svs.v1z, svs.v2x, svs.v2y, svs.v2z, svs.v3x, svs.v3y, svs.v3z,
	//		svs.w1x, svs.w1y, svs.w1z, svs.w2x, svs.w2y, svs.w2z, svs.w3x, svs.w3y, svs.w3z, svs.u1x, svs.u1y, svs.u1z, svs.u2x, svs.u2y, svs.u2z,
	//		svs.u3x, svs.u3y, svs.u3z, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
	//#else
	return Filtered_Orientation::UNCERTAIN;
	//#endif
}

int orient3D_LPI_filtered(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz)
{
	LPI_filtered_suppvars svs;

	if (!orient3D_LPI_prefilter(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, svs))
		return Filtered_Orientation::UNCERTAIN;
	return orient3D_LPI_postfilter(svs, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
}

int orient3D_TPI_filtered(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z)
{
	TPI_filtered_suppvars svs;

	if (!orient3D_TPI_prefilter(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z, u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, svs))
		return Filtered_Orientation::UNCERTAIN;
	return orient3D_TPI_postfilter(svs, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

// The following macros are fast implementations of basic expansion arithmetic due
// to Dekker, Knuth, Priest, Shewchuk, and others.

// See Y. Hida, X. S. Li,  D. H. Bailey "Algorithms for Quad-Double Precision Floating Point Arithmetic"

// Sums
#define Quick_Two_Sum(a, b, x, y) \
	x = a + b;                    \
	y = b - (x - a)
#define Two_Sum(a, b, x, y) \
	x = a + b;              \
	_bv = x - a;            \
	y = (a - (x - _bv)) + (b - _bv)
#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
	Two_Sum(a0, b, _i, x0);                \
	Two_Sum(a1, _i, x2, x1)

// Differences
#define Two_Diff(a, b, x, y) \
	x = a - b;               \
	_bv = a - x;             \
	y = (a - (x + _bv)) + (_bv - b)
#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
	Two_Diff(a0, b, _i, x0);                \
	Two_Sum(a1, _i, x2, x1)

// Products
#define Split(a, _ah, _al)            \
	_c = 1.3421772800000003e+008 * a; \
	_ah = _c - (_c - a);              \
	_al = a - _ah
#define Two_Prod_PreSplit(a, b, _bh, _bl, x, y) \
	x = a * b;                                  \
	Split(a, _ah, _al);                         \
	y = (_al * _bl) - (((x - (_ah * _bh)) - (_al * _bh)) - (_ah * _bl))
#define Two_Product_2Presplit(a, _ah, _al, b, _bh, _bl, x, y) \
	x = a * b;                                                \
	y = (_al * _bl) - (((x - _ah * _bh) - (_al * _bh)) - (_ah * _bl))

//////////////////////////////////////////////////////////////////////////////////
//
//   O R I E N T 3 D _ L P I
//
//////////////////////////////////////////////////////////////////////////////////

bool orient3D_LPI_pre_exact(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	double *a11, double *a12, double *a13,
	double *d, double *n,
	int &dl, int &nl)
{
	double a21[2], a22[2], a23[2], a31[2], a32[2], a33[2];
	double px_rx[2], py_ry[2], pz_rz[2];
	double t1[8], t2[8];
	double a2233[16], a2133[16], a2132[16];
	double tt1[64], tt2[64], tt3[64];
	double ttt1[128];
	int a2233l, a2133l, a2132l, tt1l, tt2l, tt3l, ttt1l;

	expansionObject o;

	o.two_Diff(px, qx, a11);   // a11 = px - qx;
	o.two_Diff(py, qy, a12);   // a12 = py - qy;
	o.two_Diff(pz, qz, a13);   // a13 = pz - qz;
	o.two_Diff(sx, rx, a21);   // a21 = sx - rx;
	o.two_Diff(sy, ry, a22);   // a22 = sy - ry;
	o.two_Diff(sz, rz, a23);   // a23 = sz - rz;
	o.two_Diff(tx, rx, a31);   // a31 = tx - rx;
	o.two_Diff(ty, ry, a32);   // a32 = ty - ry;
	o.two_Diff(tz, rz, a33);   // a33 = tz - rz;
	o.two_Diff(px, rx, px_rx); // px_rx = px - rx;
	o.two_Diff(py, ry, py_ry); // py_ry = py - ry;
	o.two_Diff(pz, rz, pz_rz); // pz_rz = pz - rz;

	o.Two_Two_Prod(a22[1], a22[0], a33[1], a33[0], t1); // t1 = a22 * a33;
	o.Two_Two_Prod(a23[1], a23[0], a32[1], a32[0], t2); // t2 = a23 * a32;
	o.Gen_Invert(8, t2);								// t2 = -t2;
	a2233l = o.Gen_Sum(8, t1, 8, t2, a2233);			// a2233 = t1 + t2; // = a22*a33 - a23*a32;

	o.Two_Two_Prod(a21[1], a21[0], a33[1], a33[0], t1); // t1 = a21 * a33;
	o.Two_Two_Prod(a23[1], a23[0], a31[1], a31[0], t2); // t2 = a23 * a31;
	o.Gen_Invert(8, t2);								// t2 = -t2;
	a2133l = o.Gen_Sum(8, t1, 8, t2, a2133);			// a2133 = t1 + t2; // = a21*a33 - a23*a31;

	o.Two_Two_Prod(a21[1], a21[0], a32[1], a32[0], t1); // t1 = a21 * a32;
	o.Two_Two_Prod(a22[1], a22[0], a31[1], a31[0], t2); // t2 = a22 * a31;
	o.Gen_Invert(8, t2);								// t2 = -t2;
	a2132l = o.Gen_Sum(8, t1, 8, t2, a2132);			// a2132 = t1 + t2; // = a21*a32 - a22*a31;

	tt1l = o.Gen_Product(a2233l, a2233, 2, a11, tt1); // tt1 = a2233 * a11;
	tt2l = o.Gen_Product(a2133l, a2133, 2, a12, tt2); // tt2 = a2133 * a12;
	tt3l = o.Gen_Product(a2132l, a2132, 2, a13, tt3); // tt3 = a2132 * a13;
	o.Gen_Invert(tt2l, tt2);						  // tt2 = -tt2;
	ttt1l = o.Gen_Sum(tt1l, tt1, tt2l, tt2, ttt1);	  // ttt1 = tt1 + tt2;
	dl = o.Gen_Sum(ttt1l, ttt1, tt3l, tt3, d);		  // d = ttt1 + tt3; // = tt1 + tt2 + tt3; // = a2233 * a11 - a2133 * a12 + a2132 * a13;

	if (dl == 0)
		return false;

	tt1l = o.Gen_Product(a2133l, a2133, 2, py_ry, tt1); // tt1 = a2133 * py_ry;
	tt2l = o.Gen_Product(a2233l, a2233, 2, px_rx, tt2); // tt2 = a2233 * px_rx;
	tt3l = o.Gen_Product(a2132l, a2132, 2, pz_rz, tt3); // tt3 = a2132 * pz_rz;
	o.Gen_Invert(tt2l, tt2);							// tt2 = -tt2;
	ttt1l = o.Gen_Sum(tt1l, tt1, tt2l, tt2, ttt1);		// ttt1 = tt1 + tt2;
	o.Gen_Invert(tt3l, tt3);							// tt3 = -tt3;
	nl = o.Gen_Sum(ttt1l, ttt1, tt3l, tt3, n);			// n = ttt1 + tt3; // = tt1 + tt2 + tt3; // = a2133 * py_ry - a2233 * px_rx - a2132 * pz_rz;

	if (nl == 0)
		return false;

	return true;
}

bool orient3D_LPI_pre_exact(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	LPI_exact_suppvars &s)
{
	return orient3D_LPI_pre_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, s.a11, s.a12, s.a13, s.d, s.n, s.dl, s.nl);
}

int orient3D_LPI_post_exact(
	double *a11, double *a12, double *a13,
	double *d, double *n,
	int dl, int nl,
	double px, double py, double pz,
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz)
{
	double px_cx[2], py_cy[2], pz_cz[2];
	double d21[2], d22[2], d23[2], d31[2], d32[2], d33[2];
	double a2233[16], a2133[16], a2132[16];

	double ii1[768], ii2[768];
	double d11[1536], d12[1536], d13[1536];
	double d2233[8], d2332[8], d2133[8], d2331[8], d2132[8], d2231[8];

	double s1p[256], s2p[256], s3p[256];
	double *s1 = s1p, *s2 = s2p, *s3 = s3p; // Bound is 49152 each
	int s1l = 256, s2l = 256, s3l = 256;
	double ss1p[256];
	double *ss1 = ss1p; // Bound is 98304
	int ss1l = 256;
	double detp[256];
	double *det = detp; // Bound is 147456
	int detl = 256;
	int a2233l, a2133l, a2132l, ii1l, ii2l, d11l, d12l, d13l;

	expansionObject o;

	o.two_Diff(px, cx, px_cx); // px_cx = px - cx;
	o.two_Diff(py, cy, py_cy); // py_cy = py - cy;
	o.two_Diff(pz, cz, pz_cz); // pz_cz = pz - cz;

	ii1l = o.Gen_Product(dl, d, 2, px_cx, ii1);
	ii2l = o.Gen_Product(nl, n, 2, a11, ii2);
	d11l = o.Gen_Sum(ii1l, ii1, ii2l, ii2, d11); // d11 = (d * px_cx) + (a11 * n);

	ii1l = o.Gen_Product(dl, d, 2, py_cy, ii1);
	ii2l = o.Gen_Product(nl, n, 2, a12, ii2);
	d12l = o.Gen_Sum(ii1l, ii1, ii2l, ii2, d12); // d12 = (d * py_cy) + (a12 * n);

	ii1l = o.Gen_Product(dl, d, 2, pz_cz, ii1);
	ii2l = o.Gen_Product(nl, n, 2, a13, ii2);
	d13l = o.Gen_Sum(ii1l, ii1, ii2l, ii2, d13); // d13 = (d * pz_cz) + (a13 * n);

	o.two_Diff(ax, cx, d21); // d21 = (ax - cx);
	o.two_Diff(bx, cx, d31); // d31 = (bx - cx);
	o.two_Diff(ay, cy, d22); // d22 = (ay - cy);
	o.two_Diff(by, cy, d32); // d32 = (by - cy);
	o.two_Diff(az, cz, d23); // d23 = (az - cz);
	o.two_Diff(bz, cz, d33); // d33 = (bz - cz);

	o.Two_Two_Prod(d22[1], d22[0], d33[1], d33[0], d2233);	 // d2233 = d22*d33;
	o.Two_Two_Prod(-d23[1], -d23[0], d32[1], d32[0], d2332); // d2332 = -d23*d32;
	o.Two_Two_Prod(d21[1], d21[0], d33[1], d33[0], d2133);	 // d2133 = d21*d33;
	o.Two_Two_Prod(-d23[1], -d23[0], d31[1], d31[0], d2331); // d2331 = -d23*d31;
	o.Two_Two_Prod(d21[1], d21[0], d32[1], d32[0], d2132);	 // d2132 = d21*d32;
	o.Two_Two_Prod(-d22[1], -d22[0], d31[1], d31[0], d2231); // d2231 = -d22*d31;

	a2233l = o.Gen_Sum(8, d2233, 8, d2332, a2233); // a2233 = d2233 + d2332;
	a2133l = o.Gen_Sum(8, d2133, 8, d2331, a2133); // a2133 = d2133 + d2331;
	a2132l = o.Gen_Sum(8, d2132, 8, d2231, a2132); // a2132 = d2132 + d2231;

	s1l = o.Gen_Product_With_PreAlloc(d11l, d11, a2233l, a2233, &s1, s1l); // s1 = d11 * a2233;
	s2l = o.Gen_Product_With_PreAlloc(d12l, d12, a2133l, a2133, &s2, s2l); // s2 = d12 * a2133;
	s3l = o.Gen_Product_With_PreAlloc(d13l, d13, a2132l, a2132, &s3, s3l); // s3 = d13 * a2132;

	o.Gen_Invert(s2l, s2);
	ss1l = o.Gen_Sum_With_PreAlloc(s1l, s1, s2l, s2, &ss1, ss1l);	// ss1 = s1 - s2;
	detl = o.Gen_Sum_With_PreAlloc(ss1l, ss1, s3l, s3, &det, detl); // det = ss1 + s3; // = s1 - s2 + s3;

	if (s1 != s1p)
		free(s1);
	if (s2 != s2p)
		free(s2);
	if (s3 != s3p)
		free(s3);
	if (ss1 != ss1p)
		free(ss1);

	double s = det[detl - 1];
	double sd = d[dl - 1];

	if (det != detp)
		free(det);

	if ((s > 0))
		return (sd > 0) ? (1) : (-1);
	if ((s < 0))
		return (sd > 0) ? (-1) : (1);
	return 0;
}

int orient3D_LPI_post_exact(
	LPI_exact_suppvars &s,
	double px, double py, double pz,
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz)
{
	return orient3D_LPI_post_exact(s.a11, s.a12, s.a13, s.d, s.n, s.dl, s.nl, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
}

int orient3D_LPI_exact(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz)
{
	LPI_exact_suppvars s;

	if (orient3D_LPI_pre_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, s))
		return orient3D_LPI_post_exact(s, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////
//
//   O R I E N T 3 D _ T P I
//
//////////////////////////////////////////////////////////////////////////////////

inline void o3dTPI_tf1(expansionObject &o, double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
					   double *nvx, double *nvy, double *nvz, int &nvxl, int &nvyl, int &nvzl)
{
	double v32x[2], v32y[2], v32z[2], v21x[2], v21y[2], v21z[2]; // 2
	double tp1[8], tp2[8];										 // 8

	o.two_Diff(v3x, v2x, v32x); // v32x = v3x - v2x;
	o.two_Diff(v3y, v2y, v32y); // v32y = v3y - v2y;
	o.two_Diff(v3z, v2z, v32z); // v32z = v3z - v2z;
	o.two_Diff(v2x, v1x, v21x); // v21x = v2x - v1x;
	o.two_Diff(v2y, v1y, v21y); // v21y = v2y - v1y;
	o.two_Diff(v2z, v1z, v21z); // v21z = v2z - v1z;

	o.Two_Two_Prod(v21y[1], v21y[0], v32z[1], v32z[0], tp1); // tp1 = v21y*v32z;
	o.Two_Two_Prod(v21z[1], v21z[0], v32y[1], v32y[0], tp2); // tp2 = v21z*v32y;
	o.Gen_Invert(8, tp2);
	nvxl = o.Gen_Sum(8, tp1, 8, tp2, nvx);

	o.Two_Two_Prod(v32x[1], v32x[0], v21z[1], v21z[0], tp1); // tp1 = v32x*v21z;
	o.Two_Two_Prod(v32z[1], v32z[0], v21x[1], v21x[0], tp2); // tp2 = v32z*v21x;
	o.Gen_Invert(8, tp2);
	nvyl = o.Gen_Sum(8, tp1, 8, tp2, nvy);

	o.Two_Two_Prod(v21x[1], v21x[0], v32y[1], v32y[0], tp1); // tp1 = v21x*v32y;
	o.Two_Two_Prod(v21y[1], v21y[0], v32x[1], v32x[0], tp2); // tp2 = v21y*v32x;
	o.Gen_Invert(8, tp2);
	nvzl = o.Gen_Sum(8, tp1, 8, tp2, nvz);
}

inline void o3dTPI_tf2(expansionObject &o, double *nwx, double *nwy, double *nwz, int nwxl, int nwyl, int nwzl, double *nux, double *nuy, double *nuz, int nuxl, int nuyl, int nuzl,
					   double *nwyuz, double *nwxuz, double *nwxuy, int &nwyuzl, int &nwxuzl, int &nwxuyl)
{
	double tq1[64], tq2[64]; // 64
	int tq1l, tq2l;

	tq1l = o.Gen_Product(nwyl, nwy, nuzl, nuz, tq1); // tq1 = nwy*nuz;
	tq2l = o.Gen_Product(nwzl, nwz, nuyl, nuy, tq2); // tq2 = nwz*nuy;
	o.Gen_Invert(tq2l, tq2);
	nwyuzl = o.Gen_Sum(tq1l, tq1, tq2l, tq2, nwyuz); // nwyuz = tq1 - tq2;

	tq1l = o.Gen_Product(nwxl, nwx, nuzl, nuz, tq1); // tq1 = nwx*nuz;
	tq2l = o.Gen_Product(nwzl, nwz, nuxl, nux, tq2); // tq2 = nwz*nux;
	o.Gen_Invert(tq2l, tq2);
	nwxuzl = o.Gen_Sum(tq1l, tq1, tq2l, tq2, nwxuz); // nwxuz = tq1 - tq2;

	tq1l = o.Gen_Product(nwxl, nwx, nuyl, nuy, tq1); // tq1 = nwx*nuy;
	tq2l = o.Gen_Product(nwyl, nwy, nuxl, nux, tq2); // tq2 = nwy*nux;
	o.Gen_Invert(tq2l, tq2);
	nwxuyl = o.Gen_Sum(tq1l, tq1, tq2l, tq2, nwxuy); // nwxuy = tq1 - tq2;
}

// g = a*b + c*d - e*f
inline void o3dTPI_tf3(expansionObject &o, double *a, double *b, double *c, double *d, double *e, double *f,
					   int al, int bl, int cl, int dl, int el, int fl,
					   double **g, int &gl)
{
	int s1 = al * bl;
	int s2 = el * fl;
	int s3 = cl * dl;
	if (s2 > s1)
		s1 = s2;
	int ss = s2 * s3;
	int tr1l, tr2l, tr3l, tsl;

	double tr1p[256], tr2p[256], tr3p[256], tsp[256];
	double *tr1 = tr1p, *tr2 = tr2p, *tr3 = tr3p, *ts = tsp;

	tr1l = o.Gen_Product_With_PreAlloc(al, a, bl, b, &tr1, 256);   // tr1 = a*b;
	tr2l = o.Gen_Product_With_PreAlloc(cl, c, dl, d, &tr2, 256);   // tr2 = c*d;
	tsl = o.Gen_Sum_With_PreAlloc(tr1l, tr1, tr2l, tr2, &ts, 256); // ts = tr1 + tr2;
	tr3l = o.Gen_Product_With_PreAlloc(el, e, fl, f, &tr3, 256);   // tr3 = e*f;
	o.Gen_Invert(tr3l, tr3);
	gl = o.Gen_Sum_With_PreAlloc(tsl, ts, tr3l, tr3, g, gl); // g = ts + tr3;

	if (tr1 != tr1p)
		free(tr1);
	if (tr2 != tr2p)
		free(tr2);
	if (tr3 != tr3p)
		free(tr3);
	if (ts != tsp)
		free(ts);
}

inline void o3dTPI_tf3s(expansionObject &o, double *a, double b, double *c, double d, double *e, double f,
						int al, int cl, int el,
						double *g, int &gl)
{
	double tr1[32], tr2[32]; // 32
	double ts[64];			 // 64
	int tr1l, tr2l, tsl;

	tr1l = o.Gen_Scale(al, a, b, tr1);		   // tr1 = a*b;
	tr2l = o.Gen_Scale(cl, c, d, tr2);		   // tr2 = c*d;
	tsl = o.Gen_Sum(tr1l, tr1, tr2l, tr2, ts); // ts = tr1 + tr2;
	tr1l = o.Gen_Scale(el, e, f, tr1);		   // tr1 = e*f;
	gl = o.Gen_Sum(tsl, ts, tr1l, tr1, g);	   // g = ts + tr1;
}

bool orient3D_TPI_pre_exact(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double **d, int &dl, double **n1, int &n1l, double **n2, int &n2l, double **n3, int &n3l)
{
	double nvx[16], nvy[16], nvz[16], nwx[16], nwy[16], nwz[16], nux[16], nuy[16], nuz[16]; // 16
	int nvxl, nvyl, nvzl, nwxl, nwyl, nwzl, nuxl, nuyl, nuzl;
	double p1[96], p2[96], p3[96]; // 96
	int p1l, p2l, p3l;
	double nwyuz[128], nwxuz[128], nwxuy[128], nvyuz[128], nvxuz[128], nvxuy[128], nvywz[128], nvxwz[128], nvxwy[128]; // 128
	int nwyuzl, nwxuzl, nwxuyl, nvyuzl, nvxuzl, nvxuyl, nvywzl, nvxwzl, nvxwyl;

	expansionObject o;

	o3dTPI_tf1(o, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, nvx, nvy, nvz, nvxl, nvyl, nvzl);
	o3dTPI_tf1(o, w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z, nwx, nwy, nwz, nwxl, nwyl, nwzl);
	o3dTPI_tf1(o, u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, nux, nuy, nuz, nuxl, nuyl, nuzl);

	o3dTPI_tf2(o, nwx, nwy, nwz, nwxl, nwyl, nwzl, nux, nuy, nuz, nuxl, nuyl, nuzl, nwyuz, nwxuz, nwxuy, nwyuzl, nwxuzl, nwxuyl);
	o3dTPI_tf2(o, nvx, nvy, nvz, nvxl, nvyl, nvzl, nux, nuy, nuz, nuxl, nuyl, nuzl, nvyuz, nvxuz, nvxuy, nvyuzl, nvxuzl, nvxuyl);
	o3dTPI_tf2(o, nvx, nvy, nvz, nvxl, nvyl, nvzl, nwx, nwy, nwz, nwxl, nwyl, nwzl, nvywz, nvxwz, nvxwy, nvywzl, nvxwzl, nvxwyl);

	o3dTPI_tf3s(o, nvx, v1x, nvy, v1y, nvz, v1z, nvxl, nvyl, nvzl, p1, p1l);
	o3dTPI_tf3s(o, nwx, w1x, nwy, w1y, nwz, w1z, nwxl, nwyl, nwzl, p2, p2l);
	o3dTPI_tf3s(o, nux, u1x, nuy, u1y, nuz, u1z, nuxl, nuyl, nuzl, p3, p3l);

	o3dTPI_tf3(o, p1, nwyuz, p3, nvywz, p2, nvyuz, p1l, nwyuzl, p3l, nvywzl, p2l, nvyuzl, n1, n1l); // n1 = p1*nwyuz - p2*nvyuz + p3*nvywz;
	o3dTPI_tf3(o, p3, nvxwy, p1, nwxuy, p2, nvxuy, p3l, nvxwyl, p1l, nwxuyl, p2l, nvxuyl, n3, n3l); // n3 = p3*nvxwy - p2*nvxuy + p1*nwxuy;
	o.Gen_Invert(p3l, p3);
	o3dTPI_tf3(o, p2, nvxuz, p3, nvxwz, p1, nwxuz, p2l, nvxuzl, p3l, nvxwzl, p1l, nwxuzl, n2, n2l); // n2 = p2*nvxuz - p3*nvxwz - p1*nwxuz;

	o3dTPI_tf3(o, nvx, nwyuz, nvz, nwxuy, nvy, nwxuz, nvxl, nwyuzl, nvzl, nwxuyl, nvyl, nwxuzl, d, dl);

	return (dl != 0);
}

bool orient3D_TPI_pre_exact(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	TPI_exact_suppvars &s)
{
	return orient3D_TPI_pre_exact(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
								  u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, &s.d, s.dl, &s.n1, s.n1l, &s.n2, s.n2l, &s.n3, s.n3l);
}

int orient3D_TPI_post_exact(
	double *d, int dl, double *n1, int n1l, double *n2, int n2l, double *n3, int n3l,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z)
{
	double a21[2], a22[2], a23[2], a31[2], a32[2], a33[2];			   // 2
	double a2233[8], a2332[8], a2133[8], a2331[8], a2132[8], a2231[8]; // 8
	double dd1[16], dd2[16], dd3[16];								   // 16
	int dd1l, dd2l, dd3l;

	double dq3xp[256], dq3yp[256], dq3zp[256];
	double *dq3x = dq3xp, *dq3y = dq3yp, *dq3z = dq3zp; // 24576
	int dq3xl = 256, dq3yl = 256, dq3zl = 256;

	double a11p[256], a12p[256], a13p[256];
	double *a11 = a11p, *a12 = a12p, *a13 = a13p; // 98304
	int a11l = 256, a12l = 256, a13l = 256;
	double ee1p[256], ee2p[256];
	double *ee1 = ee1p, *ee2 = ee2p; // 3145728
	int ee1l = 256, ee2l = 256;
	double ffp[256];
	double *ff = ffp; // 6291456
	int ffl = 256;
	double detp[256];
	double *det = detp; // 9437184
	int detl = 256;

	expansionObject o;

	dq3xl = o.Gen_Scale_With_PreAlloc(dl, d, -q3x, &dq3x, dq3xl); // dq3x = -d * q3x;
	dq3yl = o.Gen_Scale_With_PreAlloc(dl, d, -q3y, &dq3y, dq3yl); // dq3y = -d * q3y;
	dq3zl = o.Gen_Scale_With_PreAlloc(dl, d, -q3z, &dq3z, dq3zl); // dq3z = -d * q3z;

	a11l = o.Gen_Sum_With_PreAlloc(n1l, n1, dq3xl, dq3x, &a11, a11l); // a11 = n1 + dq3x;
	a12l = o.Gen_Sum_With_PreAlloc(n2l, n2, dq3yl, dq3y, &a12, a12l); // a12 = n2 + dq3y;
	a13l = o.Gen_Sum_With_PreAlloc(n3l, n3, dq3zl, dq3z, &a13, a13l); // a13 = n3 + dq3z;

	if (dq3x != dq3xp)
		free(dq3x);
	if (dq3y != dq3yp)
		free(dq3y);
	if (dq3z != dq3zp)
		free(dq3z);

	o.two_Diff(q1x, q3x, a21); // a21 = q1x - q3x;
	o.two_Diff(q1y, q3y, a22); // a22 = q1y - q3y;
	o.two_Diff(q1z, q3z, a23); // a23 = q1z - q3z;
	o.two_Diff(q2x, q3x, a31); // a31 = q2x - q3x;
	o.two_Diff(q2y, q3y, a32); // a32 = q2y - q3y;
	o.two_Diff(q2z, q3z, a33); // a33 = q2z - q3z;

	o.Two_Two_Prod(a22[1], a22[0], a33[1], a33[0], a2233); // a2233 = a22*a33;
	o.Two_Two_Prod(a23[1], a23[0], a32[1], a32[0], a2332); // a2332 = a23*a32;
	o.Two_Two_Prod(a21[1], a21[0], a33[1], a33[0], a2133); // a2133 = a21*a33;
	o.Two_Two_Prod(a23[1], a23[0], a31[1], a31[0], a2331); // a2331 = a23*a31;
	o.Two_Two_Prod(a21[1], a21[0], a32[1], a32[0], a2132); // a2132 = a21*a32;
	o.Two_Two_Prod(a22[1], a22[0], a31[1], a31[0], a2231); // a2231 = a22*a31;

	o.Gen_Invert(8, a2332);
	o.Gen_Invert(8, a2331);
	o.Gen_Invert(8, a2231);
	dd1l = o.Gen_Sum(8, a2233, 8, a2332, dd1); // dd1 = a2233 + a2332;
	dd2l = o.Gen_Sum(8, a2133, 8, a2331, dd2); // dd2 = a2133 + a2331;
	dd3l = o.Gen_Sum(8, a2132, 8, a2231, dd3); // dd3 = a2132 + a2231;

	ee1l = o.Gen_Product_With_PreAlloc(a11l, a11, dd1l, dd1, &ee1, ee1l); // ee1 = a11*dd1;
	ee2l = o.Gen_Product_With_PreAlloc(a13l, a13, dd3l, dd3, &ee2, ee2l); // ee2 = a13*dd3;
	ffl = o.Gen_Sum_With_PreAlloc(ee1l, ee1, ee2l, ee2, &ff, ffl);		  // ff = ee1 + ee2;
	if (ee1 != ee1p)
	{
		free(ee1);
		ee1 = ee1p;
	}
	if (ee2 != ee2p)
		free(ee2);
	ee1l = o.Gen_Product_With_PreAlloc(a12l, a12, dd2l, dd2, &ee1, ee1l); // ee1 = a12*dd2;

	if (a11 != a11p)
		free(a11);
	if (a12 != a12p)
		free(a12);
	if (a13 != a13p)
		free(a13);

	o.Gen_Invert(ee1l, ee1);
	detl = o.Gen_Sum_With_PreAlloc(ffl, ff, ee1l, ee1, &det, detl); // det = ff + ee1;
	if (ee1 != ee1p)
		free(ee1);
	if (ff != ffp)
		free(ff);

	double s = det[detl - 1];
	double sd = d[dl - 1];

	if (det != detp)
		free(det);

	if ((s > 0))
		return (sd > 0) ? (1) : (-1);
	if ((s < 0))
		return (sd > 0) ? (-1) : (1);
	return 0;
}

int orient3D_TPI_post_exact(
	TPI_exact_suppvars &s,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z)
{
	return orient3D_TPI_post_exact(s.d, s.dl, s.n1, s.n1l, s.n2, s.n2l, s.n3, s.n3l, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

TPI_exact_suppvars::TPI_exact_suppvars()
{
	d = dp;
	n1 = n1p;
	n2 = n2p;
	n3 = n3p;
	dl = n1l = n2l = n3l = 256;
}

TPI_exact_suppvars::~TPI_exact_suppvars()
{
	if (dl > 256)
		free(d);
	if (n1l > 256)
		free(n1);
	if (n2l > 256)
		free(n2);
	if (n3l > 256)
		free(n3);
}

int orient3D_TPI_exact(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z)
{
	TPI_exact_suppvars s;

	if (orient3D_TPI_pre_exact(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z, u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, s))
		return orient3D_TPI_post_exact(s, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
	return 0;
}

int orient3D_LPI(
	double px, double py, double pz,
	double qx, double qy, double qz,
	double rx, double ry, double rz,
	double sx, double sy, double sz,
	double tx, double ty, double tz,
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz)
{
	int r = orient3D_LPI_filtered(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, ax, ay, az, bx, by, bz, cx, cy, cz);
	if (r != Filtered_Orientation::UNCERTAIN)
		return r;
	return orient3D_LPI_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, ax, ay, az, bx, by, bz, cx, cy, cz);
}

int orient3D_TPI(
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z)
{
	int r = orient3D_TPI_filtered(
		v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
		w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
		u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
		q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
	if (r != Filtered_Orientation::UNCERTAIN)
		return r;
	return orient3D_TPI_exact(
		v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
		w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
		u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
		q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

int triangle_normal_filtered(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	double v3x = ov3x - ov2x;
	double v3y = ov3y - ov2y;
	double v3z = ov3z - ov2z;
	double v2x = ov2x - ov1x;
	double v2y = ov2y - ov1y;
	double v2z = ov2z - ov1z;
	double nvx1 = v2y * v3z;
	double nvx2 = v2z * v3y;
	double nvx = nvx1 - nvx2;
	double nvy1 = v3x * v2z;
	double nvy2 = v3z * v2x;
	double nvy = nvy1 - nvy2;
	double nvz1 = v2x * v3y;
	double nvz2 = v2y * v3x;
	double nvz = nvz1 - nvz2;

	double _tmp_fabs, max_var = 0;
	if ((_tmp_fabs = fabs(v3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = 8.88395e-016 * max_var * max_var;

	double nvxc = fabs(nvx);
	double nvyc = fabs(nvy);
	double nvzc = fabs(nvz);
	double nv = nvxc;
	if (nvyc > nv)
		nv = nvyc;
	if (nvzc > nv)
		nv = nvzc;

	if (nv > epsilon)
	{
		if (nv == nvx)
			return 0;
		if (nv == nvy)
			return 1;
		if (nv == nvz)
			return 2;
	}
	return -1;
}

int triangle_normal_exact(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	expansionObject o;
	double v3x[2];
	o.two_Diff(ov3x, ov2x, v3x);
	double v3y[2];
	o.two_Diff(ov3y, ov2y, v3y);
	double v3z[2];
	o.two_Diff(ov3z, ov2z, v3z);
	double v2x[2];
	o.two_Diff(ov2x, ov1x, v2x);
	double v2y[2];
	o.two_Diff(ov2y, ov1y, v2y);
	double v2z[2];
	o.two_Diff(ov2z, ov1z, v2z);
	double nvx1[8];
	o.Two_Two_Prod(v2y, v3z, nvx1);
	double nvx2[8];
	o.Two_Two_Prod(v2z, v3y, nvx2);
	double nvx[16];
	int nvx_len = o.Gen_Diff(8, nvx1, 8, nvx2, nvx);
	double nvy1[8];
	o.Two_Two_Prod(v3x, v2z, nvy1);
	double nvy2[8];
	o.Two_Two_Prod(v3z, v2x, nvy2);
	double nvy[16];
	int nvy_len = o.Gen_Diff(8, nvy1, 8, nvy2, nvy);
	double nvz1[8];
	o.Two_Two_Prod(v2x, v3y, nvz1);
	double nvz2[8];
	o.Two_Two_Prod(v2y, v3x, nvz2);
	double nvz[16];
	int nvz_len = o.Gen_Diff(8, nvz1, 8, nvz2, nvz);

	double nvxc = fabs(nvx[nvx_len - 1]);
	double nvyc = fabs(nvy[nvy_len - 1]);
	double nvzc = fabs(nvz[nvz_len - 1]);
	double nv = nvxc;
	if (nvyc > nv)
		nv = nvyc;
	if (nvzc > nv)
		nv = nvzc;

	if (nv == nvxc)
		return 0;
	if (nv == nvyc)
		return 1;
	if (nv == nvzc)
		return 2;

	assert(false);
	return -1;
}

#include "ip_filtered_ex.cpp"

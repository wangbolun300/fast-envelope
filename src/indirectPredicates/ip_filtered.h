#pragma once

namespace fastEnvelope
{
	class ip_filtered
	{

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
// orient3D_LPI returns 1 if the tetrahedron (i, a, b, c) has positive volume.
// Returns -1 if it has negative volume.
// Zero is returned if one of the following conditions holds:
// - degenerate input (L and P1 are parallel, <r, s, t> are collinear, ...)
// - i is exactly on P2
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.
	public:
		static int orient3D_LPI_filtered(
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
	// orient3D_TPI returns 1 if the tetrahedron (i, q1, q2, q3) has positive volume.
	// Returns -1 if it has negative volume.
	// Zero is returned if one of the following conditions holds:
	// - degenerate input (i is undefined)
	// - i is exactly on PQ
	// - the floating point precision is insufficient to guarantee an exact answer
	// - the input coordinates cause an under/overflow during the computation.

		static int orient3D_TPI_filtered(
			double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
			double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
			double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
			double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z);



		template<typename T>
		static int orient3D_LPI_filtered_multiprecision(
			double px, double py, double pz, double qx, double qy, double qz,
			double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz,
			double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz) {
			T a11, a12, a13, a21, a22, a23, a31, a32, a33, d21, d31, d22, d32, d23, d33;
			T px_rx, py_ry, pz_rz, px_cx, py_cy, pz_cz;
			T a2233, a2133, a2132;
			T d, n;
			T d11, d12, d13;
			T d2233, d2332, d2133, d2331, d2132, d2231;
			T det;

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

			if (d.is_class() == ZERO) {
				return -2;//not exist
			}

			if (d.is_class() == MIXED) {
				return 100;//failed in precision
			}
			if (det.is_class() == MIXED) {
				return 100;//failed in precision
			}

			if ((det.is_class() == POSITIVE)) {
				if (d.is_class() == POSITIVE) {
					return 1;
				}
				if (d.is_class() == NEGATIVE) {
					return -1;
				}
			}
			if ((det.is_class() == NEGATIVE)) {
				if (d.is_class() == POSITIVE) {
					return -1;
				}
				if (d.is_class() == NEGATIVE) {
					return 1;
				}
			}
			
			return 0;



		}

		static int orient3D_TPI_filtered_multiprecision(
			double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
			double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
			double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
			double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z);
	};
	

	
	

}
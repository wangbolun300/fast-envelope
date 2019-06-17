#pragma once
#include <fastenvelope/Types.hpp>
#include <vector>
#include <array>
#include <fenv.h>
#include <unordered_map>
namespace fastEnvelope {
	

	class FastEnvelope
	{
	private:
		static const int ORI_POSITIVE = 1;
		static const int ORI_ZERO = 0;
		static const int ORI_NEGATIVE = -1;
		static const int NOT_INTERSECTD = 2;
		static const int OUT_PRISM = 1;
		static const int IN_PRISM = 0;
		static const int CUT_COPLANAR = 4;
		static const int CUT_EMPTY = -1;
		static const int CUT_FACE = 3;
		//static const Scalar  BOX_SCALE = 1 / 10.0;
	public:
		FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar& eps, const int& spac);
		bool is_outside(const std::array<Vector3, 3> &triangle) const;
		inline int prism_size() const { return envprism.size(); }
		bool sample_triangle_outside(const std::array<Vector3, 3> &triangle, const Scalar sampleerror) const;
		void print_prisms(const std::array<Vector3, 3> &triangle) const;
	private:
		std::vector<std::array<Vector3, 12>> envprism;
		std::unordered_map<int, std::vector<int>> prismmap;
		std::vector<std::array<Vector3, 2>> cornerlist;
		Vector3 min, max;
		int subx, suby, subz;
	private:
		//static bool FastEnvelopeTest(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism);
		//static bool FastEnvelopeTestTemp(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism);
		static bool FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism);
	public:
		static void triangle_sample(const std::array<Vector3, 3> &triangle, std::vector<Vector3>& ps, const Scalar &error);
	private:
		static void get_bb_corners(const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
			min = vertices.front();
			max = vertices.front();

			for (size_t j = 0; j < vertices.size(); j++) {
				for (int i = 0; i < 3; i++) {
					min(i) = std::min(min(i), vertices[j](i));
					max(i) = std::max(max(i), vertices[j](i));
				}
			}

			const Scalar dis = (max - min).minCoeff() * 0.1;//TODO it used to be Parameters::box_scale

			for (int j = 0; j < 3; j++) {
				min[j] -= dis;
				max[j] += dis;
			}

			//cout << "min = " << min[0] << " " << min[1] << " " << min[2] << endl;
			//cout << "max = " << max[0] << " " << max[1] << " " << max[2] << endl;
			//            pausee();
		}
		static void  CornerList(const std::vector<std::array<Vector3, 12>>& prism,
			std::vector<std::array<Vector3, 2>>& list) {
			std::vector<Vector3> ver12(12);
			Vector3 min, max;
			list.resize(prism.size());//to be safer
			for (int i = 0; i < prism.size(); i++) {
				for (int j = 0; j < 12; j++) {
					ver12[j] = prism[i][j];
				}
				get_bb_corners(ver12, min, max);
				list[i] = { {min,max } };
			}
		}
		static void BoxFindCells(const Vector3& min, const Vector3& max,
			const Vector3& cellmin, const Vector3& cellmax, const int& subx, const int&suby, const int subz, std::vector<int>& intercell) {

			Vector3 delta;
			delta[0] = (cellmax - cellmin)[0] / subx;
			delta[1] = (cellmax - cellmin)[1] / suby;
			delta[2] = (cellmax - cellmin)[2] / subz;
			//intercell.reserve(int((max - min)[0] / delta[0])*int((max - min)[1] / delta[1])*int((max - min)[2] / delta[2]));
			intercell.clear();
			int location[2][3];
			for (int i = 0; i < 3; i++) {
				location[0][i] = (min[i] - cellmin[i]) / delta[i];
			}
			for (int i = 0; i < 3; i++) {
				location[1][i] = (max[i] - cellmin[i]) / delta[i];
			}
			for (int i = location[0][0]; i <= location[1][0]; i++) {
				for (int j = location[0][1]; j <= location[1][1]; j++) {
					for (int k = location[0][2]; k <= location[1][2]; k++) {
						intercell.emplace_back(k*subx*suby + j * subx + i);
					}
				}
			}


		}
		static void get_triangle_corners(const std::array<Vector3, 3> &triangle, Vector3 &mint, Vector3 &maxt) {
			mint[0] = std::min(std::min(triangle[0][0], triangle[1][0]), triangle[2][0]);
			mint[1] = std::min(std::min(triangle[0][1], triangle[1][1]), triangle[2][1]);
			mint[2] = std::min(std::min(triangle[0][2], triangle[1][2]), triangle[2][2]);
			maxt[0] = std::max(std::max(triangle[0][0], triangle[1][0]), triangle[2][0]);
			maxt[1] = std::max(std::max(triangle[0][1], triangle[1][1]), triangle[2][1]);
			maxt[2] = std::max(std::max(triangle[0][2], triangle[1][2]), triangle[2][2]);

		}

		// to check if a point is in the prisms. the jump index shows the prisms not counted in calculation, and jump is sorted from small to big
		static bool point_out_prism(const Vector3& point, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump);

		static void BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, const Scalar& epsilon);
		static int orient_3triangles(const Eigen::Matrix<Scalar, 3, 3>& A,const Eigen::Matrix<Scalar, 3, 3>& AT,
			const Eigen::Matrix<Scalar, 3, 3>& ATA,const Eigen::Matrix<Scalar, 3, 1>& B, const std::array<Vector3, 3> & triangle3);
		static int Implicit_Seg_Facet_interpoint_Out_Prism(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
			const std::vector<std::array<Vector3, 12>>& envprism,const std::vector<int>& jump);
	
		static int Implicit_Tri_Facet_Facet_interpoint_Out_Prism( const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2,
			const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump);
		static int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_M(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2,
			const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump);
		static bool is_seg_facet_intersection(const double& px, const double& py, const double& pz,
			const double& qx, const double& qy, const double& qz,
			const double& rx, const double& ry, const double& rz,
			const double& sx, const double& sy, const double& sz,
			const double& tx, const double& ty, const double& tz,
			double& a11, double& a12, double& a13,
			double& a21, double& a22, double& a23,
			double& a31, double& a32, double& a33,
			double& px_rx, double& py_ry, double& pz_rz,
			double& d, double& n);


		static int tri_cut_tri_simple(const Vector3& p1, const Vector3& p2, const Vector3& p3,const Vector3& q1, const Vector3& q2, const Vector3& q3);
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

		static	inline double det4x4(
			double a11, double a12, double a13, double a14,
			double a21, double a22, double a23, double a24,
			double a31, double a32, double a33, double a34,
			double a41, double a42, double a43, double a44)
		{
			double m12 = ((a21 * a12) - (a11 * a22));
			double m13 = ((a31 * a12) - (a11 * a32));
			double m14 = ((a41 * a12) - (a11 * a42));
			double m23 = ((a31 * a22) - (a21 * a32));
			double m24 = ((a41 * a22) - (a21 * a42));
			double m34 = ((a41 * a32) - (a31 * a42));
			double m123 = (((m23 * a13) - (m13 * a23)) + (m12 * a33));
			double m124 = (((m24 * a13) - (m14 * a23)) + (m12 * a43));
			double m134 = (((m34 * a13) - (m14 * a33)) + (m13 * a43));
			double m234 = (((m34 * a23) - (m24 * a33)) + (m23 * a43));
			return ((((m234 * a14) - (m134 * a24)) + (m124 * a34)) - (m123 * a44));
		}

		static inline double det3x3(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33) {
			return (((a11 * ((a22 * a33) - (a23 * a32))) - (a12 * ((a21 * a33) - (a23 * a31)))) + (a13 * ((a21 * a32) - (a22 * a31))));
		}

		static inline double dot(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z) { return (((v1x * v2x) + (v1y * v2y)) + (v1z * v2z)); }
		
		static inline double det2x2(double a1, double a2, double b1, double b2) { return ((a1 * b2) - (a2 * b1)); }

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

		static inline int orient3D_LPI(//right hand law
			const double& px, const double& py, const double& pz,
			const double& ax, const double& ay, const double& az,
			const double& bx, const double& by, const double& bz,
			const double& cx, const double& cy, const double& cz,
			const double& a11, const double& a12, const double& a13,
			const double& a21, const double& a22, const double& a23,
			const double& a31, const double& a32, const double& a33,
			const double& px_rx, const double& py_ry, const double& pz_rz,
			const double& d, const double& n);



		static bool is_3triangle_intersect(double& v1x,  double& v1y, double& v1z, double& v2x, double& v2y, double& v2z, double& v3x, double& v3y, double& v3z,
			double& w1x, double& w1y, double& w1z, double& w2x, double& w2y, double& w2z, double& w3x, double& w3y, double& w3z,
			double& u1x, double& u1y, double& u1z, double& u2x, double& u2y, double& u2z, double& u3x, double& u3y, double& u3z,
			double& n1, double& n2, double& n3, double& d) {
			
			::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
			
			double t2x = v2x, t2y = v2y, t2z = v2z, t3x = v3x, t3y = v3y, t3z = v3z;//still need these
			
			v3x -= v2x;
			v3y -= v2y;
			v3z -= v2z;
			v2x -= v1x;
			v2y -= v1y;
			v2z -= v1z;
			double nvx = det2x2(v2y, v2z, v3y, v3z);//v3=v2v3, v2=v1v2;v4=v3v1
			double nvy = det2x2(v3x, v3z, v2x, v2z);
			double nvz = det2x2(v2x, v2y, v3x, v3y);

			w3x -= w2x;
			w3y -= w2y;
			w3z -= w2z;
			w2x -= w1x;
			w2y -= w1y;
			w2z -= w1z;
			double nwx = det2x2(w2y, w2z, w3y, w3z);
			double nwy = det2x2(w3x, w3z, w2x, w2z);
			double nwz = det2x2(w2x, w2y, w3x, w3y);

			u3x -= u2x;
			u3y -= u2y;
			u3z -= u2z;
			u2x -= u1x;
			u2y -= u1y;
			u2z -= u1z;
			double nux = det2x2(u2y, u2z, u3y, u3z);
			double nuy = det2x2(u3x, u3z, u2x, u2z);
			double nuz = det2x2(u2x, u2y, u3x, u3y);

			double nwyuz = nwy * nuz - nwz * nuy;
			double nwxuz = nwx * nuz - nwz * nux;
			double nwxuy = nwx * nuy - nwy * nux;

			double nvyuz = nvy * nuz - nvz * nuy;
			double nvxuz = nvx * nuz - nvz * nux;
			double nvxuy = nvx * nuy - nvy * nux;

			double nvywz = nvy * nwz - nvz * nwy;
			double nvxwz = nvx * nwz - nvz * nwx;
			double nvxwy = nvx * nwy - nvy * nwx;

			d = nvx * nwyuz - nvy * nwxuz + nvz * nwxuy;

			double p1 = nvx * v1x + nvy * v1y + nvz * v1z;
			double p2 = nwx * w1x + nwy * w1y + nwz * w1z;
			double p3 = nux * u1x + nuy * u1y + nuz * u1z;

			 n1 = p1 * nwyuz - p2 * nvyuz + p3 * nvywz;
			 n2 = p2 * nvxuz - p3 * nvxwz - p1 * nwxuz;
			 n3 = p3 * nvxwy - p2 * nvxuy + p1 * nwxuy;

			

			if (d < SCALAR_ZERO_3) {// if not intersected
				return 0;
			}
		
			

			////////////////////////////////////////////////////////////////
			// this part is to predicate if the point is in the interior of triangle v
			//test : (n,v1,t2),(n,t2,t3),(n,t3,v1)
			Scalar nx, ny, nz;
			Vector3 n = Vector3(v1x + nvx, v1y + nvy, v1z + nvz);
			//n = n.normalized();//TODO maybe not needed, and this is wrong :)
			nx = n(0); ny = n(1); nz = n(2);
			int ori1, ori2, ori3;
			ori1 = orient3D_TPI(
				v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
				w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
				u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
				nx, ny, nz, v1x, v1y, v1z, t2x, t2y, t2z,
				n1,n2,n3, d);
			ori2 = orient3D_TPI(
				v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
				w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
				u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
				nx, ny, nz, t2x, t2y, t2z, t3x, t3y, t3z,
				n1,n2,n3, d);

			ori3 = orient3D_TPI(
				v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
				w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
				u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
				nx, ny, nz, t3x, t3y, t3z, v1x, v1y, v1z,
				n1,n2,n3, d);
			int mark1 = ori1 * ori2, mark2 = ori1 * ori3;
			///////////////////////////////////////////////////
			if (mark1 > 0 && mark2 > 0) {
				return 1;
			}
			return 0;
		}



		static inline int orient3D_TPI(
			const double& v1x, const double& v1y, const double& v1z, const double& v2x, const double& v2y, const double& v2z, const double& v3x, const double& v3y, const double& v3z,
			const double& w1x, const double& w1y, const double& w1z, const double& w2x, const double& w2y, const double& w2z, const double& w3x, const double& w3y, const double& w3z,
			const double& u1x, const double& u1y, const double& u1z, const double& u2x, const double& u2y, const double& u2z, const double& u3x, const double& u3y, const double& u3z,
			const double& q1x, const double& q1y, const double& q1z, const double& q2x, const double& q2y, const double& q2z, const double& q3x, const double& q3y, const double& q3z,
			const double& n1,const double& n2, const double& n3, const double& d)
		{


			// If the same intersection point must be tested against several planes,
			// code up to here can be extracted and computed only once.

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
			double fu2y = fabs(u2y);
			double fu2z = fabs(u2z);
			double fu3x = fabs(u3x);
			double fu3y = fabs(u3y);
			double fu3z = fabs(u3z);

			double fa21 = fabs(a21);
			double fa22 = fabs(a22);
			double fa23 = fabs(a23);
			double fa31 = fabs(a31);
			double fa32 = fabs(a32);
			double fa33 = fabs(a33);

			double max1, max2, max3, max4, max5, max6, max7, max8;

			max4 = fv2y;
			if (max4 < fv3y) max4 = fv3y;
			if (max4 < fw3y) max4 = fw3y;
			if (max4 < fw2y) max4 = fw2y;
			if (max4 < fu2y) max4 = fu2y;
			if (max4 < fu3y) max4 = fu3y;
			max1 = max4;
			max2 = fv3x;
			if (max2 < fv2x) max2 = fv2x;
			if (max2 < fu2x) max2 = fu2x;
			if (max2 < fu3x) max2 = fu3x;
			if (max2 < fw2x) max2 = fw2x;
			if (max2 < fw3x) max2 = fw3x;
			if (max1 < max2) max1 = max2;
			max5 = fv2z;
			if (max5 < fu2z) max5 = fu2z;
			if (max5 < fv3z) max5 = fv3z;
			if (max5 < fw3z) max5 = fw3z;
			if (max5 < fw2z) max5 = fw2z;
			if (max5 < fu3z) max5 = fu3z;
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
			max7 = fu2x;
			if (max7 < fu3x) max7 = fu3x;
			if (max7 < fw2x) max7 = fw2x;
			if (max7 < fw3x) max7 = fw3x;
			if (max7 < fa21) max7 = fa21;
			if (max7 < fa31) max7 = fa31;
			max8 = fa22;
			if (max8 < fa23) max8 = fa23;
			if (max8 < fa33) max8 = fa33;
			if (max8 < fa32) max8 = fa32;

			double eps = 3.4025182954957945e-012 * (((((((max1 * max3) * max2) * max5) * max7) * max4) * max6) * max8);

			if ((det > eps)) return (d > 0) ? (1) : (-1);
			if ((det < -eps)) return (d > 0) ? (-1) : (1);
			return 0;
		}
		

	
	};
}

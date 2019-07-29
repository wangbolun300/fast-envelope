#pragma once
#include <fastenvelope/Types.hpp>
#include <vector>
#include <array>
#include <fenv.h>
#include <unordered_map>
#include<iostream>
#include<arbitraryprecision/fprecision.h>
#include<arbitraryprecision/intervalprecision.h>
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
		static const int CUT_COPLANAR_OR_ABOVE = 4;
		static const int CUT_EMPTY = -1;
		static const int CUT_FACE = 3;

		static const int NOT_DEGENERATED = 0;
		static const int NERLY_DEGENERATED = 1;
		static const int DEGENERATED_SEGMENT = 2;
		static const int DEGENERATED_POINT = 3;
		//static const Scalar  BOX_SCALE = 1 / 10.0;
	public:
		FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps, const int spac);
		bool is_outside(const std::array<Vector3, 3> &triangle) const;
		inline int prism_size() const { return envprism.size(); }
		bool sample_triangle_outside(const std::array<Vector3, 3> &triangle, const Scalar sampleerror) const;
		void print_prisms(const std::array<Vector3, 3> &triangle) const;
		static int test_dege(const std::array<Vector3, 3>& triangle) {
			return is_triangle_degenerated(triangle);
		}
		static Vector3 accurate_normal_vector(const std::array<Vector3, 3> & triangle, const int &digit);

	private:
		

		std::vector<std::array<Vector3, 12>> envprism;
		std::unordered_map<int, std::vector<int>> prismmap;
		std::vector<std::array<Vector3, 2>> cornerlist;
		Vector3 min, max;
		int subx, suby, subz;
	private:
		//static bool FastEnvelopeTest(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism);
		//static bool FastEnvelopeTestTemp(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism);
		bool FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism) const;

		static int seg_cut_tri(const Vector3 & seg0, const Vector3 &seg1, const Vector3&t0, const Vector3&t1, const Vector3 &t2);
	public:
		static void triangle_sample(const std::array<Vector3, 3> &triangle, std::vector<Vector3>& ps, const Scalar &error);
		static void print_number();
	private:
		//static void seg_tri_cut_on_tir()()
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
		static bool point_out_prism(const Vector3& point, const std::vector<std::array<Vector3, 12>>& envprism, const int& jump);

		static void BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, const Scalar& epsilon);


		static int Implicit_Seg_Facet_interpoint_Out_Prism_redundant(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
			const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump);

		static int Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
			const std::vector<std::array<Vector3, 12>>& envprism, const int& jump);

		static int Implicit_prism_edge_triangle_interpoint_Out_Prism_multi_precision(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
			const std::vector<std::array<Vector3, 12>>& envprism, const int& jump);


		static int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_redundant(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2,
			const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump);

		static int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_multi_precision(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2,
			const std::vector<std::array<Vector3, 12>>& envprism, const int& jump1,const int &jump2);



		static int tri_cut_tri_simple(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& q1, const Vector3& q2, const Vector3& q3);





		static bool is_3_triangle_cut(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& f1, const std::array<Vector3, 3>& f2);
		public:
		static int is_triangle_degenerated(const std::array<Vector3, 3>& triangle);
		private:
		static Vector2 to_2d(const Vector3 &p, int t) {
			return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
		}

		static arbitrary_precision::interval<arbitrary_precision::float_precision> converting_Scalar_to_arbitary(const Scalar &a, const int &i);
		public:
		template<typename T>
		static bool orient3D_LPI_prefilter_multiprecision(
			const T& px, const T& py, const T& pz, const T& qx, const T& qy, const T& qz,
			const T& rx, const T& ry, const T& rz, const T& sx, const T& sy, const T& sz, const T& tx, const T& ty, const T& tz,
			T& a11, T& a12, T& a13, T& d, const std::function<int(T)> &checker) {
			
			a11 = (px - qx);// TODO is that ok?
			a12 = (py - qy);
			a13 = (pz - qz);
			T a21(sx - rx);
			T a22(sy - ry);
			T a23(sz - rz);
			T a31(tx - rx);
			T a32(ty - ry);
			T a33(tz - rz);
			T a2233((a22 * a33) - (a23 * a32));
			T a2133((a21 * a33) - (a23 * a31));
			T a2132((a21 * a32) - (a22 * a31));
			d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));//TODO maybe not safe
			int flag1 = checker(d);
			if (flag1 == -2 || flag1 == 0) {
				return false;// not enough precision
			}
			T px_rx( px - rx);
			T py_ry( py - ry);
			T pz_rz( pz - rz);

			T n ((((py_ry)* a2133) - ((px_rx)* a2233)) - ((pz_rz)* a2132));

			a11 = a11 * n;
			a12 = a12 * n;
			a13 = a13 * n;
			return true;
		}

		template<typename T>
		static bool orient3D_TPI_prefilter_multiprecision(
			const T& ov1x, const T& ov1y, const T& ov1z, const T& ov2x, const T& ov2y, const T& ov2z, const T& ov3x, const T& ov3y, const T& ov3z,
			const T& ow1x, const T& ow1y, const T& ow1z, const T& ow2x, const T& ow2y, const T& ow2z, const T& ow3x, const T& ow3y, const T& ow3z,
			const T& ou1x, const T& ou1y, const T& ou1z, const T& ou2x, const T& ou2y, const T& ou2z, const T& ou3x, const T& ou3y, const T& ou3z,
			T& d, T& n1, T& n2, T& n3, const std::function<int(T)> &checker
		)
		{
			::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

			T v3x(ov3x - ov2x);
			T v3y(ov3y - ov2y);
			T v3z(ov3z - ov2z);
			T v2x(ov2x - ov1x);
			T v2y(ov2y - ov1y);
			T v2z(ov2z - ov1z);
			T w3x(ow3x - ow2x);
			T w3y(ow3y - ow2y);
			T w3z(ow3z - ow2z);
			T w2x(ow2x - ow1x);
			T w2y(ow2y - ow1y);
			T w2z(ow2z - ow1z);
			T u3x(ou3x - ou2x);
			T u3y(ou3y - ou2y);
			T u3z(ou3z - ou2z);
			T u2x(ou2x - ou1x);
			T u2y(ou2y - ou1y);
			T u2z(ou2z - ou1z);
			
			T nvx(v2y * v3z - v2z * v3y);
			T nvy(v3x * v2z - v3z * v2x);
			T nvz(v2x * v3y - v2y * v3x);

			T nwx(w2y * w3z - w2z * w3y);
			T nwy(w3x * w2z - w3z * w2x);
			T nwz(w2x * w3y - w2y * w3x);
			
			T nux(u2y * u3z - u2z * u3y);
			T nuy(u3x * u2z - u3z * u2x);
			T nuz(u2x * u3y - u2y * u3x);
			
			T nwyuz(nwy * nuz - nwz * nuy);
			T nwxuz(nwx * nuz - nwz * nux);
			T nwxuy(nwx * nuy - nwy * nux);
			
			T nvyuz(nvy * nuz - nvz * nuy);
			T nvxuz(nvx * nuz - nvz * nux);
			T nvxuy(nvx * nuy - nvy * nux);
			
			T nvywz(nvy * nwz - nvz * nwy);
			T nvxwz(nvx * nwz - nvz * nwx);
			T nvxwy(nvx * nwy - nvy * nwx);

			d = (nvx * nwyuz - nvy * nwxuz + nvz * nwxuy);
			
			

			int flag1 = checker(d);
			if (flag1 == -2 || flag1 == 0) {
				return false;// not enough precision
			}

			T p1(nvx * ov1x + nvy * ov1y + nvz * ov1z);
			T p2(nwx * ow1x + nwy * ow1y + nwz * ow1z);
			T p3(nux * ou1x + nuy * ou1y + nuz * ou1z);

			n1 = p1 * nwyuz - p2 * nvyuz + p3 * nvywz;
			n2 = p2 * nvxuz - p3 * nvxwz - p1 * nwxuz;
			n3 = p3 * nvxwy - p2 * nvxuy + p1 * nwxuy;
			return true;
		}

		template<typename T>
		static int orient3D_TPI_postfilter_multiprecision(
			const T& d, const T& n1, const T& n2, const T& n3,
			const T& q1x, const T& q1y, const T& q1z, const T& q2x, const T& q2y, const T& q2z, const T& q3x, const T& q3y, const T& q3z, const std::function<int(T)> &checker
		)
		{
			::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

			T dq3x(d * q3x);
			T dq3y(d * q3y);
			T dq3z(d * q3z);

			T a11(n1 - dq3x);
			T a12(n2 - dq3y);
			T a13(n3 - dq3z);
			T a21(q1x - q3x);
			T a22(q1y - q3y);
			T a23(q1z - q3z);
			T a31(q2x - q3x);
			T a32(q2y - q3y);
			T a33(q2z - q3z);

			T det(a11 * (a22*a33 - a23 * a32) - a12 * (a21*a33 - a23 * a31) + a13 * (a21*a32 - a22 * a31));

			int flag2 = checker(det);
			if (flag2 == -2) {
				return 100;// not enough precision
			}
			if (flag2 == 1) {
				if (d > 0) {
					return 1;
				}
				if (d < 0) {
					return -1;
				}
			}
			if (flag2 == -1) {
				if (d > 0) {
					return -1;
				}
				if (d < 0) {
					return 1;
				}
			}
			return 0;
		}

		template<typename T>
		static int orient3D_LPI_postfilter_multiprecision(
			const T& a11, const T& a12, const T& a13, const T& d,
			const T& px, const T& py, const T& pz,
			const T& ax, const T& ay, const T& az,
			const T& bx, const T& by, const T& bz,
			const T& cx, const T& cy, const T& cz, const std::function<int(T)> &checker) {
			
			T px_cx(px - cx);
			T py_cy(py - cy);
			T pz_cz(pz - cz);

			T d11((d * px_cx) + (a11));
			T d21(ax - cx);
			T d31(bx - cx);
			T d12((d * py_cy) + (a12));
			T d22(ay - cy);
			T d32(by - cy);
			T d13((d * pz_cz) + (a13));
			T d23(az - cz);
			T d33(bz - cz);

			T d2233(d22 * d33);
			T d2332(d23 * d32);
			T d2133(d21 * d33);
			T d2331(d23 * d31);
			T d2132(d21 * d32);
			T d2231(d22 * d31);
 
			T det(d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231));

			int flag2 = checker(det);
			if (flag2 == -2) {
				return 100;// not enough precision
			}
			if (flag2 == 1) {
				if (d > 0) {
					return 1;
				}
				if (d < 0) {
					return -1;
				}
			}
			if (flag2 == -1) {
				if (d > 0) {
					return -1;
				}
				if (d < 0) {
					return 1;
				}
			}
			return 0;
		}
		public:
		template<typename T>

		static int orient3D_TPI_filtered_multiprecision(
			T v1x, T v1y, T v1z, T v2x, T v2y, T v2z, T v3x, T v3y, T v3z,
			T w1x, T w1y, T w1z, T w2x, T w2y, T w2z, T w3x, T w3y, T w3z,
			T u1x, T u1y, T u1z, T u2x, T u2y, T u2z, T u3x, T u3y, T u3z,
			T q1x, T q1y, T q1z, T q2x, T q2y, T q2z, T q3x, T q3y, T q3z, const std::function<int(T)> &checker) {

			::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

			v3x = v3x - v2x;
			v3y = v3y - v2y;
			v3z = v3z - v2z;
			v2x = v2x - v1x;
			v2y = v2y - v1y;
			v2z = v2z - v1z;
			T  nvx(v2y * v3z - v2z * v3y);
			T  nvy(v3x * v2z - v3z * v2x);
			T  nvz(v2x * v3y - v2y * v3x);

			w3x = w3x - w2x;
			w3y = w3y - w2y;
			w3z = w3z - w2z;
			w2x = w2x - w1x;
			w2y = w2y - w1y;
			w2z = w2z - w1z;
			T  nwx ( w2y * w3z - w2z * w3y);
			T  nwy ( w3x * w2z - w3z * w2x);
			T  nwz ( w2x * w3y - w2y * w3x);

			u3x = u3x - u2x;
			u3y = u3y - u2y;
			u3z = u3z - u2z;
			u2x = u2x - u1x;
			u2y = u2y - u1y;
			u2z = u2z - u1z;
			T  nux ( u2y * u3z - u2z * u3y);
			T  nuy ( u3x * u2z - u3z * u2x);
			T  nuz ( u2x * u3y - u2y * u3x);

			T  nwyuz ( nwy * nuz - nwz * nuy);
			T  nwxuz ( nwx * nuz - nwz * nux);
			T  nwxuy ( nwx * nuy - nwy * nux);

			T  nvyuz ( nvy * nuz - nvz * nuy);
			T  nvxuz ( nvx * nuz - nvz * nux);
			T  nvxuy ( nvx * nuy - nvy * nux);

			T  nvywz ( nvy * nwz - nvz * nwy);
			T  nvxwz ( nvx * nwz - nvz * nwx);
			T  nvxwy ( nvx * nwy - nvy * nwx);

			T  d ( nvx * nwyuz - nvy * nwxuz + nvz * nwxuy);

			int flag1 = checker(d);
			if (flag1 == -2) {
				return 100;// not enough precision
			}
			if (flag1 == 0) {
				return -2;// not exist
			}

			T  p1 ( nvx * v1x + nvy * v1y + nvz * v1z);
			T  p2 ( nwx * w1x + nwy * w1y + nwz * w1z);
			T  p3 ( nux * u1x + nuy * u1y + nuz * u1z);

			T  n1 ( p1 * nwyuz - p2 * nvyuz + p3 * nvywz);
			T  n2 ( p2 * nvxuz - p3 * nvxwz - p1 * nwxuz);
			T  n3 ( p3 * nvxwy - p2 * nvxuy + p1 * nwxuy);

			T  dq3x ( d * q3x);
			T  dq3y ( d * q3y);
			T  dq3z ( d * q3z);

			T  a11 ( n1 - dq3x);
			T  a12 ( n2 - dq3y);
			T  a13 ( n3 - dq3z);
			T  a21 ( q1x - q3x);
			T  a22 ( q1y - q3y);
			T  a23 ( q1z - q3z);
			T  a31 ( q2x - q3x);
			T  a32 ( q2y - q3y);
			T  a33 ( q2z - q3z);

			T  det ( a11 * (a22*a33 - a23 * a32) - a12 * (a21*a33 - a23 * a31) + a13 * (a21*a32 - a22 * a31));

			if (::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID)) return 0; // Fast reject in case of under/overflow

			int flag2 = checker(det);
			if (flag2 == -2) {
				return 100;// not enough precision
			}
			if (flag2 == 1) {
				if (flag1 == 1) {
					return 1;
				}
				if (flag1 == -1) {
					return -1;
				}
			}
			if (flag2 == -1) {
				if (flag1 == 1) {
					return -1;
				}
				if (flag1 == -1) {
					return 1;
				}
			}
			return 0;

		}



	};
}

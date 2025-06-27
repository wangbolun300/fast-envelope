#pragma once
#include <fastenvelope/Types.hpp>
#include <fastenvelope/AABB.h>

#include <vector>
#include <array>
#include <fenv.h>

#include <iostream>

namespace fastEnvelope
{

	class FastEnvelope
	{

	public:
		FastEnvelope(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, const Scalar eps);
		FastEnvelope() {}

		void init(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, const Scalar eps);
		void init(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, const std::vector<Scalar> eps);
		void init(const std::vector<std::vector<std::array<Vector3, 3>>> halfspace_input,
				  const std::vector<std::array<Vector3, 2>> cornerlist_input, const Scalar eps);
		static void printnumber();

		// check if tri or point is outside
		bool is_outside(const std::array<Vector3, 3> &triangle) const;
		bool is_outside(const Vector3 &point) const;
		bool is_outside(const Vector3 &point0, const Vector3 &point1) const;
		bool is_outside_no_optimazation(const std::array<Vector3, 3> &triangle) const;

	private:
		struct INDEX
		{
			int Pi;
			std::vector<int> FACES;
		};

	private:
		AABB tree;
		std::vector<Vector3> m_ver_p;
		std::vector<Vector3i> m_faces_p;
		std::vector<std::array<Vector3, 2>> cornerlist;
		std::vector<std::vector<std::array<Vector3, 3>>> halfspace;
		bool USE_ADJACENT_INFORMATION = true;

		bool debugcode(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const;
		bool triangle_out_simple(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const;

		bool triangle_out_of_envelope(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const;
		bool segment_out_of_envelope(const Vector3 &seg0, const Vector3 &seg1, const std::vector<unsigned int> &prismindex) const;
		bool is_two_facets_neighbouring(const int &pid, const int &i, const int &j) const;

		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const;
#ifdef ENVELOPE_WITH_GMP
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_Rational(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const;
#endif
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>> &faceorder, const int &jump, int &id) const;
#ifdef ENVELOPE_WITH_GMP
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_Rational(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>> &intersect_face, const int &jump, int &id) const;
#endif

		// this function check the adjacent polyhedrons and jump over the polyhedrons that already in the cover list
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>> &intersect_face, const std::vector<bool> &coverlist, const int &jump, int &id) const;

#ifdef ENVELOPE_WITH_GMP
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over_Rational(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>> &intersect_face, const std::vector<bool> &coverlist, const int &jump, int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over_Multiprecision(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>> &intersect_face, const std::vector<bool> &coverlist, const int &jump, int &id) const;
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_Rational(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const int &jump1, const int &jump2) const;

#endif
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const int &jump1, const int &jump2) const;
		// this function check the adjacent polyhedrons and jump over the polyhedrons that already in the cover list
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>> &intersect_face, const std::vector<bool> &coverlist, const int &jump1, const int &jump2,
			int &id) const;

		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>> &intersect_face, const int &jump1, const int &jump2,
			int &id) const;

#ifdef ENVELOPE_WITH_GMP
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_Rational(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>> &intersect_face, const int &jump1, const int &jump2,
			int &id) const;
#endif
		static bool is_3_triangle_cut_pure_multiprecision(const std::array<Vector3, 3> &triangle, TPI_exact_suppvars &s);

		static bool is_3_triangle_cut(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12,
			const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22);
		bool is_tpp_on_polyhedra(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const int &prismid, const int &faceid) const;

		bool point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const;
		bool point_out_prism_return_local_id(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const;

		// heuristics to refine the box box guess from the tree
		static int is_3_triangle_cut_float_fast(
			const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12,
			const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22);
#ifdef ENVELOPE_WITH_GMP
		static bool is_3_triangle_cut_Rational(const std::array<Vector3, 3> &triangle,
											   const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22);
#endif
		int is_triangle_cut_envelope_polyhedra(const int &cindex,
											   const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, std::vector<int> &cid) const;
		bool is_seg_cut_polyhedra(const int &cindex,
								  const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const;
#ifdef ENVELOPE_WITH_GMP
		bool is_tpp_on_polyhedra_Rational(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const int &prismid, const int &faceid) const;
#endif
		template <typename T>
		static bool orient3D_LPI_prefilter_multiprecision(
			const T &px, const T &py, const T &pz, const T &qx, const T &qy, const T &qz,
			const T &rx, const T &ry, const T &rz, const T &sx, const T &sy, const T &sz, const T &tx, const T &ty, const T &tz,
			T &a11, T &a12, T &a13, T &d, const std::function<int(T)> &checker)
		{

			a11 = (px - qx);
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
			d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132)); // TODO maybe not safe
			int flag1 = checker(d);
			if (flag1 == -2 || flag1 == 0)
			{
				return false; // not enough precision
			}
			T px_rx(px - rx);
			T py_ry(py - ry);
			T pz_rz(pz - rz);

			T n((((py_ry)*a2133) - ((px_rx)*a2233)) - ((pz_rz)*a2132));

			a11 = a11 * n;
			a12 = a12 * n;
			a13 = a13 * n;
			return true;
		}

		template <typename T>
		static bool orient3D_TPI_prefilter_multiprecision(
			const T &ov1x, const T &ov1y, const T &ov1z, const T &ov2x, const T &ov2y, const T &ov2z, const T &ov3x, const T &ov3y, const T &ov3z,
			const T &ow1x, const T &ow1y, const T &ow1z, const T &ow2x, const T &ow2y, const T &ow2z, const T &ow3x, const T &ow3y, const T &ow3z,
			const T &ou1x, const T &ou1y, const T &ou1z, const T &ou2x, const T &ou2y, const T &ou2z, const T &ou3x, const T &ou3y, const T &ou3z,
			T &d, T &n1, T &n2, T &n3, const std::function<int(T)> &checker)
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
			if (flag1 == -2 || flag1 == 0)
			{
				return false; // not enough precision
			}

			T p1(nvx * ov1x + nvy * ov1y + nvz * ov1z);
			T p2(nwx * ow1x + nwy * ow1y + nwz * ow1z);
			T p3(nux * ou1x + nuy * ou1y + nuz * ou1z);

			n1 = p1 * nwyuz - p2 * nvyuz + p3 * nvywz;
			n2 = p2 * nvxuz - p3 * nvxwz - p1 * nwxuz;
			n3 = p3 * nvxwy - p2 * nvxuy + p1 * nwxuy;
			return true;
		}

		template <typename T>
		static int orient3D_LPI_postfilter_multiprecision(
			const T &a11, const T &a12, const T &a13, const T &d,
			const T &px, const T &py, const T &pz,
			const T &ax, const T &ay, const T &az,
			const T &bx, const T &by, const T &bz,
			const T &cx, const T &cy, const T &cz, const std::function<int(T)> &checker)
		{

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
			if (flag2 == -2)
			{
				return 100; // not enough precision, only happens when using floating points
			}
			if (flag2 == 1)
			{
				if (d > 0)
				{
					return 1;
				}
				if (d < 0)
				{
					return -1;
				}
			}
			if (flag2 == -1)
			{
				if (d > 0)
				{
					return -1;
				}
				if (d < 0)
				{
					return 1;
				}
			}
			return 0;
		}

		template <typename T>
		static int orient3D_TPI_postfilter_multiprecision(
			const T &d, const T &n1, const T &n2, const T &n3,
			const T &q1x, const T &q1y, const T &q1z, const T &q2x, const T &q2y, const T &q2z, const T &q3x, const T &q3y, const T &q3z, const std::function<int(T)> &checker)
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

			T det(a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31));

			int flag2 = checker(det);
			if (flag2 == -2)
			{
				return 100; // not enough precision
			}
			if (flag2 == 1)
			{
				if (d > 0)
				{
					return 1;
				}
				if (d < 0)
				{
					return -1;
				}
			}
			if (flag2 == -1)
			{
				if (d > 0)
				{
					return -1;
				}
				if (d < 0)
				{
					return 1;
				}
			}
			return 0;
		}
	};

}
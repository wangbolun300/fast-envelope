#pragma once
#include <fastenvelope/Types.hpp>
#include <fastenvelope/AABB.h>
#include <vector>
#include <array>
#include <fenv.h>
//#include <unordered_map>
#include<iostream>
//#include<arbitraryprecision/fprecision.h>
//#include<arbitraryprecision/intervalprecision.h>
#include <fastenvelope/Multiprecision.hpp>
#include <fastenvelope/Rational.hpp>
#include <geogram/mesh/mesh.h>
#include <fastenvelope/ip_filtered.h>

namespace fastEnvelope {
	class AABB;

	class FastEnvelope
	{
	private:

		static const int NOT_INTERSECTED = 2;
		static const int INTERSECTED = 1;
		static const int OUT_PRISM = 1;
		static const int IN_PRISM = 0;
		static const int CUT_COPLANAR = 4;
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
		bool sample_triangle_outside(const std::array<Vector3, 3> &triangle, const int& pieces) const;
		void print_prisms(const std::array<Vector3, 3> &triangle) const;
		static void print_number();
		static void print_ini_number();

		static int is_triangle_degenerated(const Vector3& triangle0, const Vector3& triangle1, const Vector3& triangle2);
		static bool is_triangle_cut_bounding_box(
			const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, const Vector3 &bmin, const Vector3 &bmax);
		std::vector<std::array<Vector3, 2>> cornerlist;

	private:
		AABB tree;

		std::vector<std::array<Vector3, 12>> envprism;


		int subx, suby, subz;
		int prism_size;
		std::vector<std::array<Vector3, 8>> envcubic;

		bool FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int>& prismindex) const;

		static int seg_cut_plane(const Vector3 & seg0, const Vector3 &seg1, const Vector3&t0, const Vector3&t1, const Vector3 &t2);
		static Vector3 accurate_normal_vector(const Vector3 & p0, const Vector3 & p1,
			const Vector3 & q0, const Vector3 & q1);

		static void triangle_sample_segment(const std::array<Vector3, 3> &triangle, Vector3& ps, const int &pieces,const int & nbr);
		static void triangle_sample_point(const std::array<Vector3, 3> &triangle, Vector3& ps);
		static void triangle_sample_normal(const std::array<Vector3, 3> &triangle, Vector3& ps, const int &pieces, const int & nbr1, const int &nbr2);

		static void triangle_sample_normal_rational(const std::array<Vector3, 3> &triangle, Rational& ps0, Rational& ps1, Rational& ps2, const int &pieces, const int & nbr1, const int &nbr2);
		static void get_bb_corners(const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {//TODO why use this one
			min = vertices.front();
			max = vertices.front();

			for (size_t j = 0; j < vertices.size(); j++) {
				for (int i = 0; i < 3; i++) {
					min(i) = std::min(min(i), vertices[j](i));
					max(i) = std::max(max(i), vertices[j](i));
				}
			}

			//const Scalar dis = (max - min).minCoeff() *0;//TODO  change to 1e-5 or sth
			////const Scalar dis = 1e-4;
			//for (int j = 0; j < 3; j++) {
			//	min[j] -= dis;
			//	max[j] += dis;
			//}

		}
		static void prism_bbox(const std::array<Vector3, 12>&prism, Vector3 &min, Vector3& max);
		static void cubic_bbox(const std::array<Vector3, 8>&cubic, Vector3 &min, Vector3& max);
		static void three_facets_inter_point(const Vector3& a0, const Vector3& a1, const Vector3& a2, const Vector3& b0, const Vector3& b1,
			const Vector3& b2, const Vector3& c0, const Vector3& c1, const Vector3& c2,  Vector3& p);
		static void get_bb_corners_12(const std::array<Vector3,12> &vertices, Vector3 &min, Vector3 &max) {//TODO why use this one
			min = vertices[0];
			max = vertices[0];

			for (size_t j = 0; j < 12; j++) {
				for (int i = 0; i < 3; i++) {
					min[i] = std::min(min[i], vertices[j][i]);
					max[i] = std::max(max[i], vertices[j][i]);
				}
			}

			//const Scalar dis = (max - min).minCoeff() * 3;//TODO  change to 1e-5 or sth
			const Scalar dis = 1e-4;
			for (int j = 0; j < 3; j++) {
				min[j] -= dis;
				max[j] += dis;
			}

		}
		static void get_bb_corners_8(const std::array<Vector3, 8> &vertices, Vector3 &min, Vector3 &max) {//TODO why use this one
			min = vertices[0];
			max = vertices[0];

			for (size_t j = 0; j < 8; j++) {
				for (int i = 0; i < 3; i++) {
					min[i] = std::min(min[i], vertices[j][i]);
					max[i] = std::max(max[i], vertices[j][i]);
				}
			}

			//const Scalar dis = (max - min).minCoeff() * 3;//TODO  change to 1e-5 or sth
			const Scalar dis = 1e-4;
			for (int j = 0; j < 3; j++) {
				min[j] -= dis;
				max[j] += dis;
			}

		}
		static void  CornerList_prism(const std::vector<std::array<Vector3, 12>>& prism,
			std::vector<std::array<Vector3, 2>>& list) {

			list.resize(prism.size());//to be safer
			for (int i = 0; i < prism.size(); i++) {

				get_bb_corners_12(prism[i], list[i][0], list[i][1]);
				//prism_bbox(prism[i], list[i][0], list[i][1]);
			}
		}
		static void  CornerList_cubic(const std::vector<std::array<Vector3, 8>>& cubic,
			std::vector<std::array<Vector3, 2>>& list) {

			list.resize(cubic.size());//to be safer
			for (int i = 0; i < cubic.size(); i++) {

				get_bb_corners_8(cubic[i], list[i][0], list[i][1]);
				//cubic_bbox(cubic[i], list[i][0], list[i][1]);
			}
		}


		static void get_tri_corners(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2 , Vector3 &mint, Vector3 &maxt) {
			mint[0] = std::min(std::min(triangle0[0], triangle1[0]), triangle2[0]);
			mint[1] = std::min(std::min(triangle0[1], triangle1[1]), triangle2[1]);
			mint[2] = std::min(std::min(triangle0[2], triangle1[2]), triangle2[2]);
			maxt[0] = std::max(std::max(triangle0[0], triangle1[0]), triangle2[0]);
			maxt[1] = std::max(std::max(triangle0[1], triangle1[1]), triangle2[1]);
			maxt[2] = std::max(std::max(triangle0[2], triangle1[2]), triangle2[2]);

		}

		static bool box_box_intersection(const Vector3 &min1, const Vector3 &max1, const Vector3 &min2, const Vector3 &max2) {
			if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2]) return 0;
			if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2]) return 0;
			return 1;
		}


		// to check if a point is in the prisms. the jump index shows the prisms not counted in calculation, and jump is sorted from small to big
		bool point_out_prism(const Vector3& point, const std::vector<unsigned int>& prismindex, const int& jump) const;
		bool point_out_prism_rational(const Rational& point0, const Rational& point1, const Rational& point2, const std::vector<unsigned int>& prismindex, const int& jump) const;
		static int orient_3d_rational(const Rational& a11, const Rational& a12, const Rational& a13,
			const Vector3& a, const Vector3& b, const Vector3& c ) {
			Rational a21(a[0]), a22(a[1]), a23(a[2]), a31(b[0]), a32(b[1]), a33(b[2]), a41(c[0]), a42(c[1]), a43(c[2]);
			Rational number = a11 * a22*a33 - a12 * a23*a41 + a13 * a31*a42 - a21 * a32*a43 + a41 * a32*a23 - a42 * a33*a11
				+ a43 * a21*a12 - a31 * a22*a13 + a11 * a23*a42 - a13 * a32*a41 + a22 * a31*a43 - a12 * a21*a33
				+ a41 * a33*a12 - a43 * a22*a11 + a32 * a21*a13 - a42 * a31*a23 + a11 * a32*a43 - a22 * a33 * a41
				+ a12 * a23*a31 - a13 * a21*a42 + a41 * a22*a13 - a32 * a23*a11 + a42 * a33*a21 - a43 * a31*a12;
			return number.get_sign();
		}
		
		static void BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism,  std::vector<std::array<Vector3, 8>>& envbox, const Scalar& epsilon);
		static void seg_cube(const Vector3 &p1, const Vector3 &p2, const Scalar& width, std::array<Vector3, 8>& envbox);
		


		 struct DATA_LPI {
			 int segid;
			 int prismid;
			 int facetid;
			 int jump1;
		 };
		 
		 int Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_LPI& datalpi,const std::array<Vector3,3>&triangle, const std::vector<unsigned int>& prismindex)const;
		 
		 int Implicit_Seg_Facet_interpoint_Out_Prism_double(
			 const Scalar& a11, const Scalar&a12, const Scalar& a13, const Scalar& d, const Scalar& fa11,
			 const Scalar& fa12, const Scalar& fa13, const Scalar& max1, const Scalar&max2, const Scalar& max5,
			 const Vector3& segpoint0, const Vector3& segpoint1, const Vector3& triangle1,
			 const Vector3& triangle2, const Vector3& triangle3, const std::vector<unsigned int>& prismindex, const int& jump)const;



		struct DATA_TPI {
			int prismid1;
			int facetid1;
			int prismid2;
			int facetid2;
			int jump1;
			int jump2;
		};
		
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_TPI& datatpi, const std::array<Vector3, 3>&triangle, const std::vector<unsigned int>& prismindex)const;


		
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double(
			const Scalar& d, const Scalar& n1d, const Scalar& n2d, const Scalar& n3d,
			const Scalar& max1, const Scalar& max2, const Scalar& max3, const Scalar& max4, const Scalar& max5, const Scalar&max6, const Scalar& max7,
			const std::array<Vector3, 3>& triangle,
			const Vector3& facet10, const Vector3& facet11, const Vector3& facet12, const Vector3& facet20, const Vector3& facet21, const Vector3& facet22,
			const std::vector<unsigned int>& prismindex, const int& jump1, const int &jump2, const bool & multiflag, TPI_exact_suppvars& s) const;


		
		static bool is_3_triangle_cut_pure_multiprecision(const std::array<Vector3, 3>& triangle,  TPI_exact_suppvars& s);
		
		static bool is_3_triangle_cut_double(
			const Scalar &d, const Scalar & n1d, const Scalar &n2d, const Scalar & n3d,
			const Scalar & max1, const Scalar &max2, const Scalar &max3, const Scalar & max4, const Scalar & max5,
			const Scalar & max6, const Scalar &max7,
			const std::array<Vector3, 3>& triangle,
			const Vector3& facet10, const Vector3& facet11, const Vector3& facet12, const Vector3& facet20, const Vector3& facet21, const Vector3& facet2,
			bool & multiflag,
			TPI_exact_suppvars& s);


		static int is_3_triangle_cut_float_fast(
			const Vector3& tri0, const Vector3& tri1, const Vector3& tri2,
			const Vector3& facet10, const Vector3& facet11, const Vector3& facet12,
			const Vector3& facet20, const Vector3& facet21, const Vector3& facet22);
		//not accurate but conservative
		bool is_triangle_cut_prism(const int&pindex,
			const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, std::vector<int> &cid)const ;
		//not accurate but conservative
		bool is_triangle_cut_cube(const int&cindex,
			 const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, std::vector<int> &cid)const;
		//not accurate but conservative
		bool is_seg_cut_prism(const int&pindex,
			const Vector3& seg0, const Vector3& seg1, std::vector<int> &cid)const;
		//not accurate but conservative
		bool is_seg_cut_cube(const int&cindex,
			const Vector3& seg0, const Vector3& seg1, std::vector<int> &cid)const;

		static Vector2 to_2d(const Vector3 &p, int t) {
			return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
		}

		
		static int dot_sign(const Vector3& a, const Vector3& b) {
			Scalar t = a.dot(b);
			if (t > SCALAR_ZERO) return 1;
			if (t < -1*SCALAR_ZERO) return -1;
			Multiprecision a0(a[0]), a1(a[1]), a2(a[2]), b0(b[0]), b1(b[1]), b2(b[2]), dot;
			dot = a0 * b0 + a1 * b1 + a2 * b2;
			if (dot.get_sign() > 0) return 1;
			if (dot.get_sign() < 0) return -1;
			return 0;

		}
		static void to_geogram_mesh(const std::vector<Vector3>& V, const std::vector<Vector3i>& F, GEO::Mesh &M) {

			M.clear();

			// Setup vertices

			M.vertices.create_vertices(V.size());

			for (int i = 0; i < (int)M.vertices.nb(); ++i) {

				GEO::vec3 &p = M.vertices.point(i);

				p[0] = V[i][0];
				p[1] = V[i][1];
				p[2] = V[i][2];

			}

			// Setup faces


			M.facets.create_triangles(F.size());




			for (int c = 0; c < (int)M.facets.nb(); ++c) {

				for (int lv = 0; lv < 3; ++lv) {

					M.facets.set_vertex(c, lv, F[c][lv]);

				}

			}

		}
		static void from_geogram_mesh(const GEO::Mesh &M, std::vector<Vector3>& V, std::vector<Vector3i>& F) {

			V.resize(M.vertices.nb());

			for (int i = 0; i < (int)M.vertices.nb(); ++i) {

				GEO::vec3 p = M.vertices.point(i);

				V[i][0] = p[0];
				V[i][1] = p[1];
				V[i][2] = p[2];

			}

			assert(M.facets.are_simplices());

			F.resize(M.facets.nb());

			for (int c = 0; c < (int)M.facets.nb(); ++c) {

				for (int lv = 0; lv < 3; ++lv) {

					F[c][lv] = M.facets.vertex(c, lv);

				}

			}

			assert(M.cells.are_simplices());

		}



	};

}

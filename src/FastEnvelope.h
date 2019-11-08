#pragma once
#include <fastenvelope/Types.hpp>
#include <fastenvelope/AABB.h>

#include <fastenvelope/ip_filtered.h>


#include <fastenvelope/Rational.hpp>

#include <geogram/mesh/mesh.h>

#include <vector>
#include <array>
#include <fenv.h>

#include<iostream>


namespace fastEnvelope {

	class FastEnvelope
	{
	private:
		struct DATA_LPI
		{
			int segid;
			int prismid;
			int facetid;
			int jump1;
		};
		struct DATA_TPI
		{
			int prismid1;
			int facetid1;
			int prismid2;
			int facetid2;
			int jump1;
			int jump2;
		};

		struct INDEX
		{
			int Pi;
			std::vector<int> FACES;
		};

		static const int NOT_INTERSECTED = 2;
		static const int INTERSECTED = 1;
		static const int OUT_PRISM = 1;
		static const int IN_PRISM = 0;
		static const int CUT_COPLANAR = 4;
		static const int CUT_EMPTY = -1;
		static const int CUT_FACE = 3;
	public:
		static const int NOT_DEGENERATED = 0;
		static const int NERLY_DEGENERATED = 1;
		static const int DEGENERATED_SEGMENT = 2;
		static const int DEGENERATED_POINT = 3;


		//static const Scalar  BOX_SCALE = 1 / 10.0;
	public:
		FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps);
		static void printnumber();
		static void reset_time();
		//check if tri is outside
		bool is_outside(const std::array<Vector3, 3> &triangle) const;

		//move to another file with cpp
		static int seg_cut_plane(const Vector3 &seg0, const Vector3 &seg1, const Vector3 &t0, const Vector3 &t1, const Vector3 &t2);
		static Vector3 accurate_normal_vector(const Vector3 &p0, const Vector3 &p1, const Vector3 &q0, const Vector3 &q1);
		static int is_triangle_degenerated(const Vector3& triangle0, const Vector3& triangle1, const Vector3& triangle2);


		//sample for debugging
		bool sample_triangle_outside(const std::array<Vector3, 3> &triangle, const int &pieces) const;

		//export files for debugging
		//void print_prisms(const std::array<Vector3, 3> &triangle, const std::string &path) const;

	private:
		AABB tree;
		std::vector<std::array<Vector3, 2>> cornerlist;
		std::vector<std::vector<std::array<Vector3, 3>>> halfspace;
		



		//main pipeline
		bool debugcode(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex)const;
		bool FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int>& prismindex) const;
		bool is_two_facets_neighbouring(const int & pid, const int &i, const int &j)const;
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double_return_local_id(
			const Scalar &d, const Scalar &n1, const Scalar &n2, const Scalar &n3,
			const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5, const Scalar &max6, const Scalar &max7,
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const int &jump1, const int &jump2, const bool &multiflag,
			TPI_exact_suppvars &s, int &id) const;
		//algorithm
		int Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const;
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision_return_local_id(const DATA_TPI &datatpi,
			const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex, TPI_exact_suppvars& s, int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_double(
			const Scalar &a11, const Scalar &a12, const Scalar &a13, const Scalar &d, const Scalar &fa11,
			const Scalar &fa12, const Scalar &fa13, const Scalar &max1, const Scalar &max2, const Scalar &max5,
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle1,
			const Vector3 &triangle2, const Vector3 &triangle3, const std::vector<unsigned int> &prismindex, const int &jump) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_double_return_id(
			LPI_filtered_suppvars& s,
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_double_return_local_id(
			LPI_filtered_suppvars& sf,
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_check_id(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision_return_id(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle,
			const std::vector<unsigned int> &prismindex, int& id) const;

		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>>& faceorder, const int &jump, int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision_return_local_id(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle,
			const std::vector<unsigned int> &prismindex, int& id) const;
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_TPI &datatpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex, TPI_exact_suppvars& s) const;
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double_with_face_order(
			const Scalar &d, const Scalar &n1, const Scalar &n2, const Scalar &n3,
			const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5, const Scalar &max6, const Scalar &max7,
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::array<bool, 8>>intersect_face, const int &jump1, const int &jump2, const bool &multiflag,
			TPI_exact_suppvars &s) const;
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>>&intersect_face, const int &jump1, const int &jump2,
			int &id) const;
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double_return_id_with_face_order(
			const Scalar &d, const Scalar &n1, const Scalar &n2, const Scalar &n3,
			const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5, const Scalar &max6, const Scalar &max7,
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::array<bool, 8>>intersect_face, const int &jump1, const int &jump2, const bool &multiflag,
			TPI_exact_suppvars &s, int &id) const;
		static bool is_3_triangle_cut_pure_multiprecision(const std::array<Vector3, 3> &triangle, TPI_exact_suppvars &s);

		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double(
			const Scalar &d, const Scalar &n1d, const Scalar &n2d, const Scalar &n3d,
			const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5, const Scalar &max6, const Scalar &max7,
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const int &jump1, const int &jump2, const bool &multiflag, TPI_exact_suppvars &s) const;

		static bool is_3_triangle_cut(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12,
			const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22);
		bool is_tpp_on_polyhedra(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const int &prismid, const int &faceid)const;

		// to check if a point is in the prisms. the jump index shows the prisms not counted in calculation, and jump is sorted from small to big
		bool point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const;
		bool point_out_prism_return_id(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump, int &id)const;
		bool point_out_prism_return_id_list(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump, std::vector<int> &idlist) const;
		void halfspace_init(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
			std::vector<std::array<Vector3, 2>>& cornerlist, const Scalar &epsilon);
		void seg_cube(const Vector3 &p1, const Vector3 &p2, const Scalar &width, std::array<Vector3, 8> &envbox);
		static void resorting(const std::vector<Vector3> &V, const std::vector<Vector3i> &F, std::vector<Vector3i> &fnew);
		static bool box_box_intersection(const Vector3 &min1, const Vector3 &max1, const Vector3 &min2, const Vector3 &max2)
		{
			if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
				return 0;
			if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
				return 0;
			return 1;
		}







		//heuristics to refine the box box guess from the tree
		static int is_3_triangle_cut_float_fast(
			const Vector3& tri0, const Vector3& tri1, const Vector3& tri2,
			const Vector3& facet10, const Vector3& facet11, const Vector3& facet12,
			const Vector3& facet20, const Vector3& facet21, const Vector3& facet22);

		int is_triangle_cut_envelope_polyhedra(const int &cindex,
			const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, std::vector<int> &cid) const;
		bool is_seg_cut_polyhedra(const int &cindex,
			const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const;
	};

}
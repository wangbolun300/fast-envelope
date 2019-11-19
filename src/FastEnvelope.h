#pragma once
#include <fastenvelope/Types.hpp>
#include <fastenvelope/AABB.h>
#include <indirectpredicates/ip_filtered.h>

#include <vector>
#include <array>
#include <fenv.h>

#include<iostream>


namespace fastEnvelope {

	class FastEnvelope
	{
	
	public:
		FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps);
		static void printnumber();

		//check if tri or point is outside
		bool is_outside(const std::array<Vector3, 3> &triangle) const;
		bool is_outside(const Vector3 &point) const;
		
	private:

		struct INDEX
		{
			int Pi;
			std::vector<int> FACES;
		};

	private:
		AABB tree;
		std::vector<std::array<Vector3, 2>> cornerlist;
		std::vector<std::vector<std::array<Vector3, 3>>> halfspace;

		bool debugcode(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex)const;
		bool FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int>& prismindex) const;
		bool is_two_facets_neighbouring(const int & pid, const int &i, const int &j)const;
		
	
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const;
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>>& faceorder, const int &jump, int &id) const;
		
		int Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
			const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
			const std::vector<std::vector<int>>& intersect_face, const std::vector<bool>& coverlist, const int &jump, int &id) const;
		
		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>>&intersect_face, const std::vector<bool>& coverlist, const int &jump1, const int &jump2,
			int &id) const;

		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>>&intersect_face, const int &jump1, const int &jump2,
			int &id) const;

		static bool is_3_triangle_cut_pure_multiprecision(const std::array<Vector3, 3> &triangle, TPI_exact_suppvars &s);

		
		static bool is_3_triangle_cut(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12,
			const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22);
		bool is_tpp_on_polyhedra(
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const int &prismid, const int &faceid)const;

		bool point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const;
		bool point_out_prism_return_local_id(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump, int &id)const;

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
#pragma once
#include <fastenvelope/Types.hpp>
#include <fastenvelope/AABB.h>

#include <fastenvelope/ip_filtered.h>

#include <fastenvelope/Multiprecision.hpp>
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

		static const int NOT_DEGENERATED = 0;
		static const int NERLY_DEGENERATED = 1;
		static const int DEGENERATED_SEGMENT = 2;
		static const int DEGENERATED_POINT = 3;


		//static const Scalar  BOX_SCALE = 1 / 10.0;
	public:
		FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps);

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

		//unify


		
		std::vector<std::vector<std::array<Vector3, 3>>> halfspace;

		//main pipeline
		bool FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int>& prismindex) const;
		bool is_two_facets_neighbouring(const int & pid, const int &i, const int &j)const;

		//algorithm
		int Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const;

		int Implicit_Seg_Facet_interpoint_Out_Prism_double(
			const Scalar &a11, const Scalar &a12, const Scalar &a13, const Scalar &d, const Scalar &fa11,
			const Scalar &fa12, const Scalar &fa13, const Scalar &max1, const Scalar &max2, const Scalar &max5,
			const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle1,
			const Vector3 &triangle2, const Vector3 &triangle3, const std::vector<unsigned int> &prismindex, const int &jump) const;

		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_TPI &datatpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const;
		static bool is_3_triangle_cut_pure_multiprecision(const std::array<Vector3, 3> &triangle, TPI_exact_suppvars &s);

		int Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double(
			const Scalar &d, const Scalar &n1d, const Scalar &n2d, const Scalar &n3d,
			const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5, const Scalar &max6, const Scalar &max7,
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
			const std::vector<unsigned int> &prismindex, const int &jump1, const int &jump2, const bool &multiflag, TPI_exact_suppvars &s) const;
		static bool is_3_triangle_cut_double(
			const Scalar &d, const Scalar &n1d, const Scalar &n2d, const Scalar &n3d,
			const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5,
			const Scalar &max6, const Scalar &max7,
			const std::array<Vector3, 3> &triangle,
			const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet2,
			bool &multiflag,
			TPI_exact_suppvars &s);

		// to check if a point is in the prisms. the jump index shows the prisms not counted in calculation, and jump is sorted from small to big
		bool point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const;

		//maybe unify
		static void halfspace_generation(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, std::vector<std::array<Vector3, 12>> &envprism, std::vector<std::array<Vector3, 8>> &envbox, const Scalar &epsilon);
		/*template<typename T>*/
		//static void halfspace_generation(const std::vector<Vector3> &m_ver, const std::vector<std::vector<int>> &m_faces, const Scalar &epsilon, std::vector<std::vector<std::array<Vector3, 3>>> &hs);
		//generate boxes for segments
		static void seg_cube(const Vector3 &p1, const Vector3 &p2, const Scalar &width, std::array<Vector3, 8> &envbox);


		//build boxes for set of half planes, change maybe to unify
		static void get_bb_corners_12(const std::array<Vector3, 12> &vertices, Vector3 &min, Vector3 &max) {//TODO why use this one
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


		//heuristics to refine the box box guess from the tree
		static int is_3_triangle_cut_float_fast(
			const Vector3& tri0, const Vector3& tri1, const Vector3& tri2,
			const Vector3& facet10, const Vector3& facet11, const Vector3& facet12,
			const Vector3& facet20, const Vector3& facet21, const Vector3& facet22);
		//not accurate but conservative
		//bool is_triangle_cut_prism(const int&pindex,
		//	const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, std::vector<int> &cid)const ;
		////not accurate but conservative
		//bool is_triangle_cut_cube(const int&cindex,
		//	 const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, std::vector<int> &cid)const;
		////not accurate but conservative
		//bool is_seg_cut_prism(const int&pindex,
		//	const Vector3& seg0, const Vector3& seg1, std::vector<int> &cid)const;
		////not accurate but conservative
		//bool is_seg_cut_cube(const int&cindex,
		//	const Vector3& seg0, const Vector3& seg1, std::vector<int> &cid)const;
		bool is_triangle_cut_envelope_polyhedra(const int &cindex,
			const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, std::vector<int> &cid) const;
		bool is_seg_cut_polyhedra(const int &cindex,
			const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const;
	};

}

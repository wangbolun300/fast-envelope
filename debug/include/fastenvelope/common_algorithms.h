#pragma once
#include<fastenvelope/Types.hpp>
#include<vector>
#include<array>

namespace fastEnvelope {
	namespace algorithms {
		 int seg_cut_plane(const Vector3 &seg0, const Vector3 &seg1, const Vector3 &t0, const Vector3 &t1, const Vector3 &t2);
		
		 int is_triangle_degenerated(const Vector3& triangle0, const Vector3& triangle1, const Vector3& triangle2);
		 Vector3 accurate_normal_vector(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2);
		 void get_tri_corners(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2, Vector3 &mint, Vector3 &maxt);
		 bool box_box_intersection(const Vector3 &min1, const Vector3 &max1, const Vector3 &min2, const Vector3 &max2);//TDOO;
		 Vector2 to_2d(const Vector3 &p, int t);
		 //resort the facets using Morton's code
		 void resorting(const std::vector<Vector3> &V, const std::vector<Vector3i> &F, std::vector<Vector3i> &fnew);
		 void resorting(const std::vector<Vector3> &V, const std::vector<Vector3i> &F, std::vector<Vector3i> &fnew, std::vector<int>& new2old);
		 void seg_cube(const Vector3 &p1, const Vector3 &p2, const Scalar &width, std::array<Vector3, 8> &envbox);
		 
		 void halfspace_generation(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
			 std::vector<std::array<Vector3, 2>>& cornerlist, const Scalar &epsilon);
		 void halfspace_generation(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
			 std::vector<std::array<Vector3, 2>>& cornerlist, const std::vector<Scalar>& epsilons);
	}
}

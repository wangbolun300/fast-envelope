#pragma once

#include <vector>
#include <array>
#include<fastenvelope/Types.hpp>
#include<iostream>

namespace fastEnvelope {
	class obb
	{
	public:
		
		obb(const std::vector<std::vector<Vector3>>& envelope_vertices, const int p_face[8][3], const int c_face[6][3],
			const std::array<std::vector<int>, 8>& p_facepoint, const std::array<std::array<int, 4>, 6>& c_facepoint);
		//build obb out of a set of points, transformation matrix Trans maps the box into a box which corners are
		//[1,1,1] and [-1,-1,-1]
		//input: points
		//output: transformation matrix and it's inversion 
		static void build_obb(const std::vector<Vector3>& points,const Vector3& normal, const Scalar offset, 
			Eigen::Matrix4d &M, Eigen::Matrix4d &invM);

		static bool obb_intersection(const Eigen::Matrix4d &M1, const Eigen::Matrix4d &invM1, const Eigen::Matrix4d &M2, const Eigen::Matrix4d &invM2);
		static void test();

	private:
		std::vector<std::vector<Eigen::Matrix4d>> Trans;
		std::vector<std::vector<Eigen::Matrix4d>> invTrans;
	};
}
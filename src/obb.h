#pragma once

#include <vector>
#include <array>
#include<fastenvelope/Types.hpp>
#include<iostream>


// this flag DO_NOT_HAVE_DEGENERATED_FACES is to make code faster when it is sure we dont have any three points that
//define one halfspace is degenerated to a segment, or a point
#define DO_NOT_HAVE_DEGENERATED_FACES
namespace fastEnvelope {
	class obb
	{
	private:
		 
			Eigen::Matrix4d Trans;
			Eigen::Matrix4d invTrans;
		
	public:

		static std::vector<obb> build_obb_matrixs(const std::vector<Vector3>& prism_vertices, const int p_face[8][3], const int c_face[6][3],
			const std::array<std::vector<int>, 8>& p_facepoint, const std::array<std::array<int, 4>, 6>& c_facepoint);
		
		static obb build_triangle_obb_matrixs(const Vector3&t0, const Vector3&t1, const Vector3&t2);// can not take degenerated triangle
		bool intersects(const obb& M2)const;
		bool intersects(const obb& M2, const obb& M3)const;
		static void test();

	private:
		static const int polyhedron_point_number1 = 12;
		static const int polyhedron_point_number2 = 8;
		static const int polyhedron_face_number1 = 8;
		static const int polyhedron_face_number2 = 6;
		static const double OBB_OFFSET = 1e-4;
		
		
		
		//build obb out of a set of points, transformation matrix Trans maps the box into a box which corners are
		//[1,1,1] and [-1,-1,-1]
		//input: points, normal vector of these points(respect to the minimal eigenvalue of PCA matrix)
		//output: transformation matrix and it's inversion 
		static void build_obb(const std::vector<Vector3>& points, const Vector3& normal, const Scalar offset,
			Eigen::Matrix4d &M, Eigen::Matrix4d &invM);

		static bool obb_intersection(const Eigen::Matrix4d &M1, const Eigen::Matrix4d &invM1, const Eigen::Matrix4d &M2, const Eigen::Matrix4d &invM2);

	};
}
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
	public:
		// this function is to initialize obb for every facets of polyhedron
		void obb_init(const std::vector<std::vector<Vector3>>& envelope_vertices, const int p_face[8][3], const int c_face[6][3],
			const std::array<std::vector<int>, 8>& p_facepoint, const std::array<std::array<int, 4>, 6>& c_facepoint);
		bool intersected(const int prismid1, const int faceid1, const int prismid2, const int faceid2)const;
		static void test();

	private:
		int polyhedron_point_number1 = 12;
		int polyhedron_point_number2 = 8;
		int polyhedron_face_number1 = 8;
		int polyhedron_face_number2 = 6;
		int OBB_OFFSET = 1e-4;

		std::vector<std::vector<Eigen::Matrix4d>> Trans;
		std::vector<std::vector<Eigen::Matrix4d>> invTrans;
		//build obb out of a set of points, transformation matrix Trans maps the box into a box which corners are
		//[1,1,1] and [-1,-1,-1]
		//input: points, normal vector of these points(respect to the minimal eigenvalue of PCA matrix)
		//output: transformation matrix and it's inversion 
		static void build_obb(const std::vector<Vector3>& points, const Vector3& normal, const Scalar offset,
			Eigen::Matrix4d &M, Eigen::Matrix4d &invM);

		static bool obb_intersection(const Eigen::Matrix4d &M1, const Eigen::Matrix4d &invM1, const Eigen::Matrix4d &M2, const Eigen::Matrix4d &invM2);

	};
}
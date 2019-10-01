#pragma once

#include <fastenvelope/Types.hpp>
#include<vector>
#include<fastenvelope/Parameters.h>
namespace fastEnvelope {

	class Envelop
	{
	public:
		// generate octree, recording the cells position, input vertices list, cell list, the depth of iteration,
		// out put an orctree
		//static void Octree(const std::vector<Vector3>& vlist, const std::vector<Vector3i>& cell, const int& depth, 
		//	std::vector<std::vector<int>>& tree);
		
		// to judge if two faces have intersection,  in order to 
		//static void FaceIntersection(const std::vector<Vector3>& vlist, const std::vector<Vector3i>& cell, const int& depth,
		//	std::vector<std::vector<int>>& tree);

		//to judge if a triangle is in the envelope. input 3 vertices of the triangle, the potential intersection prisms,
		//out put a bool value to tell if the envelope is out of the envelope
		static bool EnvelopeTest(const std::array<Vector3, 3>& face, const std::vector<std::array<Vector3, 12>>& interprism);
		static bool EnvelopeTest1(const std::array<Vector3, 3>& face, const std::vector<std::array<Vector3, 12>>& interprism, const int& mark);
		static bool EnvelopeTest2(const std::array<Vector3, 3>& face, const std::vector<std::array<Vector3, 12>>& interprism);
		//to cut a triangle with a prism.
		//input a single prism, a triangle
		//output a outer triangle list
		static void OnePrismCutOneTriangle(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 12>& envprism
			, std::vector< std::array<Vector3, 3>>& olist);
		// to cut a triangle with a prism triangle face. 
		//input triangle vertices, prism triangle face, output if the triangle is intersected, and a inner triangles list and a outer triangle list
		static void PrismTriangleCutTriangle(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& p_face, bool& intersection,
			std::vector<std::array<Vector3, 3>>& innerlist, std::vector<std::array<Vector3, 3>>& outerlist);
		// to generate envelope prisms accroding to the model faces. 
		//input model vertices, model faces indexs,and bounding box diagonal.
		//out put prism vertices
		 static void PrismGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3,6>>& envprism, const Scalar& bbd);
		 
		 
		 //surface and line intersection points, input surface vertices and line vertices,out put if
		 //the line is cutted and the cut point,and the bool if the line go through the surface
		 //THIS VERSION DO NOT DO THE COLLASPING
		 static void SLIntersection(const Parameters &params, const std::array<Vector3, 3>& cutface, const Vector3& linepoints0, const Vector3& linepoints1,
			 int & cutOrnot, Vector3& interp, Scalar &t, int& ion);

		 // this is for back up of the SLIntersection, it do the collasping while calculation of intersection points
		 static void SLIntersectionBackup(const Parameters & params, const std::array<Vector3, 3>& cutface, const Vector3 & linepoints0,
			 const Vector3 & linepoints1, int & cutOrnot, Vector3 & interp, Scalar & t, int & ion);

		 //Minkowski sweeping. not considering face orientation
		 static void BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, const Scalar& bbd);
	};


		
}
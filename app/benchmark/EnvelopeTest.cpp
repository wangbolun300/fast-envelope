#include "EnvelopeTest.h"
#include<fastenvelope/Parameters.h>
#include<fastenvelope/Predicates.hpp>
namespace fastEnvelope {
	
	// EnvelopeTest = 1 means out of the envelope, =0 means in the envelope
	// this version still have degeneration problem. 
	bool Envelop::EnvelopeTest(const std::array<Vector3, 3>& face, const std::vector<std::array<Vector3, 12>>& interprism)
	{
		if (interprism.size() == 0) {
			return 1;
		}
		std::vector< std::array<Vector3, 3>> olist;
		std::vector< std::array<Vector3, 3>> testfacelist;
		testfacelist.push_back(face);
		for (int i = 0; i < interprism.size(); i++) {
			olist.clear();
			for (int j = 0; j < testfacelist.size(); j++) {
				OnePrismCutOneTriangle(testfacelist[j], interprism[i], olist);
			}
			if (olist.size() == 0) {
				return 0;
			}

			testfacelist = olist;
		}
		//std::cout << "inner prisms size " << interprism.size() << std::endl;
		//std::cout << "outer triangles size " << olist.size() << std::endl;
		/*
		for (int i = 0; i < olist.size(); i++) {
			std::cout << "outer triangles list " << olist[i][0]<<","<< olist[i][1] <<","<< olist[i][2] << std::endl;
		}
		*/
		for (int i = 0; i < olist.size(); i++) {
			if ((olist[i][0] - olist[i][1]).cross(olist[i][0] - olist[i][2]).norm()/2 > SCALAR_ZERO) {
				return 1;
			}
		}

		return 0;

	}
	
	bool Envelop::EnvelopeTest1(const std::array<Vector3, 3>& face, const std::vector<std::array<Vector3, 12>>& interprism,const int& mark)
	{
		if (interprism.size() == 0) {
			return 1;
		}
		std::vector< std::array<Vector3, 3>> olist;
		std::vector< std::array<Vector3, 3>> testfacelist;
		testfacelist.push_back(face);
		for (int i = 0; i < interprism.size(); i++) {
			olist.clear();
			for (int j = 0; j < testfacelist.size(); j++) {
				OnePrismCutOneTriangle(testfacelist[j], interprism[i], olist);
			}
			if (olist.size() == 0) {
				return 0;
			}

			testfacelist = olist;
		}
		//std::cout << "inner prisms size " << interprism.size() << std::endl;
		//std::cout << "outer triangles size " << olist.size() << std::endl;
		/*
		for (int i = 0; i < olist.size(); i++) {
			std::cout << "outer triangles list " << olist[i][0]<<","<< olist[i][1] <<","<< olist[i][2] << std::endl;
		}
		*/
		for (int i = 0; i < olist.size(); i++) {
			if ((olist[i][0] - olist[i][1]).cross(olist[i][0] - olist[i][2]).norm() / 2 > SCALAR_ZERO) {
				return 1;
			}
		}


		if (mark == 973) {
				std::cout <<"olist situation: "<< olist.size() << std::endl;
				std::cout << "olist triangle: \n" << olist[0][0]<<"\n"<<olist[0][1] << "\n" << olist[0][2] << std::endl;
			
		
		
		}
		//std::cout << "inner prisms size " << interprism.size() << std::endl;
		//std::cout << "outer triangles size " << olist.size() << std::endl;
		/*
		for (int i = 0; i < olist.size(); i++) {
			std::cout << "outer triangles list " << olist[i][0]<<","<< olist[i][1] <<","<< olist[i][2] << std::endl;
		}
		*/
		return 0;

	}
	
	bool Envelop::EnvelopeTest2(const std::array<Vector3, 3>& face, const std::vector<std::array<Vector3, 12>>& interprism) {
		if (interprism.size() == 0) {
			return 1;
		}
		std::vector< std::array<Vector3, 3>> olist;
		std::vector< std::array<Vector3, 3>> testfacelist;
		testfacelist.push_back(face);
		for (int i = 0; i < interprism.size(); i++) {
			olist.clear();
			for (int j = 0; j < testfacelist.size(); j++) {
				OnePrismCutOneTriangle(testfacelist[j], interprism[i], olist);
			}
			if (olist.size() == 0) {
				return 0;
			}

			testfacelist = olist;
		}
		//std::cout << "inner prisms size " << interprism.size() << std::endl;
		//std::cout << "outer triangles size " << olist.size() << std::endl;
		/*
		for (int i = 0; i < olist.size(); i++) {
			std::cout << "outer triangles list " << olist[i][0]<<","<< olist[i][1] <<","<< olist[i][2] << std::endl;
		}
		*/
		bool r;
		for (int i = 0; i < olist.size(); i++) {
			if ((olist[i][0] - olist[i][1]).cross(olist[i][0] - olist[i][2]).norm() / 2 < SCALAR_ZERO) {
				continue;
			}
			r = EnvelopeTest(olist[i], interprism);
			if (r == 1) {
				//std::cout << "really out side in test 2" << std::endl;
				return 1;
			}
		}
		//std::cout << " triangles inside test 2, or out triangle degenerated. and the outlist: " <<olist.size()<< std::endl;
		return 0;

	
	}
	void Envelop::OnePrismCutOneTriangle(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 12>& envprism
		, std::vector< std::array<Vector3, 3>>& olist)
	{
		//int p_face[5][3] = { {0,1,2},{3,5,4},{0,3,4},{0,2,5},{1,4,5} };//prism triangle index. all with orientation. this is for triangle prism
		int p_face[8][3] = { {0,1,2},{8,7,6},{1,0,7},{2,1,7},{3,2,8},{3,9,10},{5,4,11},{0,5,6} };//prism triangle index. all with orientation.
		std::vector<std::array<Vector3, 3>> trianglelist;
		trianglelist.push_back(triangle);
		std::array<Vector3, 3> tface;
		bool intersection;
		std::vector<std::array<Vector3, 3>> innerlist;
		//olist.resize(0);// in this way not be polluted by input
		for (int i = 0; i < 8; i++) {
			tface = { envprism[p_face[i][0]],envprism[p_face[i][1]],envprism[p_face[i][2]] };
			innerlist.clear();
			for (int j = 0; j < trianglelist.size(); j++) {

				PrismTriangleCutTriangle(trianglelist[j], tface, intersection, innerlist, olist);
				//std::cout << "i,j :" << i << "," << j << "\n olist size: " << olist.size()<<std::endl;
				//std::cout <<"\n olist size: " << olist.size() << std::endl;
			}
			trianglelist = innerlist;
		}
	}
	/*
	void Envelop::PrismTriangleCutTriangle(const std::array<Vector3, 3>& triangle,
		const std::array<Vector3, 3>& p_face, bool &intersection, std::vector<std::array<Vector3, 3>>& innerlist, std::vector<std::array<Vector3, 3>>& outerlist)
	{
		
		Vector3 A = p_face[0], B = p_face[1], C = p_face[2]; //ctr,//centroid of triangle
			//temp,//to judge if the point is in which side of the prism triangle
			//dist,//distance of the point to prism face
			//normal = (B - A).cross(C - B).normalized();
		int pointin[3], grt0=0, less0=0, eq0=0, ori;
		//std::array<Vector3, 3> vectorsnorm;
		//std::array<Vector3, 3> vectors;
		std::array<Vector3, 2> intersectionp;
		int edgepair[3][2] = { {1,2},{0,2},{0,1} };
		//ctr = (A + B + C)/3;
		for (int i = 0; i < 3; i++) {
			//vectors[i] = triangle[i] - ctr;
			//vectorsnorm[i] = vectors[i].normalized();
			//temp[i] = vectorsnorm[i].dot(normal);
			//dist[i] = abs(vectors[i].dot(normal));
			//std::cout << "temp : \n" << temp[i] << std::endl;
			//if (temp[i] > 0&&dist[i]> SCALAR_ZERO) {//TODO the parameters. TODO this version collapse points by distance
			//if (temp[i] > SCALAR_ZERO) {//TODO the parameters, and this version collapse points by angles.
			ori = Predicates::orient_3d(p_face[0], p_face[1], p_face[2], triangle[i]);
			if(ori==-1){// the point is outside
				pointin[i] =1;// to find which side is the ith triangle point in
				grt0 = grt0 + 1;
			}
			else {
				//if (temp[i] < 0&& dist[i]>SCALAR_ZERO) {
				//if (temp[i] < (-1)*SCALAR_ZERO) {
				if(ori==1){//the point is inside
					pointin[i] = -1;
					less0 = less0 + 1;
				}
				else {
					pointin[i] = 0;
					eq0 = eq0 + 1;
				}
			}
			
		}
		//std::cout << "temp[i]\n" << temp[0] << "\n" << temp[1] << "\n" << temp[2] << std::endl;
		//std::cout << "grt0,less0 and eq0\n" << grt0 << "\n" << less0 << "\n" << eq0 << std::endl;
		if (grt0 + eq0 == 3 && eq0 < 3) {// if eq0==3, it's inside the prism surface
			intersection = 0;
			outerlist.push_back(triangle);
			return;
		}
		else {
			if (less0 + eq0 == 3) {// if eq0==3, it's inside the prism surface
				intersection = 0;
				innerlist.push_back(triangle);
				//std::cout << "hhhhhhhhhhhhhhhhhhhhhhh1" << std::endl;
				return;
				
			}
			else {
				//TODO for now the triangle is cutted no matter if intersected, in this case 
				//square-triangle intersection becomes tirangle-triangle intersection.
				intersection = 1;
				if (eq0 > 0) {// one vertex on the plane
					int i;
					for (i = 0; i < 3; i++) {
						if (pointin[i] == 0) {
							intersectionp[0] = triangle[i];
							break;
						}
					}
					Parameters  params;
					Vector3 linepoints0 = triangle[edgepair[i][0]];
					Vector3 linepoints1 = triangle[edgepair[i][1]];
					int  cutOrnot;
					Vector3  interp;
					Scalar t;
					int  ion;
					SLIntersection(params, p_face, linepoints0,linepoints1, cutOrnot, interp, t, ion);// TODO may not be accurate.
					intersectionp[1] = interp;
					
	if (cutOrnot < 0) {
		std::cout << "Intersection Wrong in PrismTriangleCutTriangle" << std::endl;
		return;
	}
	else {
		if (pointin[edgepair[i][0]] > 0) {
			outerlist.push_back({ intersectionp[0],intersectionp[1],triangle[edgepair[i][0]] });
			innerlist.push_back({ intersectionp[0],intersectionp[1],triangle[edgepair[i][1]] });
		}
		else {
			outerlist.push_back({ intersectionp[0],intersectionp[1],triangle[edgepair[i][1]] });
			innerlist.push_back({ intersectionp[0],intersectionp[1],triangle[edgepair[i][0]] });
		}
	}
	return;
}// one vertex on the plane
				else {//no vertex on the plane
				if (grt0 == 1) {// 1 triangle point outside the plane
					int i;
					for (i = 0; i < 3; i++) {
						if (pointin[i] > 0) {
							break;
						}
					}
					Parameters  params;
					Vector3 linepoints0 = triangle[edgepair[i][0]];
					Vector3 linepoints1 = triangle[i];
					int  cutOrnot;
					Vector3  interp;
					Scalar t;
					int  ion;
					SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
					if (cutOrnot < 0) {
						std::cout << "Intersection Wrong in PrismTriangleCutTriangle" << std::endl;
						return;
					}
					intersectionp[0] = interp;

					linepoints0 = triangle[edgepair[i][1]];
					SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
					if (cutOrnot < 0) {
						std::cout << "Intersection Wrong in PrismTriangleCutTriangle" << std::endl;
						return;
					}
					intersectionp[1] = interp;
					outerlist.push_back({ triangle[i],intersectionp[0],intersectionp[1] });
					innerlist.push_back({ intersectionp[0],intersectionp[1] ,triangle[edgepair[i][0]] });
					innerlist.push_back({ triangle[edgepair[i][1]],intersectionp[1] ,triangle[edgepair[i][0]] });
				}// 1 triangle point outside the plane
				else {//two triangle points outside
					int i;
					for (i = 0; i < 3; i++) {
						if (pointin[i] < 0) {
							break;
						}
					}
					Parameters  params;
					Vector3 linepoints0 = triangle[edgepair[i][0]];
					Vector3 linepoints1 = triangle[i];
					int  cutOrnot;
					Vector3  interp;
					Scalar t;
					int  ion;
					SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
					if (cutOrnot < 0) {
						std::cout << "Intersection Wrong in PrismTriangleCutTriangle" << std::endl;
						return;
					}
					intersectionp[0] = interp;

					linepoints0 = triangle[edgepair[i][1]];
					SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
					if (cutOrnot < 0) {
						std::cout << "Intersection Wrong in PrismTriangleCutTriangle" << std::endl;
						return;
					}
					intersectionp[1] = interp;
					innerlist.push_back({ triangle[i],intersectionp[0],intersectionp[1] });
					outerlist.push_back({ intersectionp[0],intersectionp[1] ,triangle[edgepair[i][0]] });
					outerlist.push_back({ triangle[edgepair[i][1]],intersectionp[1] ,triangle[edgepair[i][0]] });
					//TODO generate 2 triangles from 4 points may have more than one solution

				}//two triangle points outside
				}//no vertex on the plane	
			}
		}

		//std::cout << "normal : \n" << normal <<"\n ctr:\n"<<ctr<< "\n vectors0\n" << vectors[0] << "\n vectors1\n" << vectors[1] << std::endl;
		//std::cout << "triangle in or not : \n" << grt0<<"\n"<<less0 <<"\n"<<eq0<< std::endl;
	}
	*/
	void Envelop::PrismTriangleCutTriangle(const std::array<Vector3, 3>& triangle,
		const std::array<Vector3, 3>& p_face, bool &intersection, std::vector<std::array<Vector3, 3>>& innerlist, std::vector<std::array<Vector3, 3>>& outerlist)
	{
		
		Vector3 A = p_face[0], B = p_face[1], C = p_face[2]; //ctr,//centroid of triangle
			//temp,//to judge if the point is in which side of the prism triangle
			//dist,//distance of the point to prism face
			//normal = (B - A).cross(C - B).normalized();
		int pointin[3], grt0=0, less0=0, eq0=0, ori;
		//std::array<Vector3, 3> vectorsnorm;
		//std::array<Vector3, 3> vectors;
		std::array<Vector3, 2> intersectionp;
		int edgepair[3][2] = { {1,2},{0,2},{0,1} };
		//ctr = (A + B + C)/3;
		for (int i = 0; i < 3; i++) {
			//vectors[i] = triangle[i] - ctr;
			//vectorsnorm[i] = vectors[i].normalized();
			//temp[i] = vectorsnorm[i].dot(normal);
			//dist[i] = abs(vectors[i].dot(normal));
			//std::cout << "temp : \n" << temp[i] << std::endl;
			//if (temp[i] > 0&&dist[i]> SCALAR_ZERO) {//TODO the parameters. TODO this version collapse points by distance
			//if (temp[i] > SCALAR_ZERO) {//TODO the parameters, and this version collapse points by angles.
			ori = Predicates::orient_3d(p_face[0], p_face[1], p_face[2], triangle[i]);
			if(ori==-1){// the point is outside
				pointin[i] =1;// to find which side is the ith triangle point in
				grt0 = grt0 + 1;
			}
			else {
				//if (temp[i] < 0&& dist[i]>SCALAR_ZERO) {
				//if (temp[i] < (-1)*SCALAR_ZERO) {
				if(ori==1){//the point is inside
					pointin[i] = -1;
					less0 = less0 + 1;
				}
				else {
					pointin[i] = 0;
					eq0 = eq0 + 1;
					//std::cout << "find a point on surface" << std::endl;
				}
			}
			
		}
		//std::cout << "temp[i]\n" << temp[0] << "\n" << temp[1] << "\n" << temp[2] << std::endl;
		//std::cout << "grt0,less0 and eq0\n" << grt0 << "\n" << less0 << "\n" << eq0 << std::endl;
		if (grt0 + eq0 == 3 && eq0 < 3) {// if eq0==3, it's inside the prism surface
			intersection = 0;
			outerlist.push_back(triangle);
			return;
		}
		else {
			if (less0 + eq0 == 3) {// if eq0==3, it's inside the prism surface
				intersection = 0;
				innerlist.push_back(triangle);
				//std::cout << "hhhhhhhhhhhhhhhhhhhhhhh1" << std::endl;
				return;
				
			}
			else {
				//TODO for now the triangle is cutted no matter if intersected, in this case 
				//square-triangle intersection becomes tirangle-triangle intersection.
				intersection = 1;
				if (eq0 > 0) {// one vertex on the plane
					int i;
					for (i = 0; i < 3; i++) {
						if (pointin[i] == 0) {
							intersectionp[0] = triangle[i];
							break;
						}
					}
					Parameters  params;
					Vector3 linepoints0 = triangle[edgepair[i][0]];
					Vector3 linepoints1 = triangle[edgepair[i][1]];
					int  cutOrnot;
					Vector3  interp;
					Scalar t;
					int  ion;
					SLIntersection(params, p_face, linepoints0,linepoints1, cutOrnot, interp, t, ion);// TODO may not be accurate.
					intersectionp[1] = interp;
					/*
					
					*/
					
					
						if (pointin[edgepair[i][0]] > 0) {
							outerlist.push_back({intersectionp[0],intersectionp[1],triangle[edgepair[i][0]]});
							innerlist.push_back({ intersectionp[0],intersectionp[1],triangle[edgepair[i][1]] });
						}
						else {
							outerlist.push_back({ intersectionp[0],intersectionp[1],triangle[edgepair[i][1]] });
							innerlist.push_back({ intersectionp[0],intersectionp[1],triangle[edgepair[i][0]] });
						}
					
					return;
				}// one vertex on the plane
				else {//no vertex on the plane
					if (grt0 == 1) {// 1 triangle point outside the plane
						int i;
						for (i = 0; i < 3; i++) {
							if (pointin[i] > 0) {
								break;
							}
						}
						Parameters  params;
						Vector3 linepoints0 = triangle[edgepair[i][0]];
						Vector3 linepoints1 = triangle[i];
						int  cutOrnot;
						Vector3  interp;
						Scalar t;
						int  ion;
						SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
					
						intersectionp[0] = interp;
						
						linepoints0= triangle[edgepair[i][1]];
						SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
					
						intersectionp[1] = interp;
						outerlist.push_back({ triangle[i],intersectionp[0],intersectionp[1] });
						innerlist.push_back({ intersectionp[0],intersectionp[1] ,triangle[edgepair[i][0]] });
						innerlist.push_back({ triangle[edgepair[i][1]],intersectionp[1] ,triangle[edgepair[i][0]] });
					}// 1 triangle point outside the plane
					else {//two triangle points outside
						int i;
						for (i = 0; i < 3; i++) {
							if (pointin[i] < 0) {
								break;
							}
						}
						Parameters  params;
						Vector3 linepoints0 = triangle[edgepair[i][0]];
						Vector3 linepoints1 = triangle[i];
						int  cutOrnot;
						Vector3  interp;
						Scalar t;
						int  ion;
						SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
						
						intersectionp[0] = interp;

						linepoints0 = triangle[edgepair[i][1]];
						SLIntersection(params, p_face, linepoints0, linepoints1, cutOrnot, interp, t, ion);
						
						intersectionp[1] = interp;
						innerlist.push_back({ triangle[i],intersectionp[0],intersectionp[1] });
						outerlist.push_back({ intersectionp[0],intersectionp[1] ,triangle[edgepair[i][0]] });
						outerlist.push_back({ triangle[edgepair[i][1]],intersectionp[1] ,triangle[edgepair[i][0]] });
						//TODO generate 2 triangles from 4 points may have more than one solution
						
					}//two triangle points outside
				}//no vertex on the plane	
			}
		}
		
		//std::cout << "normal : \n" << normal <<"\n ctr:\n"<<ctr<< "\n vectors0\n" << vectors[0] << "\n vectors1\n" << vectors[1] << std::endl;
		//std::cout << "triangle in or not : \n" << grt0<<"\n"<<less0 <<"\n"<<eq0<< std::endl;
	}
	
	
	void Envelop::PrismGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 6>>& envprism, const Scalar& bbd)
	{

		Vector3 AB, AC, BC, normal, ABn,ACn;
		Parameters pram;
		std::array<Vector3, 3> striangle;
		Scalar sinA, sinB, sinC, a, b, c, length, tolerance = bbd * pram.eps_rel,crsp;
		//envprism.resize(m_faces.size());// to make the results safer
		for (int i = 0; i < m_faces.size(); i++) {
			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			a = BC.norm();
			b = AC.norm();
			c = AB.norm();
			ABn = AB.normalized();
			ACn = AC.normalized();
			crsp = ABn.dot(ACn);
			normal = AB.cross(AC).normalized();
			if (abs(crsp)>1- SCALAR_ZERO) {//TODO how to judge degeneration? the sum of two edges less than third one?
				
				std::cout << "Degenerated prism: NO. " << i << std::endl;
				continue;
			}
			length = tolerance / sqrt(1 - pow((b*b + c * c - a * a) / (2 * c*b), 2));
			striangle[0] = m_ver[m_faces[i][0]] + length * ((-1)*AC / b + (-1)*AB / c);
			length = tolerance / sqrt(1 - pow((a * a + c * c - b * b) / (2 * a*c), 2));
			striangle[1] = m_ver[m_faces[i][1]] + length * (AB / c + (-1)*BC / a);
			length = tolerance / sqrt(1 - pow((a *a + b * b - c * c) / (2 * a*b), 2));
			striangle[2] = m_ver[m_faces[i][2]] + length * (AC / b + BC / a);
			envprism.push_back({ striangle[0] + tolerance * normal,striangle[1] + tolerance * normal,striangle[2] + tolerance * normal,
			striangle[0] - tolerance * normal, striangle[1] - tolerance * normal, striangle[2] - tolerance * normal });


			
		}
		//std::cout << "epsilon " << tolerance << std::endl;
		//std::cout << "triangle expanded: \n" << striangle[0] <<"\n--------\n"<< striangle[1] <<"\n---------\n"<< striangle[2] << std::endl;
		//std::cout << "prism: \n" << envbox[0][0] << "\n,\n"<<envbox[0][1] << "\n,\n" << envbox[0][2] << "\n,\n" 
		//	<< envbox[0][3] << "\n,\n" << envbox[0][4] << "\n,\n" << envbox[0][5] << "\n,\n" << std::endl;
	}

	void Envelop::SLIntersection(const Parameters & params, const std::array<Vector3, 3>& cutface, const Vector3 & linepoints0, 
		const Vector3 & linepoints1, int & cutOrnot, Vector3 & interp, Scalar & t, int & ion)
	{
			const Vector3 n = (cutface[1] - cutface[0]).cross(cutface[2] - cutface[0]);
			const Vector3 pm = linepoints1 - linepoints0;
			const Vector3 bm = linepoints0 - cutface[0];
			ion = 0;
			Vector3 A;
			Eigen::Matrix<Scalar, 3, 2> B;
			Vector2 sl;

			Vector3 nn = n.normalized();
			Vector3 npm = pm.normalized();
			//std::cout << nn << "\n" << npm << "\n"  << std::endl;
			//std::cout << nn.dot(bm) << "\n" << nn.dot(pm) << "\n" << std::endl;
			if (abs(nn.dot(npm)) < SCALAR_ZERO) {//if the line is parallel to the face. TODO check this parameter
				cutOrnot = 0;
				interp.setZero();
				return;
			}
			else
			{
				t = (-1) * nn.dot(bm) / nn.dot(pm);
				//std::cout << nn.dot(bm) << "\n" << nn.dot(pm) << "\n" << std::endl;
				if (t <= 1 && t >= 0) {
					cutOrnot = 1;
					interp = linepoints0 + t * pm;
				}
				else
				{
					cutOrnot = 0;
					interp.setZero();
					return;
				}
			}
			for (int i = 0; i < 3; i++) {
				A[i] = interp[i] - cutface[0][i];
				B(i, 0) = cutface[1][i] - cutface[0][i];
				B(i, 1) = cutface[2][i] - cutface[0][i];
			}
			sl = B.jacobiSvd(Eigen::ComputeFullV | Eigen::ComputeFullU).solve(A);
			if (sl[0] >= 0 - SCALAR_ZERO &&sl[0] <= 1 + SCALAR_ZERO &&sl[1] >= 0 - SCALAR_ZERO &&sl[1] < 1 + SCALAR_ZERO) {
				ion = 1;//the line goes through the triangle
			}
	}
	void Envelop::SLIntersectionBackup(const Parameters & params, const std::array<Vector3, 3>& cutface, const Vector3 & linepoints0,
		const Vector3 & linepoints1, int & cutOrnot, Vector3 & interp, Scalar & t, int & ion)
	{

		// Scalar a, b, c, xpm, ypm, zpm, xbm, ybm, zbm;//xpm=line[1,0]-line[0,0],xbm=line[0,0]-cutface[0,0]

		// a = (cutface[1][1] - cutface[0][1])*(cutface[2][2] - cutface[0][2]) - (cutface[2][1] - cutface[0][1])*(cutface[1][2] - cutface[0][2]);
		// b = (cutface[1][2] - cutface[0][2])*(cutface[2][0] - cutface[0][0]) - (cutface[2][2] - cutface[0][2])*(cutface[1][0] - cutface[0][0]);
		// c = (cutface[1][0] - cutface[0][0])*(cutface[2][1] - cutface[0][1]) - (cutface[2][0] - cutface[0][0])*(cutface[1][1] - cutface[0][1]);
		const Vector3 n = (cutface[1] - cutface[0]).cross(cutface[2] - cutface[0]);
		const Vector3 pm = linepoints1 - linepoints0;
		const Vector3 bm = linepoints0 - cutface[0];
		ion = 0;
		Vector3 A;
		Eigen::Matrix<Scalar, 3, 2> B;
		//Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> B(3,2);
		//int n = 20;
		//... somethin with n
		//B.resize(10, n);
		Vector2 sl;

		Vector3 nn = n.normalized();
		Vector3 npm = pm.normalized();
		//std::cout << nn << "\n" << npm << "\n"  << std::endl;
		//std::cout << nn.dot(bm) << "\n" << nn.dot(pm) << "\n" << std::endl;
		if (abs(nn.dot(npm)) < SCALAR_ZERO) {//if the line is parallel to the face. TODO check this parameter
			cutOrnot = 0;
			interp.setZero();
			return;
		}
		else
		{
			t = (-1) * nn.dot(bm) / nn.dot(pm);
			//std::cout << nn.dot(bm) << "\n" << nn.dot(pm) << "\n" << std::endl;
			if (t <= 1 && t >= 0) {
				cutOrnot = 1;
				interp = linepoints0 + t * pm;
			}
			else
			{
				cutOrnot = 0;
				interp.setZero();
				return;
			}
		}
		for (int i = 0; i < 3; i++) {
			A[i] = interp[i] - cutface[0][i];
			B(i, 0) = cutface[1][i] - cutface[0][i];
			B(i, 1) = cutface[2][i] - cutface[0][i];
		}
		sl = B.jacobiSvd(Eigen::ComputeFullV | Eigen::ComputeFullU).solve(A);
		if (sl[0] >= 0 - SCALAR_ZERO && sl[0] <= 1 + SCALAR_ZERO && sl[1] >= 0 - SCALAR_ZERO && sl[1] < 1 + SCALAR_ZERO) {
			ion = 1;//the line goes through the triangle
		}
	}
	void Envelop::BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, const Scalar& bbd)
	{

		Vector3 AB, AC, BC, normal, vector1,ABn;
		Parameters pram;
		std::array<Vector3, 6> polygon;
		std::array<Vector3, 12> polygonoff;
		Scalar  a, b, c,

			tolerance = bbd * pram.eps_rel/sqrt(3),
			//tolerance = 2*bbd * pram.eps_rel,//TODO this is not right

			area;
		//envprism.resize(m_faces.size());// to make the results safer
		for (int i = 0; i < m_faces.size(); i++) {
			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			a = BC.norm();
			b = AC.norm();
			c = AB.norm();
			//ABn = AB.normalized();
			//ACn = AC.normalized();
			//crsp = ABn.dot(ACn);
			area = 0.25*sqrt((a + b + c)*(a + b - c)*(a + c - b)*(b + c - a));
			if (area<SCALAR_ZERO) {
				continue;
			}
			normal = AB.cross(AC).normalized();
			vector1 = AB.cross(normal).normalized();
			ABn = AB.normalized();
			polygon[0] = m_ver[m_faces[i][0]] + (vector1-ABn) * tolerance;
			polygon[1] = m_ver[m_faces[i][1]] + (vector1 + ABn) * tolerance;
			if (AB.dot(BC) < 0) {
				polygon[2] = m_ver[m_faces[i][1]] + (-vector1 + ABn) * tolerance;
				polygon[3] = m_ver[m_faces[i][2]] + (-vector1 + ABn) * tolerance;
				polygon[4] = m_ver[m_faces[i][2]] + (-vector1 - ABn) * tolerance;
				if (AB.dot(AC) < 0) {
					polygon[5] = m_ver[m_faces[i][2]] + (vector1 - ABn) * tolerance;
				}
				else {
					polygon[5] = m_ver[m_faces[i][0]] + (-vector1 - ABn) * tolerance;
				}
			}
			else {
				polygon[2] = m_ver[m_faces[i][2]] + (vector1 + ABn) * tolerance;
				polygon[3] = m_ver[m_faces[i][2]] + (-vector1 + ABn) * tolerance;
				polygon[4] = m_ver[m_faces[i][2]] + (-vector1 - ABn) * tolerance;
				polygon[5] = m_ver[m_faces[i][0]] + (-vector1 - ABn) * tolerance;
			}
			//polygon[2] = m_ver[m_faces[i][1]] + (-vector1 + ABn) * tolerance;
			/*if (abs(crsp) > 1 - SCALAR_ZERO) {//TODO how to judge degeneration? the sum of two edges less than third one?

				std::cout << "Degenerated prism: NO. " << i << std::endl;
				continue;
			}
			*/ //it is correct even if it is a degenerated triangle
			for (int j = 0; j < 6; j++) {
				polygonoff[j] =  polygon[j] + normal * tolerance ;
			}
			for (int j = 6; j < 12; j++) {
				polygonoff[j] = polygon[j-6] - normal * tolerance;
			}
			envprism.push_back(polygonoff);
		
		}
			
		//std::cout << "epsilon " << tolerance << std::endl;
		//std::cout << "triangle expanded: \n" << striangle[0] <<"\n--------\n"<< striangle[1] <<"\n---------\n"<< striangle[2] << std::endl;
		//std::cout << "prism: \n" << envbox[0][0] << "\n,\n"<<envbox[0][1] << "\n,\n" << envbox[0][2] << "\n,\n" 
		//	<< envbox[0][3] << "\n,\n" << envbox[0][4] << "\n,\n" << envbox[0][5] << "\n,\n" << std::endl;
	}
}





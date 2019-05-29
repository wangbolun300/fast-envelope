#include "FastEnvelope.h"
#include<fastenvelope/Parameters.h>
#include<fastenvelope/Predicates.hpp>
#include <fenv.h>
namespace fastEnvelope {

	bool FastEnvelope::FastEnvelopeTest(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism)
	{
		int p_face[8][3] = { {0,1,2},{8,7,6},{1,0,7},{2,1,7},{3,2,8},{3,9,10},{5,4,11},{0,5,6} };//prism triangle index. all with orientation.

		std::vector<std::array<Vector3, 2>> seglist;
		std::vector<std::array<int,2>> segtoprism;//segment belongs to which prism,always 3 element shorter than seglist
		std::vector<Vector3> interp; 
		//interp.reserve(envprism.size() * 20);// TODO why this cost more time?
		
		std::vector<int> interseg,jump;
		//interseg.reserve(envprism.size() * 20);// TODO why this cost more time?
		bool out;
		seglist.push_back({ triangle[0],triangle[1] });//TODO how to use emplace back?
		seglist.push_back({ triangle[0],triangle[2] });
		seglist.push_back({ triangle[1],triangle[2] });

		for (int i = 0; i < 3; i++) {
			out = point_out_prism(triangle[i], envprism, jump);
			if (out == 1) {
				return 1;
			}
		}
		for (int i = 0; i < envprism.size(); i++) {
			for (int j = 0; j < 8; j++) {
				Segment_facet_intersection(seglist, { envprism[i][p_face[j][0]],envprism[i][p_face[j][1]], envprism[i][p_face[j][2]] },interp , interseg);
				if (interp.size() > 0) {
					segtoprism.push_back({ i,j });
					jump.clear();
					jump.emplace_back(i);
					for (int k = 0; k < 2; k++) {
						out = point_out_prism(interp[k], envprism, jump);
						if (out == 1) {
							return 1;
						}
					}
					for (int k = 2; k < interp.size(); k++) {
						jump.clear();
						jump.emplace_back(segtoprism[interseg[k] - 3][0]);
						jump.emplace_back(i);
						out = point_out_prism(interp[k], envprism, jump);
						if (out == 1) {
							return 1;
						}
					}
				}
			}	
		}
		return 0;
	}



	void FastEnvelope::SLIntersection(const Parameters & params, const std::array<Vector3, 3>& cutface, const Vector3 & linepoints0,
		const Vector3 & linepoints1, int & cutOrnot, Vector3 & interp)
	{
		Scalar t;
		const Vector3 n = (cutface[1] - cutface[0]).cross(cutface[2] - cutface[0]);
		const Vector3 pm = linepoints1 - linepoints0;
		const Vector3 bm = linepoints0 - cutface[0];


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

	}

	void FastEnvelope::Segment_facet_intersection(std::vector<std::array<Vector3, 2>>& seglist, const std::array<Vector3, 3>& facet, std::vector<Vector3>& interp, std::vector<int>& interseg)
	{
		interp.clear();
		//interp.reserve(seglist.size());
		interseg.clear();
		//interseg.reserve(seglist.size());
		Parameters params;
		int  cutOrnot;
		Vector3  ipoint;
		Scalar  t;
		int  ion,rcd=0;
		for (int i = 0; i < 3; i++) {
			SLIntersection(params, facet, seglist[i][0], seglist[i][1], cutOrnot, ipoint);
			if (cutOrnot > 0) {
				rcd = rcd + 1;
				interp.emplace_back(ipoint);
				interseg.emplace_back(i);
				if (rcd == 2) {
					break;
				}
			}
		}
		if (rcd == 0) {
			interp.clear();
			interseg.clear();
			return;
		}
		if (rcd ==1) {
			interp.clear();
			interseg.clear();
			std::cout << "ERROR ENVELOPE TEST" << std::endl;
			return;
		}
		for (int i = 3; i < seglist.size(); i++) {

			SLIntersection(params, facet, seglist[i][0], seglist[i][1], cutOrnot, ipoint);
			if (cutOrnot > 0) {
				interp.emplace_back(ipoint);
				interseg.emplace_back(i);
			}
		}
		seglist.push_back({ interp[0],interp[1] });
	}
	

	bool FastEnvelope::point_out_prism(const Vector3 & point, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump)
	{
		int jm=0,ori, p_face[8][3] = { {0,1,2},{8,7,6},{1,0,7},{2,1,7},{3,2,8},{3,9,10},{5,4,11},{0,5,6} };//prism triangle index. all with orientation.

		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm+1);
					continue;
				}
			}
			
			for (int j = 0; j < 8; j++) {

				ori = Predicates::orient_3d(envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]], point);
				if (ori == -1||ori==0) {
					break;
				}
				if (j == 7) {
					return 0;
				}
			}
			
		}
		return 1;
	}

	void FastEnvelope::BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, const Scalar& bbd)
	{

		Vector3 AB, AC, BC, normal, vector1, ABn;
		Parameters pram;
		std::array<Vector3, 6> polygon;
		std::array<Vector3, 12> polygonoff;
		Scalar  a, b, c,
			tolerance = bbd * pram.eps_rel / sqrt(3),
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
			if (area < SCALAR_ZERO) {
				continue;
			}
			normal = AB.cross(AC).normalized();
			vector1 = AB.cross(normal).normalized();
			ABn = AB.normalized();
			polygon[0] = m_ver[m_faces[i][0]] + (vector1 - ABn) * tolerance;
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
				polygonoff[j] = polygon[j] + normal * tolerance;
			}
			for (int j = 6; j < 12; j++) {
				polygonoff[j] = polygon[j - 6] - normal * tolerance;
			}
			envprism.push_back(polygonoff);

		}

		//std::cout << "epsilon " << tolerance << std::endl;
		//std::cout << "triangle expanded: \n" << striangle[0] <<"\n--------\n"<< striangle[1] <<"\n---------\n"<< striangle[2] << std::endl;
		//std::cout << "prism: \n" << envbox[0][0] << "\n,\n"<<envbox[0][1] << "\n,\n" << envbox[0][2] << "\n,\n" 
		//	<< envbox[0][3] << "\n,\n" << envbox[0][4] << "\n,\n" << envbox[0][5] << "\n,\n" << std::endl;
	}


// Indirect 3D orientation predicate with floating point filter.
//
// Input: eight points p, q, r, s, t, a, b, c
// where:
// <p, q> define a straight line L
// <r, s, t> define a plane P1
// <a, b, c> define a plane P2
//
// Output:
// Let i be the exact intersection point between L and P1.
// orient3D_LPI returns 1 if the tetrahedron (i, a, b, c) has positive volume.
// Returns -1 if it has negative volume.
// Zero is returned if one of the following conditions holds:
// - degenerate input (L and P1 are parallel, <r, s, t> are collinear, ...)
// - i is exactly on P2
// - the floating point precision is insufficient to guarantee an exact answer
// - the input coordinates cause an under/overflow during the computation.

	inline int orient3D_LPI(double px, double py, double pz,
		double qx, double qy, double qz,
		double rx, double ry, double rz,
		double sx, double sy, double sz,
		double tx, double ty, double tz,
		double ax, double ay, double az,
		double bx, double by, double bz,
		double cx, double cy, double cz)
	{
		::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

		double a11, a12, a13, a21, a22, a23, a31, a32, a33;
		double px_rx, py_ry, pz_rz;
		double a2233, a2133, a2132;
		double d, n;
		double ix, iy, iz;
		double m12, m13, m14, m23, m24, m34;
		double m123, m124, m134, m234;
		double det4x4_return_value;

		a11 = (px - qx);
		a12 = (py - qy);
		a13 = (pz - qz);
		a21 = (sx - rx);
		a22 = (sy - ry);
		a23 = (sz - rz);
		a31 = (tx - rx);
		a32 = (ty - ry);
		a33 = (tz - rz);
		px_rx = px - rx;
		py_ry = py - ry;
		pz_rz = pz - rz;
		a2233 = ((a22 * a33) - (a23 * a32));
		a2133 = ((a21 * a33) - (a23 * a31));
		a2132 = ((a21 * a32) - (a22 * a31));
		d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));
		n = ((((py_ry)* a2133) - ((px_rx)* a2233)) - ((pz_rz)* a2132));
		ix = ((d * px) + (a11 * n));
		iy = ((d * py) + (a12 * n));
		iz = ((d * pz) + (a13 * n));
		m12 = (((d * ax) * iy) - (ix * (d * ay)));
		m13 = (((d * bx) * iy) - (ix * (d * by)));
		m14 = (((d * cx) * iy) - (ix * (d * cy)));
		m23 = (((d * bx) * (d * ay)) - ((d * ax) * (d * by)));
		m24 = (((d * cx) * (d * ay)) - ((d * ax) * (d * cy)));
		m34 = (((d * cx) * (d * by)) - ((d * bx) * (d * cy)));
		m123 = (((m23 * iz) - (m13 * (d * az))) + (m12 * (d * bz)));
		m124 = (((m24 * iz) - (m14 * (d * az))) + (m12 * (d * cz)));
		m134 = (((m34 * iz) - (m14 * (d * bz))) + (m13 * (d * cz)));
		m234 = (((m34 * (d * az)) - (m24 * (d * bz))) + (m23 * (d * cz)));
		det4x4_return_value = m234 - m134 + m124 - m123;
		if (::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID)) return 0; // Fast reject in case of under/overflow


		// Almost static filter
		double fa11 = fabs(a11);
		double fa21 = fabs(a21);
		double fa31 = fabs(a31);
		double fa12 = fabs(a12);
		double fa22 = fabs(a22);
		double fa32 = fabs(a32);
		double fa13 = fabs(a13);
		double fa23 = fabs(a23);
		double fa33 = fabs(a33);
		double fax = fabs(ax);
		double fay = fabs(ay);
		double faz = fabs(az);
		double fbx = fabs(bx);
		double fby = fabs(by);
		double fbz = fabs(bz);
		double fcx = fabs(cx);
		double fcy = fabs(cy);
		double fcz = fabs(cz);
		double fpxrx = fabs(px_rx);
		double fpyry = fabs(py_ry);
		double fpzrz = fabs(pz_rz);

		double max1, max2, max3, max4, max5, max6, max7, max8;
		max4 = fa11;
		if (max4 < fa31) max4 = fa31;
		if (max4 < fa21) max4 = fa21;
		max5 = max4;
		if (max5 < fpxrx)  max5 = fpxrx;
		max1 = max5;
		if (max1 < fbx) max1 = fbx;
		if (max1 < fax) max1 = fax;
		if (max1 < fcx) max1 = fcx;
		max2 = fbz;
		if (max2 < faz) max2 = faz;
		if (max2 < fcz) max2 = fcz;
		if (max2 < fa13) max2 = fa13;
		max6 = fa12;
		if (max6 < fa22) max6 = fa22;
		if (max6 < fa32) max6 = fa32;
		max3 = max6;
		if (max3 < fay) max3 = fay;
		if (max3 < fcy) max3 = fcy;
		if (max3 < fby) max3 = fby;
		if (max3 < fpyry) max3 = fpyry;
		max7 = fa13;
		if (max7 < fa23) max7 = fa23;
		if (max7 < fa33) max7 = fa33;
		max8 = max7;
		if (max8 < fpzrz) max8 = fpzrz;

		double eps = 3.53761371545404460000e-011 * max6 * max7 * max4 * max1 * max6 * max7 * max4 * max3 * max3 * max8 * max5 * max2;
		if ((det4x4_return_value > eps)) return -1;
		else if ((det4x4_return_value < -eps)) return 1;
		else return 0;
	}
}
#include "FastEnvelopeV2.h"
#include<fastenvelope/Parameters.h>
#include<fastenvelope/Predicates.hpp>
#include <fstream>
#include <istream>
#include <igl/Timer.h>
#include<fastenvelope/intersections.h>
//#include<fastenvelope/intersections.h>

//igl::Timer ftimer1,ftimer2,ftimer3;
fastEnvelope::Scalar ftime1 = 0, ftime2 = 0, ftime3 = 0;
std::ofstream fout, fout1, fout2, fout3, fout4, fout5, fout6;
int call = 0, kmark = 0, kmark1 = 0;
int markhf = 0, recordnumber = 0, rcdn = 0;
static const int p_face[8][3] = { {0,1,2},{8,7,6},{1,0,7},{2,1,7},{3,2,8},{3,9,10},{5,4,11},{0,5,6} };//prism triangle index. all with orientation.
static const std::array<std::vector<fastEnvelope::Vector3i>, 8> p_triangle = {
		{
			{fastEnvelope::Vector3i(0,1,2),fastEnvelope::Vector3i(0,2,5),fastEnvelope::Vector3i(5,2,3),fastEnvelope::Vector3i(5,3,4)},
	{fastEnvelope::Vector3i(8,7,6),fastEnvelope::Vector3i(8,6,11),fastEnvelope::Vector3i(9,8,11),fastEnvelope::Vector3i(9,11,10)},
	{fastEnvelope::Vector3i(1,0,7),fastEnvelope::Vector3i(7,0,6)},
	{fastEnvelope::Vector3i(1,7,2),fastEnvelope::Vector3i(2,7,8)},
	{fastEnvelope::Vector3i(2,8,3),fastEnvelope::Vector3i(3,8,9)},
	{fastEnvelope::Vector3i(3,9,4),fastEnvelope::Vector3i(4,9,10)},
	{fastEnvelope::Vector3i(4,10,11),fastEnvelope::Vector3i(5,4,11)},
	{fastEnvelope::Vector3i(0,5,6),fastEnvelope::Vector3i(5,11,6)}

		}
};//prism triangle index. all with orientation.

static const std::array<std::array<int, 2>, 3> triseg = {
	
	{{{0,1}},{{0,2}},{{1,2}}}

};

extern "C++" int tri_tri_intersection_test_3d(fastEnvelope::Scalar p1[3], fastEnvelope::Scalar q1[3], fastEnvelope::Scalar r1[3],
	fastEnvelope::Scalar p2[3], fastEnvelope::Scalar q2[3], fastEnvelope::Scalar r2[3],
	int * coplanar,
	fastEnvelope::Scalar source[3], fastEnvelope::Scalar target[3]);

namespace fastEnvelope {
	using namespace std;
	bool FastEnvelope::FastEnvelopeTest(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism)
	{
		
		if (envprism.size() == 0) {
			return 1;
		}

		std::vector<std::array<Vector3, 2>> seglist;
		//seglist.reserve(envprism.size());
		std::vector<std::array<int, 2>> segtoprism;//segment belongs to which prism,always 3 element shorter than seglist
		//segtoprism.reserve(envprism.size());
		std::vector<Vector3> interp;
		//interp.reserve(envprism.size() * 20);
		std::vector<int> interseg, jump;
		//interseg.reserve(envprism.size() * 20);
		bool out;
		//seglist.push_back({ { triangle[0],triangle[1] } });
		seglist.emplace_back(std::array<Vector3, 2>{ {Vector3(triangle[0]), Vector3(triangle[1])}});
		seglist.emplace_back(std::array<Vector3, 2>{{Vector3(triangle[0]), Vector3(triangle[2])}});
		seglist.emplace_back(std::array<Vector3, 2>{ {Vector3(triangle[1]), Vector3(triangle[2])}});
		//fout.open("D:\\vs\\float project\\data\\pointlist.txt");
		//fout1.open("D:\\vs\\float project\\data\\seglist.txt");
		for (int i = 0; i < 3; i++) {
			//ftimer1.start();
			out = point_out_prism(triangle[i], envprism, jump);
			//ftime1 = ftime1 + ftimer1.getElapsedTimeInSec();
			//fout<<std::setprecision(17) << triangle[i][0] << "              " << triangle[i][1] << "             " << triangle[i][2] << std::endl;
			if (out == true) {
				return 1;
			}
		}
		for (int i = 0; i < envprism.size(); i++) {
			for (int j = 0; j < 8; j++) {
				Segment_facet_intersection(seglist, {{envprism[i][p_face[j][0]],envprism[i][p_face[j][1]], envprism[i][p_face[j][2]] }}, interp, interseg);
				if (interp.size() > 0) {
					
					segtoprism.emplace_back(std::array<int, 2>{ { i,j } });
					jump.clear();
					jump.emplace_back(i);
					
				
					for (int k = 0; k < 2; k++) {
						//ftimer1.start();
						out = point_out_prism(interp[k], envprism, jump);
						//ftime1 = ftime1 + ftimer1.getElapsedTimeInSec();
						if (out == true) {
							return 1;
						}
					}
					for (int k = 2; k < interp.size(); k++) {
						jump.clear();
						jump.emplace_back(segtoprism[interseg[k] - 3][0]);
						jump.emplace_back(i);
						//ftimer1.start();
						out = point_out_prism(interp[k], envprism, jump);
						//ftime1 = ftime1 + ftimer1.getElapsedTimeInSec();
						if (out == true) {
							return 1;
						}
					}
				}
			}
		}

		return 0;
	}

	bool FastEnvelope::FastEnvelopeTestTemp(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism)
	{

		if (envprism.size() == 0) {
			return 1;
		}

		std::vector<std::array<Vector3, 2>> seglist;
		//seglist.reserve(envprism.size()/10);//TODO it is faster but the parameter may change from time to time
		std::vector<std::array<int, 2>> segtoprism;//segment belongs to which prism,always 3 element shorter than seglist
		//segtoprism.reserve(envprism.size()/10);//TODO it is faster but the parameter may change from time to time
		std::vector<Vector3> interp;
		
		std::vector<int> interseg, jump;
		//interseg.reserve(envprism.size() * 20);
		bool out;
		int tti;//triangle-triangle intersection
		//seglist.push_back({ { triangle[0],triangle[1] } });
		seglist.emplace_back(std::array<Vector3, 2>{ {Vector3(triangle[0]), Vector3(triangle[1])}});
		seglist.emplace_back(std::array<Vector3, 2>{ {Vector3(triangle[0]), Vector3(triangle[2])}});
		seglist.emplace_back(std::array<Vector3, 2>{ {Vector3(triangle[1]), Vector3(triangle[2])}});
		//fout.open("D:\\vs\\float project\\data\\pointlist.txt");
		//fout1.open("D:\\vs\\float project\\data\\seglist.txt");
		for (int i = 0; i < 3; i++) {
			//ftimer1.start();
			out = point_out_prism(triangle[i], envprism, jump);
			//ftime1 = ftime1 + ftimer1.getElapsedTimeInSec();
			//fout<<std::setprecision(17) << triangle[i][0] << "              " << triangle[i][1] << "             " << triangle[i][2] << std::endl;
			if (out == true) {
				return 1;
			}
		}
		for (int i = 0; i < envprism.size(); i++) {
			for (int j = 0; j < 8; j++) {
				for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
					tti = tri_cut_tri_simple(triangle[0], triangle[1],triangle[2], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
					if (tti == 4 ) {
						break;
					}
					if (tti == -1) {
						continue;
					}
					Segment_facet_intersection(seglist, {{ envprism[i][p_face[j][0]],envprism[i][p_face[j][1]], envprism[i][p_face[j][2]] }}, interp, interseg);
				
					if (interp.size() > 0) {

						segtoprism.emplace_back(std::array<int, 2>{ { i, j } });
						jump.clear();
						jump.emplace_back(i);


						for (int k = 0; k < 2; k++) {
							//ftimer1.start();

							/*if (i == 31 && j == 2 && c == 1 && k == 0) {
								std::cout << "heeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeer" << std::endl;
							}*/
							out = point_out_prism(interp[k], envprism, jump);
							/*if (i == 31 && j == 2 && c == 1 && k == 0) {
								std::cout << "the point is "<<interp[0] << std::endl;
								std::cout << "heeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeer" << std::endl;
								std::cout << "out "<<out << std::endl;
							}*/
							//ftime1 = ftime1 + ftimer1.getElapsedTimeInSec();
							if (out == true) {
								return 1;
							}
						}
						for (int k = 2; k < interp.size(); k++) {
							jump.clear();
							jump.emplace_back(segtoprism[interseg[k] - 3][0]);
							jump.emplace_back(i);
							//ftimer1.start();
							out = point_out_prism(interp[k], envprism, jump);
							//ftime1 = ftime1 + ftimer1.getElapsedTimeInSec();
							if (out == true) {
								return 1;
							}
						}

					}
					break;
				}//each triangle of the facet

				
			}
		}

		return 0;
	}
	
	bool FastEnvelope::FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism)
	{

		if (envprism.size() == 0) {
			return 1;
		}
		std::vector<int> jump;
		std::vector<Vector3i> inter_ijk_list;//list of intersected triangle
		bool out;
		int inter, inter1, record1, record2, tti;//triangle-triangle intersection
		jump.clear();
		for (int i = 0; i < 3; i++) {
			out = point_out_prism(triangle[i], envprism, jump);
			if (out == true) {
				//std::cout << "111111111111111111111111111out 1" << std::endl;
				return 1;
			}
		}
		for (int i = 0; i < envprism.size(); i++) {
			for (int j = 0; j < 8; j++) {
				for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
					tti = tri_cut_tri_simple(triangle[0], triangle[1], triangle[2], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
					if (tti == 4) {
						break;
					}
					if (tti == -1) {
						continue;
					}

					record1 = 0;

					jump.clear();
					jump.emplace_back(i);
					for (int k = 0; k < 3; k++) {

						inter = Implicit_Seg_Facet_interpoint_Out_Prism(triangle[triseg[k][0]], triangle[triseg[k][1]],
							{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);

						/*if (i == 46 && j == 4 && c == 1 && k == 0) {
							std::cout << "heeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeer0" << std::endl;
						}*/
						/*inter = Implicit_Seg_Facet_interpoint_Out_Prism_Wang(triangle[triseg[k][0]], triangle[triseg[k][1]],
							{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);*/
							/*if (i == 46 && j == 4 && c == 1 && k == 0) {
								std::cout << "the seg: "<< triangle[triseg[k][0]] <<"\n,\n"<< triangle[triseg[k][1]] << std::endl;
								std::cout << "heeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeer0" << std::endl;
							}*/
						if (inter == 1) {
							/*std::cout << "!!!!!!!!!!!!!!!!!out 2, and which edge "<<k << std::endl;
							std::cout << "ijc " << i<<" "<<j<<" "<<c << std::endl;*/
							return 1;
						}
						record1 = record1 + inter;
					}
					if (record1 >= 4) {
						std::cout << "intersection predicate wrong, record " << record1 << std::endl;
						//std::cout << "inter " << inter[0]<<" "<<inter[1]<<" "<<inter[2] << std::endl;
						//return 1;
					}


					for (int e = 0; e < inter_ijk_list.size(); e++) {
						for (int f = inter_ijk_list[e][2]; f < p_triangle[inter_ijk_list[e][1]].size(); f++) {
							tti = tri_cut_tri_simple(triangle[0], triangle[1], triangle[2],
								envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]]);
							if (tti == 4) {
								break;
							}
							if (tti == -1) {
								continue;
							}
							jump.clear();
							jump.emplace_back(inter_ijk_list[e][0]);
							jump.emplace_back(i);
							inter1 = Implicit_Tri_Facet_Facet_interpoint_Out_Prism(triangle,
								{ {envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]]} },
								{ {envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]} }, envprism, jump);
							if (inter1 == 1) {
								//std::cout << "ijc " << i << " " << j << " " << c << std::endl;
								//std::cout << "inter list ijk " << inter_ijk_list[e][0] << " " << inter_ijk_list[e][1] << " " << f << std::endl;
								//std::cout << "111111111111111111111out 3" << std::endl;
								//std::ofstream fout;
								//fout.open("D:\\vs\\float project\\data\\tri.txt");
								//fout<< std::setprecision(17)<< triangle[0][0] << " " << triangle[0][1] << " " << triangle[0][2] << endl;
								//fout << std::setprecision(17) << triangle[1][0] << " " << triangle[1][1] << " " << triangle[1][2] << endl;
								//fout << std::setprecision(17) << triangle[2][0] << " " << triangle[2][1] << " " << triangle[2][2] << endl;
								//fout.close();
								//fout.open("D:\\vs\\float project\\data\\facet1.txt");
								//fout << std::setprecision(17) << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]][0] << " " << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]][1] << " " << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]][2] << endl;
								//fout << std::setprecision(17) << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]][0] << " " << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]][1] << " " << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]][2] << endl;
								//fout << std::setprecision(17) << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]][0] << " " << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]][1] << " " << envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]][2] << endl;
								//fout.close();
								//fout.open("D:\\vs\\float project\\data\\facet2.txt");
								//fout << std::setprecision(17) << envprism[i][p_triangle[j][c][0]][0] << " " << envprism[i][p_triangle[j][c][0]][1] << " " << envprism[i][p_triangle[j][c][0]][2] << endl;
								//fout << std::setprecision(17) << envprism[i][p_triangle[j][c][1]][0] << " " << envprism[i][p_triangle[j][c][1]][1] << " " << envprism[i][p_triangle[j][c][1]][2] << endl;
								//fout << std::setprecision(17) << envprism[i][p_triangle[j][c][2]][0] << " " << envprism[i][p_triangle[j][c][2]][1] << " " << envprism[i][p_triangle[j][c][2]][2] << endl;
								//fout.close();
								//
								//kmark = 1;
								////jump.clear();
								//std::cout << "jump " << jump[0]<<" "<<jump[1] << std::endl;
								//inter1=Implicit_Tri_Facet_Facet_interpoint_Out_Prism(triangle,
								//	{ {envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]]} },
								//	{ {envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]} }, envprism, jump);
								//std::cout << "here inter1 " << inter1 << std::endl;
								return 1;
							}
						}
					}
					inter_ijk_list.emplace_back(Vector3i(i, j, c));
					break;
				}//each triangle of the facet
			}
		}

		return 0;
	}



	bool FastEnvelope::is_seg_facet_intersection(const double& px, const double& py, const double& pz,
		const double& qx, const double& qy, const double& qz,
		const double& rx, const double& ry, const double& rz,
		const double& sx, const double& sy, const double& sz,
		const double& tx, const double& ty, const double& tz,
		double& a11, double& a12, double& a13,
		double& a21, double& a22, double& a23,
		double& a31, double& a32, double& a33,
		double& px_rx, double& py_ry, double& pz_rz,
		double& d, double& n) {

		double a2233, a2133, a2132;
		a11 = (px - qx);
		a12 = (py - qy);
		a13 = (pz - qz);//a1: qp
		a21 = (sx - rx);
		a22 = (sy - ry);
		a23 = (sz - rz);//a2: rs
		a31 = (tx - rx);
		a32 = (ty - ry);
		a33 = (tz - rz);//a3: rt
		px_rx = px - rx;
		py_ry = py - ry;
		pz_rz = pz - rz;//rp
		a2233 = ((a22 * a33) - (a23 * a32));
		a2133 = ((a21 * a33) - (a23 * a31));
		a2132 = ((a21 * a32) - (a22 * a31));
		d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));
		n = ((((py_ry)* a2133) - ((px_rx)* a2233)) - ((pz_rz)* a2132));

		if (d != 0) {
			if (-1 * n / d >= 0 && -1 * n / d <= 1) {
				return 1;
			}
		}
		return 0;
	}
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
		const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		double  a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n;
		int jm = 0, ori, ori1;
		bool inter = is_seg_facet_intersection(segpoint0[0], segpoint0[1], segpoint0[2],
			segpoint1[0], segpoint1[1], segpoint1[2],
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n);
		if (inter == 0) {
			return 2;//not intersected
		}
		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}
			}

			for (int j = 0; j < 8; j++) {
				//ftimer2.start();
				ori = orient3D_LPI(
					segpoint0[0], segpoint0[1], segpoint0[2],
					envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
					envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
					envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2],
					a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n);
				/////////////////////////////////////////
				ori1 = -1 * Predicates::orient_3d(envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]], segpoint0 + (segpoint0 - segpoint1)*n / d);//because n is -n
				if (ori != ori1) {
					markhf = 0;
					if (recordnumber == 0) {
						fout.open("D:\\vs\\float project\\data\\output_0522\\segpoints.txt");
						fout1.open("D:\\vs\\float project\\data\\output_0522\\segindex.txt");
						fout2.open("D:\\vs\\float project\\data\\output_0522\\triangle_points.txt");
						fout3.open("D:\\vs\\float project\\data\\output_0522\\triangle_index.txt");
						fout4.open("D:\\vs\\float project\\data\\output_0522\\facet_points.txt");
						fout5.open("D:\\vs\\float project\\data\\output_0522\\facet_index.txt");

					}
					if (recordnumber < 10000) {
						/*std::cout << "number " << recordnumber<< std::endl;
						std::cout << "ori and ori1 " << ori << " " << ori1 << std::endl;*/
						ori = orient3D_LPI(
							segpoint0[0], segpoint0[1], segpoint0[2],
							envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
							envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
							envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2],
							a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n);

						fout << std::setprecision(17) << segpoint0[0] << " " << segpoint0[1] << " " << segpoint0[2] << std::endl;
						fout << std::setprecision(17) << segpoint1[0] << " " << segpoint1[1] << " " << segpoint1[2] << std::endl;
						fout1 << 2 * recordnumber << " " << 2 * recordnumber + 1 << std::endl;
						fout2 << std::setprecision(17) << triangle[0][0] << " " << triangle[0][1] << " " << triangle[0][2] << std::endl;
						fout2 << std::setprecision(17) << triangle[1][0] << " " << triangle[1][1] << " " << triangle[1][2] << std::endl;
						fout2 << std::setprecision(17) << triangle[2][0] << " " << triangle[2][1] << " " << triangle[2][2] << std::endl;
						fout3 << 3 * recordnumber << " " << 3 * recordnumber + 1 << " " << 3 * recordnumber + 2 << endl;
						fout5 << 3 * recordnumber << " " << 3 * recordnumber + 1 << " " << 3 * recordnumber + 2 << endl;
						fout4 << std::setprecision(17) << envprism[i][p_face[j][0]][0] << " " << envprism[i][p_face[j][0]][1] << " " << envprism[i][p_face[j][0]][2] << std::endl;
						fout4 << std::setprecision(17) << envprism[i][p_face[j][1]][0] << " " << envprism[i][p_face[j][1]][1] << " " << envprism[i][p_face[j][1]][2] << std::endl;
						fout4 << std::setprecision(17) << envprism[i][p_face[j][2]][0] << " " << envprism[i][p_face[j][2]][1] << " " << envprism[i][p_face[j][2]][2] << std::endl;

						//std::cout << "show you the small d: " <<d<< std::endl;
						recordnumber++;
					}
					if (recordnumber == 10000) {
						fout.close();
						fout1.close();
						fout2.close();
						fout3.close();
						fout4.close();
						fout5.close();
					}

					markhf = 0;

				}

				///////////////////////////////////////////
				if (ori == 1 || ori == 0) {
					break;
				}
				if (j == 7) {
					//std::cout << "which prism include this point: " << i  << std::endl;
					//ftime3 = ftime3 + ftimer3.getElapsedTimeInSec();
					return 0;
				}
			}


		}
		return 1;
	}
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_Wang(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
		const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		Scalar d, n;
		int jm = 0, ori;
		bool inter = is_seg_facet_intersection_Wang(segpoint0, segpoint1, triangle, d, n);
		if (inter == false) {
			return 2;//not intersected
		}
		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}
			}

			for (int j = 0; j < 8; j++) {
				ori = Orient_LPI_Wang(segpoint0, segpoint1, d, n, { {envprism[i][p_face[j][0]],envprism[i][p_face[j][1]],envprism[i][p_face[j][2]]} });
				//ori = -1*Predicates::orient_3d(envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]], segpoint0 + (segpoint1 - segpoint0)*n / d);
				if (ori == 1 || ori == 0) {
					/*if (i == 32) {
						std::cout << "which facet is out,j " << j << std::endl;
						std::cout << "the point is " << segpoint0 + (segpoint1 - segpoint0)*n / d << std::endl;
						std::cout << "t is " << n / d << std::endl;
						std::cout << "in old way, ori is " << -1 * Predicates::orient_3d(envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]], segpoint0 + (segpoint1 - segpoint0)*n / d) << std::endl;

					}*/

					break;
				}
				if (j == 7) {
					//std::cout << "which prism include this point: " << i  << std::endl;
					//ftime3 = ftime3 + ftimer3.getElapsedTimeInSec();
					return 0;
				}
			}


		}
		return 1;
	}
	bool FastEnvelope::is_seg_facet_intersection_Wang(const Vector3& seg0, const Vector3& seg1, const std::array<Vector3, 3>& triangle, Scalar& d, Scalar&n) {
		Eigen::Matrix<Scalar, 3, 3> dm, nm;
		dm << (seg0 - seg1).transpose(),
			(triangle[1] - triangle[0]).transpose(),
			(triangle[2] - triangle[0]).transpose();
		nm << (seg0 - triangle[0]).transpose(),
			(triangle[1] - triangle[0]).transpose(),
			(triangle[2] - triangle[0]).transpose();
		d = dm.determinant();
		n = nm.determinant();
		if (fabs(d) < SCALAR_ZERO_3) {
			return 0;
		}
		Scalar t = n / d;
		if (t >= 0 && t <= 1) {
			return 1;
		}
		return 0;
	}
	int FastEnvelope::Orient_LPI_Wang(const Vector3& seg0, const Vector3& seg1, Scalar& d, Scalar&n, const std::array<Vector3, 3>& facet) {
		Eigen::Matrix<Scalar, 4, 4> M;
		M << (d*d*seg0 + n * d*(seg1 - seg0)).transpose(), 1,
			d*d*facet[0].transpose(), 1,
			d*d*facet[1].transpose(), 1,
			d*d*facet[2].transpose(), 1;
		//std::cout << M << std::endl;
		Scalar mv = M.determinant();
		if (mv > 0) {
			return 1;
		}
		if (mv < 0) {
			return -1;
		}

		return 0;
	}
	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2,
		const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		Eigen::Matrix<Scalar, 3, 3> A, AT, ATA, A1;
		Eigen::Matrix<Scalar, 3, 1> B;
		Eigen::Matrix<Scalar, 4, 4> C;
		Eigen::Matrix<Scalar, 2, 3> Prj;
		Vector3 p3d;
		Vector2 t0, t1, t2, p;
		int rcd = 0, jm = 0, ori;;
		Scalar ad, ad1, fad, fadi;
		AT << (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]),
			(facet1[0] - facet1[1]).cross(facet1[0] - facet1[2]),
			(facet2[0] - facet2[1]).cross(facet2[0] - facet2[2]);
		Scalar mv = fabs(AT(0, 0));
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				mv = mv == 0 ? fabs(AT(i, j)) : fabs(AT(i, j)) == 0 ? mv : std::min(mv, fabs(AT(i, j)));
			}
		}
		mv = mv < 1 ? 1 / mv : 1;

		A = AT.transpose();
		//////////////////////////////////////////////////////////////////
		//if (mv <= 0) {
		//	std::cout << "mv calculate wrong" << std::endl;
		//}
		//ad  = det3x3(A(0,0),A(0,1),A(0,2),
		//	A(1, 0), A(1, 1), A(1, 2),
		//	A(2, 0), A(2, 1), A(2, 2));
		//fad = fabs(ad);
		//if (fad < SCALAR_ZERO_3) {//TODO really need this?
		//	//std::cout << "           A singularity            " << std::endl;
		//	return 2;
		//}
		//fadi = fad < 1 ? mv / fad : mv;
		//A(0) = A(0)*fadi;
		//A(1) = A(1)*fadi;
		//A(2) = A(2)*fadi;
		/////////////////////////////////////////////////////////////////////////////
		A(0) = A(0)*mv;
		A(1) = A(1)*mv;
		A(2) = A(2)*mv;
		ad = det3x3(A(0, 0), A(0, 1), A(0, 2),
			A(1, 0), A(1, 1), A(1, 2),
			A(2, 0), A(2, 1), A(2, 2));
		fad = fabs(ad);
		if (fad < SCALAR_ZERO_3) {//TODO really need this?
		//	//std::cout << "           A singularity            " << std::endl;
			return 2;
		}
		fadi = fad < 1 ? 1 / fad : 1;
		A(0) = A(0)*fadi;
		A(1) = A(1)*fadi;
		A(2) = A(2)*fadi;
		////////////////////////////////////////////////////////////////////////////////
		B[0] = A(0, 0)*triangle[0][0] + A(0, 1)*triangle[0][1] + A(0, 2)*triangle[0][2];
		B[1] = A(1, 0)*facet1[0][0] + A(1, 1)*facet1[0][1] + A(1, 2)*facet1[0][2];
		B[2] = A(2, 0)*facet2[0][0] + A(2, 1)*facet2[0][1] + A(2, 2)*facet2[0][2];




		//ad1 = 1 / ad;
		//A1 = ad1 * A;
		//int t = get_t(triangle[0], triangle[1], triangle[2]);
		//t0 = to_2d(triangle[0], t);
		//t1 = to_2d(triangle[1], t);
		//t2 = to_2d(triangle[2], t);
		//////////////////////////////TODO not get p
		Vector3 f1, f2, f3;
		/*f1 = (ad1*B - A1 * triangle[1]).cross(A1*triangle[1] - A1*triangle[0]);
		f2 = (ad1*B - A1 * triangle[2]).cross(A1*triangle[2] - A1*triangle[1]);
		f3 = (ad1*B - A1 * triangle[0]).cross(A1*triangle[0] - A1*triangle[2]);*/
		f1 = (B - A * triangle[1]).cross(A*triangle[1] - A * triangle[0]);
		f2 = (B - A * triangle[2]).cross(A*triangle[2] - A * triangle[1]);
		f3 = (B - A * triangle[0]).cross(A*triangle[0] - A * triangle[2]);
		/*f1 = (A.inverse()*B -  triangle[1]).cross(triangle[1] - triangle[0]);
		f2 = (A.inverse()*B -  triangle[2]).cross(triangle[2] - triangle[1]);
		f3 = (A.inverse()*B -  triangle[0]).cross(triangle[0] - triangle[2]);*/
		int in = f1.dot(f2) > 0 && f1.dot(f3) > 0 ? 1 : 0;//TODO consider the triangle degeneration
		/*p3d = A.inverse()*B;
		for (int i = 0; i < 3; i++) {
			if (i == t) {
				continue;
			}
			p[rcd] = p3d[i];
			rcd = rcd + 1;
		}
		bool in = is_p_inside_tri_2d(p, { {to_2d(triangle[0], t),to_2d(triangle[1], t),to_2d(triangle[2], t)} });*/
		///////////////////////////////get p
		if (in == 0) {
			return 2;
		}
		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}

			}
			for (int j = 0; j < 8; j++) {
				ori = orient_3triangles(A, AT, ATA, B, { {envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]] } });
				if (ori == 1 || ori == 0) {
					/*if (i == 8) {
						kmark = 3;
						orient_3triangles(A, AT, ATA, B, { {envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]] } });
						kmark = 4;
						std::cout << "which facet did this! " << j << std::endl;
					}*/
					break;
				}
				if (j == 7) {

					return 0;
				}
			}
		}
		/*if (kmark == 1) {
			std::cout << "A determinant" << A.determinant() << std::endl;
			std::cout << "The point is " << A.inverse()*B << std::endl;
			std::cout << "f1.f2 and f1.f3 " << f1.dot(f2)<<" "<<f1.dot(f3) << std::endl;
			kmark1 = 2;
			std::cout<<"point test alone " << point_out_prism(A.inverse()*B, envprism, jump) << std::endl;
			kmark1 = 0;
			std::ofstream fout;
			fout.open("D:\\vs\\float project\\data\\singlepoint.txt");
			fout << std::setprecision(17) << (A.inverse()*B)[0] << " " << (A.inverse()*B)[1] << " " << (A.inverse()*B)[2] << endl;
			fout.close();
		}*/
		return 1;
	}




	int FastEnvelope::tri_cut_tri_simple(const Vector3& p1, const Vector3& p2, const Vector3& p3,
		const Vector3& q1, const Vector3& q2, const Vector3& q3) {
		std::array<Scalar, 3> p_1 = { {0, 0, 0} }, q_1 = { {0, 0, 0} }, r_1 = { {0, 0, 0} };
		std::array<Scalar, 3> p_2 = { {0, 0, 0} }, q_2 = { {0, 0, 0} }, r_2 = { {0, 0, 0} };
		int coplanar = 0;
		std::array<Scalar, 3> s = { {0,0,0} }, t = { {0,0,0} };
		for (int j = 0; j < 3; j++) {
			p_1[j] = p1[j];
			q_1[j] = p2[j];
			r_1[j] = p3[j];
			p_2[j] = q1[j];
			q_2[j] = q2[j];
			r_2[j] = q3[j];
		}

		if (!tri_tri_intersection_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0], &coplanar, &s[0], &t[0]))
			return -1;

		if (coplanar == 1) {
			// if(tri_tri_overlap_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0])) //TODO?
			return 4;
		}

		if (s[0] == t[0] && s[1] == t[1] && s[2] == t[2])
			return -1;

		return 3;
	}
	void FastEnvelope::Segment_facet_intersection(std::vector<std::array<Vector3, 2>>& seglist, const std::array<Vector3, 3>& facet, std::vector<Vector3>& interp, std::vector<int>& interseg)
	{
		interp.clear();
		interp.reserve(seglist.size() / 2);
		interseg.clear();
		interseg.reserve(seglist.size() / 2);
		Parameters params;
		int  cutOrnot;
		Vector3  ipoint;
		Scalar  t;
		int  ion, rcd = 0;
		//if (is_tri_tri_cutted(p1, p2, p3, q1, q2, q3)!=3)
		for (int i = 0; i < 3; i++) {
			SLIntersection(params, facet, seglist[i][0], seglist[i][1], cutOrnot, ipoint);
			if (cutOrnot == 1) {
				rcd = rcd + 1;
				interp.emplace_back(ipoint);
				interseg.emplace_back(i);
				if (rcd == 2) {// TODO possible that the segment is an edge,but may not need to fix because this may not slow much
					break;
				}
			}
		}
		if (rcd == 0) {
			interp.clear();
			interseg.clear();
			return;
		}
		if (rcd == 1) {
			interp.clear();
			interseg.clear();
			std::cout << "ERROR ENVELOPE TEST" << std::endl;
			return;
		}
		for (int i = 3; i < seglist.size(); i++) {

			SLIntersection(params, facet, seglist[i][0], seglist[i][1], cutOrnot, ipoint);
			if (cutOrnot == 1) {
				interp.emplace_back(ipoint);
				interseg.emplace_back(i);
			}
		}
		seglist.emplace_back(std::array<Vector3, 2>{ {interp[0], interp[1]} });
	}
	//void FastEnvelope::Segment_facet_intersectionV2(std::vector<std::array<Vector3, 2>>& seglist, const std::array<Vector3, 3>& facet, std::vector<Vector3>& interp, std::vector<int>& interseg) {
	//	interp.clear();
	//	interseg.clear();
	//	Parameters params;
	//	int  cutOrnot, ion, rcd = 0;
	//	Vector3  ipoint;
	//	Scalar  t;
	//	
	//	for (int i = 0; i < 3; i++) {
	//		SLIntersection(params, facet, seglist[i][0], seglist[i][1], cutOrnot, ipoint);
	//		if (cutOrnot == 1) {
	//			rcd = rcd + 1;
	//			interp.emplace_back(ipoint);
	//			interseg.emplace_back(i);
	//			if (rcd == 2) {// TODO possible that the segment is an edge,but may not need to fix because this may not slow much
	//				break;
	//			}
	//		}
	//	}
	//	if (rcd == 0) {
	//		interp.clear();
	//		interseg.clear();
	//		return;
	//	}
	//	if (rcd == 1) {
	//		interp.clear();
	//		interseg.clear();
	//		std::cout << "ERROR ENVELOPE TEST" << std::endl;
	//		return;
	//	}
	//	for (int i = 3; i < seglist.size(); i++) {

	//		SLIntersection(params, facet, seglist[i][0], seglist[i][1], cutOrnot, ipoint);
	//		if (cutOrnot == 1) {
	//			interp.emplace_back(ipoint);
	//			interseg.emplace_back(i);
	//		}
	//	}
	//	seglist.emplace_back(std::array<Vector3, 2>{ {interp[0], interp[1]} });
	//}

	int FastEnvelope::orient_3triangles(const Eigen::Matrix<Scalar, 3, 3>& A, const Eigen::Matrix<Scalar, 3, 3>& AT,
		const Eigen::Matrix<Scalar, 3, 3>& ATA, const Eigen::Matrix<Scalar, 3, 1>& B, const std::array<Vector3, 3> & triangle3) {//right hand law

		Eigen::Matrix<Scalar, 4, 4> C;
		Scalar  ad = det3x3(A(0, 0), A(0, 1), A(0, 2),
			A(1, 0), A(1, 1), A(1, 2),
			A(2, 0), A(2, 1), A(2, 2)), ad1 = 1 / ad, cd;
		/*if (fabs(ad) > 1e8) {
			std::cout << "not right in enlarge A" << std::endl;
		}*/
		Eigen::Matrix<Scalar, 3, 3> A1;
		/*C << AT * B, ATA*triangle3[0], ATA*triangle3[1], ATA*triangle3[2],
			1, 1, 1, 1;*/

			//if (kmark == 3) {
			//	
			//	std::cout << "if  accurate calculate C: " << C.determinant() << std::endl;
			//	C << A.inverse()*B, triangle3[0], triangle3[1], triangle3[2],
			//		1, 1, 1, 1;
			//	std::cout << "if not accurate calculate C: " << C.determinant() << std::endl;
			//	C << B, A*triangle3[0], A*triangle3[1], A*triangle3[2],
			//		1, 1, 1, 1;
			//	std::cout << "if half accurate calculate C: " <<A.determinant()* C.determinant() << std::endl;
			//	std::cout << "A: " << A << std::endl;
			//	std::cout << "A.inverse: " << A.inverse() << std::endl;
			//	std::cout << "A*A-1: " << A*A.inverse() << std::endl;
			//	std::cout << "A.det and A-1.det " << A.determinant()<<" "<<A.inverse().determinant()<< endl;
			//	Eigen::Matrix<Scalar, 3, 1> r;
			//	Eigen::Matrix<Scalar, 1, 3> d;
			//	r << 0, 0, 0;
			//	d << 0, 0, 0;
			//	C << ATA, r,
			//		d, 1;
			//	std::cout << "see the para matrix " << C << endl;
			//	std::cout << "see the para  " << C.determinant() << endl;
			//}//TODO potential fix is to let c.det>eps
			/*C << A.inverse()*B, triangle3[0], triangle3[1], triangle3[2],
				1, 1, 1, 1;
			if (C.determinant() > 0) {
				return 1;
			}
			if (C.determinant() < 0) {
				return -1;
			}
			return 0; */
		C << B, A*triangle3[0], A*triangle3[1], A*triangle3[2],
			1, 1, 1, 1;
		cd = det4x4(C(0, 0), C(0, 1), C(0, 2), C(0, 3),
			C(1, 0), C(1, 1), C(1, 2), C(1, 3),
			C(2, 0), C(2, 1), C(2, 2), C(2, 3),
			C(3, 0), C(3, 1), C(3, 2), C(3, 3));
		if (ad*cd > 0) {
			return 1;
		}
		if (ad*cd < 0) {
			return -1;
		}
		return 0;
		/*A1 = ad1* A;
		C << ad1*B, A1*triangle3[0], A1*triangle3[1], A1*triangle3[2],
			1, 1, 1, 1;
		if (C.determinant() > 0) {
			return 1;
		}
		if (C.determinant() < 0) {
			return -1;
		}
		return 0;*/
		/*Vector3 p;
		p = A.inverse()*B;
		int in = Predicates::orient_3d(triangle3[0], triangle3[1], triangle3[2], p);
		if (in > 0) {
			return -1;
		}
		if (in < 0) {
			return 1;
		}
		return 0;*/

	}


	bool FastEnvelope::point_out_prism(const Vector3 & point, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump)
	{
		//ftimer3.start();
		int jm = 0, ori;

		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}
			}

			for (int j = 0; j < 8; j++) {
				//ftimer2.start();
				ori = Predicates::orient_3d(envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]], point);
				//ftime2 = ftime2 + ftimer2.getElapsedTimeInSec();
				if (ori == -1 || ori == 0) {
					break;
				}
				if (j == 7) {
					/*if (kmark1 == 2) {
						std::cout << "which prism include,i " << i<< "/"<<envprism.size()-1 << std::endl;
					}*/
					//std::cout << "which prism include this point: " << i  << std::endl;
					//ftime3 = ftime3 + ftimer3.getElapsedTimeInSec();
					return false;
				}
			}


		}
		//ftime3 = ftime3 + ftimer3.getElapsedTimeInSec();
		return true;
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
		if (nn.dot(npm) == 0) {//if the line is parallel to the face. TODO check this parameter
			cutOrnot = 0;
			//interp.setZero();
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
				//interp.setZero();
				return;
			}
		}

	}

	void FastEnvelope::BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, const Scalar& bbd)
	{
		envprism.reserve(m_faces.size());
		Vector3 AB, AC, BC, normal, vector1, ABn;
		Parameters pram;
		std::array<Vector3, 6> polygon;
		std::array<Vector3, 12> polygonoff;
		Scalar  a, b, c,
			tolerance = bbd * pram.eps_rel / sqrt(3),
			//tolerance = 2*bbd * pram.epsilon_rel,//TODO this is not right
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
			envprism.emplace_back(polygonoff);

		}

		//std::cout << "epsilon " << tolerance << std::endl;
		//std::cout << "triangle expanded: \n" << striangle[0] <<"\n--------\n"<< striangle[1] <<"\n---------\n"<< striangle[2] << std::endl;
		//std::cout << "prism: \n" << envbox[0][0] << "\n,\n"<<envbox[0][1] << "\n,\n" << envbox[0][2] << "\n,\n" 
		//	<< envbox[0][3] << "\n,\n" << envbox[0][4] << "\n,\n" << envbox[0][5] << "\n,\n" << std::endl;
	}

	void FastEnvelope::timerecord()
	{
		//std::cout << "time look into function: " << ftime1 << std::endl;
		//std::cout << "orientation time " << ftime2 <<"  call times "<< call<< std::endl;
		//std::cout << "point out prism time in total " << ftime3 << std::endl;
		//std::cout << "p_triangle " << p_triangle[0][0][2] << std::endl;
		//std::cout << "p_triangle " << p_triangle[0][1][0] << std::endl;
		//std::cout << "p_triangle " << p_triangle[1][1][2] << std::endl;
	}
	/*
	int FastEnvelope::tri_cut_tri(const Vector3 & p1, const Vector3 & p2, const Vector3 & p3, const Vector3 & q1, const Vector3 & q2, const Vector3 & q3)
		{

			std::array<Scalar, 3> p_1 = { {0, 0, 0} }, q_1 = { {0, 0, 0} }, r_1 = { {0, 0, 0} };
			std::array<Scalar, 3> p_2 = { {0, 0, 0} }, q_2 = { {0, 0, 0} }, r_2 = { {0, 0, 0} };
			int coplanar = 0;
			std::array<Scalar, 3> s = { {0,0,0} }, t = { {0,0,0} };
			for (int j = 0; j < 3; j++) {
				p_1[j] = p1[j];
				q_1[j] = p2[j];
				r_1[j] = p3[j];
				p_2[j] = q1[j];
				q_2[j] = q2[j];
				r_2[j] = q3[j];
			}

			if (!tri_tri_intersection_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0], &coplanar, &s[0], &t[0]))
				return CUT_EMPTY;

			if (coplanar == 1) {
				// if(tri_tri_overlap_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0])) //TODO?
				return CUT_COPLANAR;
			}

			if (s[0] == t[0] && s[1] == t[1] && s[2] == t[2])
				return CUT_EMPTY;

			auto is_collinear = [&](const Vector3 &st, const Vector3 &q12) {
				Vector3 cross = st.cross(q12);
				if (cross[0] <= SCALAR_ZERO && cross[1] <= SCALAR_ZERO && cross[2] <= SCALAR_ZERO)
					return true;
				return false;
			};

			Vector3 st(t[0] - s[0], t[1] - s[1], t[2] - s[2]);
			if (is_collinear(st, q1 - q2))
				return CUT_EDGE_0;
			if (is_collinear(st, q2 - q3))
				return CUT_EDGE_1;
			if (is_collinear(st, q3 - q1))
				return CUT_EDGE_2;

			return CUT_FACE;
		}
	*/
	int FastEnvelope::orient3D_LPI(//right hand law
		const double& px, const double& py, const double& pz,
		const double& ax, const double& ay, const double& az,
		const double& bx, const double& by, const double& bz,
		const double& cx, const double& cy, const double& cz,
		const double& a11, const double& a12, const double& a13,
		const double& a21, const double& a22, const double& a23,
		const double& a31, const double& a32, const double& a33,
		const double& px_rx, const double& py_ry, const double& pz_rz,
		const double& d, const double& n)
	{
		::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);


		double dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
		double ix, iy, iz;
		double m12, m13, m14, m23, m24, m34;
		double m123, m124, m134, m234;
		double det4x4_return_value;


		ix = ((d * px) + (a11 * n));
		iy = ((d * py) + (a12 * n));
		iz = ((d * pz) + (a13 * n));
		dax = d * ax;
		day = d * ay;
		daz = d * az;
		dbx = d * bx;
		dby = d * by;
		dbz = d * bz;
		dcx = d * cx;
		dcy = d * cy;
		dcz = d * cz;
		m12 = ((dax * iy) - (ix * day));
		m13 = ((dbx * iy) - (ix * dby));
		m14 = ((dcx * iy) - (ix * dcy));
		m23 = ((dbx * day) - (dax * dby));
		m24 = ((dcx * day) - (dax * dcy));
		m34 = ((dcx * dby) - (dbx * dcy));
		m123 = (((m23 * iz) - (m13 * daz)) + (m12 * dbz));
		m124 = (((m24 * iz) - (m14 * daz)) + (m12 * dcz));
		m134 = (((m34 * iz) - (m14 * dbz)) + (m13 * dcz));
		m234 = (((m34 * daz) - (m24 * dbz)) + (m23 * dcz));
		det4x4_return_value = (m234 - m134 + m124 - m123);
		if (::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID))
			return 0; // Fast reject in case of under/overflow


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

		double eps = 3.5376137154540446e-011 * max6 * max7 * max4 * max1 *
			max6 * max7 * max4 * max3 * max3 * max8 * max5 * max2;
		/////////////////////////////////
		if (markhf == 0) {
			if (rcdn == 0) {
				fout6.open("D:\\vs\\float project\\data\\output_0522\\eps_and_det.txt");
			}
			if (det4x4_return_value == 0) {
				/*std::cout << "zero matrix: " << d * px + n * a11 << " " << d * py + n * a12 << " " << d * pz + n * a13 << "\n"
					<< ax << " " << ay << " " << az << "\n"
					<< bx << " " << by << " " << bz << "\n"
					<< cx << " " << cy << " " << cz << "\n";*/
					/*std::cout << "every eps: " << eps << std::endl;
					std::cout << "det_4x4: " << det4x4_return_value << std::endl;*/
			}
			fout6 << eps << " " << det4x4_return_value << std::endl;
			std::cout << "every eps: " << eps << std::endl;
			std::cout << "det_4x4: " << det4x4_return_value << std::endl;
			rcdn = rcdn + 1;
			if (rcdn == 10000) {
				fout6.close();
			}

		}
		////////////////////////////////
		if ((det4x4_return_value > eps)) return (d > 0) ? (1) : (-1);
		if ((det4x4_return_value < -eps)) return (d > 0) ? (-1) : (1);
		//if ((det4x4_return_value > 0)) return (d > 0) ? (1) : (-1);//changed
		//if ((det4x4_return_value < 0)) return (d > 0) ? (-1) : (1);//changed
		return 0;
	}
	
	
	
}



	
	

	
		
	



	
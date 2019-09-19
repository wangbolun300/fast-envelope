#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Parameters.h>
#include <fastenvelope/Predicates.hpp>
#include <fstream>
#include <istream>
#include <igl/Timer.h>
#include <fastenvelope/ip_filtered.h>
#include <arbitraryprecision/fprecision.h>
#include<fastenvelope/mesh_AABB.h>
#include <geogram/mesh/mesh_reorder.h>
#include <igl/Timer.h>

static igl::Timer timer, timer_bigpart,timer_s,timerc,timerdetail,timer_u,timer_a;
static double time_multi = 0.0, time_p1 = 0.0, time_p2 = 0.0, time_p3 = 0.0, time_pm = 0.0, time_searching = 0.0, time_checking = 0.0,
time_p3d = 0,time_p3m=0,timein1=0,timein2=0,timein3=0,timein4=0,timetpp1=0,timetpp2 = 0,timetpp3 = 0,timetpp4 = 0,
timeinit1 = 0,timeinit2 = 0,timeinit3 = 0,timeinit4 = 0,timeinit5=0;
static const int p_face[8][3] = { {0,1,3},{7,6,9},{1,0,7},{2,1,7},{3,2,8},{3,9,10},{5,4,11},{0,5,6} };//prism triangle index. all with orientation.
static const int c_face[6][3] = { {0,1,2},{4,7,6},{0,3,4},{1,0,4},{1,5,2},{2,6,3} };
static const fastEnvelope::Vector3 origin = fastEnvelope::Vector3(0, 0, 0);
static const std::array<std::vector<int>, 8> p_facepoint = {
		{
			{0,1,2,3,4,5},
	{8,7,6,11,10,9},
	{7,1,0,6},
	{2,1,7,8},
	{3,2,8,9},
	{4,3,9,10},
	{4,10,11,5},
	{6,0,5,11}

		}
};
static const std::array<std::array<int, 4>, 6> c_facepoint = {
		{
			{{0,1,2,3}},
	{{4,7,6,5}},
	{{4,0,3,7}},
	{{1,0,4,5}},
	{{2,1,5,6}},
	{{3,2,6,7}}


		}
};
static const int p_facenumber = 8;
static const int c_facenumber = 6;
static const double conserve_number = 0.5;
static const int prism_map[64][2] = {
{-1, -1},
 {-1, -1},
 {0, 1},
 {1, 2},
 {2, 3},
 {3, 4},
 {4, 5},
 {0, 5},
 {-1, -1},
 {-1, -1},
 {6, 7},
 {7, 8},
 {8, 9},
 {9, 10},
 {10, 11},
 {6, 11},
 {0, 1},
 {6, 7},
 {-1, -1},
 {1, 7},
 {-1, -1},
 {-1, -1},
 {-1, -1},
 {0, 6},
 {1, 2},
 {7, 8},
 {1, 7},
 {-1, -1},
 {2, 8},
 {-1, -1},
 {-1, -1},
 {-1, -1},
 {2, 3},
 {8, 9},
 {-1, -1},
 {2, 8},
 {-1, -1},
 {3, 9},
 {-1, -1},
 {-1, -1},
 {3, 4},
 {9, 10},
 {-1, -1},
 {-1, -1},
 {3, 9},
 {-1, -1},
 {4, 10},
 {-1, -1},
 {4, 5},
 {10, 11},
 {-1, -1},
 {-1, -1},
 {-1, -1},
 {4, 10},
 {-1, -1},
 {5, 11},
 {0, 5},
 {6, 11},
 {0, 6},
 {-1, -1},
 {-1, -1},
 {-1, -1},
 {5, 11},
 {-1, -1}
};

static const int cubic_map[36][2] = {
	{-1, -1},
 {-1, -1},
 {0, 3},
 {0, 1},
 {1, 2},
 {2, 3},
 {-1, -1},
 {-1, -1},
 {4, 7},
 {4, 5},
 {5, 6},
 {6, 7},
 {0, 3},
 {4, 7},
 {-1, -1},
 {0, 4},
 {-1, -1},
 {3, 7},
 {0, 1},
 {4, 5},
 {0, 4},
 {-1, -1},
 {1, 5},
 {-1, -1},
 {1, 2},
 {5, 6},
 {-1, -1},
 {1, 5},
 {-1, -1},
 {2, 6},
 {2, 3},
 {6, 7},
 {3, 7},
 {-1, -1},
 {2, 6},
 {-1, -1}
};


static const   std::function<int(double)> check_double = [](double v) {
	if (fabs(v) < 1e-10)
		return -2;

	if (v > 0)
		return 1;

	if (v < 0)
		return -1;

	return 0;
};

static const   std::function<int(fastEnvelope::Rational)> check_Rational = [](fastEnvelope::Rational v) {

	return v.get_sign();
	// if (v > 0)
	// 	return 1;

	// if (v < 0)
	// 	return -1;

	// return 0;
};

static const   std::function<int(fastEnvelope::Multiprecision)> check_Multiprecision = [](fastEnvelope::Multiprecision v) {

	return v.get_sign();
	// if (v > 0)
	// 	return 1;

	// if (v < 0)
	// 	return -1;

	// return 0;
};




static const std::array<std::array<int, 2>, 3> triseg = {

	{{{0,1}},{{0,2}},{{1,2}}}

};



namespace fastEnvelope {
	//using namespace std;

	void FastEnvelope::print_number() {
		//std::cout << "lpi filter number " << filternumberlpi << " lpi total number " << totalnumberlpi << " percentage " << float(filternumberlpi )/ float(totalnumberlpi) << std::endl;
		//std::cout << "tpi filter number " << filternumbertpi << " tpi total number " << totalnumbertpi << " percentage " << float(filternumbertpi) / float(totalnumbertpi) << std::endl;
		//std::cout << "triangle_intersection filter number " << filternumber1 << " tpi total number " << totalnumber1 << " percentage " << float(filternumber1) / float(totalnumber1) << std::endl;
		//std::cout << "triangle_intersection filter number lpi -2 " << filternumberlpi2 << " percentage " << float(filternumberlpi2) / float(totalnumberlpi) << std::endl;
		//std::cout << "triangle_intersection filter number tpi -2 " << filternumbertpi2 << " percentage " << float(filternumbertpi2) / float(totalnumbertpi+ totalnumber1) << std::endl;
		/*std::cout << "lpi 1 " << float(after11) / float(after11 + after12 + after10) << " lpi -1 " << after12 / float(after11 + after12 + after10) << " lpi 0 " << after10 / float(after11 + after12 + after10) << " tot  " << after11 + after12 + after10 << std::endl;
		std::cout << "tpi 1 " << after21 / float(after21 + after22 + after20) << " tpi -1 " << after22 / float(after21 + after22 + after20) << " tpi 0 " << after20 / float(after21 + after22 + after20) << " tot  " << after21 + after22 + after20 << std::endl;
		std::cout << "go1 " << go1 << " go2 " << go2 << std::endl;*/
		//std::cout << "same " << float(diff1) / float(diff1 + diff2 + diff3) << " diff " << float(diff2) / float(diff1 + diff2 + diff3) << " wrong " << float(diff3) / float(diff1 + diff2 + diff3) << std::endl;
	//	std::cout << "cut tri number original " << ct1 << " conservative " << ct2  <<" rate "<<float(ct1)/float(ct2)<< std::endl;
		//std::cout << "total " <<diff1+diff2+diff3 << "   " << ct1 << "  " << ct2 << std::endl;
		std::cout << "time in multi, " << time_multi << std::endl;
		std::cout << "time part 1, " << time_p1<<"\ntime part 2, "<<time_p2<<"\ntime part 3, "<<time_p3 << std::endl;
		std::cout << "time searching, " << time_searching << std::endl;
		std::cout << "time detail in part 3_1, " << time_p3d << "\ntime detail in part 3_2, " << time_p3m <<  std::endl;
		std::cout << "time in part 3 which part? 1, " << timein1 << "\ntime in part 3 which part? 2, " << timein2 << "\ntime in part 3 which part? 3, "
			<< timein3 << "\ntime in part 3 which part? 4, " << timein4 << " " << std::endl;
		std::cout << "time in part 3 doubt double_1, " << timetpp1 <<"\ntime in part 3 doubt double_2, "<< timetpp2 << " " << std::endl;

	}
	void FastEnvelope::print_ini_number() {
		std::cout << "init time 1, " << timeinit1 << "\ninit time 2, " << timeinit2 << "\ninit time 3, "<<
			timeinit3<< "\ninit time 4, " << timeinit4 << "\ninit time 5, " << timeinit5<< std::endl;
	}
	FastEnvelope::FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps, const int spac)
	{

		//Multiprecision::set_precision(256);
		timer.start();
		Vector3 min, max;
		get_bb_corners(m_ver, min, max);
		timeinit1 += timer.getElapsedTimeInSec();
		Scalar bbd = (max - min).norm();
		Scalar epsilon = bbd * eps; //eps*bounding box diagnal
		timer.start();
		GEO::Mesh M;
		
		to_geogram_mesh(m_ver, m_faces, M);
		GEO::mesh_reorder(M, GEO::MESH_ORDER_MORTON);
		
		std::vector<Vector3> ver_new;
		ver_new.resize(m_ver.size());

		std::vector<Vector3i> faces_new;
		faces_new.resize(m_faces.size());

		from_geogram_mesh(M, ver_new, faces_new);
		timeinit2 += timer.getElapsedTimeInSec();
		timer.start();
		BoxGeneration(ver_new, faces_new, envprism, envcubic, epsilon);
		timeinit3 += timer.getElapsedTimeInSec();
		//build a  hash function
		prism_size = envprism.size();
		timer.start();
		CornerList_prism(envprism, cornerlist);
		std::vector<std::array<Vector3, 2>> cubiconors;
		CornerList_cubic(envcubic, cubiconors);
		cornerlist.insert(cornerlist.end(), cubiconors.begin(), cubiconors.end());
		timeinit4 += timer.getElapsedTimeInSec();
		timer.start();
		tree.init_envelope(cornerlist, false);
		timeinit5 += timer.getElapsedTimeInSec();
		std::cout << "prism size " << prism_size << std::endl;
		std::cout << "cubic size " << envcubic.size() << std::endl;
	}
	bool FastEnvelope::is_outside(const std::array<Vector3, 3> &triangle) const {
		timer_s.start();

		 std::vector<unsigned int> querylist;
		 //tree.facet_in_envelope(triangle[0], triangle[1], triangle[2], querylist);
		 Vector3 tmin, tmax;
		 get_tri_corners(triangle[0], triangle[1], triangle[2], tmin, tmax);
		 for (int i = 0; i < cornerlist.size(); i++) {
			 if (box_box_intersection(tmin, tmax, cornerlist[i][0], cornerlist[i][1])) querylist.emplace_back(i);
		 }
		time_searching += timer_s.getElapsedTimeInSec();

		return FastEnvelopeTestImplicit(triangle, querylist);
	}

	void FastEnvelope::print_prisms(const std::array<Vector3, 3> &triangle) const {
		
		std::vector<unsigned int> querylist;
		tree.facet_in_envelope(triangle[0], triangle[1], triangle[2], querylist);
		std::ofstream fout;
		fout.open("D:\\vs\\fast_envelope_csv\\problems\\visualtriangle.txt");
		
		for (int i = 0; i < 3; i++) {

			fout << std::setprecision(17) << triangle[i][0] << " " << triangle[i][1] << " " << triangle[i][2] << std::endl;

		}
		fout.close();
		fout.open("D:\\vs\\fast_envelope_csv\\problems\\prisms.txt");
		for (int i = 0; i < querylist.size(); i++) {
			if (querylist[i] >= prism_size)   continue;
			for (int j = 0; j < 12; j++) {
				fout << std::setprecision(17) << envprism[querylist[i]][j][0] << " " << envprism[querylist[i]][j][1] << " " << envprism[querylist[i]][j][2] << std::endl;

			}
		}
		fout.close();
		fout.open("D:\\vs\\fast_envelope_csv\\problems\\cubes.txt");
		for (int i = 0; i < querylist.size(); i++) {
			if (querylist[i] < prism_size) continue;
			for (int j = 0; j < 8; j++) {
				fout << std::setprecision(17) << envcubic[querylist[i]][j][0] << " " << envcubic[querylist[i]][j][1] << " " << envcubic[querylist[i]][j][2] << std::endl;

			}
		}
		fout.close();
	}
	bool FastEnvelope::sample_triangle_outside(const std::array<Vector3, 3> &triangle, const int& pieces) const {

		bool flagr = 0;
		bool out;
		Vector3  point;
		std::vector<unsigned int> querylist;
		tree.facet_in_envelope(triangle[0], triangle[1], triangle[2], querylist);
		int jump = -1;
		if (querylist.size() == 0) return 1;
		if (flagr == 0) {
			int deg = is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

			
			if (deg == DEGENERATED_POINT) {
				out = point_out_prism(triangle[0], querylist, jump);
				if (out == true) {

					return 1;

				}
				return 0;
			}
			if (deg == DEGENERATED_SEGMENT) {
				for (int i = 0; i < pieces; i++) {
					triangle_sample_segment(triangle, point, pieces, i);
					out = point_out_prism(point, querylist, jump);
					if (out == true) {

						return 1;

					}

				}
				return 0;
			}


			for (int i = 0; i < pieces; i++) {
				for (int j = 0; j <= i; j++) {
					triangle_sample_normal(triangle, point, pieces, i, j);
					out = point_out_prism(point, querylist, jump);
					if (out == true) {

						return 1;

					}


				}

			}
		}
		else {
			//std::cout << "*     using rational" << std::endl;
			Rational ps0, ps1, ps2;
			for (int i = 0; i < pieces; i++) {
				for (int j = 0; j <= i; j++) {
					triangle_sample_normal_rational(triangle, ps0, ps1, ps2,
						pieces, i, j);
					out = point_out_prism_rational( ps0, ps1, ps2, querylist, jump);
					if (out == true) {

						return 1;

					}


				}

			}


			
		}
		return 0;
	}

	void FastEnvelope::triangle_sample_segment(const std::array<Vector3, 3> &triangle, Vector3& ps, const int &pieces, const int & nbr) {

		int t = pieces-1;
		if (triangle[1] - triangle[0] == Vector3(0, 0, 0)) {

			ps = (triangle[0] + (triangle[2] - triangle[0])*nbr / t);

			return;
		}
		if (triangle[2] - triangle[0] == Vector3(0, 0, 0)) {


			ps = (triangle[0] + (triangle[1] - triangle[0])*nbr / t);

			return;
		}
		if (triangle[2] - triangle[1] == Vector3(0, 0, 0)) {

			ps = (triangle[0] + (triangle[1] - triangle[0])*nbr / t);

			return;
		}

		Scalar d1 = (triangle[1] - triangle[0]).norm(), d2 = (triangle[2] - triangle[0]).norm(), d3 = (triangle[1] - triangle[2]).norm();
		if (d1 >= d2 && d1 >= d3) {
			ps = (triangle[0] + (triangle[1] - triangle[0])*nbr / t);

			return;
		}
		if (d2 >= d1 && d2 >= d3) {
			ps = (triangle[0] + (triangle[2] - triangle[0])*nbr / t);

			return;
		}
		if (d3 >= d1 && d3 >= d2) {
			ps = (triangle[1] + (triangle[2] - triangle[1])*nbr / t);

			return;
		}
	}
	void  FastEnvelope::triangle_sample_point(const std::array<Vector3, 3> &triangle, Vector3& ps) {
		ps = triangle[0];
	}
	void FastEnvelope::triangle_sample_normal(const std::array<Vector3, 3> &triangle,Vector3& ps, const int &pieces, const int & nbr1, const int &nbr2) {
		int l1s = pieces - 1;//
		Vector3 p1 = triangle[0] + (triangle[1] - triangle[0])*nbr1 / l1s, d = (triangle[2] - triangle[1]) / l1s;
		ps = p1 + d * nbr2;

	}
	void FastEnvelope::triangle_sample_normal_rational(const std::array<Vector3, 3> &triangle, Rational& ps0, Rational& ps1, Rational& ps2, const int &pieces, const int & nbr1, const int &nbr2) {
		int l1s = pieces - 1;//
		Rational t00(triangle[0][0]), t01(triangle[0][1]), t02(triangle[0][2]), t10(triangle[1][0]), t11(triangle[1][1]), 
			t12(triangle[1][2]), t20(triangle[2][0]), t21(triangle[2][1]), t22(triangle[2][2]), nbr1r(nbr1), nbr2r(nbr2), l1sr(l1s);
		
		Rational  p0 = t00 + (t10 - t00)*nbr1r / l1sr, d0 = (t20 - t10) / l1sr;
		Rational  p1 = t01 + (t11 - t01)*nbr1r / l1sr, d1 = (t21 - t11) / l1sr;
		Rational  p2 = t02 + (t12 - t02)*nbr1r / l1sr, d2 = (t22 - t12) / l1sr;
		ps0 = p0 + d0 * nbr2;
		ps1 = p1 + d1 * nbr2;
		ps2 = p2 + d2 * nbr2;

	}
	struct DATA_LPI {
		int segid;
		int prismid;
		int facetid;
		int jump1;
	};
	struct DATA_TPI {
		int prismid1;
		int facetid1;
		int prismid2;
		int facetid2;
		int jump1;
		int jump2;
	};
	bool FastEnvelope::FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int>& prismindex)const

	{

		static const std::function<int(fastEnvelope::Multiprecision)> checker = check_Multiprecision;
		static const std::function<int(fastEnvelope::Multiprecision)> checker1 = check_Multiprecision;
		if (prismindex.size() == 0) {

			return 1;

		}

		int jump1, jump2;

		std::vector<std::array<int,2>> inter_ijk_list;//list of intersected triangle

		bool out,cut;

		int inter, inter1, record1, record2,

			tti;//triangle-triangle intersection


		jump1 = -1;
		for (int i = 0; i < 3; i++) {

			out = point_out_prism(triangle[i], prismindex, jump1);

			if (out == true) {

				return 1;

			}

		}

		if (prismindex.size() == 1) return 0;



		////////////////////degeneration fix
		timer_bigpart.start();
		int degeneration = is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

		if (degeneration == DEGENERATED_POINT) {//case 1 degenerate to a point

			return 0;

		}//case 1 degenerate to a point
		DATA_LPI datalpi;
		std::vector<DATA_LPI> lpi_list;
		if (degeneration == DEGENERATED_SEGMENT) {

			for (int we = 0; we < 3; we++) {//case 2 degenerated as a segment, at most test 2 segments,but still we need to test 3, because

											// of the endpoint-triangle intersection will be ignored

											// the segment is {triangle[triseg[we][0]], triangle[triseg[we][1]]}



				for (int i = 0; i < prismindex.size(); i++) {
					jump1 = prismindex[i];
					if (prismindex[i] < prism_size) {
						std::vector<int>cid;
						cut = is_seg_cut_prism(jump1, triangle[triseg[we][0]], triangle[triseg[we][1]], cid);
						if (cut == true) {
							for (int j = 0; j < cid.size(); j++) {
								Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
								bool precom = ip_filtered::orient3D_LPI_prefilter(//
									triangle[triseg[we][0]][0], triangle[triseg[we][0]][1], triangle[triseg[we][0]][2],
									triangle[triseg[we][1]][0], triangle[triseg[we][1]][1], triangle[triseg[we][1]][2],
									envprism[prismindex[i]][p_face[cid[j]][0]][0], envprism[prismindex[i]][p_face[cid[j]][0]][1], envprism[prismindex[i]][p_face[cid[j]][0]][2],
									envprism[prismindex[i]][p_face[cid[j]][1]][0], envprism[prismindex[i]][p_face[cid[j]][1]][1], envprism[prismindex[i]][p_face[cid[j]][1]][2],
									envprism[prismindex[i]][p_face[cid[j]][2]][0], envprism[prismindex[i]][p_face[cid[j]][2]][1], envprism[prismindex[i]][p_face[cid[j]][2]][2],
									a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
								if (precom == true) {
									inter = Implicit_Seg_Facet_interpoint_Out_Prism_double(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
										triangle[triseg[we][0]], triangle[triseg[we][1]],
										envprism[prismindex[i]][p_face[cid[j]][0]], envprism[prismindex[i]][p_face[cid[j]][1]], envprism[prismindex[i]][p_face[cid[j]][2]],
										prismindex, jump1, checker);
									if (inter == 1) {

										return 1;
									}
								}
								else {
									datalpi.segid = we;
									datalpi.prismid = prismindex[i];
									datalpi.facetid = cid[j];
									datalpi.jump1 = jump1;
									lpi_list.emplace_back(datalpi);
								}


							}
						}
					}
					else {
						int cindex = prismindex[i] - prism_size;
						std::vector<int>cid;
						cut = is_seg_cut_cube(cindex, triangle[triseg[we][0]], triangle[triseg[we][1]], cid);
						if (cut == true) {
							for (int j = 0; j < cid.size(); j++) {
								Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
								bool precom = ip_filtered::orient3D_LPI_prefilter(//
									triangle[triseg[we][0]][0], triangle[triseg[we][0]][1], triangle[triseg[we][0]][2],
									triangle[triseg[we][1]][0], triangle[triseg[we][1]][1], triangle[triseg[we][1]][2],
									envcubic[cindex][c_face[cid[j]][0]][0], envcubic[cindex][c_face[cid[j]][0]][1], envcubic[cindex][c_face[cid[j]][0]][2],
									envcubic[cindex][c_face[cid[j]][1]][0], envcubic[cindex][c_face[cid[j]][1]][1], envcubic[cindex][c_face[cid[j]][1]][2],
									envcubic[cindex][c_face[cid[j]][2]][0], envcubic[cindex][c_face[cid[j]][2]][1], envcubic[cindex][c_face[cid[j]][2]][2],
									a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
								if (precom == true) {
									inter = Implicit_Seg_Facet_interpoint_Out_Prism_double(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
										triangle[triseg[we][0]], triangle[triseg[we][1]],
										envcubic[cindex][c_face[cid[j]][0]], envcubic[cindex][c_face[cid[j]][1]], envcubic[cindex][c_face[cid[j]][2]],
										prismindex, jump1, checker);
									if (inter == 1) {

										return 1;
									}
								}
								else {
									datalpi.segid = we;
									datalpi.prismid = prismindex[i];
									datalpi.facetid = cid[j];
									datalpi.jump1 = jump1;
									lpi_list.emplace_back(datalpi);
								}


							}
						}


					}


				}

			}//case 2 case 2 degenerated as a segment

			for (int i = 0; i < lpi_list.size(); i++) {
				inter = Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(lpi_list[i], triangle, prismindex,checker);
				if (inter == 1) return 1;
			}
			return 0;

		}
		//
		////////////////////////////////degeneration fix over
		time_p1 += timer_bigpart.getElapsedTimeInSec();
		timer_bigpart.start();
		for (int i = 0; i < prismindex.size(); i++) {
			jump1 = prismindex[i];

			if (prismindex[i] < prism_size) {

				std::vector<int> cidl;
				cut = is_triangle_cut_prism(prismindex[i],
					triangle[0], triangle[1], triangle[2], cidl);
				//if (cut) std::cout << "cut_id from 1: "<<i+1<<" out of "<< prismindex.size() << std::endl;
				if (cut == false) continue;

				for (int j = 0; j < cidl.size(); j++) {
					for (int k = 0; k < 3; k++) {
						tti = seg_cut_plane(triangle[triseg[k][0]], triangle[triseg[k][1]],
							envprism[prismindex[i]][p_face[cidl[j]][0]],
							envprism[prismindex[i]][p_face[cidl[j]][1]],
							envprism[prismindex[i]][p_face[cidl[j]][2]]);
						if (tti != CUT_FACE) continue;
						Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
						bool precom = ip_filtered::orient3D_LPI_prefilter(//
							triangle[triseg[k][0]][0], triangle[triseg[k][0]][1], triangle[triseg[k][0]][2],
							triangle[triseg[k][1]][0], triangle[triseg[k][1]][1], triangle[triseg[k][1]][2],
							envprism[prismindex[i]][p_face[cidl[j]][0]][0], envprism[prismindex[i]][p_face[cidl[j]][0]][1], envprism[prismindex[i]][p_face[cidl[j]][0]][2],
							envprism[prismindex[i]][p_face[cidl[j]][1]][0], envprism[prismindex[i]][p_face[cidl[j]][1]][1], envprism[prismindex[i]][p_face[cidl[j]][1]][2],
							envprism[prismindex[i]][p_face[cidl[j]][2]][0], envprism[prismindex[i]][p_face[cidl[j]][2]][1], envprism[prismindex[i]][p_face[cidl[j]][2]][2],
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
						if (precom == true) {
							inter = Implicit_Seg_Facet_interpoint_Out_Prism_double(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
								triangle[triseg[k][0]], triangle[triseg[k][1]],
								envprism[prismindex[i]][p_face[cidl[j]][0]], envprism[prismindex[i]][p_face[cidl[j]][1]], envprism[prismindex[i]][p_face[cidl[j]][2]],
								prismindex, jump1, checker);
							if (inter == 1) {

								return 1;
							}
						}
						else {
							datalpi.segid = k;
							datalpi.prismid = prismindex[i];
							datalpi.facetid = cidl[j];
							datalpi.jump1 = jump1;
							lpi_list.emplace_back(datalpi);
						}

					}
					inter_ijk_list.push_back({ {int(prismindex[i]), cidl[j]} });
				}

			}
			else {
				int cindex = prismindex[i] - prism_size;
				std::vector<int> cidl;
				cut = is_triangle_cut_cube(cindex,
					triangle[0], triangle[1], triangle[2], cidl);
				if (cut == false) continue;
				for (int j = 0; j < cidl.size(); j++) {
					for (int k = 0; k < 3; k++) {
						tti = seg_cut_plane(triangle[triseg[k][0]], triangle[triseg[k][1]],
							envcubic[cindex][c_face[cidl[j]][0]],
							envcubic[cindex][c_face[cidl[j]][1]],
							envcubic[cindex][c_face[cidl[j]][2]]);
						if (tti != CUT_FACE) continue;
						Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
						bool precom = ip_filtered::orient3D_LPI_prefilter(//
							triangle[triseg[k][0]][0], triangle[triseg[k][0]][1], triangle[triseg[k][0]][2],
							triangle[triseg[k][1]][0], triangle[triseg[k][1]][1], triangle[triseg[k][1]][2],
							envcubic[cindex][c_face[cidl[j]][0]][0], envcubic[cindex][c_face[cidl[j]][0]][1], envcubic[cindex][c_face[cidl[j]][0]][2],
							envcubic[cindex][c_face[cidl[j]][1]][0], envcubic[cindex][c_face[cidl[j]][1]][1], envcubic[cindex][c_face[cidl[j]][1]][2],
							envcubic[cindex][c_face[cidl[j]][2]][0], envcubic[cindex][c_face[cidl[j]][2]][1], envcubic[cindex][c_face[cidl[j]][2]][2],
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
						if (precom == true) {
							inter = Implicit_Seg_Facet_interpoint_Out_Prism_double(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
								triangle[triseg[k][0]], triangle[triseg[k][1]],
								envcubic[cindex][c_face[cidl[j]][0]], envcubic[cindex][c_face[cidl[j]][1]], envcubic[cindex][c_face[cidl[j]][2]],
								prismindex, jump1, checker);
							if (inter == 1) {

								return 1;
							}
						}
						else {
							datalpi.segid = k;
							datalpi.prismid = prismindex[i];
							datalpi.facetid = cidl[j];
							datalpi.jump1 = jump1;
							lpi_list.emplace_back(datalpi);
						}

					}
					inter_ijk_list.push_back({ {int(prismindex[i]), cidl[j]} });
				}

			}


		}
		for (int i = 0; i < lpi_list.size(); i++) {
			inter = Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(lpi_list[i], triangle, prismindex, checker);
			if (inter == 1) return 1;
		}//TODO consider this part put here or the end of the algorithm

		time_p2 += timer_bigpart.getElapsedTimeInSec();
		timer_bigpart.start();
		timerdetail.start();
		int listsize = inter_ijk_list.size();
		DATA_TPI datatpi;
		std::vector<DATA_TPI> tpilist;
		tpilist.reserve(listsize/5);


		int id, id0 = 0;
		for (int i = 0; i < listsize; i++) {
			jump1 = inter_ijk_list[i][0];
			for (int j = i + 1; j < listsize; j++) {

				//check triangle{ { envprism[list[i][0]][p_triangle[list[i][1]][list[i][2]][0]], ...[1],...[2] } } and triangle{ { envprism[list[j][0]][p_triangle[list[j][1]][list[j][2]][0]], ...[1],...[2] } }

				//and T
				if (inter_ijk_list[i][0] == inter_ijk_list[j][0]) {
					/*if (inter_ijk_list[i][0] < prism_size) {
						id = inter_ijk_list[i][1] * 8 + inter_ijk_list[j][1];
						id0 = prism_map[id][0];

					}
					else {
						id = inter_ijk_list[i][1] * 6 + inter_ijk_list[j][1];
						id0 = cubic_map[id][0];
					}
					if (id0 == -1) continue;*/
					//TODO temp change
					continue;

				}


					//find prism_map[list[i][1]*8+list[j][1]][0],prism_map[list[i][1]*8+list[j][1]][1]
					timer_u.start();
					Vector3 t00, t01, t02, t10, t11, t12;
					int n1, n2;
					if (inter_ijk_list[i][0] < prism_size) {
						n1 = p_facenumber;
					}
					else {
						n1 = c_facenumber;
					}
					if (inter_ijk_list[j][0] < prism_size) {
						n2 = p_facenumber;
					}
					else {
						n2 = c_facenumber;
					}




					jump2 = inter_ijk_list[j][0];

					if (n1 == p_facenumber) {
						t00 = envprism[inter_ijk_list[i][0]][p_face[inter_ijk_list[i][1]][0]];
						t01 = envprism[inter_ijk_list[i][0]][p_face[inter_ijk_list[i][1]][1]];
						t02 = envprism[inter_ijk_list[i][0]][p_face[inter_ijk_list[i][1]][2]];
					}
					else {
						t00 = envcubic[inter_ijk_list[i][0] - prism_size][c_face[inter_ijk_list[i][1]][0]];
						t01 = envcubic[inter_ijk_list[i][0] - prism_size][c_face[inter_ijk_list[i][1]][1]];
						t02 = envcubic[inter_ijk_list[i][0] - prism_size][c_face[inter_ijk_list[i][1]][2]];
					}
					if (n2 == p_facenumber) {
						t10 = envprism[inter_ijk_list[j][0]][p_face[inter_ijk_list[j][1]][0]];
						t11 = envprism[inter_ijk_list[j][0]][p_face[inter_ijk_list[j][1]][1]];
						t12 = envprism[inter_ijk_list[j][0]][p_face[inter_ijk_list[j][1]][2]];
					}
					else {
						t10 = envcubic[inter_ijk_list[j][0] - prism_size][c_face[inter_ijk_list[j][1]][0]];
						t11 = envcubic[inter_ijk_list[j][0] - prism_size][c_face[inter_ijk_list[j][1]][1]];
						t12 = envcubic[inter_ijk_list[j][0] - prism_size][c_face[inter_ijk_list[j][1]][2]];
					}

					timein1 += timer_u.getElapsedTimeInSec();
					timer_u.start();
					Scalar d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7;
					bool multiflag;

					bool pre = ip_filtered::orient3D_TPI_prefilter(triangle[0][0], triangle[0][1], triangle[0][2],
						triangle[1][0], triangle[1][1], triangle[1][2], triangle[2][0], triangle[2][1], triangle[2][2],
						t00[0], t00[1], t00[2],
						t01[0], t01[1], t01[2],
						t02[0], t02[1], t02[2],
						t10[0], t10[1], t10[2],
						t11[0], t11[1], t11[2],
						t12[0], t12[1], t12[2],
						d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7);
					timein2 += timer_u.getElapsedTimeInSec();
					if (pre == true) {

						static Multiprecision dr, n1r, n2r, n3r;
						cut =is_3_triangle_cut_double(d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7,
							triangle, t00, t01, t02, t10, t11, t12, multiflag, check_Multiprecision, dr, n1r, n2r, n3r);
						if (cut == false) continue;
						timer_u.start();
						inter= Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double(//TODO takes most of time
							d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7,
							triangle, t00, t01, t02, t10, t11, t12,
							prismindex, jump1, jump2, multiflag, check_Multiprecision, dr, n1r, n2r, n3r);
						timein3 += timer_u.getElapsedTimeInSec();
						if (inter == 1) {

							return 1;
						}

					}
					else {
						timer_u.start();
						datatpi.prismid1 = inter_ijk_list[i][0];
						datatpi.facetid1 = inter_ijk_list[i][1];
						datatpi.prismid2 = inter_ijk_list[j][0];
						datatpi.facetid2 = inter_ijk_list[j][1];
						datatpi.jump1 = jump1;
						datatpi.jump2 = jump2;
						tpilist.emplace_back(datatpi);
						timein4 += timer_u.getElapsedTimeInSec();
					}

			}

		}
		time_p3d += timerdetail.getElapsedTimeInSec();
		timerdetail.start();
		for (int i = 0; i < tpilist.size(); i++) {
			inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(tpilist[i], triangle, prismindex, checker);//is_3_intersection is already in it
			if (inter == 1) return 1;
		}
		time_p3m += timerdetail.getElapsedTimeInSec();
		time_p3 += timer_bigpart.getElapsedTimeInSec();

		return 0;

	}




	struct INDEX {
		int Pi;
		std::vector<int> FACES;
	};

	template<typename T>
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_LPI& datalpi, const std::array<Vector3, 3>&triangle, const std::vector<unsigned int>& prismindex, const std::function<int(T)> &checker)const {
		int tot,ori;
		static T
			s00, s01, s02, s10, s11, s12,
			t00, t01, t02,
			t10, t11, t12,
			t20, t21, t22,
			a11r, a12r, a13r, dr,

			e00, e01, e02,
			e10, e11, e12,
			e20, e21, e22;
		if (datalpi.prismid < prism_size) {
			t00 = envprism[datalpi.prismid][p_face[datalpi.facetid][0]][0];
			t01 = envprism[datalpi.prismid][p_face[datalpi.facetid][0]][1];
			t02 = envprism[datalpi.prismid][p_face[datalpi.facetid][0]][2];

			t10 = envprism[datalpi.prismid][p_face[datalpi.facetid][1]][0];
			t11 = envprism[datalpi.prismid][p_face[datalpi.facetid][1]][1];
			t12 = envprism[datalpi.prismid][p_face[datalpi.facetid][1]][2];

			t20 = envprism[datalpi.prismid][p_face[datalpi.facetid][2]][0];
			t21 = envprism[datalpi.prismid][p_face[datalpi.facetid][2]][1];
			t22 = envprism[datalpi.prismid][p_face[datalpi.facetid][2]][2];
		}
		else {
			t00 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][0]][0];
			t01 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][0]][1];
			t02 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][0]][2];
			t10 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][1]][0];
			t11 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][1]][1];
			t12 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][1]][2];
			t20 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][2]][0];
			t21 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][2]][1];
			t22 = envcubic[datalpi.prismid-prism_size][c_face[datalpi.facetid][2]][2];

		}


		s00 = triangle[triseg[datalpi.segid][0]][0];
		s01 = triangle[triseg[datalpi.segid][0]][1];
		s02 = triangle[triseg[datalpi.segid][0]][2];

		s10 = triangle[triseg[datalpi.segid][1]][0];
		s11 = triangle[triseg[datalpi.segid][1]][1];
		s12 = triangle[triseg[datalpi.segid][1]][2];
		timer.start();
		bool premulti = orient3D_LPI_prefilter_multiprecision(s00, s01, s02, s10, s11, s12,
			t00, t01, t02, t10, t11, t12, t20, t21, t22, a11r, a12r, a13r, dr, checker);
		time_multi += timer.getElapsedTimeInSec();

		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == datalpi.jump1) {

				continue;
			}

			if (prismindex[i] < prism_size) {

				tot = 0;
				for (int j = 0; j < p_facenumber; j++) {


					e00 = envprism[prismindex[i]][p_face[j][0]][0]; e01 = envprism[prismindex[i]][p_face[j][0]][1]; e02 = envprism[prismindex[i]][p_face[j][0]][2];
					e10 = envprism[prismindex[i]][p_face[j][1]][0]; e11 = envprism[prismindex[i]][p_face[j][1]][1]; e12 = envprism[prismindex[i]][p_face[j][1]][2];
					e20 = envprism[prismindex[i]][p_face[j][2]][0]; e21 = envprism[prismindex[i]][p_face[j][2]][1]; e22 = envprism[prismindex[i]][p_face[j][2]][2];
					timer.start();
					ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22, checker);
					time_multi += timer.getElapsedTimeInSec();
					if (ori == 1 || ori == 0) {
						break;
					}

					if (ori == -1) {
						tot++;
					}

				}
				if (tot == p_facenumber) {

					return IN_PRISM;
				}
			}
			else {
				tot = 0;
				for (int j = 0; j < c_facenumber; j++) {


					e00 = envcubic[prismindex[i] - prism_size][c_face[j][0]][0]; e01 = envcubic[prismindex[i] - prism_size][c_face[j][0]][1]; e02 = envcubic[prismindex[i] - prism_size][c_face[j][0]][2];
					e10 = envcubic[prismindex[i] - prism_size][c_face[j][1]][0]; e11 = envcubic[prismindex[i] - prism_size][c_face[j][1]][1]; e12 = envcubic[prismindex[i] - prism_size][c_face[j][1]][2];
					e20 = envcubic[prismindex[i] - prism_size][c_face[j][2]][0]; e21 = envcubic[prismindex[i] - prism_size][c_face[j][2]][1]; e22 = envcubic[prismindex[i] - prism_size][c_face[j][2]][2];
					timer.start();
					ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22, checker);
					time_multi += timer.getElapsedTimeInSec();

					if (ori == 1 || ori == 0) {
						break;
					}

					if (ori == -1) {
						tot++;
					}

				}
				if (tot == c_facenumber) {

					return IN_PRISM;
				}
			}


		}
		return OUT_PRISM;
	}

	template<typename T>
	int  FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_double(
		const Scalar& a11, const Scalar&a12, const Scalar& a13, const Scalar& d, const Scalar& fa11,
		const Scalar& fa12, const Scalar& fa13, const Scalar& max1, const Scalar&max2, const Scalar& max5,
		const Vector3& segpoint0, const Vector3& segpoint1, const Vector3& triangle0,
		const Vector3& triangle1, const Vector3& triangle2, const std::vector<unsigned int>& prismindex, const int& jump, const std::function<int(T)> &checker) const {

		int tot;
		int  ori, ori1;
		INDEX index;
		std::vector<INDEX> recompute;
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump) {

				continue;
			}
			if (prismindex[i] < prism_size) {
				index.FACES.clear();
				tot = 0;
				for (int j = 0; j < p_facenumber; j++) {
					//ftimer2.start();
					ori = ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							segpoint0[0], segpoint0[1], segpoint0[2],
							envprism[prismindex[i]][p_face[j][0]][0], envprism[prismindex[i]][p_face[j][0]][1], envprism[prismindex[i]][p_face[j][0]][2],
							envprism[prismindex[i]][p_face[j][1]][0], envprism[prismindex[i]][p_face[j][1]][1], envprism[prismindex[i]][p_face[j][1]][2],
							envprism[prismindex[i]][p_face[j][2]][0], envprism[prismindex[i]][p_face[j][2]][1], envprism[prismindex[i]][p_face[j][2]][2]);


					if (ori == 1) {
						break;
					}
					if (ori == 0) {
						index.FACES.emplace_back(j);
					}

					else if (ori == -1) {
						tot++;
					}

				}
				if (tot == p_facenumber) {

					return IN_PRISM;
				}

				if (ori != 1) {
					assert(!index.FACES.empty());
					index.Pi = prismindex[i];
					recompute.emplace_back(index);
				}
			}
			else {
				index.FACES.clear();
				tot = 0;
				for (int j = 0; j < c_facenumber; j++) {
					//ftimer2.start();
					ori = ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							segpoint0[0], segpoint0[1], segpoint0[2],
							envcubic[prismindex[i] - prism_size][c_face[j][0]][0], envcubic[prismindex[i] - prism_size][c_face[j][0]][1], envcubic[prismindex[i] - prism_size][c_face[j][0]][2],
							envcubic[prismindex[i] - prism_size][c_face[j][1]][0], envcubic[prismindex[i] - prism_size][c_face[j][1]][1], envcubic[prismindex[i] - prism_size][c_face[j][1]][2],
							envcubic[prismindex[i] - prism_size][c_face[j][2]][0], envcubic[prismindex[i] - prism_size][c_face[j][2]][1], envcubic[prismindex[i] - prism_size][c_face[j][2]][2]);

					if (ori == 1) {
						break;
					}
					if (ori == 0) {
						index.FACES.emplace_back(j);
					}

					else if (ori == -1) {
						tot++;
					}

				}
				if (tot == c_facenumber) {

					return IN_PRISM;
				}

				if (ori != 1) {
					assert(!index.FACES.empty());
					index.Pi = prismindex[i];
					recompute.emplace_back(index);
				}
			}

		}

		if (!recompute.empty()) {
			static T
				s00, s01, s02, s10, s11, s12,
				t00, t01, t02,
				t10, t11, t12,
				t20, t21, t22,
				a11r, a12r, a13r, dr,

				e00, e01, e02,
				e10, e11, e12,
				e20, e21, e22;
			s00 = segpoint0[0]; s01 = segpoint0[1]; s02 = segpoint0[2]; s10 = segpoint1[0]; s11 = segpoint1[1]; s12 = segpoint1[2];
			t00 = triangle0[0]; t01 = triangle0[1]; t02 = triangle0[2];
			t10 = triangle1[0]; t11 = triangle1[1]; t12 = triangle1[2];
			t20 = triangle2[0]; t21 = triangle2[1]; t22 = triangle2[2];
			timer.start();
			bool premulti = orient3D_LPI_prefilter_multiprecision(s00, s01, s02, s10, s11, s12,
				t00, t01, t02, t10, t11, t12, t20, t21, t22, a11r, a12r, a13r, dr, checker);
			time_multi += timer.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++) {
				int in1 = recompute[k].Pi;
				if (in1 < prism_size) {
					for (int j = 0; j < recompute[k].FACES.size(); j++) {
						int in2 = recompute[k].FACES[j];


						e00 = envprism[in1][p_face[in2][0]][0]; e01 = envprism[in1][p_face[in2][0]][1]; e02 = envprism[in1][p_face[in2][0]][2];
						e10 = envprism[in1][p_face[in2][1]][0]; e11 = envprism[in1][p_face[in2][1]][1]; e12 = envprism[in1][p_face[in2][1]][2];
						e20 = envprism[in1][p_face[in2][2]][0]; e21 = envprism[in1][p_face[in2][2]][1]; e22 = envprism[in1][p_face[in2][2]][2];
						timer.start();
						ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
							e00, e01, e02, e10, e11, e12,
							e20, e21, e22, checker);
						time_multi += timer.getElapsedTimeInSec();
						if (ori == 1 || ori == 0) break;
					}

					if (ori == -1) return IN_PRISM;
				}
				else {
					for (int j = 0; j < recompute[k].FACES.size(); j++) {
						int in2 = recompute[k].FACES[j];


						e00 = envcubic[in1 - prism_size][c_face[in2][0]][0]; e01 = envcubic[in1 - prism_size][c_face[in2][0]][1]; e02 = envcubic[in1 - prism_size][c_face[in2][0]][2];
						e10 = envcubic[in1 - prism_size][c_face[in2][1]][0]; e11 = envcubic[in1 - prism_size][c_face[in2][1]][1]; e12 = envcubic[in1 - prism_size][c_face[in2][1]][2];
						e20 = envcubic[in1 - prism_size][c_face[in2][2]][0]; e21 = envcubic[in1 - prism_size][c_face[in2][2]][1]; e22 = envcubic[in1 - prism_size][c_face[in2][2]][2];
						timer.start();
						ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
							e00, e01, e02, e10, e11, e12,
							e20, e21, e22, checker);
						time_multi += timer.getElapsedTimeInSec();


						if (ori == 1 || ori == 0) break;
					}

					if (ori == -1) return IN_PRISM;
				}

			}


		}

		return OUT_PRISM;
	}


	template<typename T>
	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_TPI& datatpi, const std::array<Vector3, 3>&triangle, const std::vector<unsigned int>& prismindex, const std::function<int(T)> &checker)const {


		int tot,ori;
		int jump1 = datatpi.jump1, jump2 = datatpi.jump2;

		static T
			t00, t01, t02,
			t10, t11, t12,
			t20, t21, t22,

			f100, f101, f102,
			f110, f111, f112,
			f120, f121, f122,

			f200, f201, f202,
			f210, f211, f212,
			f220, f221, f222,
			dr, n1r, n2r, n3r,

			e00, e01, e02,
			e10, e11, e12,
			e20, e21, e22;
		t00 = (triangle[0][0]); t01 = (triangle[0][1]); t02 = (triangle[0][2]);
		t10 = (triangle[1][0]); t11 = (triangle[1][1]); t12 = (triangle[1][2]);
		t20 = (triangle[2][0]); t21 = (triangle[2][1]); t22 = (triangle[2][2]);

		if (datatpi.prismid1 < prism_size) {
			f100 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][0]][0];
			f101 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][0]][1];
			f102 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][0]][2];

			f110 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][1]][0];
			f111 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][1]][1];
			f112 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][1]][2];

			f120 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][2]][0];
			f121 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][2]][1];
			f122 = envprism[datatpi.prismid1][p_face[datatpi.facetid1][2]][2];
		}
		else {
			f100 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][0]][0];
			f101 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][0]][1];
			f102 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][0]][2];

			f110 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][1]][0];
			f111 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][1]][1];
			f112 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][1]][2];

			f120 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][2]][0];
			f121 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][2]][1];
			f122 = envcubic[datatpi.prismid1-prism_size][c_face[datatpi.facetid1][2]][2];
		}
		if (datatpi.prismid2 < prism_size) {
			f200 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][0]][0];
			f201 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][0]][1];
			f202 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][0]][2];
			f210 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][1]][0];
			f211 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][1]][1];
			f212 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][1]][2];
			f220 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][2]][0];
			f221 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][2]][1];
			f222 = envprism[datatpi.prismid2][p_face[datatpi.facetid2][2]][2];
		}
		else {
			f200 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][0]][0];
			f201 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][0]][1];
			f202 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][0]][2];
			f210 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][1]][0];
			f211 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][1]][1];
			f212 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][1]][2];
			f220 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][2]][0];
			f221 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][2]][1];
			f222 = envcubic[datatpi.prismid2-prism_size][c_face[datatpi.facetid2][2]][2];
		}

		timer.start();
		bool premulti = orient3D_TPI_prefilter_multiprecision(t00, t01, t02, t10, t11, t12, t20, t21, t22,
			f100, f101, f102, f110, f111, f112, f120, f121, f122,
			f200, f201, f202, f210, f211, f212, f220, f221, f222,
			dr, n1r, n2r, n3r, checker);
		time_multi += timer.getElapsedTimeInSec();
		if (premulti == false) return 2;//means have parallel facets
		bool cut = is_3_triangle_cut_pure_multiprecision(triangle, dr, n1r, n2r, n3r, checker);
		if (cut == false) return 2;
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
			if (prismindex[i] < prism_size) {
				tot = 0;
				for (int j = 0; j < p_facenumber; j++) {

					e00 = (envprism[prismindex[i]][p_face[j][0]][0]); e01 = (envprism[prismindex[i]][p_face[j][0]][1]); e02 = (envprism[prismindex[i]][p_face[j][0]][2]);
					e10 = (envprism[prismindex[i]][p_face[j][1]][0]); e11 = (envprism[prismindex[i]][p_face[j][1]][1]); e12 = (envprism[prismindex[i]][p_face[j][1]][2]);
					e20 = (envprism[prismindex[i]][p_face[j][2]][0]); e21 = (envprism[prismindex[i]][p_face[j][2]][1]); e22 = (envprism[prismindex[i]][p_face[j][2]][2]);
					timer.start();
					ori = orient3D_TPI_postfilter_multiprecision(dr, n1r, n2r, n3r, e00, e01, e02, e10, e11, e12, e20, e21, e22, checker);
					time_multi += timer.getElapsedTimeInSec();
					if (ori == 1 || ori == 0) {
						break;
					}

					if (ori == -1) {
						tot++;
					}

				}
				if (tot == p_facenumber) {

					return IN_PRISM;
				}
			}
			else {
				tot = 0;
				for (int j = 0; j < c_facenumber; j++) {

					e00 = (envcubic[prismindex[i] - prism_size][c_face[j][0]][0]); e01 = (envcubic[prismindex[i] - prism_size][c_face[j][0]][1]); e02 = (envcubic[prismindex[i] - prism_size][c_face[j][0]][2]);
					e10 = (envcubic[prismindex[i] - prism_size][c_face[j][1]][0]); e11 = (envcubic[prismindex[i] - prism_size][c_face[j][1]][1]); e12 = (envcubic[prismindex[i] - prism_size][c_face[j][1]][2]);
					e20 = (envcubic[prismindex[i] - prism_size][c_face[j][2]][0]); e21 = (envcubic[prismindex[i] - prism_size][c_face[j][2]][1]); e22 = (envcubic[prismindex[i] - prism_size][c_face[j][2]][2]);
					timer.start();
					ori = orient3D_TPI_postfilter_multiprecision(dr, n1r, n2r, n3r, e00, e01, e02, e10, e11, e12, e20, e21, e22, checker);
					time_multi += timer.getElapsedTimeInSec();
					if (ori == 1 || ori == 0) {
						break;
					}

					if (ori == -1) {
						tot++;
					}

				}
				if (tot == c_facenumber) {

					return IN_PRISM;
				}
			}


		}
		return OUT_PRISM;
	}

	template<typename T>
	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double(
		const Scalar& d, const Scalar& n1, const Scalar& n2, const Scalar& n3,
		const Scalar& max1, const Scalar& max2, const Scalar& max3, const Scalar& max4, const Scalar& max5, const Scalar&max6, const Scalar& max7,
		const std::array<Vector3, 3>& triangle,
		const Vector3& facet10, const Vector3& facet11, const Vector3& facet12, const Vector3& facet20, const Vector3& facet21, const Vector3& facet22,
		const std::vector<unsigned int>& prismindex, const int& jump1, const int &jump2, const bool & multiflag, const std::function<int(T)> &checker,
		 T& dr,  T & n1r,  T & n2r,  T& n3r) const {

		int ori;
		int tot;

		INDEX index;
		std::vector<INDEX> recompute;
		
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
			if (prismindex[i] < prism_size) {
				index.FACES.clear();
				tot = 0;
				for (int j = 0; j < p_facenumber; j++) {

					timer_a.start();



					ori = ip_filtered::
						orient3D_TPI_postfilter(
							d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							envprism[prismindex[i]][p_face[j][0]][0], envprism[prismindex[i]][p_face[j][0]][1], envprism[prismindex[i]][p_face[j][0]][2],
							envprism[prismindex[i]][p_face[j][1]][0], envprism[prismindex[i]][p_face[j][1]][1], envprism[prismindex[i]][p_face[j][1]][2],
							envprism[prismindex[i]][p_face[j][2]][0], envprism[prismindex[i]][p_face[j][2]][1], envprism[prismindex[i]][p_face[j][2]][2]);
					timetpp1 += timer_a.getElapsedTimeInSec();

					if (ori == 1) {
						break;
					}
					if (ori == 0) {
						index.FACES.emplace_back(j);
					}

					else if (ori == -1) {
						tot++;
					}
					

				}
				if (tot == p_facenumber) {

					return IN_PRISM;
				}

				if (ori != 1) {
					timer_a.start();
					index.Pi = prismindex[i];
					recompute.emplace_back(index);
					timetpp1 += timer_a.getElapsedTimeInSec();
				}
			}
			else {
				index.FACES.clear();
				tot = 0;
				for (int j = 0; j < c_facenumber; j++) {
					
					timer_a.start();
					ori = ip_filtered::
						orient3D_TPI_postfilter(
							d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							envcubic[prismindex[i] - prism_size][c_face[j][0]][0], envcubic[prismindex[i] - prism_size][c_face[j][0]][1], envcubic[prismindex[i] - prism_size][c_face[j][0]][2],
							envcubic[prismindex[i] - prism_size][c_face[j][1]][0], envcubic[prismindex[i] - prism_size][c_face[j][1]][1], envcubic[prismindex[i] - prism_size][c_face[j][1]][2],
							envcubic[prismindex[i] - prism_size][c_face[j][2]][0], envcubic[prismindex[i] - prism_size][c_face[j][2]][1], envcubic[prismindex[i] - prism_size][c_face[j][2]][2]);
					timetpp1 += timer_a.getElapsedTimeInSec();


					if (ori == 1) {
						break;
					}
					if (ori == 0) {
						index.FACES.emplace_back(j);
					}

					else if (ori == -1) {
						tot++;
					}

				}
				if (tot == c_facenumber) {

					return IN_PRISM;
				}

				if (ori != 1) {
					timer_a.start();
					index.Pi = prismindex[i];
					recompute.emplace_back(index);
					timetpp1 += timer_a.getElapsedTimeInSec();
				}
			}


		}
		
		
		if (recompute.size() > 0) {
			timer_a.start();
			static T
				t00, t01, t02,
				t10, t11, t12,
				t20, t21, t22,

				f100, f101, f102,
				f110, f111, f112,
				f120, f121, f122,

				f200, f201, f202,
				f210, f211, f212,
				f220, f221, f222,

				e00, e01, e02,
				e10, e11, e12,
				e20, e21, e22;;
			t00 = (triangle[0][0]); t01 = (triangle[0][1]); t02 = (triangle[0][2]);
			t10 = (triangle[1][0]); t11 = (triangle[1][1]); t12 = (triangle[1][2]);
			t20 = (triangle[2][0]); t21 = (triangle[2][1]); t22 = (triangle[2][2]);

			f100 = (facet10[0]); f101 = (facet10[1]); f102 = (facet10[2]);
			f110 = (facet11[0]); f111 = (facet11[1]); f112 = (facet11[2]);
			f120 = (facet12[0]); f121 = (facet12[1]); f122 = (facet12[2]);

			f200 = (facet20[0]); f201 = (facet20[1]); f202 = (facet20[2]);
			f210 = (facet21[0]); f211 = (facet21[1]); f212 = (facet21[2]);
			f220 = (facet22[0]); f221 = (facet22[1]); f222 = (facet22[2]);

			if (multiflag == false) {
				timer.start();
				bool premulti = orient3D_TPI_prefilter_multiprecision(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					dr, n1r, n2r, n3r, checker);
				time_multi += timer.getElapsedTimeInSec();
			}

			timetpp2 += timer_a.getElapsedTimeInSec();




			for (int k = 0; k < recompute.size(); k++) {
				int in1 = recompute[k].Pi;

				if (in1 < prism_size) {
					for (int j = 0; j < recompute[k].FACES.size(); j++) {
						int in2 = recompute[k].FACES[j];
						timer_a.start();
						e00 = (envprism[in1][p_face[in2][0]][0]); e01 = (envprism[in1][p_face[in2][0]][1]); e02 = (envprism[in1][p_face[in2][0]][2]);
						e10 = (envprism[in1][p_face[in2][1]][0]); e11 = (envprism[in1][p_face[in2][1]][1]); e12 = (envprism[in1][p_face[in2][1]][2]);
						e20 = (envprism[in1][p_face[in2][2]][0]); e21 = (envprism[in1][p_face[in2][2]][1]); e22 = (envprism[in1][p_face[in2][2]][2]);
						timer.start();
						ori = orient3D_TPI_postfilter_multiprecision(dr, n1r, n2r, n3r,
							e00, e01, e02, e10, e11, e12,
							e20, e21, e22, checker);
						time_multi += timer.getElapsedTimeInSec();
						timetpp2 += timer_a.getElapsedTimeInSec();
						if (ori == 1 || ori == 0) break;
					}

					if (ori == -1) return IN_PRISM;
				}
				else {
					for (int j = 0; j < recompute[k].FACES.size(); j++) {
						int in2 = recompute[k].FACES[j];
						timer_a.start();

						e00 = (envcubic[in1 - prism_size][c_face[in2][0]][0]); e01 = (envcubic[in1 - prism_size][c_face[in2][0]][1]); e02 = (envcubic[in1 - prism_size][c_face[in2][0]][2]);
						e10 = (envcubic[in1 - prism_size][c_face[in2][1]][0]); e11 = (envcubic[in1 - prism_size][c_face[in2][1]][1]); e12 = (envcubic[in1 - prism_size][c_face[in2][1]][2]);
						e20 = (envcubic[in1 - prism_size][c_face[in2][2]][0]); e21 = (envcubic[in1 - prism_size][c_face[in2][2]][1]); e22 = (envcubic[in1 - prism_size][c_face[in2][2]][2]);
						timer.start();
						ori = orient3D_TPI_postfilter_multiprecision(dr, n1r, n2r, n3r,
							e00, e01, e02, e10, e11, e12,
							e20, e21, e22, checker);
						time_multi += timer.getElapsedTimeInSec();
						timetpp2 += timer_a.getElapsedTimeInSec();

						//if (ori == -2) std::cout << "impossible thing happens in lpi" << std::endl;
						if (ori == 1 || ori == 0) break;
					}

					if (ori == -1) return IN_PRISM;
				}

			}

		}

		

		return OUT_PRISM;
	}
#include<ctime>
	template<typename T>
	bool FastEnvelope::is_3_triangle_cut_pure_multiprecision(const std::array<Vector3, 3>& triangle, const T& dr, const T& n1r, const T& n2r, const T& n3r, const std::function<int(T)> &checker) {
		int o1, o2, o3;
		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];

		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
			std::cout << "Degeneration happens" << std::endl;
			//move this guy in constructor and use fixed seed
			srand(int(time(0)));
			n = { {Vector3(rand(),rand(),rand()) } };
		}
		static T nr0, nr1, nr2,
			t00, t01, t02,
			t10, t11, t12,
			t20, t21, t22;
		nr0 = n[0];
		nr1 = n[1];
		nr2 = n[2];
		t00 = triangle[0][0];
		t01 = triangle[0][1];
		t02 = triangle[0][2];
		t10 = triangle[1][0];
		t11 = triangle[1][1];
		t12 = triangle[1][2];
		t20 = triangle[2][0];
		t21 = triangle[2][1];
		t22 = triangle[2][2];
		timer.start();
		o1 = orient3D_TPI_postfilter_multiprecision(
			dr, n1r, n2r, n3r,
			nr0, nr1, nr2,
			t00, t01, t02,
			t10, t11, t12, checker);
		time_multi += timer.getElapsedTimeInSec();
		if (o1 == 0) return false;
		timer.start();
		o2 = orient3D_TPI_postfilter_multiprecision(
			dr, n1r, n2r, n3r,
			nr0, nr1, nr2,
			t10, t11, t12,
			t20, t21, t22, checker);
		time_multi += timer.getElapsedTimeInSec();
		if (o2 == 0 || o1 + o2 == 0) return false;
		timer.start();
		o3 = orient3D_TPI_postfilter_multiprecision(
			dr, n1r, n2r, n3r,
			nr0, nr1, nr2,
			t20, t21, t22,
			t00, t01, t02, checker);
		time_multi += timer.getElapsedTimeInSec();
		if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0) return false;
		return true;
	}

	template<typename T>
	bool FastEnvelope::is_3_triangle_cut_double(
		const Scalar &d, const Scalar & n1, const Scalar &n2, const Scalar & n3,
		const Scalar & max1, const Scalar &max2, const Scalar &max3, const Scalar & max4, const Scalar & max5,
		const Scalar & max6, const Scalar &max7,
		const std::array<Vector3, 3>& triangle,
		const Vector3& facet10, const Vector3& facet11, const Vector3& facet12, const Vector3& facet20, const Vector3& facet21, const Vector3& facet22,
		bool & multiflag,
		const std::function<int(T)> &checker, T &dr, T & n1r, T &n2r, T & n3r) {
		multiflag = false;
		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];

		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
			std::cout << "Degeneration happens" << std::endl;
			//move this guy in constructor and use fixed seed
			srand(int(time(0)));
			n = { {Vector3(rand(),rand(),rand()) } };
		}



		static T
			t00, t01, t02,
			t10, t11, t12,
			t20, t21, t22,

			f100, f101, f102,
			f110, f111, f112,
			f120, f121, f122,

			f200, f201, f202,
			f210, f211, f212,
			f220, f221, f222,

			nr0, nr1, nr2;

		bool premulti = false;
		int o1 = ip_filtered::orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2]);
		if (o1 == 0) {


			t00 = (triangle[0][0]); t01 = (triangle[0][1]); t02 = (triangle[0][2]);
			t10 = (triangle[1][0]); t11 = (triangle[1][1]); t12 = (triangle[1][2]);
			t20 = (triangle[2][0]); t21 = (triangle[2][1]); t22 = (triangle[2][2]);

			f100 = (facet10[0]); f101 = (facet10[1]); f102 = (facet10[2]);
			f110 = (facet11[0]); f111 = (facet11[1]); f112 = (facet11[2]);
			f120 = (facet12[0]); f121 = (facet12[1]); f122 = (facet12[2]);

			f200 = (facet20[0]); f201 = (facet20[1]); f202 = (facet20[2]);
			f210 = (facet21[0]); f211 = (facet21[1]); f212 = (facet21[2]);
			f220 = (facet22[0]); f221 = (facet22[1]); f222 = (facet22[2]);

			nr0 = (n[0]); nr1 = (n[1]); nr2 = (n[2]);
			timer.start();
			premulti = orient3D_TPI_prefilter_multiprecision(
				t00, t01, t02, t10, t11, t12, t20, t21, t22,
				f100, f101, f102, f110, f111, f112, f120, f121, f122,
				f200, f201, f202, f210, f211, f212, f220, f221, f222,
				dr, n1r, n2r, n3r, checker);
			time_multi += timer.getElapsedTimeInSec();
			timer.start();
			o1 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t00, t01, t02,
				t10, t11, t12, checker);
			time_multi += timer.getElapsedTimeInSec();
		}

		if (o1 == 0) return false;

		int o2 = ip_filtered::orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2]);
		if (o2 == 0) {
			if (premulti == false) {
				t00 = (triangle[0][0]); t01 = (triangle[0][1]); t02 = (triangle[0][2]);
				t10 = (triangle[1][0]); t11 = (triangle[1][1]); t12 = (triangle[1][2]);
				t20 = (triangle[2][0]); t21 = (triangle[2][1]); t22 = (triangle[2][2]);

				f100 = (facet10[0]); f101 = (facet10[1]); f102 = (facet10[2]);
				f110 = (facet11[0]); f111 = (facet11[1]); f112 = (facet11[2]);
				f120 = (facet12[0]); f121 = (facet12[1]); f122 = (facet12[2]);

				f200 = (facet20[0]); f201 = (facet20[1]); f202 = (facet20[2]);
				f210 = (facet21[0]); f211 = (facet21[1]); f212 = (facet21[2]);
				f220 = (facet22[0]); f221 = (facet22[1]); f222 = (facet22[2]);

				nr0 = (n[0]); nr1 = (n[1]); nr2 = (n[2]);
				timer.start();
				premulti = orient3D_TPI_prefilter_multiprecision(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					dr, n1r, n2r, n3r, checker);
				time_multi += timer.getElapsedTimeInSec();
			}
			timer.start();
			o2 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t10, t11, t12,
				t20, t21, t22, checker);
			time_multi += timer.getElapsedTimeInSec();
			/*if (o2 == 1) after21++;
			if (o2 == -1) after22++;
			if (o2 == 0) after20++;*/
		}
		if (o2 == 0 || o1 + o2 == 0) return false;

		int o3 = ip_filtered::orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			triangle[0][0], triangle[0][1], triangle[0][2]);
		if (o3 == 0) {
			if (premulti == false) {
				t00 = (triangle[0][0]); t01 = (triangle[0][1]); t02 = (triangle[0][2]);
				t10 = (triangle[1][0]); t11 = (triangle[1][1]); t12 = (triangle[1][2]);
				t20 = (triangle[2][0]); t21 = (triangle[2][1]); t22 = (triangle[2][2]);

				f100 = (facet10[0]); f101 = (facet10[1]); f102 = (facet10[2]);
				f110 = (facet11[0]); f111 = (facet11[1]); f112 = (facet11[2]);
				f120 = (facet12[0]); f121 = (facet12[1]); f122 = (facet12[2]);

				f200 = (facet20[0]); f201 = (facet20[1]); f202 = (facet20[2]);
				f210 = (facet21[0]); f211 = (facet21[1]); f212 = (facet21[2]);
				f220 = (facet22[0]); f221 = (facet22[1]); f222 = (facet22[2]);

				nr0 = (n[0]); nr1 = (n[1]); nr2 = (n[2]);
				timer.start();
				premulti = orient3D_TPI_prefilter_multiprecision(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					dr, n1r, n2r, n3r, checker);
				time_multi += timer.getElapsedTimeInSec();
			}
			timer.start();
			o3 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t20, t21, t22,
				t00, t01, t02, checker);
			time_multi += timer.getElapsedTimeInSec();
			/*if (o3 == 1) after21++;
			if (o3 == -1) after22++;
			if (o3 == 0) after20++;*/
		}
		if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0) return false;
		if (premulti == true) multiflag = true;
		return true;
	}
	int FastEnvelope::is_3_triangle_cut_float_fast(
		const Vector3& tri0, const Vector3& tri1, const Vector3& tri2,
		const Vector3& facet10, const Vector3& facet11, const Vector3& facet12,
		const Vector3& facet20, const Vector3& facet21, const Vector3& facet22) {

		Vector3 n = (tri0 - tri1).cross(tri0 - tri2) + tri0;

		if (Predicates::orient_3d(n, tri0, tri1, tri2) == 0) {
			std::cout << "Degeneration happens !" << std::endl;

			srand(int(time(0)));
			n = { {Vector3(rand(),rand(),rand()) } };
		}
		Scalar d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7;
		bool pre = ip_filtered::
			orient3D_TPI_prefilter(
				tri0[0], tri0[1], tri0[2],
				tri1[0], tri1[1], tri1[2],
				tri2[0], tri2[1], tri2[2],
				facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
				d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

		if (pre == false) return 2;// means we dont know
		int o1 = ip_filtered::orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			tri0[0], tri0[1], tri0[2],
			tri1[0], tri1[1], tri1[2]);
		int o2 = ip_filtered::orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			tri1[0], tri1[1], tri1[2],
			tri2[0], tri2[1], tri2[2]);
		if (o1*o2 == -1) return 0;
		int o3 = ip_filtered::orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			tri2[0], tri2[1], tri2[2],
			tri0[0], tri0[1], tri0[2]);
		if (o1*o3 == -1 || o2 * o3 == -1) return 0;
		if (o1*o2*o3 == 0) return 2;// means we dont know
		return 1;

	}





	int FastEnvelope::seg_cut_plane(const Vector3 & seg0, const Vector3 &seg1, const Vector3&t0, const Vector3&t1, const Vector3 &t2) {
		int o1, o2;
		o1 = Predicates::orient_3d(seg0, t0, t1, t2);
		o2 = Predicates::orient_3d(seg1, t0, t1, t2);
		int op = o1 * o2;
		if (op >= 0) {
			return CUT_COPLANAR;//in fact, coplanar and not on this plane
		}
		return CUT_FACE;
	}


	bool FastEnvelope::is_triangle_cut_prism(const int&pindex,
		const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, std::vector<int> &cid)const{


		bool cut[8];
		for (int i = 0; i < 8; i++) {
			cut[i] = false;
		}
		int o1[8], o2[8], o3[8],ori=0;
		std::vector<int> cutp;

		for (int i = 0; i < 8; i++) {

			o1[i] = Predicates::orient_3d(tri0, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			o2[i] = Predicates::orient_3d(tri1, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			o3[i] = Predicates::orient_3d(tri2, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			if (o1[i] + o2[i] + o3[i] >= 3) {
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1) {
				return false;
			}
			if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0) {
				return false;
			}
			if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0) {
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 0) {
				return false;
			}


			if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1) cutp.emplace_back(i);
		}
		if (cutp.size() ==0) {
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++) {
			if (o1[cutp[i]] * o2[cutp[i]] == -1) {

				bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri1[0], tri1[1], tri1[2],

					envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
					envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
					envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false) {
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++) {
					if (i == j) continue;
					ori= ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
							envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
							envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

					if (ori == 1) break;

				}
				if (ori != 1) {
					cut[cutp[i]] = true;
				}

			}
			if (cut[cutp[i]] == true) continue;
			ori = 0;
			if (o1[cutp[i]] * o3[cutp[i]] == -1) {

				bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri2[0], tri2[1], tri2[2],
					envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
					envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
					envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false) {
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++) {
					if (i == j) continue;
					ori = ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
							envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
							envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

					if (ori == 1) break;

				}
				if (ori != 1) {
					cut[cutp[i]] = true;
				}
			}


			if (cut[cutp[i]] == true) continue;
			ori = 0;
			if (o2[cutp[i]] * o3[cutp[i]] == -1) {

				bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
					tri1[0], tri1[1], tri1[2],
					tri2[0], tri2[1], tri2[2],
					envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
					envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
					envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false) {
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++) {
					if (i == j) continue;
					ori = ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri1[0], tri1[1], tri1[2],
							envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
							envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
							envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

					if (ori == 1) break;

				}
				if (ori != 1) {
					cut[cutp[i]] = true;
				}
			}

		}

		if (cutp.size() <= 2) {
			for (int i = 0; i < 8; i++) {
				if (cut[i] == true) cid.emplace_back(i);
			}
			return true;
		}
		// triangle-facet-facet intersection
		Scalar  n1, n2, n3, max3, max4, max6, max7;
		for (int i = 0; i < cutp.size(); i++) {
			for (int j = i + 1; j < cutp.size(); j++) {
				if (cut[cutp[i]] == true && cut[cutp[j]] == true) continue;

				int id = cutp[i] * 8 + cutp[j];
				int id0 = prism_map[id][0];
				if (id0 == -1) continue;
				int inter = is_3_triangle_cut_float_fast(
					tri0, tri1, tri2,
					envprism[pindex][p_face[cutp[i]][0]],
					envprism[pindex][p_face[cutp[i]][1]],
					envprism[pindex][p_face[cutp[i]][2]],
					envprism[pindex][p_face[cutp[j]][0]],
					envprism[pindex][p_face[cutp[j]][1]],
					envprism[pindex][p_face[cutp[j]][2]]);
				if (inter == 2) {//we dont know if point exist or if inside of triangle
					cut[cutp[i]] == true;
					cut[cutp[j]] == true;
					continue;
				}
				if (inter == 0) continue;// sure not inside

				bool pre = ip_filtered::
					orient3D_TPI_prefilter(
						tri0[0], tri0[1], tri0[2],
						tri1[0], tri1[1], tri1[2],
						tri2[0], tri2[1], tri2[2],
						envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
						envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
						envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
						envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
						envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
						envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2],
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

				for (int k = 0; k < cutp.size(); k++) {

					if (k == i || k == j) continue;

					ori = ip_filtered::
						orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							envprism[pindex][p_face[cutp[k]][0]][0], envprism[pindex][p_face[cutp[k]][0]][1], envprism[pindex][p_face[cutp[k]][0]][2],
							envprism[pindex][p_face[cutp[k]][1]][0], envprism[pindex][p_face[cutp[k]][1]][1], envprism[pindex][p_face[cutp[k]][1]][2],
							envprism[pindex][p_face[cutp[k]][2]][0], envprism[pindex][p_face[cutp[k]][2]][1], envprism[pindex][p_face[cutp[k]][2]][2]);

					if (ori == 1) break;

				}

				if (ori != 1) {
					cut[cutp[i]] = true;
					cut[cutp[j]] = true;
				}
			}
		}

		for (int i = 0; i < 8; i++) {
			if (cut[i] == true) cid.emplace_back(i);
		}

		return true;

	}

	bool FastEnvelope::is_seg_cut_prism(const int&pindex,
		const Vector3& seg0, const Vector3& seg1, std::vector<int> &cid)const {
		bool cut[8];
		for (int i = 0; i < 8; i++) {
			cut[i] = false;
		}
		int o1[8], o2[8], ori = 0;
		std::vector<int> cutp;

		for (int i = 0; i < 8; i++) {

			o1[i] = Predicates::orient_3d(seg0, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			o2[i] = Predicates::orient_3d(seg1, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);

			if (o1[i] + o2[i] >= 1) {
				return false;
			}

			if (o1[i] == 0 && o2[i] == 0) {
				return false;
			}

			if (o1[i] * o2[i] == -1) cutp.emplace_back(i);
		}
		if (cutp.size() == 0) {
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++) {


			bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
				seg0[0], seg0[1], seg0[2],
				seg1[0], seg1[1], seg1[2],
				envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
				envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
				envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
				a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
			if (precom == false) {
				cut[cutp[i]] = true;
				continue;
			}
			for (int j = 0; j < cutp.size(); j++) {
				if (i == j) continue;
				ori = ip_filtered::
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						seg0[0], seg0[1], seg0[2],
						envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
						envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
						envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

				if (ori == 1) break;

			}
			if (ori != 1) {
				cut[cutp[i]] = true;
			}

		}

		for (int i = 0; i < 8; i++) {
			if (cut[i] == true) cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::is_triangle_cut_cube(const int&cindex,
		const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, std::vector<int> &cid)const {

		bool cut[6];
		for (int i = 0; i < 6; i++) {
			cut[i] = false;
		}
		int o1[6], o2[6], o3[6], ori = 0;
		std::vector<int> cutp;

		for (int i = 0; i < 6; i++) {

			o1[i] = Predicates::orient_3d(tri0, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			o2[i] = Predicates::orient_3d(tri1, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			o3[i] = Predicates::orient_3d(tri2, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			if (o1[i] + o2[i] + o3[i] >= 3) {
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1) {
				return false;
			}
			if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0) {
				return false;
			}
			if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0) {
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 0) {
				return false;
			}


			if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1) cutp.emplace_back(i);
		}
		if (cutp.size() == 0) {
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++) {
			if (o1[cutp[i]] * o2[cutp[i]] == -1) {

				bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri1[0], tri1[1], tri1[2],

					envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
					envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
					envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false) {
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++) {
					if (i == j) continue;
					ori = ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
							envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
							envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]
						);

					if (ori == 1) break;

				}
				if (ori != 1) {
					cut[cutp[i]] = true;
				}

			}
			if (cut[cutp[i]] == true) continue;
			ori = 0;
			if (o1[cutp[i]] * o3[cutp[i]] == -1) {

				bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri2[0], tri2[1], tri2[2],
					envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
					envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
					envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false) {
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++) {
					if (i == j) continue;
					ori = ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
							envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
							envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]);

					if (ori == 1) break;

				}
				if (ori != 1) {
					cut[cutp[i]] = true;
				}
			}


			if (cut[cutp[i]] == true) continue;
			ori = 0;
			if (o2[cutp[i]] * o3[cutp[i]] == -1) {

				bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
					tri1[0], tri1[1], tri1[2],
					tri2[0], tri2[1], tri2[2],
					envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
					envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
					envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false) {
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++) {
					if (i == j) continue;
					ori = ip_filtered::
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri1[0], tri1[1], tri1[2],
							envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
							envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
							envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]);

					if (ori == 1) break;

				}
				if (ori != 1) {
					cut[cutp[i]] = true;
				}
			}

		}

		if (cutp.size() <= 2) {
			for (int i = 0; i < 6; i++) {
				if (cut[i] == true) cid.emplace_back(i);
			}
			return true;
		}
		// triangle-facet-facet intersection
		Scalar  n1, n2, n3, max3, max4, max6, max7;
		for (int i = 0; i < cutp.size(); i++) {
			for (int j = i + 1; j < cutp.size(); j++) {
				if (cut[cutp[i]] == true && cut[cutp[j]] == true) continue;

				int id = cutp[i] * 6 + cutp[j];
				int id0 = prism_map[id][0];
				if (id0 == -1) continue;
				int inter = is_3_triangle_cut_float_fast(
					tri0, tri1, tri2,
					envcubic[cindex][c_face[cutp[i]][0]],
					envcubic[cindex][c_face[cutp[i]][1]],
					envcubic[cindex][c_face[cutp[i]][2]],
					envcubic[cindex][c_face[cutp[j]][0]],
					envcubic[cindex][c_face[cutp[j]][1]],
					envcubic[cindex][c_face[cutp[j]][2]]);
				if (inter == 2) {//we dont know if point exist or if inside of triangle
					cut[cutp[i]] == true;
					cut[cutp[j]] == true;
					continue;
				}
				if (inter == 0) continue;// sure not inside

				bool pre = ip_filtered::
					orient3D_TPI_prefilter(
						tri0[0], tri0[1], tri0[2],
						tri1[0], tri1[1], tri1[2],
						tri2[0], tri2[1], tri2[2],
						envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
						envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
						envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
						envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
						envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
						envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2],
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

				for (int k = 0; k < cutp.size(); k++) {

					if (k == i || k == j) continue;

					ori = ip_filtered::
						orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							envcubic[cindex][c_face[cutp[k]][0]][0], envcubic[cindex][c_face[cutp[k]][0]][1], envcubic[cindex][c_face[cutp[k]][0]][2],
							envcubic[cindex][c_face[cutp[k]][1]][0], envcubic[cindex][c_face[cutp[k]][1]][1], envcubic[cindex][c_face[cutp[k]][1]][2],
							envcubic[cindex][c_face[cutp[k]][2]][0], envcubic[cindex][c_face[cutp[k]][2]][1], envcubic[cindex][c_face[cutp[k]][2]][2]);

					if (ori == 1) break;

				}

				if (ori != 1) {
					cut[cutp[i]] = true;
					cut[cutp[j]] = true;
				}
			}
		}

		for (int i = 0; i < 6; i++) {
			if (cut[i] == true) cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::is_seg_cut_cube(const int&cindex,
		const Vector3& seg0, const Vector3& seg1, std::vector<int> &cid) const{
		bool cut[6];
		for (int i = 0; i < 6; i++) {
			cut[i] = false;
		}
		int o1[6], o2[6], ori = 0;
		std::vector<int> cutp;

		for (int i = 0; i < 6; i++) {

			o1[i] = Predicates::orient_3d(seg0, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			o2[i] = Predicates::orient_3d(seg1, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);

			if (o1[i] + o2[i] >= 1) {
				return false;
			}

			if (o1[i] == 0 && o2[i] == 0) {
				return false;
			}

			if (o1[i] * o2[i] == -1) cutp.emplace_back(i);
		}
		if (cutp.size() == 0) {
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++) {


			bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
				seg0[0], seg0[1], seg0[2],
				seg1[0], seg1[1], seg1[2],
				envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
				envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
				envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
				a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
			if (precom == false) {
				cut[cutp[i]] = true;
				continue;
			}
			for (int j = 0; j < cutp.size(); j++) {
				if (i == j) continue;
				ori = ip_filtered::
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						seg0[0], seg0[1], seg0[2],
						envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
						envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
						envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]);

				if (ori == 1) break;

			}
			if (ori != 1) {
				cut[cutp[i]] = true;
			}

		}

		for (int i = 0; i < 6; i++) {
			if (cut[i] == true) cid.emplace_back(i);
		}

		return true;
	}

	bool FastEnvelope::is_triangle_cut_bounding_box(
		const Vector3& tri0, const Vector3& tri1, const Vector3& tri2, const Vector3 &bmin, const Vector3 &bmax) {
		Vector3 tmin, tmax;
		get_tri_corners(tri0, tri1, tri2, tmin, tmax);
		return box_box_intersection(tmin, tmax, bmin, bmax);
	}


	bool FastEnvelope::point_out_prism(const Vector3 & point, const std::vector<unsigned int>& prismindex, const int& jump)const
	{

		int  ori;
		int psize = envprism.size();
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump) continue;
			if (prismindex[i] < psize) {
				for (int j = 0; j < p_facenumber; j++) {

					ori = Predicates::orient_3d(envprism[prismindex[i]][p_face[j][0]], envprism[prismindex[i]][p_face[j][1]], envprism[prismindex[i]][p_face[j][2]], point);
					if (ori == -1 || ori == 0) {
						break;
					}
					if (j == 7) {

						return false;
					}
				}
			}
			else {
				int boxid = prismindex[i] - psize;
				for (int j = 0; j < c_facenumber; j++) {

					ori = Predicates::orient_3d(envcubic[boxid][c_face[j][0]], envcubic[boxid][c_face[j][1]], envcubic[boxid][c_face[j][2]], point);

					if (ori == -1 || ori == 0) {
						break;
					}
					if (j == 5) {

						return false;
					}
				}
			}

		}

		return true;
	}
	bool FastEnvelope::point_out_prism_rational(const Rational& point0, const Rational& point1, const Rational& point2,  const std::vector<unsigned int>& prismindex, const int& jump)const
	{

		int  ori;
		int psize = envprism.size();
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump) continue;
			if (prismindex[i] < psize) {
				for (int j = 0; j < p_facenumber; j++) {

					ori = orient_3d_rational(point0,point1,point2,envprism[prismindex[i]][p_face[j][0]], envprism[prismindex[i]][p_face[j][1]], envprism[prismindex[i]][p_face[j][2]]);
					if (ori == 1 || ori == 0) {
						break;
					}
					if (j == 7) {

						return false;
					}
				}
			}
			else {
				int boxid = prismindex[i] - psize;
				for (int j = 0; j < c_facenumber; j++) {

					ori = orient_3d_rational(point0, point1, point2, envcubic[boxid][c_face[j][0]], envcubic[boxid][c_face[j][1]], envcubic[boxid][c_face[j][2]]);

					if (ori == 1 || ori == 0) {
						break;
					}
					if (j == 5) {

						return false;
					}
				}
			}

		}

		return true;
	}
	int FastEnvelope::is_triangle_degenerated(const Vector3& triangle0, const Vector3& triangle1, const Vector3& triangle2) {

		Vector3 a = triangle0 - triangle1, b = triangle0 - triangle2;
		Vector3 normal = a.cross(b);
		Scalar nbr = normal.norm();

		if (nbr > SCALAR_ZERO) {
			return NOT_DEGENERATED;
		}
		int ori;
		std::array < Vector2, 3> p;
		for (int j = 0; j < 3; j++) {

			p[0] = to_2d(triangle0, j);
			p[1] = to_2d(triangle1, j);
			p[2] = to_2d(triangle2, j);

			ori = Predicates::orient_2d(p[0], p[1], p[2]);
			if (ori != 0) {
				return NERLY_DEGENERATED;
			}
		}

		if (triangle0[0] != triangle1[0] || triangle0[1] != triangle1[1] || triangle0[2] != triangle1[2]) {
			return DEGENERATED_SEGMENT;
		}
		if (triangle0[0] != triangle2[0] || triangle0[1] != triangle2[1] || triangle0[2] != triangle2[2]) {
			return DEGENERATED_SEGMENT;
		}
		return DEGENERATED_POINT;

	}
	void FastEnvelope::BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, std::vector<std::array<Vector3, 8>>& envbox, const Scalar& epsilon)
	{
		envprism.reserve(m_faces.size());
		Vector3 AB, AC, BC, normal, vector1, ABn;
		Parameters pram;
		std::array<Vector3, 6> polygon;
		std::array<Vector3, 12> polygonoff;
		std::array<Vector3, 8> box;
		std::array<Vector3, 8> boxorder = {
			{
				{1,1,1},
		{-1,1,1},
		{-1,-1,1},
		{1,-1,1},
		{1,1,-1},
		{-1,1,-1},
		{-1,-1,-1},
		{1,-1,-1},
		}
		};

		Scalar
			tolerance = epsilon / sqrt(3),

			de;

		for (int i = 0; i < m_faces.size(); i++) {
			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			de = is_triangle_degenerated(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);



			if (de == DEGENERATED_POINT) {
				std::cout << "Envelope Triangle Degeneration- Point" << std::endl;
				for (int j = 0; j < 8; j++) {
					box[j] = m_ver[m_faces[i][0]] + boxorder[j] * tolerance;
				}
				envbox.emplace_back(box);
				continue;
			}
			if (de == DEGENERATED_SEGMENT) {
				std::cout << "Envelope Triangle Degeneration- Segment" << std::endl;
				Scalar length1 = AB.norm(), length2 = AC.norm(), length3 = BC.norm();
				if (length1 >= length2 && length1 >= length3) {
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], tolerance, box);
					envbox.emplace_back(box);
				}
				if (length2 >= length1 && length2 >= length3) {
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][2]], tolerance, box);
					envbox.emplace_back(box);
				}
				if (length3 >= length1 && length3 >= length2) {
					seg_cube(m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tolerance, box);
					envbox.emplace_back(box);
				}
				continue;
			}
			if (de == NERLY_DEGENERATED) {
				std::cout << "Envelope Triangle Degeneration- Nearly" << std::endl;
				//std::cout << std::setprecision(17) <<m_ver[m_faces[i][0]][0] << " " << m_ver[m_faces[i][0]][1] << " " << m_ver[m_faces[i][0]][2] << std::endl;
				//std::cout << std::setprecision(17) <<m_ver[m_faces[i][1]][0] << " " << m_ver[m_faces[i][1]][1] << " " << m_ver[m_faces[i][1]][2] << std::endl;
				//std::cout << std::setprecision(17) <<m_ver[m_faces[i][2]][0] << " " << m_ver[m_faces[i][2]][1] << " " << m_ver[m_faces[i][2]][2] << std::endl;

				normal = accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][0]], m_ver[m_faces[i][2]]);
				//std::cout << "pass1" << std::endl;
				vector1 = accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], origin, normal);
				//std::cout << "pass2" << std::endl;

			}
			else {
				normal = AB.cross(AC).normalized();
				vector1 = AB.cross(normal).normalized();
			}

			ABn = AB.normalized();
			polygon[0] = m_ver[m_faces[i][0]] + (vector1 - ABn) * tolerance;
			polygon[1] = m_ver[m_faces[i][1]] + (vector1 + ABn) * tolerance;
			if (dot_sign(AB,BC) < 0) {
				polygon[2] = m_ver[m_faces[i][1]] + (-vector1 + ABn) * tolerance;
				polygon[3] = m_ver[m_faces[i][2]] + (-vector1 + ABn) * tolerance;
				polygon[4] = m_ver[m_faces[i][2]] + (-vector1 - ABn) * tolerance;
				if (dot_sign(AB, AC) < 0) {
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

			for (int j = 0; j < 6; j++) {
				polygonoff[j] = polygon[j] + normal * tolerance;
			}
			for (int j = 6; j < 12; j++) {
				polygonoff[j] = polygon[j - 6] - normal * tolerance;
			}
			envprism.emplace_back(polygonoff);

		}

	}
	void FastEnvelope::seg_cube(const Vector3 &p1, const Vector3 &p2, const Scalar& width, std::array<Vector3, 8>& envbox) {
		Vector3 v1, v2, v = p2 - p1;//p12
		if (v[0] != 0) {
			v1 = Vector3((0 - v[1] - v[2]) / v[0], 1, 1);
		}
		else {
			if (v[1] != 0) {
				v1 = Vector3(1, (0 - v[2]) / v[1], 1);

			}
			else {
				v1 = Vector3(1, 1, 0);

			}
		}
		v2 = v.cross(v1);
		v = v.normalized();
		v1 = v1.normalized();
		v2 = v2.normalized();
		envbox[0] = p2 + width * (v + v1 + v2);
		envbox[1] = p2 + width * (v - v1 + v2);
		envbox[2] = p2 + width * (v - v1 - v2);
		envbox[3] = p2 + width * (v + v1 - v2);//right hand out direction
		envbox[4] = p1 + width * (-v + v1 + v2);
		envbox[5] = p1 + width * (-v - v1 + v2);
		envbox[6] = p1 + width * (-v - v1 - v2);
		envbox[7] = p1 + width * (-v + v1 - v2);//right hand in direction
	}

	Vector3 FastEnvelope::accurate_normal_vector(const Vector3 & p0, const Vector3 & p1,
		const Vector3 & q0, const Vector3 & q1) {

		const Multiprecision ax = p1[0] - p0[0];
		const Multiprecision ay = p1[1] - p0[1];
		const Multiprecision az = p1[2] - p0[2];

		const Multiprecision bx = q1[0] - q0[0];
		const Multiprecision by = q1[1] - q0[1];
		const Multiprecision bz = q1[2] - q0[2];

		Multiprecision x = ay * bz - az * by;
		Multiprecision y = az * bx - ax * bz;
		Multiprecision z = ax * by - ay * bx;
		Multiprecision ssum = x * x + y * y + z * z;
		if (ssum == 0) {
			std::cout << "divided by zero in accuratexxx" << std::endl;
			/*std::cout << std::setprecision(17) << p[0] << " " << p[1] << " " << p[2] << std::endl;
			std::cout << std::setprecision(17) << q[0] << " " << q[1] << " " << q[2] << std::endl;*/
			Rational p00r(p0[0]), p01r(p0[1]), p02r(p0[2]),
				p10r(p1[0]), p11r(p1[1]), p12r(p1[2]),
				q00r(q0[0]), q01r(q0[1]), q02r(q0[2]),
				q10r(q1[0]), q11r(q1[1]), q12r(q1[2]);
			Rational axr(p10r - p00r), ayr(p11r - p01r), azr(p12r - p02r),
				bxr(q10r - q00r), byr(q11r - q01r), bzr(q12r - q02r);
			Rational xr = ayr * bzr - azr * byr;
			Rational yr = azr * bxr - axr * bzr;
			Rational zr = axr * byr - ayr * bxr;
			Rational ssumr = xr * xr + yr * yr + zr * zr;
			Scalar sum = ssumr.to_double();
			Scalar l = sqrt(sum);
			Scalar xd = xr.to_double(), yd = yr.to_double(), zd = zr.to_double();
			Scalar fx = xd / l, fy = yd / l, fz = zd / l;
			return Vector3(fx, fy, fz);
		}
		const Multiprecision length = ssum.sqrt(ssum);

		x = x / length; y = y / length; z = z / length;

		Scalar fx = x.to_double(), fy = y.to_double(), fz = z.to_double();
		return Vector3(fx, fy, fz);

	}

	void FastEnvelope::prism_bbox(const std::array<Vector3, 12>&prism, Vector3 &min, Vector3& max) {
		int id ;
		int id0;
		bool flag = 0;
		Vector3 p;
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				if (i == j) continue;
				id = i * 8 + j;
				id0 = prism_map[id][0];
				if (id0 == -1) continue;

				for (int k = 0; k < 8; k++) {
					if (k == i || k == j) continue;
					id = i * 8 + k;
					id0 = prism_map[id][0];
					if (id0 == -1) continue;
					id = j * 8 + k;
					id0 = prism_map[id][0];
					if (id0 == -1) continue;
					three_facets_inter_point(
						prism[p_face[i][0]], prism[p_face[i][1]], prism[p_face[i][2]],
						prism[p_face[j][0]], prism[p_face[j][1]], prism[p_face[j][2]],
						prism[p_face[k][0]], prism[p_face[k][1]], prism[p_face[k][2]],
						p);
					if (flag == 0) {
						min = p;
						max = p;
					}
					else {
						for (int t = 0; t < 3; t++) {
							min[t] = std::min(min[t], p[t]);
							max[t] = std::max(max[t], p[t]);
						}
					}
					flag = 1;
				}
			}
		}
		min[0] -= conserve_number;
		min[1] -= conserve_number;
		min[2] -= conserve_number;
		max[0] += conserve_number;
		max[1] += conserve_number;
		max[2] += conserve_number;
	}

	void FastEnvelope::cubic_bbox(const std::array<Vector3, 8>&cubic, Vector3 &min, Vector3& max) {
		int id;
		int id0;
		bool flag = 0;
		Vector3 p;
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				if (i == j) continue;
				id = i * 6 + j;
				id0 = cubic_map[id][0];
				if (id0 == -1) continue;

				for (int k = 0; k < 6; k++) {
					if (k == i || k == j) continue;
					id = i * 6 + k;
					id0 = cubic_map[id][0];
					if (id0 == -1) continue;
					id = j * 6 + k;
					id0 = cubic_map[id][0];
					if (id0 == -1) continue;
					three_facets_inter_point(
						cubic[c_face[i][0]], cubic[c_face[i][1]], cubic[c_face[i][2]],
						cubic[c_face[j][0]], cubic[c_face[j][1]], cubic[c_face[j][2]],
						cubic[c_face[k][0]], cubic[c_face[k][1]], cubic[c_face[k][2]],
						p);
					if (flag == 0) {
						min = p;
						max = p;
					}
					else {
						for (int t = 0; t < 3; t++) {
							min[t] = std::min(min[t], p[t]);
							max[t] = std::max(max[t], p[t]);
						}
					}
					flag = 1;
				}
			}
		}
		min[0] -= conserve_number;
		min[1] -= conserve_number;
		min[2] -= conserve_number;
		max[0] += conserve_number;
		max[1] += conserve_number;
		max[2] += conserve_number;
	}
	void  FastEnvelope::three_facets_inter_point(const Vector3& a0, const Vector3& a1, const Vector3& a2, const Vector3& b0, const Vector3& b1,
		const Vector3& b2, const Vector3& c0, const Vector3& c1, const Vector3& c2, Vector3& p) {
		Matrix3 A;
		Vector3 B;
		Vector3 ae1, ae2, be1, be2, ce1, ce2;
		ae1 = a1 - a0;
		ae2 = a2 - a0;
		be1 = b1 - b0;
		be2 = b2 - b0;
		ce1 = c1 - c0;
		ce2 = c2 - c0;
		A(0,0) = ae1[1] * ae2[2] - ae1[2] * ae2[1];
		A(0,1) = -ae1[0] * ae2[2] + ae1[2] * ae2[0];
		A(0,2) = ae1[0] * ae2[1] - ae1[1] * ae2[0];

		A(1, 0) = be1[1] * be2[2] - be1[2] * be2[1];
		A(1, 1) = -be1[0] * be2[2] + be1[2] * be2[0];
		A(1, 2) = be1[0] * be2[1] - be1[1] * be2[0];

		A(2, 0) = ce1[1] * ce2[2] - ce1[2] * ce2[1];
		A(2, 1) = -ce1[0] * ce2[2] + ce1[2] * ce2[0];
		A(2, 2) = ce1[0] * ce2[1] - ce1[1] * ce2[0];
		B[0] = A(0, 0)*a0[0] + A(0, 1)*a0[1] + A(0, 2)*a0[2];
		B[1] = A(1, 0)*b0[0] + A(1, 1)*b0[1] + A(1, 2)*b0[2];
		B[2] = A(2, 0)*c0[0] + A(2, 1)*c0[1] + A(2, 2)*c0[2];
		p = A.inverse()*B;
	}

}













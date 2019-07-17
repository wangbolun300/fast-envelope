#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Parameters.h>
#include <fastenvelope/Predicates.hpp>
#include <fstream>
#include <istream>
#include <igl/Timer.h>
#include <fastenvelope/ip_filtered.h>
#include <arbitraryprecision/fprecision.h>
#include <fastenvelope/Rational.hpp>


int markhf = 0, markhf1 = 0, i_time = 10, after11 = 0, after12 = 0, after10 = 0, after21 = 0, after22 = 0, after20 = 0;
int recordnumber = 0, recordnumber1 = 0, recordnumber2 = 0, recordnumber3 = 0, recordnumber4 = 0;
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

static const   std::function<int(fastEnvelope::Rational)> check_rational = [](fastEnvelope::Rational v) {
	
	if (v > 0)
		return 1;

	if (v < 0)
		return -1;

	return 0;
};

static const std::function<int(arbitrary_precision::interval<arbitrary_precision::float_precision>)> check_interval =
[](arbitrary_precision::interval<arbitrary_precision::float_precision> v) {
	const auto clazz = v.is_class();
	if (clazz == arbitrary_precision::MIXED || clazz == arbitrary_precision::NO_CLASS)
		return -2;

	if (clazz == arbitrary_precision::POSITIVE)
		return 1;

	if (clazz == arbitrary_precision::NEGATIVE)
		return -1;

	assert(clazz == arbitrary_precision::ZERO);
	return 0;
};




static const std::array<std::array<int, 2>, 3> triseg = {

	{{{0,1}},{{0,2}},{{1,2}}}

};

extern "C++" int tri_tri_intersection_test_3d(fastEnvelope::Scalar p1[3], fastEnvelope::Scalar q1[3], fastEnvelope::Scalar r1[3],
	fastEnvelope::Scalar p2[3], fastEnvelope::Scalar q2[3], fastEnvelope::Scalar r2[3],
	int * coplanar,
	fastEnvelope::Scalar source[3], fastEnvelope::Scalar target[3]);

namespace fastEnvelope {
	//using namespace std;

	void FastEnvelope::print_number() {
		//std::cout << "lpi filter number " << filternumberlpi << " lpi total number " << totalnumberlpi << " percentage " << float(filternumberlpi )/ float(totalnumberlpi) << std::endl;
		//std::cout << "tpi filter number " << filternumbertpi << " tpi total number " << totalnumbertpi << " percentage " << float(filternumbertpi) / float(totalnumbertpi) << std::endl;
		//std::cout << "triangle_intersection filter number " << filternumber1 << " tpi total number " << totalnumber1 << " percentage " << float(filternumber1) / float(totalnumber1) << std::endl;
		//std::cout << "triangle_intersection filter number lpi -2 " << filternumberlpi2 << " percentage " << float(filternumberlpi2) / float(totalnumberlpi) << std::endl;
		//std::cout << "triangle_intersection filter number tpi -2 " << filternumbertpi2 << " percentage " << float(filternumbertpi2) / float(totalnumbertpi+ totalnumber1) << std::endl;
		std::cout << "lpi 1 " << float(after11)/float(after11+ after12+ after10) << " lpi -1 " << after12 / float(after11 + after12 + after10) << " lpi 0 " << after10 / float(after11 + after12 + after10) << " tot  " << after11 + after12 + after10 << std::endl;
		std::cout << "tpi 1 " << after21 / float(after21 + after22 + after20) << " tpi -1 " << after22 / float(after21 + after22 + after20) << " tpi 0 " << after20 / float(after21 + after22 + after20) << " tot  " << after21 + after22 + after20<< std::endl;
	}
	FastEnvelope::FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar& eps, const int& spac)
	{
		get_bb_corners(m_ver, min, max);
		Scalar bbd = (max - min).norm();
		Scalar epsilon = bbd * eps; //eps*bounding box diagnal
		FastEnvelope::BoxGeneration(m_ver, m_faces, envprism, epsilon);
		//build a  hash function



		const Scalar boxlength = std::min(std::min(max[0] - min[0], max[1] - min[1]), max[2] - min[2]) / spac;//TODO a better strategy?
		subx = (max[0] - min[0]) / boxlength, suby = (max[1] - min[1]) / boxlength, subz = (max[2] - min[2]) / boxlength;


		CornerList(envprism, cornerlist);
		std::vector<int> intercell;
		int ct = 0, prismsize = envprism.size();

		prismmap.reserve(spac*spac*spac / 10);
		for (int i = 0; i < cornerlist.size(); i++) {
			BoxFindCells(cornerlist[i][0], cornerlist[i][1], min, max, subx, suby, subz, intercell);
			for (int j = 0; j < intercell.size(); j++) {
				prismmap[intercell[j]].emplace_back(i);
			}
		}
		std::cout << "map size " << prismmap.size() << std::endl;
	}
	bool FastEnvelope::is_outside(const std::array<Vector3, 3> &triangle) const {
		Vector3 tmin, tmax;
		std::vector<int> inumber;
		std::vector<int> intercell;
		std::vector<std::array<Vector3, 12>> interenvprism;
		get_triangle_corners(triangle, tmin, tmax);
		BoxFindCells(tmin, tmax, min, max, subx, suby, subz, intercell);
		inumber.clear();
		for (int j = 0; j < intercell.size(); j++) {
			auto search = prismmap.find(intercell[j]);
			if (search != prismmap.end()) {
				inumber.insert(inumber.end(), search->second.begin(), search->second.end());
			}
		}
		sort(inumber.begin(), inumber.end());
		inumber.erase(unique(inumber.begin(), inumber.end()), inumber.end());
		interenvprism.reserve(inumber.size());
		for (int j = 0; j < inumber.size(); j++) {
			interenvprism.emplace_back(envprism[inumber[j]]);
		}
		return FastEnvelope::FastEnvelopeTestImplicit(triangle, interenvprism);
	}

	/*bool FastEnvelope::is_outside_signal(const std::array<Vector3, 3> &triangle, int &signal) const {
		Vector3 tmin, tmax;
		std::vector<int> inumber;
		std::vector<int> intercell;
		std::vector<std::array<Vector3, 12>> interenvprism;
		get_triangle_corners(triangle, tmin, tmax);
		BoxFindCells(tmin, tmax, min, max, subx, suby, subz, intercell);
		inumber.clear();
		for (int j = 0; j < intercell.size(); j++) {
			auto search = prismmap.find(intercell[j]);
			if (search != prismmap.end()) {
				inumber.insert(inumber.end(), search->second.begin(), search->second.end());
			}
		}
		sort(inumber.begin(), inumber.end());
		inumber.erase(unique(inumber.begin(), inumber.end()), inumber.end());
		interenvprism.reserve(inumber.size());
		for (int j = 0; j < inumber.size(); j++) {
			interenvprism.emplace_back(envprism[inumber[j]]);
		}
		return FastEnvelope::FastEnvelopeTestImplicit_signal(triangle, interenvprism, signal);
	}*/
	void FastEnvelope::print_prisms(const std::array<Vector3, 3> &triangle) const {

		Vector3 tmin, tmax;
		std::vector<int> inumber;
		std::vector<int> intercell;
		std::vector<std::array<Vector3, 12>> interenvprism;
		get_triangle_corners(triangle, tmin, tmax);
		BoxFindCells(tmin, tmax, min, max, subx, suby, subz, intercell);
		inumber.clear();
		for (int j = 0; j < intercell.size(); j++) {
			auto search = prismmap.find(intercell[j]);
			if (search != prismmap.end()) {
				inumber.insert(inumber.end(), search->second.begin(), search->second.end());
			}
		}
		sort(inumber.begin(), inumber.end());
		inumber.erase(unique(inumber.begin(), inumber.end()), inumber.end());
		interenvprism.reserve(inumber.size());
		for (int j = 0; j < inumber.size(); j++) {
			interenvprism.emplace_back(envprism[inumber[j]]);
		}

		std::ofstream fout;
		fout.open("D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\visualprism.txt");
		for (int i = 0; i < interenvprism.size(); i++) {
			for (int j = 0; j < 12; j++) {

				fout << std::setprecision(17) << interenvprism[i][j][0] << " " << interenvprism[i][j][1] << " " << interenvprism[i][j][2] <<std:: endl;

			}
		}

		fout.close();
	}
	bool FastEnvelope::sample_triangle_outside(const std::array<Vector3, 3> &triangle, const Scalar sampleerror) const {
		std::vector<Vector3> ps;
		triangle_sample(triangle, ps, sampleerror);//dd is used for sapmling
		bool out;
		Vector3 tmin, tmax;
		std::vector<int> inumber;
		std::vector<int> intercell;
		std::vector<std::array<Vector3, 12>> interenvprism;
		get_triangle_corners(triangle, tmin, tmax);
		BoxFindCells(tmin, tmax, min, max, subx, suby, subz, intercell);
		inumber.clear();
		for (int j = 0; j < intercell.size(); j++) {
			auto search = prismmap.find(intercell[j]);
			if (search != prismmap.end()) {
				inumber.insert(inumber.end(), search->second.begin(), search->second.end());
			}
		}
		sort(inumber.begin(), inumber.end());
		inumber.erase(unique(inumber.begin(), inumber.end()), inumber.end());
		interenvprism.reserve(inumber.size());
		for (int j = 0; j < inumber.size(); j++) {
			interenvprism.emplace_back(envprism[inumber[j]]);
		}
		std::vector<int> jump;
		for (int i = 0; i < ps.size(); i++) {
			out = point_out_prism(ps[i], interenvprism, jump);

			if (out == true) {

				return 1;

			}

		}

		return 0;

	}
	void FastEnvelope::triangle_sample(const std::array<Vector3, 3> &triangle, std::vector<Vector3>& ps, const Scalar &error) {
		ps.clear();
		Scalar l1 = (triangle[1] - triangle[0]).norm(), l2 = (triangle[2] - triangle[0]).norm(), l3 = (triangle[2] - triangle[1]).norm();//length
		int de = is_triangle_degenerated(triangle);
		if (de==DEGENERATED_POINT) {
			ps.push_back(triangle[0]);
			return;
		}
		if (de == DEGENERATED_SEGMENT) {
			if (triangle[1] - triangle[0] == Vector3(0, 0, 0)) {
				//std::cout << "here1 " << std::endl;
				int t = l2 / error + 1;

				for (int i = 0; i <= t; i++) {
					ps.push_back(triangle[0] + (triangle[2] - triangle[0])*i / t);
				}
				return;
			}
			if (triangle[2] - triangle[0] == Vector3(0, 0, 0)) {
				//std::cout << "here2 " << std::endl;
				int t = l1 / error + 1;

				for (int i = 0; i <= t; i++) {
					ps.push_back(triangle[0] + (triangle[1] - triangle[0])*i / t);
				}
				return;
			}
			if (triangle[2] - triangle[1] == Vector3(0, 0, 0)) {
				//std::cout << "here3 " << std::endl;
				int t = l1 / error + 1;

				for (int i = 0; i <= t; i++) {
					ps.push_back(triangle[0] + (triangle[1] - triangle[0])*i / t);
				}
				return;
			}
		}
		//std::cout << "here " << std::endl;
		/*if (l1 == 0 || l2 == 0 || l3 == 0) {



				return;

		}*/
		int l1s = l1 / error + 1, l2s = l2 / error + 1, l2sn;//subdivided

		Scalar e1 = l1 / l1s, e2 = l2 / l2s, e3;//length of every piece
		Vector3 subl1 = (triangle[1] - triangle[0]) / l1s, subl2 = (triangle[2] - triangle[0]) / l2s,//vector of piece
			temp;

		for (int i = 0; i <= l1s; i++) {
			//Scalar length = (l1s - i)*e1*l2 / l1;// length of this line
			//l2sn = length / e2 + 1;// subdivided of this line
			//Vector3 subl3 = (triangle[2] - triangle[0])*length / l2sn / l2;//vector of each piece
			//for (int j = 0; j <= l2sn; j++) {
			//	temp = subl1 * i + triangle[0] + subl3 * j;
			//	ps.push_back(temp);
			//}


			//segment 0-1 is subdivided into l1s sub-pieces, so is segment 1-2
			Vector3 p1 = triangle[0] + (triangle[1] - triangle[0])*i / l1s, p2 = triangle[0]+(triangle[2]-triangle[0])*i/l1s;
			ps.push_back(p1 );
			ps.push_back(p2);
			Scalar length = (p1 - p2).norm();
			l2sn = length / e2 + 1;
			for (int j = 0; j <= l2sn; j++) {
				ps.push_back(p1 + (p2 - p1)* j / l2sn);
			}
			
		}
		return;
	}




	bool FastEnvelope::FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<std::array<Vector3, 12>>& envprism)
	{

		if (envprism.size() == 0) {
			return 1;
		}
		std::vector<int> jump;
		std::vector<Vector3i> inter_ijk_list;//list of intersected triangle
		bool out;
		int inter, inter1, record1, record2,
			tti;//triangle-triangle intersection
		jump.clear();
		for (int i = 0; i < 3; i++) {
			out = point_out_prism(triangle[i], envprism, jump);
			if (out == true) {
				return 1;
			}
		}


		////////////////////degeneration fix
		int degeneration= is_triangle_degenerated(triangle);
		if (degeneration== DEGENERATED_POINT) {//case 1 degenerate to a point
			return 0;
		}//case 1 degenerate to a point
		if (degeneration == DEGENERATED_SEGMENT) {
			for (int we = 0; we < 3; we++) {//case 2 degenerated as a segment, at most test 2 segments,but still we need to test 3, because
											// of the endpoint-triangle intersection will be ignored
											// the segment is {triangle[triseg[we][0]], triangle[triseg[we][1]]}
					for (int i = 0; i < envprism.size(); i++) {
						for (int j = 0; j < 8; j++) {
							for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
								tti = seg_cut_tri(triangle[triseg[we][0]], triangle[triseg[we][1]], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
								if (tti == CUT_COPLANAR) {
									break;
								}
								if (tti == CUT_EMPTY) {//this is not redundant
									continue;
								}
								jump.clear();
								jump.emplace_back(i);


								inter = Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(triangle[triseg[we][0]], triangle[triseg[we][1]],
									{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);
								

								if (inter == 1) {
									return 1;
								}
								break;

							}
						}
					}
			}//case 2 case 2 degenerated as a segment
			return 0;
		}


		////////////////////////////////degeneration fix over



		for (int i = 0; i < envprism.size(); i++) {
			for (int j = 0; j < 8; j++) {
				for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
					tti = tri_cut_tri_simple(triangle[0], triangle[1], triangle[2], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
					if (tti == CUT_COPLANAR) {
						break;
					}
					if (tti == CUT_EMPTY) {//TODO maybe redundant because we want "float above" case leading to break
						continue;
					}

					record1 = 0;

					jump.clear();
					jump.emplace_back(i);
					for (int k = 0; k < 3; k++) {

						inter = Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(triangle[triseg[k][0]], triangle[triseg[k][1]],
							{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);


						if (inter == 1) {
							return 1;
						}
						record1 = record1 + inter;
					}
					if (record1 >= 4) {
						std::cout << "intersection predicate wrong1, record " << record1 << std::endl;

					}


					for (int e = 0; e < inter_ijk_list.size(); e++) {
						for (int f = inter_ijk_list[e][2]; f < p_triangle[inter_ijk_list[e][1]].size(); f++) {
							tti = tri_cut_tri_simple(triangle[0], triangle[1], triangle[2],
								envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]]);
							if (tti == CUT_COPLANAR) {
								break;
							}
							if (tti == CUT_EMPTY) {
								continue;
							}
							jump.clear();
							jump.emplace_back(inter_ijk_list[e][0]);
							jump.emplace_back(i);

							int inter2 = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_multi_precision(triangle,
								{ {envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]]} },
								{ {envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]} }, envprism, jump);

							if (inter2 == 1) {

								return 1;//out
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


	




	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
		const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		int jm = 0, ori;
		int inter = seg_cut_tri(segpoint0, segpoint1, triangle[0], triangle[1], triangle[2]);
		std::vector<int> INDEX;
		int tot;
		if (inter == CUT_COPLANAR) {// we can not add "CUT_EMPTY" to this, because we use tri-tri intersection, not tri-facet intersection
									//so even if seg cut tri or next tri, seg_cut_tri may returns cut_empty
			return NOT_INTERSECTD;//not intersected
		}
		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}
			}
			INDEX.clear();
			tot = 0;
			for (int j = 0; j < 8; j++) {
				//ftimer2.start();
				ori = ip_filtered::
					orient3D_LPI_filtered(
						segpoint0[0], segpoint0[1], segpoint0[2], segpoint1[0], segpoint1[1], segpoint1[2],
						triangle[0][0], triangle[0][1], triangle[0][2], triangle[1][0], triangle[1][1], triangle[1][2], triangle[2][0], triangle[2][1], triangle[2][2],
						envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2], envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2], envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2]);
				
				/////////////////////////////////////////
				
				
				
				/*if (ori == 0) {
					if (markhf == 0) {
						Rational s00(segpoint0[0]), s01(segpoint0[1]), s02(segpoint0[2]), s10(segpoint1[0]), s11(segpoint1[1]), s12(segpoint1[2]),
							t00(triangle[0][0]), t01(triangle[0][1]), t02(triangle[0][2]), 
							t10(triangle[1][0]), t11(triangle[1][1]), t12(triangle[1][2]), 
							t20(triangle[2][0]), t21(triangle[2][1]), t22(triangle[2][2]),
							e00(envprism[i][p_face[j][0]][0]), e01(envprism[i][p_face[j][0]][1]), e02(envprism[i][p_face[j][0]][2]), 
							e10(envprism[i][p_face[j][1]][0]), e11(envprism[i][p_face[j][1]][1]), e12(envprism[i][p_face[j][1]][2]),
							e20(envprism[i][p_face[j][2]][0]), e21(envprism[i][p_face[j][2]][1]), e22(envprism[i][p_face[j][2]][2]);

						ori = orient3D_LPI_filtered_multiprecision(
							s00, s01, s02, s10, s11, s12,
							t00, t01, t02, t10, t11, t12, t20, t21, t22,
							e00, e01, e02, e10, e11, e12, e20, e21, e22, check_rational);

						



					}
				}*/
				
				///////////////////////////////////////////
				if (ori == 1) {
					break;
				}
				if (ori == -2 || ori == 0) {
					INDEX.push_back(j);
				}
				
				if (ori == -1) {
					tot++;
				}
				if (tot == 8) {

					return IN_PRISM;
				}
			}
			if (ori != 1) {
				for (int k = 0; k < INDEX.size(); k++) {
					Rational s00(segpoint0[0]), s01(segpoint0[1]), s02(segpoint0[2]), s10(segpoint1[0]), s11(segpoint1[1]), s12(segpoint1[2]),
						t00(triangle[0][0]), t01(triangle[0][1]), t02(triangle[0][2]),
						t10(triangle[1][0]), t11(triangle[1][1]), t12(triangle[1][2]),
						t20(triangle[2][0]), t21(triangle[2][1]), t22(triangle[2][2]),
						e00(envprism[i][p_face[INDEX[k]][0]][0]), e01(envprism[i][p_face[INDEX[k]][0]][1]), e02(envprism[i][p_face[INDEX[k]][0]][2]),
						e10(envprism[i][p_face[INDEX[k]][1]][0]), e11(envprism[i][p_face[INDEX[k]][1]][1]), e12(envprism[i][p_face[INDEX[k]][1]][2]),
						e20(envprism[i][p_face[INDEX[k]][2]][0]), e21(envprism[i][p_face[INDEX[k]][2]][1]), e22(envprism[i][p_face[INDEX[k]][2]][2]);

					ori = orient3D_LPI_filtered_multiprecision(
						s00, s01, s02, s10, s11, s12,
						t00, t01, t02, t10, t11, t12, t20, t21, t22,
						e00, e01, e02, e10, e11, e12, e20, e21, e22, check_rational);
					if (ori == 1) after11++;
					if (ori == -1) after12++;
					if (ori == 0) after10++;
					if (ori == -2) std::cout << "impossible thing happens in lpi" << std::endl;
					if (ori == 1 || ori == 0) break;
				}
				if (ori == -1) return IN_PRISM;
			}

		}
		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_multi_precision(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		int jm = 0, ori;
		std::vector<int> INDEX;
		int tot;
		bool in = is_3_triangle_cut(triangle, facet1, facet2);

		if (in == 0) {
			return NOT_INTERSECTD;
		}


		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}

			}
			INDEX.clear();
			tot = 0;
			for (int j = 0; j < 8; j++) {
				ori = ip_filtered::orient3D_TPI_filtered(triangle[0][0], triangle[0][1], triangle[0][2],
					triangle[1][0], triangle[1][1], triangle[1][2],
					triangle[2][0], triangle[2][1], triangle[2][2],
					facet1[0][0], facet1[0][1], facet1[0][2],
					facet1[1][0], facet1[1][1], facet1[1][2],
					facet1[2][0], facet1[2][1], facet1[2][2],
					facet2[0][0], facet2[0][1], facet2[0][2],
					facet2[1][0], facet2[1][1], facet2[1][2],
					facet2[2][0], facet2[2][1], facet2[2][2],
					envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
					envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
					envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2]);

				////////////////////////////////////////////////////////////////////////
				/*if (ori == 0) {
					ori = orient3D_TPI_filtered_multiprecision(
						Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]),
						Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]),
						Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]),
						Rational(facet1[0][0]), Rational(facet1[0][1]), Rational(facet1[0][2]),
						Rational(facet1[1][0]), Rational(facet1[1][1]), Rational(facet1[1][2]),
						Rational(facet1[2][0]), Rational(facet1[2][1]), Rational(facet1[2][2]),
						Rational(facet2[0][0]), Rational(facet2[0][1]), Rational(facet2[0][2]),
						Rational(facet2[1][0]), Rational(facet2[1][1]), Rational(facet2[1][2]),
						Rational(facet2[2][0]), Rational(facet2[2][1]), Rational(facet2[2][2]),
						Rational(envprism[i][p_face[j][0]][0]), Rational(envprism[i][p_face[j][0]][1]), Rational(envprism[i][p_face[j][0]][2]),
						Rational(envprism[i][p_face[j][1]][0]), Rational(envprism[i][p_face[j][1]][1]), Rational(envprism[i][p_face[j][1]][2]),
						Rational(envprism[i][p_face[j][2]][0]), Rational(envprism[i][p_face[j][2]][1]), Rational(envprism[i][p_face[j][2]][2]), check_rational);
				}*/


				////////////////////////////////////////////////////////////////////////

				if (ori == 1) {
					break;
				}
				if (ori == -2 || ori == 0) {
					INDEX.push_back(j);
				}

				if (ori == -1) {
					tot++;
				}
				if (tot == 8) {

					return IN_PRISM;
				}
			}
			if (ori != 1) {
				for (int k = 0; k < INDEX.size(); k++) {
					ori = orient3D_TPI_filtered_multiprecision(
						Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]),
						Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]),
						Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]),
						Rational(facet1[0][0]), Rational(facet1[0][1]), Rational(facet1[0][2]),
						Rational(facet1[1][0]), Rational(facet1[1][1]), Rational(facet1[1][2]),
						Rational(facet1[2][0]), Rational(facet1[2][1]), Rational(facet1[2][2]),
						Rational(facet2[0][0]), Rational(facet2[0][1]), Rational(facet2[0][2]),
						Rational(facet2[1][0]), Rational(facet2[1][1]), Rational(facet2[1][2]),
						Rational(facet2[2][0]), Rational(facet2[2][1]), Rational(facet2[2][2]),
						Rational(envprism[i][p_face[INDEX[k]][0]][0]), Rational(envprism[i][p_face[INDEX[k]][0]][1]), Rational(envprism[i][p_face[INDEX[k]][0]][2]),
						Rational(envprism[i][p_face[INDEX[k]][1]][0]), Rational(envprism[i][p_face[INDEX[k]][1]][1]), Rational(envprism[i][p_face[INDEX[k]][1]][2]),
						Rational(envprism[i][p_face[INDEX[k]][2]][0]), Rational(envprism[i][p_face[INDEX[k]][2]][1]), Rational(envprism[i][p_face[INDEX[k]][2]][2]), check_rational);
					if (ori == 1) after21++;
					if (ori == -1) after22++;
					if (ori == 0) after20++;
					if (ori == 1 || ori == 0) break;
					if (ori == -2) std::cout << "impossible thing happens in tpi" << std::endl;
				}
				if (ori == -1) return IN_PRISM;
			}

		}

		return OUT_PRISM;
	}
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_redundant(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
		const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		int jm = 0, ori;
		int inter = seg_cut_tri(segpoint0, segpoint1, triangle[0], triangle[1], triangle[2]);

		if (inter == CUT_COPLANAR) {// we can not add "CUT_EMPTY" to this, because we use tri-tri intersection, not tri-facet intersection
									//so even if seg cut tri or next tri, seg_cut_tri may returns cut_empty
			return NOT_INTERSECTD;//not intersected
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
				ori = ip_filtered::
					orient3D_LPI_filtered(
						segpoint0[0], segpoint0[1], segpoint0[2], segpoint1[0], segpoint1[1], segpoint1[2],
						triangle[0][0], triangle[0][1], triangle[0][2], triangle[1][0], triangle[1][1], triangle[1][2], triangle[2][0], triangle[2][1], triangle[2][2],
						envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2], envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2], envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2]);
				
				/////////////////////////////////////////



				/*if (ori == 0) {
					if (markhf == 0) {
						Rational s00(segpoint0[0]), s01(segpoint0[1]), s02(segpoint0[2]), s10(segpoint1[0]), s11(segpoint1[1]), s12(segpoint1[2]),
							t00(triangle[0][0]), t01(triangle[0][1]), t02(triangle[0][2]),
							t10(triangle[1][0]), t11(triangle[1][1]), t12(triangle[1][2]),
							t20(triangle[2][0]), t21(triangle[2][1]), t22(triangle[2][2]),
							e00(envprism[i][p_face[j][0]][0]), e01(envprism[i][p_face[j][0]][1]), e02(envprism[i][p_face[j][0]][2]),
							e10(envprism[i][p_face[j][1]][0]), e11(envprism[i][p_face[j][1]][1]), e12(envprism[i][p_face[j][1]][2]),
							e20(envprism[i][p_face[j][2]][0]), e21(envprism[i][p_face[j][2]][1]), e22(envprism[i][p_face[j][2]][2]);

						ori = orient3D_LPI_filtered_multiprecision(
							s00, s01, s02, s10, s11, s12,
							t00, t01, t02, t10, t11, t12, t20, t21, t22,
							e00, e01, e02, e10, e11, e12, e20, e21, e22, check_rational);





					}
				}*/

				///////////////////////////////////////////
				if (ori == -2) {
					return NOT_INTERSECTD;
				}
				if (ori == 1 || ori == 0) {
					break;
				}

				if (j == 7) {

					return IN_PRISM;
				}
			}


		}
		return OUT_PRISM;
	}



	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_redundant(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump)
	{
		int jm = 0, ori;

		bool in = is_3_triangle_cut(triangle, facet1, facet2);
		
		if (in == 0) {
			return NOT_INTERSECTD;
		}


		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}

			}
			for (int j = 0; j < 8; j++) {
				ori = ip_filtered::orient3D_TPI_filtered(triangle[0][0], triangle[0][1], triangle[0][2],
					triangle[1][0], triangle[1][1], triangle[1][2],
					triangle[2][0], triangle[2][1], triangle[2][2],
					facet1[0][0], facet1[0][1], facet1[0][2],
					facet1[1][0], facet1[1][1], facet1[1][2],
					facet1[2][0], facet1[2][1], facet1[2][2],
					facet2[0][0], facet2[0][1], facet2[0][2],
					facet2[1][0], facet2[1][1], facet2[1][2],
					facet2[2][0], facet2[2][1], facet2[2][2],
					envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
					envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
					envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2]);
				
				////////////////////////////////////////////////////////////////////////
				/*if (ori == 0) {
					ori = orient3D_TPI_filtered_multiprecision(
						Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]),
						Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]),
						Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]),
						Rational(facet1[0][0]), Rational(facet1[0][1]), Rational(facet1[0][2]),
						Rational(facet1[1][0]), Rational(facet1[1][1]), Rational(facet1[1][2]),
						Rational(facet1[2][0]), Rational(facet1[2][1]), Rational(facet1[2][2]),
						Rational(facet2[0][0]), Rational(facet2[0][1]), Rational(facet2[0][2]),
						Rational(facet2[1][0]), Rational(facet2[1][1]), Rational(facet2[1][2]),
						Rational(facet2[2][0]), Rational(facet2[2][1]), Rational(facet2[2][2]),
						Rational(envprism[i][p_face[j][0]][0]), Rational(envprism[i][p_face[j][0]][1]), Rational(envprism[i][p_face[j][0]][2]),
						Rational(envprism[i][p_face[j][1]][0]), Rational(envprism[i][p_face[j][1]][1]), Rational(envprism[i][p_face[j][1]][2]),
						Rational(envprism[i][p_face[j][2]][0]), Rational(envprism[i][p_face[j][2]][1]), Rational(envprism[i][p_face[j][2]][2]), check_rational);
				}*/


				////////////////////////////////////////////////////////////////////////

				if (ori == 1 || ori == 0) {

					break;
				}
				if (ori == -2) {
					return NOT_INTERSECTD;
				}
				if (j == 7) {

					return IN_PRISM;
				}
			}
		}

		return OUT_PRISM;
	}
#include<ctime>
	bool FastEnvelope::is_3_triangle_cut(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& f1, const std::array<Vector3, 3>& f2) {
		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];
		//Vector3 n = max;
		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
			std::cout << "Degeneration happens" << std::endl;
			srand(int(time(0)));
			n = { {Vector3(rand(),rand(),rand()) } };
		}
		int o1 = ip_filtered::orient3D_TPI_filtered(
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			f1[0][0], f1[0][1], f1[0][2], f1[1][0], f1[1][1], f1[1][2], f1[2][0], f1[2][1], f1[2][2],
			f2[0][0], f2[0][1], f2[0][2], f2[1][0], f2[1][1], f2[1][2], f2[2][0], f2[2][1], f2[2][2],
			n[0], n[1], n[2],
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2]);

		int o2 = ip_filtered::orient3D_TPI_filtered(
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			f1[0][0], f1[0][1], f1[0][2], f1[1][0], f1[1][1], f1[1][2], f1[2][0], f1[2][1], f1[2][2],
			f2[0][0], f2[0][1], f2[0][2], f2[1][0], f2[1][1], f2[1][2], f2[2][0], f2[2][1], f2[2][2],
			n[0], n[1], n[2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2]);

		int o3 = ip_filtered::orient3D_TPI_filtered(
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			f1[0][0], f1[0][1], f1[0][2], f1[1][0], f1[1][1], f1[1][2], f1[2][0], f1[2][1], f1[2][2],
			f2[0][0], f2[0][1], f2[0][2], f2[1][0], f2[1][1], f2[1][2], f2[2][0], f2[2][1], f2[2][2],
			n[0], n[1], n[2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			triangle[0][0], triangle[0][1], triangle[0][2]);
		

		////////////////////////////////////////////////////////////////////////////
		if (o1 == 0) {
			o1 = orient3D_TPI_filtered_multiprecision(Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]),
				Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]),
				Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]),
				Rational(f1[0][0]), Rational(f1[0][1]), Rational(f1[0][2]), Rational(f1[1][0]), Rational(f1[1][1]), Rational(f1[1][2]), Rational(f1[2][0]), Rational(f1[2][1]), Rational(f1[2][2]),
				Rational(f2[0][0]), Rational(f2[0][1]), Rational(f2[0][2]), Rational(f2[1][0]), Rational(f2[1][1]), Rational(f2[1][2]), Rational(f2[2][0]), Rational(f2[2][1]), Rational(f2[2][2]),
				Rational(n[0]), Rational(n[1]), Rational(n[2]),
				Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]),
				Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]), check_rational);
		}
		if (o2 == 0) {
			o2= orient3D_TPI_filtered_multiprecision(Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]),
				Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]),
				Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]),
				Rational(f1[0][0]), Rational(f1[0][1]), Rational(f1[0][2]), Rational(f1[1][0]), Rational(f1[1][1]), Rational(f1[1][2]), Rational(f1[2][0]), Rational(f1[2][1]), Rational(f1[2][2]),
				Rational(f2[0][0]), Rational(f2[0][1]), Rational(f2[0][2]), Rational(f2[1][0]), Rational(f2[1][1]), Rational(f2[1][2]), Rational(f2[2][0]), Rational(f2[2][1]), Rational(f2[2][2]),
				Rational(n[0]), Rational(n[1]), Rational(n[2]),
				Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]),
				Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]), check_rational);
		}
		if (o3 == 0) {
			o3 = orient3D_TPI_filtered_multiprecision(Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]),
				Rational(triangle[1][0]), Rational(triangle[1][1]), Rational(triangle[1][2]),
				Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]),
				Rational(f1[0][0]), Rational(f1[0][1]), Rational(f1[0][2]), Rational(f1[1][0]), Rational(f1[1][1]), Rational(f1[1][2]), Rational(f1[2][0]), Rational(f1[2][1]), Rational(f1[2][2]),
				Rational(f2[0][0]), Rational(f2[0][1]), Rational(f2[0][2]), Rational(f2[1][0]), Rational(f2[1][1]), Rational(f2[1][2]), Rational(f2[2][0]), Rational(f2[2][1]), Rational(f2[2][2]),
				Rational(n[0]), Rational(n[1]), Rational(n[2]),
				Rational(triangle[2][0]), Rational(triangle[2][1]), Rational(triangle[2][2]),
				Rational(triangle[0][0]), Rational(triangle[0][1]), Rational(triangle[0][2]), check_rational);
		}

		//////////////////////////////////////////////////////////////////////////////



		if (o1 == 1 && o2 == 1 && o3 == 1) {// interpoint in the interior of triangle
			return 1;
		}
		if (o1 == -1 && o2 == -1 && o3 == -1) {// interpoint in the interior of triangle
			return 1;
		}
		return 0;
	}

	int FastEnvelope::tri_cut_tri_simple(const Vector3& p1, const Vector3& p2, const Vector3& p3,//even if only edge connected, regarded as intersected
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
			return CUT_EMPTY;

		if (coplanar == 1) {
			return CUT_COPLANAR;
		}

		if (s[0] == t[0] && s[1] == t[1] && s[2] == t[2])
			return CUT_EMPTY;

		return CUT_FACE;
	}
	
	int FastEnvelope::seg_cut_tri(const Vector3 & seg0, const Vector3 &seg1, const Vector3&t0, const Vector3&t1, const Vector3 &t2) {
		int o1, o2, o3, o4, o5;
		o1 = Predicates::orient_3d(seg0, t0, t1, t2);
		o2 = Predicates::orient_3d(seg1, t0, t1, t2);
		int op = o1 * o2;
		if (op >= 0) {
			return CUT_COPLANAR;//in fact, coplanar and not on this plane
		}

		//s0,t0,t1; s0,t1,t2;s0,t2,t0;
		o3 = Predicates::orient_3d(seg1, seg0, t0, t1);
		o4 = Predicates::orient_3d(seg1, seg0, t1, t2);
		o5 = Predicates::orient_3d(seg1, seg0, t2, t0);
		/*if (o3*o4 == 1 && o3*o5 == 1) {
			return CUT_FACE;
		}*/
		if (o3+o4+o5 >= 2 || o3+o4+o5 <= -2) {// in fact, cut through triangle or segment cut triangle edge
			return CUT_FACE;
		}
		return CUT_EMPTY;
	}





	bool FastEnvelope::point_out_prism(const Vector3 & point, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump)
	{

		int jm = 0, ori;

		for (int i = 0; i < envprism.size(); i++) {
			if (jump.size() > 0) {
				if (i == jump[jm]) {//TODO jump avoid vector
					jm = (jm + 1) >= jump.size() ? 0 : (jm + 1);
					continue;
				}
			}

			for (int j = 0; j < 8; j++) {

				ori = Predicates::orient_3d(envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]], point);
				if (ori == -1 || ori == 0) {
					break;
				}
				if (j == 7) {

					return false;
				}
			}


		}
		return true;
	}

	int FastEnvelope::is_triangle_degenerated(const std::array<Vector3, 3>& triangle) {//TODO temporary version of degeneration detection

		Vector3 a = triangle[0] - triangle[1], b = triangle[0] - triangle[2];
		Vector3 normal = a.cross(b);
		Scalar nbr = normal.norm();

		if (nbr > SCALAR_ZERO) {
			return NOT_DEGENERATED;
		}
		int ori;
		std::array < Vector2, 3> p;
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
				p[i] = to_2d(triangle[i], j);
			}
			ori = Predicates::orient_2d(p[0], p[1], p[2]);
			if (ori != 0) {
				return NERLY_DEGENERATED;
			}
		}

		if (triangle[0][0] != triangle[1][0] || triangle[0][1] != triangle[1][1] || triangle[0][2] != triangle[1][2]) {
			return DEGENERATED_SEGMENT;
		}
		if (triangle[0][0] != triangle[2][0] || triangle[0][1] != triangle[2][1] || triangle[0][2] != triangle[2][2]) {
			return DEGENERATED_SEGMENT;
		}
		return DEGENERATED_POINT;
		//TODO not finished
	}
	void FastEnvelope::BoxGeneration(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, std::vector<std::array<Vector3, 12>>& envprism, const Scalar& epsilon)
	{
		envprism.reserve(m_faces.size());
		Vector3 AB, AC, BC, normal, vector1, ABn;
		Parameters pram;
		std::array<Vector3, 6> polygon;
		std::array<Vector3, 12> polygonoff;
		Scalar  a, b, c,
			tolerance = epsilon/ sqrt(3),

			area;

		for (int i = 0; i < m_faces.size(); i++) {
			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			a = BC.norm();
			b = AC.norm();
			c = AB.norm();
			area = 0.25*sqrt((a + b + c)*(a + b - c)*(a + c - b)*(b + c - a));
			if (area < SCALAR_ZERO) {//TODO fix this with degeneration detection function
				std::cout << "Envelope Triangle Degeneration" << std::endl;
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

			for (int j = 0; j < 6; j++) {
				polygonoff[j] = polygon[j] + normal * tolerance;
			}
			for (int j = 6; j < 12; j++) {
				polygonoff[j] = polygon[j - 6] - normal * tolerance;
			}
			envprism.emplace_back(polygonoff);

		}


	}


	Vector3 FastEnvelope::accurate_normal_vector(const std::array<Vector3, 3> & triangle, const int &digit) {
		using namespace arbitrary_precision;
		const float_precision ax = float_precision(triangle[0][0], digit) - float_precision(triangle[1][0], digit);
		const float_precision ay = float_precision(triangle[0][1], digit) - float_precision(triangle[1][1], digit);
		const float_precision az = float_precision(triangle[0][2], digit) - float_precision(triangle[1][2], digit);

		const float_precision bx = float_precision(triangle[0][0], digit) - float_precision(triangle[2][0], digit);
		const float_precision by = float_precision(triangle[0][1], digit) - float_precision(triangle[2][1], digit);
		const float_precision bz = float_precision(triangle[0][2], digit) - float_precision(triangle[2][2], digit);

		float_precision x = ay * bz - az * by;
		float_precision y = az * bx - ax * bz;
		float_precision z = ax * by - ay * bx;
		std::cout << "precision " << x.precision() << " vs " << ax.precision() << " vs " << digit << std::endl;
		const float_precision length = sqrt(x * x + y * y + z * z);
		x = x / length; y = y / length; z = z / length;

		
		std::cout << "value " << x << " vs " << ax << " vs " << digit << std::endl;
		Scalar fx = x, fy = y, fz = z;
		return Vector3(fx, fy, fz);

	}

	arbitrary_precision::interval<arbitrary_precision::float_precision> FastEnvelope::converting_Scalar_to_arbitary(const Scalar &a, const int &i) {
		arbitrary_precision::interval<arbitrary_precision::float_precision> b;
		b.ref_lower()->precision(i);
		b.ref_upper()->precision(i);
		b = arbitrary_precision::float_precision(a, 16);
		//std::cout << "pre " << b.ref_lower()->precision() << std::endl;
		return b;
	}


}













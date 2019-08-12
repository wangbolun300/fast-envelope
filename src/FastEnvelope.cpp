#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Parameters.h>
#include <fastenvelope/Predicates.hpp>
#include <fstream>
#include <istream>
#include <igl/Timer.h>
#include <fastenvelope/ip_filtered.h>
#include <arbitraryprecision/fprecision.h>
#include <fastenvelope/Rational.hpp>
#include <igl/Timer.h>



int markhf = 0, markhf1 = 0, i_time = 10, after11 = 0, after12 = 0, after10 = 0, after21 = 0, after22 = 0, after20 = 0;
int recordnumber = 0, recordnumber1 = 0, recordnumber2 = 0, recordnumber3 = 0, recordnumber4 = 0;
int go1 = 0, go2 = 0;
igl::Timer timer;
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
static const int p_facenumber = 8;


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

	if (v > 0)
		return 1;

	if (v < 0)
		return -1;

	return 0;
};

static const   std::function<int(fastEnvelope::Multiprecision)> check_Multiprecision = [](fastEnvelope::Multiprecision v) {

	if (v > 0)
		return 1;

	if (v < 0)
		return -1;

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
		std::cout << "lpi 1 " << float(after11) / float(after11 + after12 + after10) << " lpi -1 " << after12 / float(after11 + after12 + after10) << " lpi 0 " << after10 / float(after11 + after12 + after10) << " tot  " << after11 + after12 + after10 << std::endl;
		std::cout << "tpi 1 " << after21 / float(after21 + after22 + after20) << " tpi -1 " << after22 / float(after21 + after22 + after20) << " tpi 0 " << after20 / float(after21 + after22 + after20) << " tot  " << after21 + after22 + after20 << std::endl;
		std::cout << "go1 " << go1 << " go2 " << go2 << std::endl;

	}

	FastEnvelope::FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps, const int spac)
	{
		get_bb_corners(m_ver, min, max);
		Scalar bbd = (max - min).norm();
		Scalar epsilon = bbd * eps; //eps*bounding box diagnal
		BoxGeneration(m_ver, m_faces, envprism, envcubic, epsilon);
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


		get_triangle_corners(triangle, tmin, tmax);
		BoxFindCells(tmin, tmax, min, max, subx, suby, subz, intercell);

		for (int j = 0; j < intercell.size(); j++) {
			auto search = prismmap.find(intercell[j]);
			if (search != prismmap.end()) {
				inumber.insert(inumber.end(), search->second.begin(), search->second.end());
			}
		}
		sort(inumber.begin(), inumber.end());
		inumber.erase(unique(inumber.begin(), inumber.end()), inumber.end());


		return FastEnvelopeTestImplicit(triangle, inumber);
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

				fout << std::setprecision(17) << interenvprism[i][j][0] << " " << interenvprism[i][j][1] << " " << interenvprism[i][j][2] << std::endl;

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

		int jump = -1;
		for (int i = 0; i < ps.size(); i++) {
			out = point_out_prism(ps[i], inumber, jump);

			if (out == true) {

				return 1;

			}

		}

		return 0;

	}
	void FastEnvelope::triangle_sample(const std::array<Vector3, 3> &triangle, std::vector<Vector3>& ps, const Scalar &error) {
		ps.clear();
		Scalar l1 = (triangle[1] - triangle[0]).norm(), l2 = (triangle[2] - triangle[0]).norm(), l3 = (triangle[2] - triangle[1]).norm();//length
		int de = is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);
		if (de == DEGENERATED_POINT) {
			ps.emplace_back(triangle[0]);
			return;
		}
		if (de == DEGENERATED_SEGMENT) {
			if (triangle[1] - triangle[0] == Vector3(0, 0, 0)) {
				//std::cout << "here1 " << std::endl;
				int t = l2 / error + 1;

				for (int i = 0; i <= t; i++) {
					ps.emplace_back(triangle[0] + (triangle[2] - triangle[0])*i / t);
				}
				return;
			}
			if (triangle[2] - triangle[0] == Vector3(0, 0, 0)) {
				//std::cout << "here2 " << std::endl;
				int t = l1 / error + 1;

				for (int i = 0; i <= t; i++) {
					ps.emplace_back(triangle[0] + (triangle[1] - triangle[0])*i / t);
				}
				return;
			}
			if (triangle[2] - triangle[1] == Vector3(0, 0, 0)) {
				//std::cout << "here3 " << std::endl;
				int t = l1 / error + 1;

				for (int i = 0; i <= t; i++) {
					ps.emplace_back(triangle[0] + (triangle[1] - triangle[0])*i / t);
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
			Vector3 p1 = triangle[0] + (triangle[1] - triangle[0])*i / l1s, p2 = triangle[0] + (triangle[2] - triangle[0])*i / l1s;
			ps.emplace_back(p1);
			ps.emplace_back(p2);
			Scalar length = (p1 - p2).norm();
			l2sn = length / e2 + 1;
			for (int j = 0; j <= l2sn; j++) {
				ps.emplace_back(p1 + (p2 - p1)* j / l2sn);
			}

		}
		return;
	}

	
	bool FastEnvelope::FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<int>& prismindex)const

	{

		static const std::function<int(fastEnvelope::Rational)> checker = check_Rational;
		static const std::function<int(fastEnvelope::Multiprecision)> checker1 = check_Multiprecision;
		if (prismindex.size() == 0) {

			return 1;

		}

		int jump1, jump2;

		std::vector<Vector3i> inter_ijk_list;//list of intersected triangle

		bool out;

		int inter, inter1, record1, record2,

			tti;//triangle-triangle intersection


		jump1 = -1;
		for (int i = 0; i < 3; i++) {

			out = point_out_prism(triangle[i], prismindex, jump1);

			if (out == true) {

				return 1;

			}

		}





		////////////////////degeneration fix

		int degeneration = is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

		if (degeneration == DEGENERATED_POINT) {//case 1 degenerate to a point

			return 0;

		}//case 1 degenerate to a point

		if (degeneration == DEGENERATED_SEGMENT) {

			for (int we = 0; we < 3; we++) {//case 2 degenerated as a segment, at most test 2 segments,but still we need to test 3, because

											// of the endpoint-triangle intersection will be ignored

											// the segment is {triangle[triseg[we][0]], triangle[triseg[we][1]]}

				for (int i = 0; i < prismindex.size(); i++) {
					jump1 = prismindex[i];

					for (int j = 0; j < p_facenumber; j++) {

						
						if (j == 0) {
							tti = seg_cut_polygon_6(triangle[triseg[we][0]], triangle[triseg[we][1]],
								envprism[prismindex[i]][0], envprism[prismindex[i]][1], envprism[prismindex[i]][2],
								envprism[prismindex[i]][3], envprism[prismindex[i]][4], envprism[prismindex[i]][5]);
						}
						else {
							if (j == 1) {
								tti = seg_cut_polygon_6(triangle[triseg[we][0]], triangle[triseg[we][1]],
									envprism[prismindex[i]][6], envprism[prismindex[i]][7], envprism[prismindex[i]][8],
									envprism[prismindex[i]][9], envprism[prismindex[i]][10], envprism[prismindex[i]][11]);
							}
							else {
								tti = seg_cut_polygon_4(triangle[triseg[we][0]], triangle[triseg[we][1]],
									envprism[prismindex[i]][p_facepoint[j][0]], envprism[prismindex[i]][p_facepoint[j][1]],
									envprism[prismindex[i]][p_facepoint[j][2]], envprism[prismindex[i]][p_facepoint[j][3]]);
							}
							
						}
						
					
						if (tti == CUT_COPLANAR) {

							continue;

						}

						if (tti == CUT_EMPTY) {//this is not redundant

							continue;

						}

						inter = Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(triangle[triseg[we][0]], triangle[triseg[we][1]],

							envprism[prismindex[i]][p_triangle[j][0][0]], envprism[prismindex[i]][p_triangle[j][0][1]], envprism[prismindex[i]][p_triangle[j][0][2]],
							prismindex, jump1,checker);//rational
						//inter1 = Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(triangle[triseg[we][0]], triangle[triseg[we][1]],

						//	envprism[prismindex[i]][p_triangle[j][0][0]], envprism[prismindex[i]][p_triangle[j][0][1]], envprism[prismindex[i]][p_triangle[j][0][2]],
						//	prismindex, jump1, checker1);//multi
						//if (inter != inter1) {
						//	std::cout << "different happens between rational and multi" << inter << " " << inter1 << std::endl;
						//}
						//TODO can be replaced by box box intersection

						if (inter == 1) {

							return 1;

						}

						break;

					}

				}

			}//case 2 case 2 degenerated as a segment

			return 0;

		}





		////////////////////////////////degeneration fix over







		for (int i = 0; i < prismindex.size(); i++) {
			jump1 = prismindex[i];

			for (int j = 0; j < p_facenumber; j++) {

				for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet

					tti = tri_cut_tri_simple(triangle[0], triangle[1], triangle[2], envprism[prismindex[i]][p_triangle[j][c][0]], envprism[prismindex[i]][p_triangle[j][c][1]], envprism[prismindex[i]][p_triangle[j][c][2]]);
					//TODO can be replaced by box-box intersection
					if (tti == CUT_COPLANAR) {
						break;
					}

					if (tti == CUT_EMPTY) {
						continue;
					}

					record1 = 0;




					for (int k = 0; k < 3; k++) {
						//TODO change to seg-facet cut maybe. we need triangle-facet cut, only way is triangulation.
						tti = seg_cut_tri(triangle[triseg[k][0]], triangle[triseg[k][1]], envprism[prismindex[i]][p_triangle[j][0][0]], envprism[prismindex[i]][p_triangle[j][0][1]], envprism[prismindex[i]][p_triangle[j][0][2]]);
						if (tti == CUT_COPLANAR) {
							continue;
						}
						inter = Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(triangle[triseg[k][0]], triangle[triseg[k][1]],
							envprism[prismindex[i]][p_triangle[j][0][0]], envprism[prismindex[i]][p_triangle[j][0][1]], envprism[prismindex[i]][p_triangle[j][0][2]], 
							prismindex, jump1, checker);

						go1++;

						if (inter == 1) {

							return 1;

						}

						record1 = record1 + inter;

					}

					if (record1 >= 4) {

						std::cout << "intersection predicate wrong1, record " << record1 << std::endl;



					}

					inter_ijk_list.emplace_back(Vector3i(prismindex[i], j, c));
					break;
				}

			}

		}



		int listsize = inter_ijk_list.size();



		for (int i = 1; i < listsize; i++) {
			jump1 = inter_ijk_list[i][0];
			for (int j = 0; j < i; j++) {

				//check triangle{ { envprism[list[i][0]][p_triangle[list[i][1]][list[i][2]][0]], ...[1],...[2] } } and triangle{ { envprism[list[j][0]][p_triangle[list[j][1]][list[j][2]][0]], ...[1],...[2] } }

				//and T

				if (inter_ijk_list[i][0] != inter_ijk_list[j][0]) {//belong to two different prisms
					for (int c1 = 0; c1 < p_triangle[inter_ijk_list[i][1]].size(); c1++) {
						for (int c2 = 0; c2 < p_triangle[inter_ijk_list[j][1]].size(); c2++) {
							tti = tri_cut_tri_simple(
								envprism[inter_ijk_list[i][0]][p_triangle[inter_ijk_list[i][1]][c1][0]],
								envprism[inter_ijk_list[i][0]][p_triangle[inter_ijk_list[i][1]][c1][1]],
								envprism[inter_ijk_list[i][0]][p_triangle[inter_ijk_list[i][1]][c1][2]],
								envprism[inter_ijk_list[j][0]][p_triangle[inter_ijk_list[j][1]][c2][0]],
								envprism[inter_ijk_list[j][0]][p_triangle[inter_ijk_list[j][1]][c2][1]],
								envprism[inter_ijk_list[j][0]][p_triangle[inter_ijk_list[j][1]][c2][2]]);
							if (tti == CUT_FACE || tti == CUT_COPLANAR) break;
						}
						if (tti == CUT_FACE || tti == CUT_COPLANAR) break;
					}
					if (tti == CUT_COPLANAR || tti == CUT_EMPTY) continue;
					//TODO can use box box intersection


					jump2 = inter_ijk_list[j][0];

					int inter2 = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_multi_precision(triangle,

						envprism[inter_ijk_list[i][0]][p_triangle[inter_ijk_list[i][1]][0][0]], envprism[inter_ijk_list[i][0]][p_triangle[inter_ijk_list[i][1]][0][1]], envprism[inter_ijk_list[i][0]][p_triangle[inter_ijk_list[i][1]][0][2]],
																																																										
						envprism[inter_ijk_list[j][0]][p_triangle[inter_ijk_list[j][1]][0][0]], envprism[inter_ijk_list[j][0]][p_triangle[inter_ijk_list[j][1]][0][1]], envprism[inter_ijk_list[j][0]][p_triangle[inter_ijk_list[j][1]][0][2]],

						prismindex, jump1, jump2, checker);

					go2++;

					if (inter2 == 1) {



						return 1;//out

					}

				}

				else {//belong to one same prism



					//find prism_map[list[i][1]*8+list[j][1]][0],prism_map[list[i][1]*8+list[j][1]][1]

					int id = inter_ijk_list[i][1] * 8 + inter_ijk_list[j][1];

					int id0 = prism_map[id][0], id1 = prism_map[id][1];

					if (id0 != -1) {//find map

						tti = seg_cut_tri(envprism[inter_ijk_list[i][0]][id0], envprism[inter_ijk_list[i][0]][id1], triangle[0], triangle[1], triangle[2]);

						if (tti == CUT_COPLANAR || tti == CUT_EMPTY) {
							continue;//not intersected
						}
						//the segment is envprism[inter_ijk_list[i][0]][prism_map[list[i][1]*8+list[j][1]][0]],envprism[inter_ijk_list[i][0]][prism_map[list[i][1]*8+list[j][1]][1]]

						int inter2 = Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(envprism[inter_ijk_list[i][0]][id0], envprism[inter_ijk_list[i][0]][id1], triangle[0], triangle[1], triangle[2], 
							prismindex, jump1, checker);

						go2++;
						if (inter2 == 1) {



							return 1;//out

						}

					}



				}





			}

		}





		return 0;

	}




	struct INDEX {
		int Pi;
		std::vector<int> FACES;
	};
	template<typename T>
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_multi_precision(const Vector3& segpoint0, const Vector3& segpoint1, const Vector3& triangle0,
		const Vector3& triangle1, const Vector3& triangle2, const std::vector<int>& prismindex, const int& jump, const std::function<int(T)>& checker) const {
		int  ori,ori1;
		//int inter = seg_cut_tri(segpoint0, segpoint1, triangle0, triangle1, triangle2);

		//if (inter == CUT_COPLANAR) {// we can not add "CUT_EMPTY" to this, because we use tri-tri intersection, not tri-facet intersection
		//							//so even if seg cut tri or next tri, seg_cut_tri may returns cut_empty
		//	return NOT_INTERSECTD;//not intersected
		//}

		int tot;
		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
			segpoint0[0], segpoint0[1], segpoint0[2],
			segpoint1[0], segpoint1[1], segpoint1[2],
			triangle0[0], triangle0[1], triangle0[2],
			triangle1[0], triangle1[1], triangle1[2],
			triangle2[0], triangle2[1], triangle2[2],
			a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);

		if (precom == false) {
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

			bool premulti = orient3D_LPI_prefilter_multiprecision(s00, s01, s02, s10, s11, s12,
				t00, t01, t02, t10, t11, t12, t20, t21, t22, a11r, a12r, a13r, dr, checker);
			for (int i = 0; i < prismindex.size(); i++) {

				if (prismindex[i] == jump) {

					continue;
				}

				tot = 0;
				for (int j = 0; j < p_facenumber; j++) {


					e00 = envprism[prismindex[i]][p_face[j][0]][0]; e01 = envprism[prismindex[i]][p_face[j][0]][1]; e02 = envprism[prismindex[i]][p_face[j][0]][2];
					e10 = envprism[prismindex[i]][p_face[j][1]][0]; e11 = envprism[prismindex[i]][p_face[j][1]][1]; e12 = envprism[prismindex[i]][p_face[j][1]][2];
					e20 = envprism[prismindex[i]][p_face[j][2]][0]; e21 = envprism[prismindex[i]][p_face[j][2]][1]; e22 = envprism[prismindex[i]][p_face[j][2]][2];
					ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22, checker);
					ori1 = orient3D_LPI_filtered_multiprecision(
						Multiprecision(segpoint0[0]), Multiprecision(segpoint0[1]), Multiprecision(segpoint0[2]),
						Multiprecision(segpoint1[0]), Multiprecision(segpoint1[1]), Multiprecision(segpoint1[2]),
						Multiprecision(triangle0[0]), Multiprecision(triangle0[1]), Multiprecision(triangle0[2]),
						Multiprecision(triangle1[0]), Multiprecision(triangle1[1]), Multiprecision(triangle1[2]),
						Multiprecision(triangle2[0]), Multiprecision(triangle2[1]), Multiprecision(triangle2[2]),
						Multiprecision(envprism[prismindex[i]][p_face[j][0]][0]), Multiprecision(envprism[prismindex[i]][p_face[j][0]][1]), Multiprecision(envprism[prismindex[i]][p_face[j][0]][2]),
						Multiprecision(envprism[prismindex[i]][p_face[j][1]][0]), Multiprecision(envprism[prismindex[i]][p_face[j][1]][1]), Multiprecision(envprism[prismindex[i]][p_face[j][1]][2]),
						Multiprecision(envprism[prismindex[i]][p_face[j][2]][0]), Multiprecision(envprism[prismindex[i]][p_face[j][2]][1]), Multiprecision(envprism[prismindex[i]][p_face[j][2]][2]),
						check_Multiprecision);
					if (ori != ori1) {
						std::cout << "result diff in rat and mul " << ori << " " << ori1 << std::endl;
					}
					if (ori == 1) after11++;
					if (ori == -1) after12++;
					if (ori == 0) after10++;
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
			return OUT_PRISM;
		}


		INDEX index;
		std::vector<INDEX> recompute;
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump) {

				continue;
			}

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
			bool premulti = orient3D_LPI_prefilter_multiprecision(s00, s01, s02, s10, s11, s12,
				t00, t01, t02, t10, t11, t12, t20, t21, t22, a11r, a12r, a13r, dr, checker);


			for (int k = 0; k < recompute.size(); k++) {
				int in1 = recompute[k].Pi;
				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];


					e00 = envprism[in1][p_face[in2][0]][0]; e01 = envprism[in1][p_face[in2][0]][1]; e02 = envprism[in1][p_face[in2][0]][2];
					e10 = envprism[in1][p_face[in2][1]][0]; e11 = envprism[in1][p_face[in2][1]][1]; e12 = envprism[in1][p_face[in2][1]][2];
					e20 = envprism[in1][p_face[in2][2]][0]; e21 = envprism[in1][p_face[in2][2]][1]; e22 = envprism[in1][p_face[in2][2]][2];

					ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22, checker);
					ori1 = orient3D_LPI_filtered_multiprecision(
						Multiprecision(segpoint0[0]), Multiprecision(segpoint0[1]), Multiprecision(segpoint0[2]),
						Multiprecision(segpoint1[0]), Multiprecision(segpoint1[1]), Multiprecision(segpoint1[2]),
						Multiprecision(triangle0[0]), Multiprecision(triangle0[1]), Multiprecision(triangle0[2]),
						Multiprecision(triangle1[0]), Multiprecision(triangle1[1]), Multiprecision(triangle1[2]),
						Multiprecision(triangle2[0]), Multiprecision(triangle2[1]), Multiprecision(triangle2[2]),
					Multiprecision(envprism[in1][p_face[in2][0]][0]), Multiprecision(envprism[in1][p_face[in2][0]][1]),Multiprecision( envprism[in1][p_face[in2][0]][2]),
					Multiprecision(envprism[in1][p_face[in2][1]][0]), Multiprecision(envprism[in1][p_face[in2][1]][1]),Multiprecision( envprism[in1][p_face[in2][1]][2]),
					Multiprecision(envprism[in1][p_face[in2][2]][0]), Multiprecision(envprism[in1][p_face[in2][2]][1]),Multiprecision( envprism[in1][p_face[in2][2]][2]),
						check_Multiprecision);
					if (ori != ori1) {
						std::cout << "result diff in rat and mul " << ori << " " << ori1 << std::endl;
					}
					if (ori == 1) after11++;
					if (ori == -1) after12++;
					if (ori == 0) after10++;

					if (ori == 1 || ori == 0) break;
				}

				if (ori == -1) return IN_PRISM;
			}


		}

		return OUT_PRISM;
	}
	/*
	int FastEnvelope::Implicit_prism_edge_triangle_interpoint_Out_Prism_multi_precision(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
		const std::vector<int>& prismindex, const int& jump) const {
		
		int  ori;
		int tot;
		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		bool precom = ip_filtered::orient3D_LPI_prefilter(// it is boolean maybe need considering
			segpoint0[0], segpoint0[1], segpoint0[2],
			segpoint1[0], segpoint1[1], segpoint1[2],
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);

		if (precom == false) {
			static typeprec
				s00, s01, s02, s10, s11, s12,
				t00, t01, t02,
				t10, t11, t12,
				t20, t21, t22,
				a11r, a12r, a13r, dr,

				e00, e01, e02,
				e10, e11, e12,
				e20, e21, e22;
			s00 = segpoint0[0]; s01 = segpoint0[1]; s02 = segpoint0[2]; s10 = segpoint1[0]; s11 = segpoint1[1]; s12 = segpoint1[2];
			t00 = triangle[0][0]; t01 = triangle[0][1]; t02 = triangle[0][2];
			t10 = triangle[1][0]; t11 = triangle[1][1]; t12 = triangle[1][2];
			t20 = triangle[2][0]; t21 = triangle[2][1]; t22 = triangle[2][2];

			bool premulti = orient3D_LPI_prefilter_multiprecision(s00, s01, s02, s10, s11, s12,
				t00, t01, t02, t10, t11, t12, t20, t21, t22, a11r, a12r, a13r, dr, check_Rational);

			for (int i = 0; i < prismindex.size(); i++) {

				if (prismindex[i] == jump) {

					continue;
				}

				tot = 0;
				for (int j = 0; j < p_facenumber; j++) {
					//ftimer2.start();


					e00 = (envprism[prismindex[i]][p_face[j][0]][0]); e01 = (envprism[prismindex[i]][p_face[j][0]][1]); e02 = (envprism[prismindex[i]][p_face[j][0]][2]);
					e10 = (envprism[prismindex[i]][p_face[j][1]][0]); e11 = (envprism[prismindex[i]][p_face[j][1]][1]); e12 = (envprism[prismindex[i]][p_face[j][1]][2]);
					e20 = (envprism[prismindex[i]][p_face[j][2]][0]); e21 = (envprism[prismindex[i]][p_face[j][2]][1]); e22 = (envprism[prismindex[i]][p_face[j][2]][2]);
					ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22, check_Rational);
					if (ori == 1) after11++;
					if (ori == -1) after12++;
					if (ori == 0) after10++;
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
			return OUT_PRISM;
		}


		INDEX index;
		std::vector<INDEX> recompute;
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump) {

				continue;
			}

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

		if (!recompute.empty()) {
			static typeprec
				s00, s01, s02, s10, s11, s12,
				t00, t01, t02,
				t10, t11, t12,
				t20, t21, t22,
				a11r, a12r, a13r, dr,

				e00, e01, e02,
				e10, e11, e12,
				e20, e21, e22;
			s00 = segpoint0[0]; s01 = segpoint0[1]; s02 = segpoint0[2]; s10 = segpoint1[0]; s11 = segpoint1[1]; s12 = segpoint1[2];
			t00 = triangle[0][0]; t01 = triangle[0][1]; t02 = triangle[0][2];
			t10 = triangle[1][0]; t11 = triangle[1][1]; t12 = triangle[1][2];
			t20 = triangle[2][0]; t21 = triangle[2][1]; t22 = triangle[2][2];
			bool premulti = orient3D_LPI_prefilter_multiprecision(s00, s01, s02, s10, s11, s12,
				t00, t01, t02, t10, t11, t12, t20, t21, t22, a11r, a12r, a13r, dr, check_Rational);


			for (int k = 0; k < recompute.size(); k++) {
				int in1 = recompute[k].Pi;
				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];

					e00 = envprism[in1][p_face[in2][0]][0]; e01 = envprism[in1][p_face[in2][0]][1]; e02 = envprism[in1][p_face[in2][0]][2];
					e10 = envprism[in1][p_face[in2][1]][0]; e11 = envprism[in1][p_face[in2][1]][1]; e12 = envprism[in1][p_face[in2][1]][2];
					e20 = envprism[in1][p_face[in2][2]][0]; e21 = envprism[in1][p_face[in2][2]][1]; e22 = envprism[in1][p_face[in2][2]][2];

					ori = orient3D_LPI_postfilter_multiprecision(a11r, a12r, a13r, dr, s00, s01, s02,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22, check_Rational);
					if (ori == 1) after11++;
					if (ori == -1) after12++;
					if (ori == 0) after10++;
					//if (ori == -2) std::cout << "impossible thing happens in lpi" << std::endl;
					if (ori == 1 || ori == 0) break;
				}

				if (ori == -1) return IN_PRISM;
			}


		}

		return OUT_PRISM;
	}
	*/
template<typename T>
	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_multi_precision(const std::array<Vector3, 3>& triangle,
		const Vector3& facet10, const Vector3& facet11, const Vector3& facet12, const Vector3& facet20, const Vector3& facet21, const Vector3& facet22,
		const std::vector<int>& prismindex, const int& jump1, const int &jump2, const std::function<int(T)>& checker) const {
		int ori;
		int tot;
		bool in = is_3_triangle_cut(triangle, facet10, facet11, facet12, facet20, facet21, facet22,checker);

		if (in == 0) {
			return NOT_INTERSECTD;
		}

		Scalar n1, n2, n3, d, max1, max2, max3, max4, max5, max6, max7;
		bool precom = ip_filtered::orient3D_TPI_prefilter(
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			facet10[0], facet10[1], facet10[2],
			facet11[0], facet11[1], facet11[2],
			facet12[0], facet12[1], facet12[2],
			facet20[0], facet20[1], facet20[2],
			facet21[0], facet21[1], facet21[2],
			facet22[0], facet22[1], facet22[2], d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

		if (precom == false) {
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
			bool premulti = orient3D_TPI_prefilter_multiprecision(t00, t01, t02, t10, t11, t12, t20, t21, t22,
				f100, f101, f102, f110, f111, f112, f120, f121, f122,
				f200, f201, f202, f210, f211, f212, f220, f221, f222,
				dr, n1r, n2r, n3r, checker);

			for (int i = 0; i < prismindex.size(); i++) {
				if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
				tot = 0;
				for (int j = 0; j < p_facenumber; j++) {

					e00 = (envprism[prismindex[i]][p_face[j][0]][0]); e01 = (envprism[prismindex[i]][p_face[j][0]][1]); e02 = (envprism[prismindex[i]][p_face[j][0]][2]);
					e10 = (envprism[prismindex[i]][p_face[j][1]][0]); e11 = (envprism[prismindex[i]][p_face[j][1]][1]); e12 = (envprism[prismindex[i]][p_face[j][1]][2]);
					e20 = (envprism[prismindex[i]][p_face[j][2]][0]); e21 = (envprism[prismindex[i]][p_face[j][2]][1]); e22 = (envprism[prismindex[i]][p_face[j][2]][2]);
					ori = orient3D_TPI_postfilter_multiprecision(dr, n1r, n2r, n3r, e00, e01, e02, e10, e11, e12, e20, e21, e22, checker);
					if (ori == 1) after21++;
					if (ori == -1) after22++;
					if (ori == 0) after20++;
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
			return OUT_PRISM;
		}

		INDEX index;
		std::vector<INDEX> recompute;
		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;


			index.FACES.clear();
			tot = 0;
			for (int j = 0; j < p_facenumber; j++) {
				//ftimer2.start();
				ori = ip_filtered::
					orient3D_TPI_postfilter(
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
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
				index.Pi = prismindex[i];
				recompute.emplace_back(index);
			}
		}

		if (recompute.size() > 0) {
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
			bool premulti = orient3D_TPI_prefilter_multiprecision(t00, t01, t02, t10, t11, t12, t20, t21, t22,
				f100, f101, f102, f110, f111, f112, f120, f121, f122,
				f200, f201, f202, f210, f211, f212, f220, f221, f222,
				dr, n1r, n2r, n3r, checker);



			for (int k = 0; k < recompute.size(); k++) {
				int in1 = recompute[k].Pi;
				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];
					//static Rational e00, e01, ...;
					//e00 = ...; e01 = ...;

					e00 = (envprism[in1][p_face[in2][0]][0]); e01 = (envprism[in1][p_face[in2][0]][1]); e02 = (envprism[in1][p_face[in2][0]][2]);
					e10 = (envprism[in1][p_face[in2][1]][0]); e11 = (envprism[in1][p_face[in2][1]][1]); e12 = (envprism[in1][p_face[in2][1]][2]);
					e20 = (envprism[in1][p_face[in2][2]][0]); e21 = (envprism[in1][p_face[in2][2]][1]); e22 = (envprism[in1][p_face[in2][2]][2]);
					ori = orient3D_TPI_postfilter_multiprecision(dr, n1r, n2r, n3r,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22, checker);

					if (ori == 1) after21++;
					if (ori == -1) after22++;
					if (ori == 0) after20++;

					//if (ori == -2) std::cout << "impossible thing happens in lpi" << std::endl;
					if (ori == 1 || ori == 0) break;
				}

				if (ori == -1) return IN_PRISM;
			}

		}



		return OUT_PRISM;
	}


#include<ctime>
	template<typename T>

	bool FastEnvelope::is_3_triangle_cut(const std::array<Vector3, 3>& triangle,
		const Vector3& facet10, const Vector3& facet11, const Vector3& facet12, const Vector3& facet20, const Vector3& facet21, const Vector3& facet22, const std::function<int(T)> &checker) {
		//make this guy static
		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];

		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
			std::cout << "Degeneration happens" << std::endl;
			//move this guy in constructor and use fixed seed
			srand(int(time(0)));
			n = { {Vector3(rand(),rand(),rand()) } };
		}
		Scalar d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7;
		bool pre = ip_filtered::
			orient3D_TPI_prefilter(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],
				facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
				d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

		if (pre == false) {
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

				nr0, nr1, nr2,

				dr, n1r, n2r, n3r;

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
			bool premulti = orient3D_TPI_prefilter_multiprecision(
				t00, t01, t02, t10, t11, t12, t20, t21, t22,
				f100, f101, f102, f110, f111, f112, f120, f121, f122,
				f200, f201, f202, f210, f211, f212, f220, f221, f222,
				dr, n1r, n2r, n3r, checker);
			if (premulti == false) return false;

			int o1 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t00, t01, t02,
				t10, t11, t12, checker);
			/*if (o1 == 1) after21++;
			if (o1 == -1) after22++;
			if (o1 == 0) after20++;*/
			if (o1 == 0) return false;

			int o2 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t10, t11, t12,
				t20, t21, t22, checker);
			/*if (o2 == 1) after21++;
			if (o2 == -1) after22++;
			if (o2 == 0) after20++;*/
			if (o2 == 0 || o1 + o2 == 0) return false;

			int o3 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t20, t21, t22,
				t00, t01, t02, checker);
			/*if (o3 == 1) after21++;
			if (o3 == -1) after22++;
			if (o3 == 0) after20++;*/
			if (o3 == 0 || o1 + o3 == 0) {
				return false;
			}

			return true;
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

			nr0, nr1, nr2,

			dr, n1r, n2r, n3r;

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

			premulti = orient3D_TPI_prefilter_multiprecision(
				t00, t01, t02, t10, t11, t12, t20, t21, t22,
				f100, f101, f102, f110, f111, f112, f120, f121, f122,
				f200, f201, f202, f210, f211, f212, f220, f221, f222,
				dr, n1r, n2r, n3r, checker);
			o1 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t00, t01, t02,
				t10, t11, t12, checker);
			/*if (o1 == 1) after21++;
			if (o1 == -1) after22++;
			if (o1 == 0) after20++;*/
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
				premulti = orient3D_TPI_prefilter_multiprecision(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					dr, n1r, n2r, n3r, checker);
			}
			o2 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t10, t11, t12,
				t20, t21, t22, checker);
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
				premulti = orient3D_TPI_prefilter_multiprecision(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					dr, n1r, n2r, n3r, checker);
			}
			o3 = orient3D_TPI_postfilter_multiprecision(
				dr, n1r, n2r, n3r,
				nr0, nr1, nr2,
				t20, t21, t22,
				t00, t01, t02, checker);
			/*if (o3 == 1) after21++;
			if (o3 == -1) after22++;
			if (o3 == 0) after20++;*/
		}
		if (o3 == 0 || o1 + o3 == 0) return false;

		return true;
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
		if (o3 + o4 + o5 >= 2 || o3 + o4 + o5 <= -2) {// in fact, cut through triangle or segment cut triangle edge
			return CUT_FACE;
		}
		return CUT_EMPTY;
	}

	int FastEnvelope::seg_cut_polygon_4(const Vector3 & seg0, const Vector3 &seg1, const Vector3&t0, const Vector3&t1, const Vector3 &t2, const Vector3 &t3) {
		
		int o1, o2, o3, o4, o5,o6;
		o1 = Predicates::orient_3d(seg0, t0, t1, t2);
		o2 = Predicates::orient_3d(seg1, t0, t1, t2);
		int op = o1 * o2;
		if (op >= 0) {
			return CUT_COPLANAR;//in fact, coplanar and not on this plane
		}
		//now cutted
		//s0,t0,t1; s0,t1,t2;s0,t2,t0;
		o3 = Predicates::orient_3d(seg1, seg0, t0, t1);
		o4 = Predicates::orient_3d(seg1, seg0, t1, t2);
		o5 = Predicates::orient_3d(seg1, seg0, t2, t3);
		o6 = Predicates::orient_3d(seg1, seg0, t3, t0);
		/*if (o3*o4 == 1 && o3*o5 == 1) {
			return CUT_FACE;
		}*/
		if (o3 + o4 + o5 + o6 >= 3 || o3 + o4 + o5 + o6 <= -3) {// in fact, cut through triangle or segment cut triangle edge
			return CUT_FACE;
		}
		return CUT_EMPTY;
	}


	int FastEnvelope::seg_cut_polygon_6(const Vector3 & seg0, const Vector3 &seg1, const Vector3&t0, const Vector3&t1, const Vector3 &t2, const Vector3 &t3, const Vector3 &t4, const Vector3 &t5) {

		int o1, o2, o3, o4, o5, o6,o7,o8;
		o1 = Predicates::orient_3d(seg0, t0, t1, t2);
		o2 = Predicates::orient_3d(seg1, t0, t1, t2);
		int op = o1 * o2;
		if (op >= 0) {
			return CUT_COPLANAR;//in fact, coplanar and not on this plane
		}
		//now cutted
		//s0,t0,t1; s0,t1,t2;s0,t2,t0;
		o3 = Predicates::orient_3d(seg1, seg0, t0, t1);
		o4 = Predicates::orient_3d(seg1, seg0, t1, t2);
		o5 = Predicates::orient_3d(seg1, seg0, t2, t3);
		o6 = Predicates::orient_3d(seg1, seg0, t3, t4);
		o7 = Predicates::orient_3d(seg1, seg0, t4, t5);
		o8 = Predicates::orient_3d(seg1, seg0, t5, t0);


		/*if (o3*o4 == 1 && o3*o5 == 1) {
			return CUT_FACE;
		}*/
		if (o3 + o4 + o5 + o6 + o7 + o8 >= 5 || o3 + o4 + o5 + o6 + o7 + o8 <= -5) {// in fact, cut through triangle or segment cut triangle edge
			return CUT_FACE;
		}
		return CUT_EMPTY;
	}


	bool FastEnvelope::point_out_prism(const Vector3 & point, const std::vector<int>& prismindex, const int& jump)const
	{

		int  ori;

		for (int i = 0; i < prismindex.size(); i++) {
			if (prismindex[i] == jump) continue;

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
				envbox.push_back(box);
				continue;
			}
			if (de == DEGENERATED_SEGMENT) {
				std::cout << "Envelope Triangle Degeneration- Segment" << std::endl;
				Scalar length1 = AB.norm(), length2 = AC.norm(), length3 = BC.norm();
				if (length1 >= length2 && length1 >= length3) {
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], tolerance, box);
					envbox.push_back(box);
				}
				if (length2 >= length1 && length2 >= length3) {
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][2]], tolerance, box);
					envbox.push_back(box);
				}
				if (length3 >= length1 && length3 >= length2) {
					seg_cube(m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tolerance, box);
					envbox.push_back(box);
				}
				continue;
			}
			if (de == NERLY_DEGENERATED) {
				std::cout << "Envelope Triangle Degeneration- Nearly" << std::endl;

				//normal = accurate_normal_vector(AB, AC);
				//vector1= accurate_normal_vector(AB, normal);
				continue;
			}
			else {
				normal = AB.cross(AC).normalized();
				vector1 = AB.cross(normal).normalized();
			}

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
		envbox[3] = p2 + width * (v + v1 - v2);
		envbox[4] = p1 + width * (-v + v1 + v2);
		envbox[5] = p1 + width * (-v - v1 + v2);
		envbox[6] = p1 + width * (-v - v1 - v2);
		envbox[7] = p1 + width * (-v + v1 - v2);
	}

	Vector3 FastEnvelope::accurate_normal_vector(const Vector3 & p, const Vector3 & q) {
		
		const Multiprecision ax = p[0];
		const Multiprecision ay = p[1];
		const Multiprecision az = p[2];
								   
		const Multiprecision bx = q[0];
		const Multiprecision by = q[1];
		const Multiprecision bz = q[2];

		Multiprecision x = ay * bz - az * by;
		Multiprecision y = az * bx - ax * bz;
		Multiprecision z = ax * by - ay * bx;
		Multiprecision ssum = x * x + y * y + z * z;
		const Multiprecision length = ssum.sqrt(ssum);
		x = x / length; y = y / length; z = z / length;
		
		Scalar fx = x.to_double(), fy = y.to_double(), fz = z.to_double();
		return Vector3(fx, fy, fz);

	}

	


}













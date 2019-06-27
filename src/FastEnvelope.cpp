#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Parameters.h>
#include <fastenvelope/Predicates.hpp>
#include <fstream>
#include <istream>
#include <igl/Timer.h>
#include <fastenvelope/ip_filtered.h>
#include <arbitraryprecision/fprecision.h>



int markhf = 0, markhf1=0;
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

static const std::array<std::array<int, 2>, 3> triseg = {

	{{{0,1}},{{0,2}},{{1,2}}}

};

extern "C++" int tri_tri_intersection_test_3d(fastEnvelope::Scalar p1[3], fastEnvelope::Scalar q1[3], fastEnvelope::Scalar r1[3],
	fastEnvelope::Scalar p2[3], fastEnvelope::Scalar q2[3], fastEnvelope::Scalar r2[3],
	int * coplanar,
	fastEnvelope::Scalar source[3], fastEnvelope::Scalar target[3]);

namespace fastEnvelope {
	//using namespace std;


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

		if (l1 == 0 && l2 == 0 && l3 == 0) {
			ps.push_back(triangle[0]);
			return;
		}
		if (l1 == 0) {

			int t = l2 / error + 1;
			Vector3 vct = (triangle[2] - triangle[0]) / t;
			for (int i = 0; i <= t; i++) {
				ps.push_back(triangle[0] + vct * i);
			}
			return;
		}
		if (l2 == 0) {

			int t = l1 / error + 1;
			Vector3 vct = (triangle[1] - triangle[0]) / t;
			for (int i = 0; i <= t; i++) {
				ps.push_back(triangle[0] + vct * i);
			}
			return;
		}
		if (l3 == 0) {

			int t = l1 / error + 1;
			Vector3 vct = (triangle[1] - triangle[0]) / t;
			for (int i = 0; i <= t; i++) {
				ps.push_back(triangle[0] + vct * i);
			}
			return;
		}
		/*if (l1 == 0 || l2 == 0 || l3 == 0) {



				return;

		}*/
		int l1s = l1 / error + 1, l2s = l2 / error + 1, l2sn;//subdivided

		Scalar e1 = l1 / l1s, e2 = l2 / l2s, e3;//length of every piece
		Vector3 subl1 = (triangle[1] - triangle[0]) / l1s, subl2 = (triangle[2] - triangle[0]) / l2s,//vector of piece
			temp;

		for (int i = 0; i <= l1s; i++) {
			Scalar length = (l1s - i)*e1*l2 / l1;// length of this line
			l2sn = length / e2 + 1;// subdivided of this line
			Vector3 subl3 = (triangle[2] - triangle[0])*length / l2sn / l2;//vector of each piece
			for (int j = 0; j <= l2sn; j++) {
				temp = subl1 * i + triangle[0] + subl3 * j;
				ps.push_back(temp);
			}
		}
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
			for (int we = 0; we < 2; we++) {//case 2 degenerated as a segment, at most test 2 segments
											// the segment is {triangle[triseg[we][0]], triangle[triseg[we][1]]}
					for (int i = 0; i < envprism.size(); i++) {
						for (int j = 0; j < 8; j++) {

							for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
								tti = seg_cut_tri(triangle[triseg[we][0]], triangle[triseg[we][1]], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
								if (tti == CUT_COPLANAR) {
									break;
								}
								if (tti == CUT_EMPTY) {
									continue;
								}
								jump.clear();
								jump.emplace_back(i);


								inter = Implicit_Seg_Facet_interpoint_Out_Prism_redundant(triangle[triseg[we][0]], triangle[triseg[we][1]],
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
					if (tti == CUT_EMPTY) {
						continue;
					}

					record1 = 0;

					jump.clear();
					jump.emplace_back(i);
					for (int k = 0; k < 3; k++) {

						inter = Implicit_Seg_Facet_interpoint_Out_Prism_redundant(triangle[triseg[k][0]], triangle[triseg[k][1]],
							{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);


						if (inter == 1) {
							return 1;
						}
						record1 = record1 + inter;
					}
					if (record1 >= 4) {
						std::cout << "intersection predicate wrong, record " << record1 << std::endl;

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

							int inter2 = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_redundant(triangle,
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
	/*
	bool FastEnvelope::FastEnvelopeTestImplicit_signal(const std::array<Vector3, 3> &triangle,
		const std::vector<std::array<Vector3, 12>>& envprism, int &signal)
	{

		if (envprism.size() == 0) {
			return 1;
		}

		std::vector<int> jump;
		std::vector<Vector3i> inter_ijk_list;//list of intersected triangle
		bool out;
		int inter, inter1, record1, record2,
			tti;//triangle-triangle intersection
		Scalar de[3];
		jump.clear();

		for (int i = 0; i < 3; i++) {
			out = point_out_prism(triangle[i], envprism, jump);
			if (out == true) {
				return 1;
			}
		}


		////////////////////degeneration fix
		for (int i = 0; i < 3; i++) {
			de[i] = (triangle[triseg[i][0]] - triangle[triseg[i][1]]).norm();

		}
		//if (0.25*sqrt((de[0] + de[1] + de[2])*(de[0] + de[1] - de[2])*(de[0] + de[2] - de[1])*(de[1] + de[2] - de[0])) < SCALAR_ZERO) {
		if ((de[0] + de[1] + de[2])*(de[0] + de[1] - de[2])*(de[0] + de[2] - de[1])*(de[1] + de[2] - de[0]) == 0) {
			if (signal == 1) {
				std::cout << "de " << de[0] << " " << de[1] << " " << de[2] << std::endl;

				std::cout << "degeneration happens" << std::endl;
			}
			if (de[0] == 0 && de[1] == 0) {//case 1 degenerate to a point
				return out = point_out_prism(triangle[0], envprism, jump);
			}//case 1 degenerate to a point

			for (int we = 0; we < 3; we++) {//case 2 two points are the same point
				if (de[we] == 0) {
					int  wt = (we + 1) % 3;// the segment is triangle[triseg[wt][0]], triangle[triseg[wt][1]]
					for (int i = 0; i < envprism.size(); i++) {
						for (int j = 0; j < 8; j++) {

							for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
								tti = seg_cut_tri(triangle[triseg[wt][0]], triangle[triseg[wt][1]], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
								if (tti == CUT_COPLANAR) {//TODO better add a parallel detection. for now not need because
									break;
								}
								if (tti == CUT_EMPTY) {
									continue;
								}
								jump.clear();
								jump.emplace_back(i);


								inter = Implicit_Seg_Facet_interpoint_Out_Prism_redundant(triangle[triseg[wt][0]], triangle[triseg[wt][1]],
									{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);

								if (signal == 1) {
									std::ofstream fout;
									fout.open("D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\seg.txt");
									fout << std::setprecision(17) << triangle[triseg[wt][0]][0] << " " << triangle[triseg[wt][0]][1] << " " << triangle[triseg[wt][0]][2]
										<< " " << triangle[triseg[wt][1]][0] << " " << triangle[triseg[wt][1]][1]
										<< " " << triangle[triseg[wt][1]][2] << " " << i << " " << j << " " << c << std::endl;
									fout.close();
								}

								if (inter == 1) {
									return 1;
								}

								break;

							}
						}
					}
					return 0;
				}
			}//case 2 two points are the same point

			//TODO maybe we can delete case 2, but maybe having case 2 can be faster
			//case 3 three points on one line
			Scalar maxi = de[0];
			int wt = 0;
			for (int i = 1; i < 3; i++) {
				if (de[i] > maxi) {
					maxi = de[i];
					wt = i;
				}
			}

			// the segment is triangle[triseg[wt][0]], triangle[triseg[wt][1]]
			for (int i = 0; i < envprism.size(); i++) {
				for (int j = 0; j < 8; j++) {

					for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
						tti = seg_cut_tri(triangle[triseg[wt][0]], triangle[triseg[wt][1]], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
						if (tti == CUT_COPLANAR) {//TODO better add a parallel detection. for now not need because
							break;
						}
						if (tti == CUT_EMPTY) {
							continue;
						}
						jump.clear();
						jump.emplace_back(i);


						inter = Implicit_Seg_Facet_interpoint_Out_Prism_redundant(triangle[triseg[wt][0]], triangle[triseg[wt][1]],
							{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);


						if (inter == 1) {
							return 1;
						}

						break;

					}
				}
			}
			return 0;

			//case 3 three points on one line

		}
		////////////////////////////////degeneration fix over



		for (int i = 0; i < envprism.size(); i++) {
			for (int j = 0; j < 8; j++) {
				for (int c = 0; c < p_triangle[j].size(); c++) {//each triangle of the facet
					tti = tri_cut_tri_simple(triangle[0], triangle[1], triangle[2], envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]);
					if (tti == CUT_COPLANAR) {
						break;
					}
					if (tti == CUT_EMPTY) {
						continue;
					}

					record1 = 0;

					jump.clear();
					jump.emplace_back(i);
					for (int k = 0; k < 3; k++) {

						inter = Implicit_Seg_Facet_interpoint_Out_Prism_redundant(triangle[triseg[k][0]], triangle[triseg[k][1]],
							{ { envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]] } }, envprism, jump);


						if (inter == 1) {
							return 1;
						}
						record1 = record1 + inter;
					}
					if (record1 >= 4) {
						std::cout << "intersection predicate wrong, record " << record1 << std::endl;

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
							inter1 = Implicit_Tri_Facet_Facet_interpoint_Out_Prism(triangle,
								{ {envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]]} },
								{ {envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]} }, envprism, jump);
							//////////////////////////////////////////////////////////////////////////////
							int inter2 = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_M(triangle,
								{ {envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][0]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][1]], envprism[inter_ijk_list[e][0]][p_triangle[inter_ijk_list[e][1]][f][2]]} },
								{ {envprism[i][p_triangle[j][c][0]], envprism[i][p_triangle[j][c][1]], envprism[i][p_triangle[j][c][2]]} }, envprism, jump);
							if (inter1 != inter2) {

								//cout << "difference in 3 triangle in-prism test, number"<<recordnumber2<<" " << inter1 << " " << inter2 << endl;
								recordnumber2++;
								if (inter2 == 1) {
									recordnumber4++;
									//cout << "test on the new facets " << recordnumber4 << endl;
								}
							}


							recordnumber3++;
							//cout << "see the total test number, number" << recordnumber3 << endl;

						//////////////////////////////////////////////////////////////////////////////
							if (signal == 1) {
								std::cout << "here in 3 triangle intersection, in or not "<< inter2<< std::endl;

							}
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
	*/
/*
bool FastEnvelope::is_seg_facet_intersection(const double& px, const double& py, const double& pz,
		const double& qx, const double& qy, const double& qz,
		const double& rx, const double& ry, const double& rz,
		const double& sx, const double& sy, const double& sz,
		const double& tx, const double& ty, const double& tz,
		double& a11, double& a12, double& a13,
		double& a21, double& a22, double& a23,
		double& a31, double& a32, double& a33,
		double& px_rx, double& py_ry, double& pz_rz,
		double& d, double& n) {//TODO do we have redundant calculation between this and tri_cut_tri() ?
		::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
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
		if (::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID)) return 0; // Fast reject in case of under/overflow
		if (d != 0) {
			if (-1 * n / d >= 0 && -1 * n / d <= 1) {
				return 1;
			}
		}
		return 0;
	}
*/

/*
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
					if (recordnumber < 1000) {
						//std::cout << "number " << recordnumber << std::endl;

						//std::cout << "ori and ori1 " << ori << " " << ori1 << std::endl;
						ori = orient3D_LPI(
							segpoint0[0], segpoint0[1], segpoint0[2],
							envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
							envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
							envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2],
							a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n);
						recordnumber++;
					}
					markhf = 0;
				}
				///////////////////////////////////////////
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
*/


	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_redundant(const Vector3& segpoint0, const Vector3& segpoint1, const std::array<Vector3, 3>& triangle,
		const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		int jm = 0, ori;
		int inter = seg_cut_tri(segpoint0, segpoint1, triangle[0], triangle[1], triangle[2]);

		if (inter == CUT_COPLANAR) {
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
				/*recordnumber1++;
				if (recordnumber1>8000000&&markhf1 == 0) {
					std::cout << "report you the  ssssssssssssssssss" << std::endl;
					markhf1 = 1;
				}*/

				//ori1 = -1 * Predicates::orient_3d(envprism[i][p_face[j][0]], envprism[i][p_face[j][1]], envprism[i][p_face[j][2]], segpoint0 + (segpoint0 - segpoint1)*n / d);//because n is -n
				/*if (ori != ori1) {
					markhf = 0;
					if (recordnumber < 1000) {
						//std::cout << "number " << recordnumber << std::endl;

						//std::cout << "ori and ori1 " << ori << " " << ori1 << std::endl;
						ori = orient3D_LPI(
							segpoint0[0], segpoint0[1], segpoint0[2],
							envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
							envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
							envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2],
							a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n);
						recordnumber++;
					}
					markhf = 0;
				}
				*/
				///////////////////////////////////////////
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




	/*
	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2,
		const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump) {
		Eigen::Matrix<Scalar, 3, 3> A, AT, ATA, A1;
		Eigen::Matrix<Scalar, 3, 1> B;
		Eigen::Matrix<Scalar, 4, 4> C;
		Eigen::Matrix<Scalar, 2, 3> Prj;
		Vector3 p3d;
		Vector2 t0, t1, t2, p;
		int rcd = 0, jm = 0, ori;
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

		A(0) = A(0)*mv;
		A(1) = A(1)*mv;
		A(2) = A(2)*mv;
		ad = det3x3(A(0, 0), A(0, 1), A(0, 2),
			A(1, 0), A(1, 1), A(1, 2),
			A(2, 0), A(2, 1), A(2, 2));
		fad = fabs(ad);
		if (fad < SCALAR_ZERO_3) {
			return 2;
		}
		fadi = fad < 1 ? 1 / fad : 1;
		A(0) = A(0)*fadi;
		A(1) = A(1)*fadi;
		A(2) = A(2)*fadi;
		B[0] = A(0, 0)*triangle[0][0] + A(0, 1)*triangle[0][1] + A(0, 2)*triangle[0][2];
		B[1] = A(1, 0)*facet1[0][0] + A(1, 1)*facet1[0][1] + A(1, 2)*facet1[0][2];
		B[2] = A(2, 0)*facet2[0][0] + A(2, 1)*facet2[0][1] + A(2, 2)*facet2[0][2];



		//////////////////////////////////////////////////////////////////////////

		Scalar m11, m12, m13, d;
		std::array<Vector3, 3>tri = triangle;//TODO really need to do this?
		std::array<Vector3, 3>fct1 = facet1;
		std::array<Vector3, 3>fct2 = facet2;


		bool inter = is_3triangle_intersect(
			tri[0][0], tri[0][1], tri[0][2],
			tri[1][0], tri[1][1], tri[1][2],
			tri[2][0], tri[2][1], tri[2][2],
			fct1[0][0], fct1[0][1], fct1[0][2],
			fct1[1][0], fct1[1][1], fct1[1][2],
			fct1[2][0], fct1[2][1], fct1[2][2],
			fct2[0][0], fct2[0][1], fct2[0][2],
			fct2[1][0], fct2[1][1], fct2[1][2],
			fct2[2][0], fct2[2][1], fct2[2][2],
			m11, m12, m13, d);

		////////////////////////////////////////////////////////////////////////////

		Vector3 f1, f2, f3;

		f1 = (B - A * triangle[1]).cross(A*triangle[1] - A * triangle[0]);
		f2 = (B - A * triangle[2]).cross(A*triangle[2] - A * triangle[1]);
		f3 = (B - A * triangle[0]).cross(A*triangle[0] - A * triangle[2]);

		bool in = f1.dot(f2) > 0 && f1.dot(f3) > 0 ? 1 : 0;//TODO consider the triangle degeneration
		/////////////////////////////////////////////////
		if (inter != in) {
			//std::cout << "intersect test diff in 3 triangles "<< in<<" "<<inter << endl;
		}
		//////////////////////////////////////////////////
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
				////////////////////////////////////////////////////////////////////////////////
				int ori1 = orient3D_TPI(
					tri[0][0], tri[0][1], tri[0][2],
					tri[1][0], tri[1][1], tri[1][2],
					tri[2][0], tri[2][1], tri[2][2],
					fct1[0][0], fct1[0][1], fct1[0][2],
					fct1[1][0], fct1[1][1], fct1[1][2],
					fct1[2][0], fct1[2][1], fct1[2][2],
					fct2[0][0], fct2[0][1], fct2[0][2],
					fct2[1][0], fct2[1][1], fct2[1][2],
					fct2[2][0], fct2[2][1], fct2[2][2],
					envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
					envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
					envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2],
					m11, m12, m13, d);
				if (ori != ori1) {
					//cout << "ori in 3 triangles different, " << ori << " " << ori1 << endl;
				}
				////////////////////////////////////////////////////////////////////////////////

				if (ori == 1 || ori == 0) {

					break;
				}
				if (j == 7) {

					return 0;
				}
			}
		}

		return 1;
	}
	*/

/*
int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_M(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump)
	{
		int jm = 0, ori;
		Scalar m11, m12, m13, d;
		std::array<Vector3, 3>tri = triangle;//TODO really need to do this?
		std::array<Vector3, 3>fct1 = facet1;
		std::array<Vector3, 3>fct2 = facet2;


		bool in = is_3triangle_intersect(
			tri[0][0], tri[0][1], tri[0][2],
			tri[1][0], tri[1][1], tri[1][2],
			tri[2][0], tri[2][1], tri[2][2],
			fct1[0][0], fct1[0][1], fct1[0][2],
			fct1[1][0], fct1[1][1], fct1[1][2],
			fct1[2][0], fct1[2][1], fct1[2][2],
			fct2[0][0], fct2[0][1], fct2[0][2],
			fct2[1][0], fct2[1][1], fct2[1][2],
			fct2[2][0], fct2[2][1], fct2[2][2],
			m11, m12, m13, d);
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
				ori = orient3D_TPI(
					tri[0][0], tri[0][1], tri[0][2],
					tri[1][0], tri[1][1], tri[1][2],
					tri[2][0], tri[2][1], tri[2][2],
					fct1[0][0], fct1[0][1], fct1[0][2],
					fct1[1][0], fct1[1][1], fct1[1][2],
					fct1[2][0], fct1[2][1], fct1[2][2],
					fct2[0][0], fct2[0][1], fct2[0][2],
					fct2[1][0], fct2[1][1], fct2[1][2],
					fct2[2][0], fct2[2][1], fct2[2][2],
					envprism[i][p_face[j][0]][0], envprism[i][p_face[j][0]][1], envprism[i][p_face[j][0]][2],
					envprism[i][p_face[j][1]][0], envprism[i][p_face[j][1]][1], envprism[i][p_face[j][1]][2],
					envprism[i][p_face[j][2]][0], envprism[i][p_face[j][2]][1], envprism[i][p_face[j][2]][2],
					m11, m12, m13, d);
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

*/

	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_redundant(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& facet1, const std::array<Vector3, 3>& facet2, const std::vector<std::array<Vector3, 12>>& envprism, const std::vector<int>& jump)
	{
		int jm = 0, ori;

		bool in = is_3_triangle_cut(triangle, facet1, facet2);
		/*is_3triangle_intersect(
			tri[0][0], tri[0][1], tri[0][2],
			tri[1][0], tri[1][1], tri[1][2],
			tri[2][0], tri[2][1], tri[2][2],
			fct1[0][0], fct1[0][1], fct1[0][2],
			fct1[1][0], fct1[1][1], fct1[1][2],
			fct1[2][0], fct1[2][1], fct1[2][2],
			fct2[0][0], fct2[0][1], fct2[0][2],
			fct2[1][0], fct2[1][1], fct2[1][2],
			fct2[2][0], fct2[2][1], fct2[2][2],
			m11, m12, m13, d);*/
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

	bool FastEnvelope::is_3_triangle_cut(const std::array<Vector3, 3>& triangle, const std::array<Vector3, 3>& f1, const std::array<Vector3, 3>& f2) {
		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];
		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
			std::cout << "Degeneration happens" << std::endl;
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
		int o = o1 + o2 + o3;
		if (o == 3 || o == -3) {
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
	int FastEnvelope::tri_cut_tri_Wang(const Vector3& p1, const Vector3& p2, const Vector3& p3,//even if only edge connected, regarded as intersected
		const Vector3& q1, const Vector3& q2, const Vector3& q3) {
		int o1, o2, o3, o4, o5, o6;
		o1 = Predicates::orient_3d(p1, q1, q2, q3);
		o2 = Predicates::orient_3d(p2, q1, q2, q3);
		o3 = Predicates::orient_3d(p3, q1, q2, q3);
		if (o1 == o2 && o2 == o3 && o3 == 0) {
			return CUT_COPLANAR;

		}
		o4 = Predicates::orient_3d(q1, p1, p2, p3);
		o5 = Predicates::orient_3d(q2, p1, p2, p3);
		o6 = Predicates::orient_3d(q3, p1, p2, p3);
		if (o1*o2 == 1 && o1*o3 == 1) {
			return CUT_EMPTY;
		}
		if (o3*o4 == 1 && o3*o5 == 1) {
			return CUT_EMPTY;
		}

		//Not_Finished
		assert(false);
		return 0;
	}
	int FastEnvelope::seg_cut_tri(const Vector3 & seg0, const Vector3 &seg1, const Vector3&t0, const Vector3&t1, const Vector3 &t2) {
		int o1, o2, o3, o4, o5;
		o1 = Predicates::orient_3d(seg0, t0, t1, t2);
		o2 = Predicates::orient_3d(seg1, t0, t1, t2);
		int op = o1 * o2;
		if (op >= 0) {
			return CUT_COPLANAR;//infact, coplanar and not on this plane
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
			if (area < SCALAR_ZERO) {
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
		float_precision ax = float_precision(triangle[0][0] - triangle[1][0], digit), ay = float_precision(triangle[0][1] - triangle[1][1], digit), az = float_precision(triangle[0][2] - triangle[1][2], digit),
			bx = float_precision(triangle[0][0] - triangle[2][0], digit), by = float_precision(triangle[0][1] - triangle[2][1], digit), bz = float_precision(triangle[0][2] - triangle[2][2], digit);
		float_precision x = ay * bz - az * by, y = az * bx - ax * bz, z = ax * by - ay * bx;
		float_precision length = sqrt(x * x + y * y + z * z);
		x = x / length; y = y / length; z = z / length;
		Scalar fx = x, fy = y, fz = z;
		return Vector3(fx, fy, fz);

	}




}













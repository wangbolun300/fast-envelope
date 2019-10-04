#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Predicates.hpp>
#include <fastenvelope/Logger.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_reorder.h>

#include <igl/Timer.h>

#include <fstream>
#include <istream>


static const int p_face[8][3] = { {0, 1, 3}, {7, 6, 9}, {1, 0, 7}, {2, 1, 7}, {3, 2, 8}, {3, 9, 10}, {5, 4, 11}, {0, 5, 6} }; //prism triangle index. all with orientation.
static const int c_face[6][3] = { {0, 1, 2}, {4, 7, 6}, {0, 3, 4}, {1, 0, 4}, {1, 5, 2}, {2, 6, 3} };

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
	{-1, -1} };
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


static const int p_facenumber = 8;
static const int c_facenumber = 6;

static const std::array<std::array<int, 2>, 3> triseg = {
	{{{0, 1}}, {{0, 2}}, {{1, 2}}}
};

namespace fastEnvelope
{
	namespace {
		const std::function<int(double)> check_double = [](double v) {
			if (fabs(v) < 1e-10)
				return -2;

			if (v > 0)
				return 1;

			if (v < 0)
				return -1;

			return 0;
		};

		const std::function<int(fastEnvelope::Rational)> check_Rational = [](fastEnvelope::Rational v) {
			return v.get_sign();
		};

		const std::function<int(fastEnvelope::Multiprecision)> check_Multiprecision = [](fastEnvelope::Multiprecision v) {
			return v.get_sign();
		};

		void to_geogram_mesh(const std::vector<Vector3> &V, const std::vector<Vector3i> &F, GEO::Mesh &M)
		{
			M.clear();

			// Setup vertices
			M.vertices.create_vertices(V.size());

			for (int i = 0; i < (int)M.vertices.nb(); ++i)
			{
				GEO::vec3 &p = M.vertices.point(i);

				p[0] = V[i][0];
				p[1] = V[i][1];
				p[2] = V[i][2];
			}

			// Setup faces

			M.facets.create_triangles(F.size());

			for (int c = 0; c < (int)M.facets.nb(); ++c)
			{
				for (int lv = 0; lv < 3; ++lv)
				{
					M.facets.set_vertex(c, lv, F[c][lv]);
				}
			}
		}

		void from_geogram_mesh(const GEO::Mesh &M, std::vector<Vector3> &V, std::vector<Vector3i> &F)
		{
			V.resize(M.vertices.nb());

			for (int i = 0; i < (int)M.vertices.nb(); ++i)
			{
				GEO::vec3 p = M.vertices.point(i);

				V[i][0] = p[0];
				V[i][1] = p[1];
				V[i][2] = p[2];
			}

			assert(M.facets.are_simplices());

			F.resize(M.facets.nb());

			for (int c = 0; c < (int)M.facets.nb(); ++c)
			{
				for (int lv = 0; lv < 3; ++lv)
				{
					F[c][lv] = M.facets.vertex(c, lv);
				}
			}

			assert(M.cells.are_simplices());
		}

		void triangle_sample_segment(const std::array<Vector3, 3> &triangle, Vector3 &ps, const int &pieces, const int &nbr)
		{
			int t = pieces - 1;
			if (triangle[1] - triangle[0] == Vector3(0, 0, 0))
			{

				ps = (triangle[0] + (triangle[2] - triangle[0]) * nbr / t);

				return;
			}
			if (triangle[2] - triangle[0] == Vector3(0, 0, 0))
			{

				ps = (triangle[0] + (triangle[1] - triangle[0]) * nbr / t);

				return;
			}
			if (triangle[2] - triangle[1] == Vector3(0, 0, 0))
			{

				ps = (triangle[0] + (triangle[1] - triangle[0]) * nbr / t);

				return;
			}

			Scalar d1 = (triangle[1] - triangle[0]).norm(), d2 = (triangle[2] - triangle[0]).norm(), d3 = (triangle[1] - triangle[2]).norm();
			if (d1 >= d2 && d1 >= d3)
			{
				ps = (triangle[0] + (triangle[1] - triangle[0]) * nbr / t);

				return;
			}
			if (d2 >= d1 && d2 >= d3)
			{
				ps = (triangle[0] + (triangle[2] - triangle[0]) * nbr / t);

				return;
			}
			if (d3 >= d1 && d3 >= d2)
			{
				ps = (triangle[1] + (triangle[2] - triangle[1]) * nbr / t);

				return;
			}
		}
	}



	FastEnvelope::FastEnvelope(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, const Scalar eps, const int spac)
	{
		igl::Timer timer;
		//Multiprecision::set_precision(256);
		timer.start();
		Vector3 min, max;
		min = m_ver.front();
		max = m_ver.front();

		for (size_t j = 0; j < m_ver.size(); j++)
		{
			for (int i = 0; i < 3; i++)
			{
				min(i) = std::min(min(i), m_ver[j](i));
				max(i) = std::max(max(i), m_ver[j](i));
			}
		}
		timer.stop();
		logger().info("Get bb corner time {}s", timer.getElapsedTimeInSec());

		const Scalar bbd = (max - min).norm();
		const Scalar epsilon = bbd * eps; //eps*bounding box diagnal


		timer.start();
		GEO::Mesh M;

		std::vector<Vector3> ver_new;
		std::vector<Vector3i> faces_new;

		to_geogram_mesh(m_ver, m_faces, M);
		GEO::mesh_reorder(M, GEO::MESH_ORDER_MORTON);
		from_geogram_mesh(M, ver_new, faces_new);
		timer.stop();
		logger().info("Resorting mesh time {}s", timer.getElapsedTimeInSec());

		timer.start();
		halfspace_generation(ver_new, faces_new, envprism, envcubic, epsilon);//TODO take lambda function to generate any shape envelope convex polyhedra
		timer.stop();
		logger().info("Box generation time {}s", timer.getElapsedTimeInSec());

		//build a hash function
		std::vector<std::array<Vector3, 2>> cornerlist;
		timer.start();
		CornerList_prism(envprism, cornerlist);
		std::vector<std::array<Vector3, 2>> cubiconors;
		CornerList_cubic(envcubic, cubiconors);

		//////////////////////////////////////////
		envelope.resize(envprism.size() + envcubic.size());

		for (int i = 0; i < envprism.size(); i++) {
			envelope[i].resize(12);
			for (int j = 0; j < 12; j++) {
				envelope[i][j] = envprism[i][j];
			}
		}
		for (int i = 0; i < envcubic.size(); i++) {
			envelope[i + envprism.size()].resize(8);
			for (int j = 0; j < 8; j++) {
				envelope[i + envprism.size()][j] = envcubic[i][j];
			}
		}
		//////////////////////////////////////////////////////////
		cornerlist.insert(cornerlist.end(), cubiconors.begin(), cubiconors.end());
		timer.stop();
		logger().info("Corner generation time {}s", timer.getElapsedTimeInSec());

		timer.start();
		tree.init_envelope(cornerlist);
		timer.stop();
		logger().info("Tree init time {}s", timer.getElapsedTimeInSec());

		//initializing types
		initFPU();

		//logger().debug("envelope size {}", envelope.size());
	}

	bool FastEnvelope::is_outside(const std::array<Vector3, 3> &triangle) const
	{
		igl::Timer timer;
		timer.start();
		std::vector<unsigned int> querylist;
		tree.facet_in_envelope(triangle[0], triangle[1], triangle[2], querylist);

		const auto res = FastEnvelopeTestImplicit(triangle, querylist);
		timer.stop();
		//logger().info("Query time time {}s", timer.getElapsedTimeInSec());

		return res;
	}

	void FastEnvelope::print_prisms(const std::array<Vector3, 3> &triangle, const std::string &path) const
	{
		const int prism_size = envprism.size();
		std::vector<unsigned int> querylist;
		tree.facet_in_envelope(triangle[0], triangle[1], triangle[2], querylist);
		std::ofstream fout;
		fout.open(path + "visualtriangle.txt");

		for (int i = 0; i < 3; i++)
		{
			fout << std::setprecision(17) << triangle[i][0] << " " << triangle[i][1] << " " << triangle[i][2] << std::endl;
		}
		fout.close();
		fout.open(path + "prisms.txt");
		for (int i = 0; i < querylist.size(); i++)
		{
			if (querylist[i] >= prism_size)
				continue;
			for (int j = 0; j < 12; j++)
			{
				fout << std::setprecision(17) << envprism[querylist[i]][j][0] << " " << envprism[querylist[i]][j][1] << " " << envprism[querylist[i]][j][2] << std::endl;
			}
		}
		fout.close();
		fout.open(path + "cubes.txt");
		for (int i = 0; i < querylist.size(); i++)
		{
			if (querylist[i] < prism_size)
				continue;
			for (int j = 0; j < 8; j++)
			{
				fout << std::setprecision(17) << envcubic[querylist[i]][j][0] << " " << envcubic[querylist[i]][j][1] << " " << envcubic[querylist[i]][j][2] << std::endl;
			}
		}
		fout.close();
	}

	bool FastEnvelope::sample_triangle_outside(const std::array<Vector3, 3> &triangle, const int &pieces) const
	{

		const auto triangle_sample_normal = [](const std::array<Vector3, 3> &triangle, Vector3 &ps, const int &pieces, const int &nbr1, const int &nbr2)
		{
			int l1s = pieces - 1; //
			Vector3 p1 = triangle[0] + (triangle[1] - triangle[0]) * nbr1 / l1s, d = (triangle[2] - triangle[1]) / l1s;
			ps = p1 + d * nbr2;
		};

		const auto triangle_sample_normal_rational = [](const std::array<Vector3, 3> &triangle, Rational &ps0, Rational &ps1, Rational &ps2, const int &pieces, const int &nbr1, const int &nbr2)
		{
			int l1s = pieces - 1; //
			Rational t00(triangle[0][0]), t01(triangle[0][1]), t02(triangle[0][2]), t10(triangle[1][0]), t11(triangle[1][1]),
				t12(triangle[1][2]), t20(triangle[2][0]), t21(triangle[2][1]), t22(triangle[2][2]), nbr1r(nbr1), nbr2r(nbr2), l1sr(l1s);

			Rational p0 = t00 + (t10 - t00) * nbr1r / l1sr, d0 = (t20 - t10) / l1sr;
			Rational p1 = t01 + (t11 - t01) * nbr1r / l1sr, d1 = (t21 - t11) / l1sr;
			Rational p2 = t02 + (t12 - t02) * nbr1r / l1sr, d2 = (t22 - t12) / l1sr;
			ps0 = p0 + d0 * nbr2;
			ps1 = p1 + d1 * nbr2;
			ps2 = p2 + d2 * nbr2;
		};

		bool out;
		Vector3 point;
		std::vector<unsigned int> querylist;
		tree.facet_in_envelope(triangle[0], triangle[1], triangle[2], querylist);
		int jump = -1;
		if (querylist.size() == 0)
			return 1;

		int deg = is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

		if (deg == DEGENERATED_POINT)
		{
			out = point_out_prism(triangle[0], querylist, jump);
			if (out == true)
			{

				return 1;
			}
			return 0;
		}
		if (deg == DEGENERATED_SEGMENT)
		{
			for (int i = 0; i < pieces; i++)
			{
				triangle_sample_segment(triangle, point, pieces, i);
				out = point_out_prism(point, querylist, jump);
				if (out == true)
				{

					return 1;
				}
			}
			return 0;
		}

		for (int i = 0; i < pieces; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				triangle_sample_normal(triangle, point, pieces, i, j);
				out = point_out_prism(point, querylist, jump);
				if (out == true)
				{

					return 1;
				}
			}
		}

		return 0;
	}


	bool FastEnvelope::FastEnvelopeTestImplicit(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const
	{
		if (prismindex.size() == 0)
		{
			return true;
		}

		const int prism_size = envprism.size();



		int jump1, jump2;

		std::vector<std::array<int, 2>> inter_ijk_list; //list of intersected triangle

		bool out, cut;

		int inter, inter1, record1, record2,

			tti; //triangle-triangle intersection

		jump1 = -1;
		for (int i = 0; i < 3; i++)
		{

			out = point_out_prism(triangle[i], prismindex, jump1);

			if (out == true)
			{

				return 1;
			}
		}

		if (prismindex.size() == 1)
			return 0;

		////////////////////degeneration fix
		igl::Timer timer;
		timer.start();
		int degeneration = is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

		if (degeneration == DEGENERATED_POINT)
		{ //case 1 degenerate to a point

			return 0;

		} //case 1 degenerate to a point
		DATA_LPI datalpi;
		std::vector<DATA_LPI> lpi_list;
		if (degeneration == DEGENERATED_SEGMENT)
		{

			for (int we = 0; we < 3; we++)
			{ //case 2 degenerated as a segment, at most test 2 segments,but still we need to test 3, because

				// of the endpoint-triangle intersection will be ignored

				// the segment is {triangle[triseg[we][0]], triangle[triseg[we][1]]}

				for (int i = 0; i < prismindex.size(); i++)
				{
					jump1 = prismindex[i];
					int face[3];
					std::vector<int> cid;
					/*if (envelope[prismindex[i]].size() == 12) {
						cut = is_seg_cut_prism(prismindex[i], triangle[triseg[we][0]], triangle[triseg[we][1]], cid);
					}
					if (envelope[prismindex[i]].size() == 8) {
						cut = is_seg_cut_cube(prismindex[i] - prism_size, triangle[triseg[we][0]], triangle[triseg[we][1]], cid);
					}*/
					cut = is_seg_cut_polyhedra(prismindex[i], triangle[triseg[we][0]], triangle[triseg[we][1]], cid);
					if (cut == false) continue;

					for (int j = 0; j < cid.size(); j++)
					{
						if (envelope[prismindex[i]].size() == 12) {
							face[0] = p_face[cid[j]][0];
							face[1] = p_face[cid[j]][1];
							face[2] = p_face[cid[j]][2];
						}
						if (envelope[prismindex[i]].size() == 8) {
							face[0] = c_face[cid[j]][0];
							face[1] = c_face[cid[j]][1];
							face[2] = c_face[cid[j]][2];
						}
						Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
						bool precom = orient3D_LPI_prefilter( //
							triangle[triseg[we][0]][0], triangle[triseg[we][0]][1], triangle[triseg[we][0]][2],
							triangle[triseg[we][1]][0], triangle[triseg[we][1]][1], triangle[triseg[we][1]][2],
							envelope[prismindex[i]][face[0]][0], envelope[prismindex[i]][face[0]][1], envelope[prismindex[i]][face[0]][2],
							envelope[prismindex[i]][face[1]][0], envelope[prismindex[i]][face[1]][1], envelope[prismindex[i]][face[1]][2],
							envelope[prismindex[i]][face[2]][0], envelope[prismindex[i]][face[2]][1], envelope[prismindex[i]][face[2]][2],
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
						if (precom == true)
						{
							inter = Implicit_Seg_Facet_interpoint_Out_Prism_double(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
								triangle[triseg[we][0]], triangle[triseg[we][1]],
								envelope[prismindex[i]][face[0]], envelope[prismindex[i]][face[1]], envelope[prismindex[i]][face[2]],
								prismindex, jump1);
							if (inter == 1)
							{

								return 1;
							}
						}
						else
						{
							datalpi.segid = we;
							datalpi.prismid = prismindex[i];
							datalpi.facetid = cid[j];
							datalpi.jump1 = jump1;
							lpi_list.emplace_back(datalpi);
						}
					}

				}

			} //case 2 case 2 degenerated as a segment

			for (int i = 0; i < lpi_list.size(); i++)
			{
				inter = Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(lpi_list[i], triangle, prismindex);
				if (inter == 1)
					return 1;
			}
			return 0;
		}
		//
		////////////////////////////////degeneration fix over
		timer.stop();
		// logger().info("Bit time? time {}s", timer.getElapsedTimeInSec())

		// timer_bigpart.start();
		for (int i = 0; i < prismindex.size(); i++)
		{
			std::vector<int> cidl;
			jump1 = prismindex[i];
			if (envelope[prismindex[i]].size() == 12) {
				cut = is_triangle_cut_prism(prismindex[i],
					triangle[0], triangle[1], triangle[2], cidl);
			}
			if (envelope[prismindex[i]].size() == 8) {
				cut = is_triangle_cut_cube(prismindex[i]-prism_size,
					triangle[0], triangle[1], triangle[2], cidl);
			}
			std::vector<int> cidd;
			bool cut1 = is_triangle_cut_envelope_polyhedra(prismindex[i],
				triangle[0], triangle[1], triangle[2], cidd);
			if (cut != cut1) std::cout << "diff in cut " << cut - cut1 << std::endl;
			if (cut == cut1) {
				if (cidl != cidd) {
					//std::cout << "size of list " << cidl.size() << " " << cidd.size()<<" "<< envelope[prismindex[i]].size() << std::endl;
				}
				//if (envelope[prismindex[i]].size() == 8) std::cout << "!!!" << std::endl;
			}
			if (cut == false) continue;
			int face[3];
			for (int j = 0; j < cidl.size(); j++) {
				if (envelope[prismindex[i]].size() == 12) {
					face[0] = p_face[cidl[j]][0];
					face[1] = p_face[cidl[j]][1];
					face[2] = p_face[cidl[j]][2];
				}
				if (envelope[prismindex[i]].size() == 8) {
					face[0] = c_face[cidl[j]][0];
					face[1] = c_face[cidl[j]][1];
					face[2] = c_face[cidl[j]][2];
				}
				for (int k = 0; k < 3; k++) {
					tti = seg_cut_plane(triangle[triseg[k][0]], triangle[triseg[k][1]],
						envelope[prismindex[i]][face[0]],
						envelope[prismindex[i]][face[1]],
						envelope[prismindex[i]][face[2]]); 
					if (tti != CUT_FACE) continue;
					Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
					bool precom = orient3D_LPI_prefilter( //
						triangle[triseg[k][0]][0], triangle[triseg[k][0]][1], triangle[triseg[k][0]][2],
						triangle[triseg[k][1]][0], triangle[triseg[k][1]][1], triangle[triseg[k][1]][2],
						envelope[prismindex[i]][face[0]][0], envelope[prismindex[i]][face[0]][1], envelope[prismindex[i]][face[0]][2],
						envelope[prismindex[i]][face[1]][0], envelope[prismindex[i]][face[1]][1], envelope[prismindex[i]][face[1]][2],
						envelope[prismindex[i]][face[2]][0], envelope[prismindex[i]][face[2]][1], envelope[prismindex[i]][face[2]][2],
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
					if (precom == true)
					{
						inter = Implicit_Seg_Facet_interpoint_Out_Prism_double(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							triangle[triseg[k][0]], triangle[triseg[k][1]],
							envelope[prismindex[i]][face[0]], envelope[prismindex[i]][face[1]], envelope[prismindex[i]][face[2]],
							prismindex, jump1);
						if (inter == 1)
						{

							return 1;
						}
					}
					else
					{
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
		for (int i = 0; i < lpi_list.size(); i++)
		{
			inter = Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(lpi_list[i], triangle, prismindex);
			if (inter == 1)
				return 1;
		} //TODO consider this part put here or the end of the algorithm

		int listsize = inter_ijk_list.size();
		DATA_TPI datatpi;
		std::vector<DATA_TPI> tpilist;
		tpilist.reserve(listsize / 5);

		int id, id0 = 0;
		for (int i = 0; i < listsize; i++)
		{
			jump1 = inter_ijk_list[i][0];
			for (int j = i + 1; j < listsize; j++)
			{

				//check triangle{ { envprism[list[i][0]][p_triangle[list[i][1]][list[i][2]][0]], ...[1],...[2] } } and triangle{ { envprism[list[j][0]][p_triangle[list[j][1]][list[j][2]][0]], ...[1],...[2] } }
				//and T
				if (inter_ijk_list[i][0] == inter_ijk_list[j][0]) continue;


				//find prism_map[list[i][1]*8+list[j][1]][0],prism_map[list[i][1]*8+list[j][1]][1]
				// timer_u.start();
				Vector3 t00, t01, t02, t10, t11, t12;
				
				if (envelope[inter_ijk_list[i][0]].size() == 12)
				{
					t00 = envelope[inter_ijk_list[i][0]][p_face[inter_ijk_list[i][1]][0]];
					t01 = envelope[inter_ijk_list[i][0]][p_face[inter_ijk_list[i][1]][1]];
					t02 = envelope[inter_ijk_list[i][0]][p_face[inter_ijk_list[i][1]][2]];
				}
				if (envelope[inter_ijk_list[i][0]].size() == 8)
				{
					t00 = envelope[inter_ijk_list[i][0]][c_face[inter_ijk_list[i][1]][0]];
					t01 = envelope[inter_ijk_list[i][0]][c_face[inter_ijk_list[i][1]][1]];
					t02 = envelope[inter_ijk_list[i][0]][c_face[inter_ijk_list[i][1]][2]];
				}
				if (envelope[inter_ijk_list[j][0]].size() == 12)
				{
					t10 = envelope[inter_ijk_list[j][0]][p_face[inter_ijk_list[j][1]][0]];
					t11 = envelope[inter_ijk_list[j][0]][p_face[inter_ijk_list[j][1]][1]];
					t12 = envelope[inter_ijk_list[j][0]][p_face[inter_ijk_list[j][1]][2]];
				}
				if (envelope[inter_ijk_list[j][0]].size() == 8)
				{
					t10 = envelope[inter_ijk_list[j][0]][c_face[inter_ijk_list[j][1]][0]];
					t11 = envelope[inter_ijk_list[j][0]][c_face[inter_ijk_list[j][1]][1]];
					t12 = envelope[inter_ijk_list[j][0]][c_face[inter_ijk_list[j][1]][2]];
				}
				

				jump2 = inter_ijk_list[j][0];

				
				
				Scalar d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7;
				bool multiflag;

				bool pre = orient3D_TPI_prefilter(triangle[0][0], triangle[0][1], triangle[0][2],
					triangle[1][0], triangle[1][1], triangle[1][2], triangle[2][0], triangle[2][1], triangle[2][2],
					t00[0], t00[1], t00[2],
					t01[0], t01[1], t01[2],
					t02[0], t02[1], t02[2],
					t10[0], t10[1], t10[2],
					t11[0], t11[1], t11[2],
					t12[0], t12[1], t12[2],
					d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7);
				// timein2 += timer_u.getElapsedTimeInSec();
				if (pre == true)
				{

					TPI_exact_suppvars s;
					cut = is_3_triangle_cut_double(d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7,
						triangle, t00, t01, t02, t10, t11, t12, multiflag, s);
					if (cut == false)
						continue;
					// timer_u.start();
					inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double( //TODO takes most of time
						d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7,
						triangle, t00, t01, t02, t10, t11, t12,
						prismindex, jump1, jump2, multiflag, s);
					// timein3 += timer_u.getElapsedTimeInSec();
					if (inter == 1)
					{

						return 1;
					}
				}
				else
				{
					// timer_u.start();
					datatpi.prismid1 = inter_ijk_list[i][0];
					datatpi.facetid1 = inter_ijk_list[i][1];
					datatpi.prismid2 = inter_ijk_list[j][0];
					datatpi.facetid2 = inter_ijk_list[j][1];
					datatpi.jump1 = jump1;
					datatpi.jump2 = jump2;
					tpilist.emplace_back(datatpi);
					// timein4 += timer_u.getElapsedTimeInSec();
				}
			}
		}
		// time_p3d += timerdetail.getElapsedTimeInSec();
		// timerdetail.start();
		for (int i = 0; i < tpilist.size(); i++)
		{
			inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(tpilist[i], triangle, prismindex); //is_3_intersection is already in it
			if (inter == 1)
				return 1;
		}
		// time_p3m += timerdetail.getElapsedTimeInSec();
		// time_p3 += timer_bigpart.getElapsedTimeInSec();

		return 0;
	}


	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const
	{
		const int prism_size = envprism.size();
		int tot, ori;
		static Scalar
			s00,
			s01, s02, s10, s11, s12,
			t00, t01, t02,
			t10, t11, t12,
			t20, t21, t22,

			e00, e01, e02,
			e10, e11, e12,
			e20, e21, e22;
		int face[3];
		if (envelope[datalpi.prismid].size() ==12)
		{
			face[0] = p_face[datalpi.facetid][0];
			face[1] = p_face[datalpi.facetid][1];
			face[2] = p_face[datalpi.facetid][2];
	
		}
		if (envelope[datalpi.prismid].size() == 8)
		{
			face[0] = c_face[datalpi.facetid][0];
			face[1] = c_face[datalpi.facetid][1];
			face[2] = c_face[datalpi.facetid][2];
			
		}

		t00 = envelope[datalpi.prismid][face[0]][0];
		t01 = envelope[datalpi.prismid][face[0]][1];
		t02 = envelope[datalpi.prismid][face[0]][2];

		t10 = envelope[datalpi.prismid][face[1]][0];
		t11 = envelope[datalpi.prismid][face[1]][1];
		t12 = envelope[datalpi.prismid][face[1]][2];

		t20 = envelope[datalpi.prismid][face[2]][0];
		t21 = envelope[datalpi.prismid][face[2]][1];
		t22 = envelope[datalpi.prismid][face[2]][2];
	

		s00 = triangle[triseg[datalpi.segid][0]][0];
		s01 = triangle[triseg[datalpi.segid][0]][1];
		s02 = triangle[triseg[datalpi.segid][0]][2];

		s10 = triangle[triseg[datalpi.segid][1]][0];
		s11 = triangle[triseg[datalpi.segid][1]][1];
		s12 = triangle[triseg[datalpi.segid][1]][2];
		// timer.start();
		LPI_exact_suppvars s;
		bool premulti = orient3D_LPI_pre_exact(s00, s01, s02, s10, s11, s12,
			t00, t01, t02, t10, t11, t12, t20, t21, t22, s);

		// time_multi += timer.getElapsedTimeInSec();

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == datalpi.jump1)
			{
				continue;
			}
			tot = 0;
			int number;
			if (envelope[prismindex[i]].size() == 12)
			{
				number = p_facenumber;
			}
			if (envelope[prismindex[i]].size() == 8)
			{
				number = c_facenumber;

			}
			for (int j = 0; j < number; j++) {
				if (envelope[prismindex[i]].size() == 12)
				{
					face[0] = p_face[j][0];
					face[1] = p_face[j][1];
					face[2] = p_face[j][2];

				}
				if (envelope[prismindex[i]].size() == 8)
				{
					face[0] = c_face[j][0];
					face[1] = c_face[j][1];
					face[2] = c_face[j][2];

				}
				e00 = envelope[prismindex[i]][face[0]][0];
				e01 = envelope[prismindex[i]][face[0]][1];
				e02 = envelope[prismindex[i]][face[0]][2];
				e10 = envelope[prismindex[i]][face[1]][0];
				e11 = envelope[prismindex[i]][face[1]][1];
				e12 = envelope[prismindex[i]][face[1]][2];
				e20 = envelope[prismindex[i]][face[2]][0];
				e21 = envelope[prismindex[i]][face[2]][1];
				e22 = envelope[prismindex[i]][face[2]][2];
				// timer.start();
				ori = orient3D_LPI_post_exact(s, s00, s01, s02,
					e00, e01, e02, e10, e11, e12,
					e20, e21, e22);
				if (ori == 1 || ori == 0)
				{
					break;
				}

				if (ori == -1)
				{
					tot++;
				}
			}
			if (tot == number)
			{
				return IN_PRISM;
			}
		}
			
			
		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_double(
		const Scalar &a11, const Scalar &a12, const Scalar &a13, const Scalar &d, const Scalar &fa11,
		const Scalar &fa12, const Scalar &fa13, const Scalar &max1, const Scalar &max2, const Scalar &max5,
		const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
		const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump) const
	{
		const int prism_size = envprism.size();
		int tot;
		int ori, ori1;
		INDEX index;
		std::vector<INDEX> recompute;
		int face[3];
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump) continue;

			index.FACES.clear();
			tot = 0;
			
			int number;
			if (envelope[prismindex[i]].size() == 12)
			{
				number = p_facenumber;
			}
			if (envelope[prismindex[i]].size() == 8)
			{
				number = c_facenumber;

			}
			for (int j = 0; j < number; j++) {
				if (envelope[prismindex[i]].size() == 12) {
					face[0] = p_face[j][0];
					face[1] = p_face[j][1];
					face[2] = p_face[j][2];
				}
				if (envelope[prismindex[i]].size() == 8) {
					face[0] = c_face[j][0];
					face[1] = c_face[j][1];
					face[2] = c_face[j][2];
				}
				ori =
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						segpoint0[0], segpoint0[1], segpoint0[2],
						envelope[prismindex[i]][face[0]][0], envelope[prismindex[i]][face[0]][1], envelope[prismindex[i]][face[0]][2],
						envelope[prismindex[i]][face[1]][0], envelope[prismindex[i]][face[1]][1], envelope[prismindex[i]][face[1]][2],
						envelope[prismindex[i]][face[2]][0], envelope[prismindex[i]][face[2]][1], envelope[prismindex[i]][face[2]][2]);

				if (ori == 1)
				{
					break;
				}
				if (ori == 0)
				{
					index.FACES.emplace_back(j);
				}

				else if (ori == -1)
				{
					tot++;
				}

			}
			if (tot == number)
			{

				return IN_PRISM;
			}

			if (ori != 1)
			{
				assert(!index.FACES.empty());
				index.Pi = prismindex[i];
				recompute.emplace_back(index);
			}
		}

		if (!recompute.empty())
		{
			static Scalar
				s00, s01, s02, s10, s11, s12,
				t00, t01, t02,
				t10, t11, t12,
				t20, t21, t22,

				e00, e01, e02,
				e10, e11, e12,
				e20, e21, e22;
			s00 = segpoint0[0];
			s01 = segpoint0[1];
			s02 = segpoint0[2];
			s10 = segpoint1[0];
			s11 = segpoint1[1];
			s12 = segpoint1[2];
			t00 = triangle0[0];
			t01 = triangle0[1];
			t02 = triangle0[2];
			t10 = triangle1[0];
			t11 = triangle1[1];
			t12 = triangle1[2];
			t20 = triangle2[0];
			t21 = triangle2[1];
			t22 = triangle2[2];
			// timer.start();
			LPI_exact_suppvars s;
			bool premulti = orient3D_LPI_pre_exact(s00, s01, s02, s10, s11, s12,
				t00, t01, t02, t10, t11, t12, t20, t21, t22, s);
			// time_multi += timer.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = recompute[k].Pi;

				
				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];
					if (envelope[in1].size() == 12) {
						face[0] = p_face[in2][0];
						face[1] = p_face[in2][1];
						face[2] = p_face[in2][2];
					}
					if (envelope[in1].size() == 8) {
						face[0] = c_face[in2][0];
						face[1] = c_face[in2][1];
						face[2] = c_face[in2][2];
					}
					e00 = envelope[in1][face[0]][0];
					e01 = envelope[in1][face[0]][1];
					e02 = envelope[in1][face[0]][2];
					e10 = envelope[in1][face[1]][0];
					e11 = envelope[in1][face[1]][1];
					e12 = envelope[in1][face[1]][2];
					e20 = envelope[in1][face[2]][0];
					e21 = envelope[in1][face[2]][1];
					e22 = envelope[in1][face[2]][2];
					// timer.start();
					ori = orient3D_LPI_post_exact(s, s00, s01, s02,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22);
					// time_multi += timer.getElapsedTimeInSec();
					if (ori == 1 || ori == 0)
						break;
					
				}
				if (ori == -1) return IN_PRISM;

			}
		}

		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_TPI &datatpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const
	{
		const int prism_size = envprism.size();
		int tot, ori;
		int jump1 = datatpi.jump1, jump2 = datatpi.jump2;

		static Scalar
			t00,
			t01, t02,
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
			e20, e21, e22;
		t00 = (triangle[0][0]);
		t01 = (triangle[0][1]);
		t02 = (triangle[0][2]);
		t10 = (triangle[1][0]);
		t11 = (triangle[1][1]);
		t12 = (triangle[1][2]);
		t20 = (triangle[2][0]);
		t21 = (triangle[2][1]);
		t22 = (triangle[2][2]);
		int face1[3],face2[3];
		if (envelope[datatpi.prismid1].size() == 12) {
			face1[0] = p_face[datatpi.facetid1][0];
			face1[1] = p_face[datatpi.facetid1][1];
			face1[2] = p_face[datatpi.facetid1][2];
		}
		if (envelope[datatpi.prismid1].size() == 8) {
			face1[0] = c_face[datatpi.facetid1][0];
			face1[1] = c_face[datatpi.facetid1][1];
			face1[2] = c_face[datatpi.facetid1][2];
		}

		if (envelope[datatpi.prismid2].size() == 12) {
			face2[0] = p_face[datatpi.facetid2][0];
			face2[1] = p_face[datatpi.facetid2][1];
			face2[2] = p_face[datatpi.facetid2][2];
		}
		if (envelope[datatpi.prismid2].size() == 8) {
			face2[0] = c_face[datatpi.facetid2][0];
			face2[1] = c_face[datatpi.facetid2][1];
			face2[2] = c_face[datatpi.facetid2][2];
		}
		f100 = envelope[datatpi.prismid1][face1[0]][0];
		f101 = envelope[datatpi.prismid1][face1[0]][1];
		f102 = envelope[datatpi.prismid1][face1[0]][2];
		f110 = envelope[datatpi.prismid1][face1[1]][0];
		f111 = envelope[datatpi.prismid1][face1[1]][1];
		f112 = envelope[datatpi.prismid1][face1[1]][2];
		f120 = envelope[datatpi.prismid1][face1[2]][0];
		f121 = envelope[datatpi.prismid1][face1[2]][1];
		f122 = envelope[datatpi.prismid1][face1[2]][2];
	
		f200 = envelope[datatpi.prismid2][face2[0]][0];
		f201 = envelope[datatpi.prismid2][face2[0]][1];
		f202 = envelope[datatpi.prismid2][face2[0]][2];
		f210 = envelope[datatpi.prismid2][face2[1]][0];
		f211 = envelope[datatpi.prismid2][face2[1]][1];
		f212 = envelope[datatpi.prismid2][face2[1]][2];
		f220 = envelope[datatpi.prismid2][face2[2]][0];
		f221 = envelope[datatpi.prismid2][face2[2]][1];
		f222 = envelope[datatpi.prismid2][face2[2]][2];

		TPI_exact_suppvars s;
		bool premulti = orient3D_TPI_pre_exact(t00, t01, t02, t10, t11, t12, t20, t21, t22,
			f100, f101, f102, f110, f111, f112, f120, f121, f122,
			f200, f201, f202, f210, f211, f212, f220, f221, f222,
			s);
		// time_multi += timer.getElapsedTimeInSec();
		if (premulti == false) return 2; //means have parallel facets
		bool cut = is_3_triangle_cut_pure_multiprecision(triangle, s);
		if (cut == false) return 2;
		int number;
		int face[3];
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump1 || prismindex[i] == jump2)
				continue;
			if (envelope[prismindex[i]].size() == 12) {
				number = p_facenumber;
			}
			if (envelope[prismindex[i]].size() == 8) {
				number = c_facenumber;
			}
			tot = 0;
			for (int j = 0; j < number; j++) {
				if (envelope[prismindex[i]].size() == 12) {
					face[0] = p_face[j][0];
					face[1] = p_face[j][1];
					face[2] = p_face[j][2];
				}
				if (envelope[prismindex[i]].size() == 8) {
					face[0] = c_face[j][0];
					face[1] = c_face[j][1];
					face[2] = c_face[j][2];
				}
				e00 = (envelope[prismindex[i]][face[0]][0]);
				e01 = (envelope[prismindex[i]][face[0]][1]);
				e02 = (envelope[prismindex[i]][face[0]][2]);
				e10 = (envelope[prismindex[i]][face[1]][0]);
				e11 = (envelope[prismindex[i]][face[1]][1]);
				e12 = (envelope[prismindex[i]][face[1]][2]);
				e20 = (envelope[prismindex[i]][face[2]][0]);
				e21 = (envelope[prismindex[i]][face[2]][1]);
				e22 = (envelope[prismindex[i]][face[2]][2]);
				ori = orient3D_TPI_post_exact(s, e00, e01, e02, e10, e11, e12, e20, e21, e22);
				// time_multi += timer.getElapsedTimeInSec();
				if (ori == 1 || ori == 0)
				{
					break;
				}

				if (ori == -1)
				{
					tot++;
				}
			}
			if (tot == number)
			{

				return IN_PRISM;
			}
			
		}
		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double(
		const Scalar &d, const Scalar &n1, const Scalar &n2, const Scalar &n3,
		const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5, const Scalar &max6, const Scalar &max7,
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
		const std::vector<unsigned int> &prismindex, const int &jump1, const int &jump2, const bool &multiflag,
		TPI_exact_suppvars &s) const
	{
		const int prism_size = envprism.size();
		int ori;
		int tot;
		int face[3];
		int number;
		INDEX index;
		std::vector<INDEX> recompute;

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump1 || prismindex[i] == jump2)
				continue;
			if (envelope[prismindex[i]].size() == 12) {
				number = p_facenumber;
			}
			if (envelope[prismindex[i]].size() == 8) {
				number = c_facenumber;
			}
			index.FACES.clear();
			tot = 0;
			for (int j = 0; j < number; j++) {
				if (envelope[prismindex[i]].size() == 12) {
					face[0] = p_face[j][0];
					face[1] = p_face[j][1];
					face[2] = p_face[j][2];
				}
				if (envelope[prismindex[i]].size() == 8) {
					face[0] = c_face[j][0];
					face[1] = c_face[j][1];
					face[2] = c_face[j][2];
				}
				
				ori =
					orient3D_TPI_postfilter(
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
						envelope[prismindex[i]][face[0]][0], envelope[prismindex[i]][face[0]][1], envelope[prismindex[i]][face[0]][2],
						envelope[prismindex[i]][face[1]][0], envelope[prismindex[i]][face[1]][1], envelope[prismindex[i]][face[1]][2],
						envelope[prismindex[i]][face[2]][0], envelope[prismindex[i]][face[2]][1], envelope[prismindex[i]][face[2]][2]);
				// timetpp1 += timer_a.getElapsedTimeInSec();

				if (ori == 1)
				{
					break;
				}
				if (ori == 0)
				{
					index.FACES.emplace_back(j);
				}

				else if (ori == -1)
				{
					tot++;
				}
			}
			if (tot == number)
			{

				return IN_PRISM;
			}

			if (ori != 1)
			{
				index.Pi = prismindex[i];
				recompute.emplace_back(index);

			}
		}

		if (recompute.size() > 0)
		{
			// timer_a.start();
			static Scalar
				t00,
				t01, t02,
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
				e20, e21, e22;
		
			t00 = (triangle[0][0]);
			t01 = (triangle[0][1]);
			t02 = (triangle[0][2]);
			t10 = (triangle[1][0]);
			t11 = (triangle[1][1]);
			t12 = (triangle[1][2]);
			t20 = (triangle[2][0]);
			t21 = (triangle[2][1]);
			t22 = (triangle[2][2]);

			f100 = (facet10[0]);
			f101 = (facet10[1]);
			f102 = (facet10[2]);
			f110 = (facet11[0]);
			f111 = (facet11[1]);
			f112 = (facet11[2]);
			f120 = (facet12[0]);
			f121 = (facet12[1]);
			f122 = (facet12[2]);

			f200 = (facet20[0]);
			f201 = (facet20[1]);
			f202 = (facet20[2]);
			f210 = (facet21[0]);
			f211 = (facet21[1]);
			f212 = (facet21[2]);
			f220 = (facet22[0]);
			f221 = (facet22[1]);
			f222 = (facet22[2]);

			if (multiflag == false)
			{
				// timer.start();
				bool premulti = orient3D_TPI_pre_exact(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					s);
				// time_multi += timer.getElapsedTimeInSec();
			}

			// timetpp2 += timer_a.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = recompute[k].Pi;

				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];
					if (envelope[in1].size() == 12) {
						face[0] = p_face[in2][0];
						face[1] = p_face[in2][1];
						face[2] = p_face[in2][2];
					}
					if (envelope[in1].size() == 8) {
						face[0] = c_face[in2][0];
						face[1] = c_face[in2][1];
						face[2] = c_face[in2][2];
					}
					e00 = (envelope[in1][face[0]][0]);
					e01 = (envelope[in1][face[0]][1]);
					e02 = (envelope[in1][face[0]][2]);
					e10 = (envelope[in1][face[1]][0]);
					e11 = (envelope[in1][face[1]][1]);
					e12 = (envelope[in1][face[1]][2]);
					e20 = (envelope[in1][face[2]][0]);
					e21 = (envelope[in1][face[2]][1]);
					e22 = (envelope[in1][face[2]][2]);
					ori = orient3D_TPI_post_exact(s,
						e00, e01, e02, e10, e11, e12,
						e20, e21, e22);

					if (ori == 1 || ori == 0)	break;
				}
				if (ori == -1) return IN_PRISM;
			}
		}

		return OUT_PRISM;
	}


	bool FastEnvelope::is_3_triangle_cut_pure_multiprecision(const std::array<Vector3, 3> &triangle, TPI_exact_suppvars &s)
	{
		int o1, o2, o3;
		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];

		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0)
		{
			logger().debug("Degeneration happens");
			//move this guy in constructor and use fixed seed
			n = { {Vector3(rand(), rand(), rand())} };
		}
		static Scalar
			t00,
			t01, t02,
			t10, t11, t12,
			t20, t21, t22;

		t00 = triangle[0][0];
		t01 = triangle[0][1];
		t02 = triangle[0][2];
		t10 = triangle[1][0];
		t11 = triangle[1][1];
		t12 = triangle[1][2];
		t20 = triangle[2][0];
		t21 = triangle[2][1];
		t22 = triangle[2][2];
		// timer.start();
		o1 = orient3D_TPI_post_exact(
			s,
			n[0], n[1], n[2],
			t00, t01, t02,
			t10, t11, t12);
		// time_multi += timer.getElapsedTimeInSec();
		if (o1 == 0)
			return false;
		// timer.start();
		o2 = orient3D_TPI_post_exact(
			s,
			n[0], n[1], n[2],
			t10, t11, t12,
			t20, t21, t22);
		// time_multi += timer.getElapsedTimeInSec();
		if (o2 == 0 || o1 + o2 == 0)
			return false;
		// timer.start();
		o3 = orient3D_TPI_post_exact(
			s,
			n[0], n[1], n[2],
			t20, t21, t22,
			t00, t01, t02);
		// time_multi += timer.getElapsedTimeInSec();
		if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0)
			return false;
		return true;
	}

	bool FastEnvelope::is_3_triangle_cut_double(
		const Scalar &d, const Scalar &n1, const Scalar &n2, const Scalar &n3,
		const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5,
		const Scalar &max6, const Scalar &max7,
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
		bool &multiflag,
		TPI_exact_suppvars &s)
	{
		multiflag = false;
		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];

		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0)
		{
			logger().debug("Degeneration happens");
			//move this guy in constructor and use fixed seed
			n = { {Vector3(rand(), rand(), rand())} };
		}

		static Scalar
			t00,
			t01, t02,
			t10, t11, t12,
			t20, t21, t22,

			f100, f101, f102,
			f110, f111, f112,
			f120, f121, f122,

			f200, f201, f202,
			f210, f211, f212,
			f220, f221, f222;

		bool premulti = false;
		int o1 = orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2]);
		if (o1 == 0)
		{

			t00 = (triangle[0][0]);
			t01 = (triangle[0][1]);
			t02 = (triangle[0][2]);
			t10 = (triangle[1][0]);
			t11 = (triangle[1][1]);
			t12 = (triangle[1][2]);
			t20 = (triangle[2][0]);
			t21 = (triangle[2][1]);
			t22 = (triangle[2][2]);

			f100 = (facet10[0]);
			f101 = (facet10[1]);
			f102 = (facet10[2]);
			f110 = (facet11[0]);
			f111 = (facet11[1]);
			f112 = (facet11[2]);
			f120 = (facet12[0]);
			f121 = (facet12[1]);
			f122 = (facet12[2]);

			f200 = (facet20[0]);
			f201 = (facet20[1]);
			f202 = (facet20[2]);
			f210 = (facet21[0]);
			f211 = (facet21[1]);
			f212 = (facet21[2]);
			f220 = (facet22[0]);
			f221 = (facet22[1]);
			f222 = (facet22[2]);

			// timer.start();
			premulti = orient3D_TPI_pre_exact(
				t00, t01, t02, t10, t11, t12, t20, t21, t22,
				f100, f101, f102, f110, f111, f112, f120, f121, f122,
				f200, f201, f202, f210, f211, f212, f220, f221, f222,
				s);
			// time_multi += timer.getElapsedTimeInSec();
			// timer.start();
			o1 = orient3D_TPI_post_exact(
				s,
				n[0], n[1], n[2],
				t00, t01, t02,
				t10, t11, t12);
			// time_multi += timer.getElapsedTimeInSec();
		}

		if (o1 == 0)
			return false;

		int o2 = orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			triangle[1][0], triangle[1][1], triangle[1][2],
			triangle[2][0], triangle[2][1], triangle[2][2]);
		if (o2 == 0)
		{
			if (premulti == false)
			{
				t00 = (triangle[0][0]);
				t01 = (triangle[0][1]);
				t02 = (triangle[0][2]);
				t10 = (triangle[1][0]);
				t11 = (triangle[1][1]);
				t12 = (triangle[1][2]);
				t20 = (triangle[2][0]);
				t21 = (triangle[2][1]);
				t22 = (triangle[2][2]);

				f100 = (facet10[0]);
				f101 = (facet10[1]);
				f102 = (facet10[2]);
				f110 = (facet11[0]);
				f111 = (facet11[1]);
				f112 = (facet11[2]);
				f120 = (facet12[0]);
				f121 = (facet12[1]);
				f122 = (facet12[2]);

				f200 = (facet20[0]);
				f201 = (facet20[1]);
				f202 = (facet20[2]);
				f210 = (facet21[0]);
				f211 = (facet21[1]);
				f212 = (facet21[2]);
				f220 = (facet22[0]);
				f221 = (facet22[1]);
				f222 = (facet22[2]);

				// timer.start();
				premulti = orient3D_TPI_pre_exact(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					s);
				// time_multi += timer.getElapsedTimeInSec();
			}
			// timer.start();
			o2 = orient3D_TPI_post_exact(
				s,
				n[0], n[1], n[2],
				t10, t11, t12,
				t20, t21, t22);
			// time_multi += timer.getElapsedTimeInSec();
			/*if (o2 == 1) after21++;
				if (o2 == -1) after22++;
				if (o2 == 0) after20++;*/
		}
		if (o2 == 0 || o1 + o2 == 0)
			return false;

		int o3 = orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			triangle[2][0], triangle[2][1], triangle[2][2],
			triangle[0][0], triangle[0][1], triangle[0][2]);
		if (o3 == 0)
		{
			if (premulti == false)
			{
				t00 = (triangle[0][0]);
				t01 = (triangle[0][1]);
				t02 = (triangle[0][2]);
				t10 = (triangle[1][0]);
				t11 = (triangle[1][1]);
				t12 = (triangle[1][2]);
				t20 = (triangle[2][0]);
				t21 = (triangle[2][1]);
				t22 = (triangle[2][2]);

				f100 = (facet10[0]);
				f101 = (facet10[1]);
				f102 = (facet10[2]);
				f110 = (facet11[0]);
				f111 = (facet11[1]);
				f112 = (facet11[2]);
				f120 = (facet12[0]);
				f121 = (facet12[1]);
				f122 = (facet12[2]);

				f200 = (facet20[0]);
				f201 = (facet20[1]);
				f202 = (facet20[2]);
				f210 = (facet21[0]);
				f211 = (facet21[1]);
				f212 = (facet21[2]);
				f220 = (facet22[0]);
				f221 = (facet22[1]);
				f222 = (facet22[2]);

				// timer.start();
				premulti = orient3D_TPI_pre_exact(
					t00, t01, t02, t10, t11, t12, t20, t21, t22,
					f100, f101, f102, f110, f111, f112, f120, f121, f122,
					f200, f201, f202, f210, f211, f212, f220, f221, f222,
					s);
				// time_multi += timer.getElapsedTimeInSec();
			}
			// timer.start();
			o3 = orient3D_TPI_post_exact(
				s,
				n[0], n[1], n[2],
				t20, t21, t22,
				t00, t01, t02);
			// time_multi += timer.getElapsedTimeInSec();
			/*if (o3 == 1) after21++;
				if (o3 == -1) after22++;
				if (o3 == 0) after20++;*/
		}
		if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0)
			return false;
		if (premulti == true)
			multiflag = true;
		return true;
	}
	int FastEnvelope::is_3_triangle_cut_float_fast(
		const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12,
		const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22)
	{

		Vector3 n = (tri0 - tri1).cross(tri0 - tri2) + tri0;

		if (Predicates::orient_3d(n, tri0, tri1, tri2) == 0)
		{
			logger().debug("Degeneration happens !");

			// srand(int(time(0)));
			n = { {Vector3(rand(), rand(), rand())} };
		}
		Scalar d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7;
		bool pre =
			orient3D_TPI_prefilter(
				tri0[0], tri0[1], tri0[2],
				tri1[0], tri1[1], tri1[2],
				tri2[0], tri2[1], tri2[2],
				facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
				d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

		if (pre == false)
			return 2; // means we dont know
		int o1 = orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			tri0[0], tri0[1], tri0[2],
			tri1[0], tri1[1], tri1[2]);
		int o2 = orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			tri1[0], tri1[1], tri1[2],
			tri2[0], tri2[1], tri2[2]);
		if (o1 * o2 == -1)
			return 0;
		int o3 = orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7, n[0], n[1], n[2],
			tri2[0], tri2[1], tri2[2],
			tri0[0], tri0[1], tri0[2]);
		if (o1 * o3 == -1 || o2 * o3 == -1)
			return 0;
		if (o1 * o2 * o3 == 0)
			return 2; // means we dont know
		return 1;
	}

	int FastEnvelope::seg_cut_plane(const Vector3 &seg0, const Vector3 &seg1, const Vector3 &t0, const Vector3 &t1, const Vector3 &t2)
	{
		int o1, o2;
		o1 = Predicates::orient_3d(seg0, t0, t1, t2);
		o2 = Predicates::orient_3d(seg1, t0, t1, t2);
		int op = o1 * o2;
		if (op >= 0)
		{
			return CUT_COPLANAR; //in fact, coplanar and not on this plane
		}
		return CUT_FACE;
	}

	bool FastEnvelope::is_triangle_cut_prism(const int &pindex,
		const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, std::vector<int> &cid) const
	{



		bool cut[8];
		for (int i = 0; i < 8; i++)
		{
			cut[i] = false;
		}
		int o1[8], o2[8], o3[8], ori = 0;
		std::vector<int> cutp;

		for (int i = 0; i < 8; i++)
		{

			o1[i] = Predicates::orient_3d(tri0, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			o2[i] = Predicates::orient_3d(tri1, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			o3[i] = Predicates::orient_3d(tri2, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			if (o1[i] + o2[i] + o3[i] >= 3)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1)
			{
				return false;
			}
			if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 0)
			{
				return false;
			}

			if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1)
				cutp.emplace_back(i);
		}
		if (cutp.size() == 0)
		{
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++)
		{
			if (o1[cutp[i]] * o2[cutp[i]] == -1)
			{

				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri1[0], tri1[1], tri1[2],

					envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
					envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
					envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false)
				{
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++)
				{
					if (i == j)
						continue;
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
							envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
							envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

					if (ori == 1)
						break;
				}
				if (ori != 1)
				{
					cut[cutp[i]] = true;
				}
			}
			if (cut[cutp[i]] == true)
				continue;
			ori = 0;
			if (o1[cutp[i]] * o3[cutp[i]] == -1)
			{

				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri2[0], tri2[1], tri2[2],
					envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
					envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
					envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false)
				{
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++)
				{
					if (i == j)
						continue;
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
							envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
							envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

					if (ori == 1)
						break;
				}
				if (ori != 1)
				{
					cut[cutp[i]] = true;
				}
			}

			if (cut[cutp[i]] == true)
				continue;
			ori = 0;
			if (o2[cutp[i]] * o3[cutp[i]] == -1)
			{

				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					tri1[0], tri1[1], tri1[2],
					tri2[0], tri2[1], tri2[2],
					envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
					envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
					envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false)
				{
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++)
				{
					if (i == j)
						continue;
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri1[0], tri1[1], tri1[2],
							envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
							envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
							envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

					if (ori == 1)
						break;
				}
				if (ori != 1)
				{
					cut[cutp[i]] = true;
				}
			}
		}

		if (cutp.size() <= 2)
		{
			for (int i = 0; i < 8; i++)
			{
				if (cut[i] == true)
					cid.emplace_back(i);
			}
			return true;
		}
		// triangle-facet-facet intersection
		Scalar n1, n2, n3, max3, max4, max6, max7;
		for (int i = 0; i < cutp.size(); i++)
		{
			for (int j = i + 1; j < cutp.size(); j++)
			{
				if (cut[cutp[i]] == true && cut[cutp[j]] == true)
					continue;

				int id = cutp[i] * 8 + cutp[j];
				int id0 = prism_map[id][0];
				if (id0 == -1)
					continue;
				int inter = is_3_triangle_cut_float_fast(
					tri0, tri1, tri2,
					envprism[pindex][p_face[cutp[i]][0]],
					envprism[pindex][p_face[cutp[i]][1]],
					envprism[pindex][p_face[cutp[i]][2]],
					envprism[pindex][p_face[cutp[j]][0]],
					envprism[pindex][p_face[cutp[j]][1]],
					envprism[pindex][p_face[cutp[j]][2]]);
				if (inter == 2)
				{ //we dont know if point exist or if inside of triangle
					cut[cutp[i]] == true;
					cut[cutp[j]] == true;
					continue;
				}
				if (inter == 0)
					continue; // sure not inside

				bool pre =
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

				for (int k = 0; k < cutp.size(); k++)
				{

					if (k == i || k == j)
						continue;

					ori =
						orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							envprism[pindex][p_face[cutp[k]][0]][0], envprism[pindex][p_face[cutp[k]][0]][1], envprism[pindex][p_face[cutp[k]][0]][2],
							envprism[pindex][p_face[cutp[k]][1]][0], envprism[pindex][p_face[cutp[k]][1]][1], envprism[pindex][p_face[cutp[k]][1]][2],
							envprism[pindex][p_face[cutp[k]][2]][0], envprism[pindex][p_face[cutp[k]][2]][1], envprism[pindex][p_face[cutp[k]][2]][2]);

					if (ori == 1)
						break;
				}

				if (ori != 1)
				{
					cut[cutp[i]] = true;
					cut[cutp[j]] = true;
				}
			}
		}

		for (int i = 0; i < 8; i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return true;
	}

	bool FastEnvelope::is_seg_cut_prism(const int &pindex,
		const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const
	{
		bool cut[8];
		for (int i = 0; i < 8; i++)
		{
			cut[i] = false;
		}
		int o1[8], o2[8], ori = 0;
		std::vector<int> cutp;

		for (int i = 0; i < 8; i++)
		{

			o1[i] = Predicates::orient_3d(seg0, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);
			o2[i] = Predicates::orient_3d(seg1, envprism[pindex][p_face[i][0]], envprism[pindex][p_face[i][1]], envprism[pindex][p_face[i][2]]);

			if (o1[i] + o2[i] >= 1)
			{
				return false;
			}

			if (o1[i] == 0 && o2[i] == 0)
			{
				return false;
			}

			if (o1[i] * o2[i] == -1)
				cutp.emplace_back(i);
		}
		if (cutp.size() == 0)
		{
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++)
		{

			bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
				seg0[0], seg0[1], seg0[2],
				seg1[0], seg1[1], seg1[2],
				envprism[pindex][p_face[cutp[i]][0]][0], envprism[pindex][p_face[cutp[i]][0]][1], envprism[pindex][p_face[cutp[i]][0]][2],
				envprism[pindex][p_face[cutp[i]][1]][0], envprism[pindex][p_face[cutp[i]][1]][1], envprism[pindex][p_face[cutp[i]][1]][2],
				envprism[pindex][p_face[cutp[i]][2]][0], envprism[pindex][p_face[cutp[i]][2]][1], envprism[pindex][p_face[cutp[i]][2]][2],
				a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
			if (precom == false)
			{
				cut[cutp[i]] = true;
				continue;
			}
			for (int j = 0; j < cutp.size(); j++)
			{
				if (i == j)
					continue;
				ori =
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						seg0[0], seg0[1], seg0[2],
						envprism[pindex][p_face[cutp[j]][0]][0], envprism[pindex][p_face[cutp[j]][0]][1], envprism[pindex][p_face[cutp[j]][0]][2],
						envprism[pindex][p_face[cutp[j]][1]][0], envprism[pindex][p_face[cutp[j]][1]][1], envprism[pindex][p_face[cutp[j]][1]][2],
						envprism[pindex][p_face[cutp[j]][2]][0], envprism[pindex][p_face[cutp[j]][2]][1], envprism[pindex][p_face[cutp[j]][2]][2]);

				if (ori == 1)
					break;
			}
			if (ori != 1)
			{
				cut[cutp[i]] = true;
			}
		}

		for (int i = 0; i < 8; i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::is_triangle_cut_cube(const int &cindex,
		const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, std::vector<int> &cid) const
	{

		bool cut[6];
		for (int i = 0; i < 6; i++)
		{
			cut[i] = false;
		}
		int o1[6], o2[6], o3[6], ori = 0;
		std::vector<int> cutp;

		for (int i = 0; i < 6; i++)
		{

			o1[i] = Predicates::orient_3d(tri0, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			o2[i] = Predicates::orient_3d(tri1, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			o3[i] = Predicates::orient_3d(tri2, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			if (o1[i] + o2[i] + o3[i] >= 3)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1)
			{
				return false;
			}
			if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 0)
			{
				return false;
			}

			if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1)
				cutp.emplace_back(i);
		}
		if (cutp.size() == 0)
		{
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++)
		{
			if (o1[cutp[i]] * o2[cutp[i]] == -1)
			{

				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri1[0], tri1[1], tri1[2],

					envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
					envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
					envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false)
				{
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++)
				{
					if (i == j)
						continue;
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
							envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
							envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]);

					if (ori == 1)
						break;
				}
				if (ori != 1)
				{
					cut[cutp[i]] = true;
				}
			}
			if (cut[cutp[i]] == true)
				continue;
			ori = 0;
			if (o1[cutp[i]] * o3[cutp[i]] == -1)
			{

				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					tri0[0], tri0[1], tri0[2],
					tri2[0], tri2[1], tri2[2],
					envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
					envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
					envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false)
				{
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++)
				{
					if (i == j)
						continue;
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri0[0], tri0[1], tri0[2],
							envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
							envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
							envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]);

					if (ori == 1)
						break;
				}
				if (ori != 1)
				{
					cut[cutp[i]] = true;
				}
			}

			if (cut[cutp[i]] == true)
				continue;
			ori = 0;
			if (o2[cutp[i]] * o3[cutp[i]] == -1)
			{

				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					tri1[0], tri1[1], tri1[2],
					tri2[0], tri2[1], tri2[2],
					envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
					envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
					envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false)
				{
					cut[cutp[i]] = true;
					continue;
				}
				for (int j = 0; j < cutp.size(); j++)
				{
					if (i == j)
						continue;
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							tri1[0], tri1[1], tri1[2],
							envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
							envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
							envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]);

					if (ori == 1)
						break;
				}
				if (ori != 1)
				{
					cut[cutp[i]] = true;
				}
			}
		}

		if (cutp.size() <= 2)
		{
			for (int i = 0; i < 6; i++)
			{
				if (cut[i] == true)
					cid.emplace_back(i);
			}
			return true;
		}
		// triangle-facet-facet intersection
		Scalar n1, n2, n3, max3, max4, max6, max7;
		for (int i = 0; i < cutp.size(); i++)
		{
			for (int j = i + 1; j < cutp.size(); j++)
			{
				if (cut[cutp[i]] == true && cut[cutp[j]] == true)
					continue;

				int id = cutp[i] * 6 + cutp[j];
				int id0 = cubic_map[id][0];
				if (id0 == -1)
					continue;
				int inter = is_3_triangle_cut_float_fast(
					tri0, tri1, tri2,
					envcubic[cindex][c_face[cutp[i]][0]],
					envcubic[cindex][c_face[cutp[i]][1]],
					envcubic[cindex][c_face[cutp[i]][2]],
					envcubic[cindex][c_face[cutp[j]][0]],
					envcubic[cindex][c_face[cutp[j]][1]],
					envcubic[cindex][c_face[cutp[j]][2]]);
				if (inter == 2)
				{ //we dont know if point exist or if inside of triangle
					cut[cutp[i]] == true;
					cut[cutp[j]] == true;
					continue;
				}
				if (inter == 0)
					continue; // sure not inside

				bool pre =
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

				for (int k = 0; k < cutp.size(); k++)
				{

					if (k == i || k == j)
						continue;

					ori =
						orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							envcubic[cindex][c_face[cutp[k]][0]][0], envcubic[cindex][c_face[cutp[k]][0]][1], envcubic[cindex][c_face[cutp[k]][0]][2],
							envcubic[cindex][c_face[cutp[k]][1]][0], envcubic[cindex][c_face[cutp[k]][1]][1], envcubic[cindex][c_face[cutp[k]][1]][2],
							envcubic[cindex][c_face[cutp[k]][2]][0], envcubic[cindex][c_face[cutp[k]][2]][1], envcubic[cindex][c_face[cutp[k]][2]][2]);

					if (ori == 1)
						break;
				}

				if (ori != 1)
				{
					cut[cutp[i]] = true;
					cut[cutp[j]] = true;
				}
			}
		}

		for (int i = 0; i < 6; i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::is_triangle_cut_envelope_polyhedra(const int &cindex,
		const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, std::vector<int> &cid) const
	{
		int number;
		if (envelope[cindex].size() == 12) number = 8;
		if (envelope[cindex].size() == 8) number = 6;
		std::vector<bool> cut;
		cut.resize(number);
		for (int i = 0; i < number; i++)
		{
			cut[i] = false;
		}
		std::vector<int> o1, o2, o3, cutp;
		o1.resize(number);
		o2.resize(number);
		o3.resize(number);
		int  ori = 0,face[3];
		

		for (int i = 0; i < number; i++)
		{
			if (envelope[cindex].size() == 12) {
				face[0] = p_face[i][0];
				face[1] = p_face[i][1];
				face[2] = p_face[i][2];
			}
			if (envelope[cindex].size() == 8) {
				face[0] = c_face[i][0];
				face[1] = c_face[i][1];
				face[2] = c_face[i][2];
			}
			o1[i] = Predicates::orient_3d(tri0, envelope[cindex][face[0]], envelope[cindex][face[1]], envelope[cindex][face[2]]);
			o2[i] = Predicates::orient_3d(tri1, envelope[cindex][face[0]], envelope[cindex][face[1]], envelope[cindex][face[2]]);
			o3[i] = Predicates::orient_3d(tri2, envelope[cindex][face[0]], envelope[cindex][face[1]], envelope[cindex][face[2]]);
			if (o1[i] + o2[i] + o3[i] >= 3)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1)
			{
				return false;
			}
			if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0)
			{
				return false;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 0)
			{
				return false;
			}

			if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1)
				cutp.emplace_back(i);
		}
		if (cutp.size() == 0)
		{
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		std::array<Scalar, 2> seg00, seg01, seg02, seg10, seg11, seg12;
		for (int i = 0; i < cutp.size(); i++)
		{
			int temp = 0;
			if (o1[cutp[i]] * o2[cutp[i]] == -1) {
				seg00[temp] = tri0[0];
				seg01[temp] = tri0[1];
				seg02[temp] = tri0[2];
				seg10[temp] = tri1[0];
				seg11[temp] = tri1[1];
				seg12[temp] = tri1[2];
				temp++;
			}
			if (o1[cutp[i]] * o3[cutp[i]] == -1) {
				seg00[temp] = tri0[0];
				seg01[temp] = tri0[1];
				seg02[temp] = tri0[2];
				seg10[temp] = tri2[0];
				seg11[temp] = tri2[1];
				seg12[temp] = tri2[2];
				temp++;
			}
			if (o2[cutp[i]] * o3[cutp[i]] == -1) {
				seg00[temp] = tri1[0];
				seg01[temp] = tri1[1];
				seg02[temp] = tri1[2];
				seg10[temp] = tri2[0];
				seg11[temp] = tri2[1];
				seg12[temp] = tri2[2];
				temp++;
			}
			if (temp == 3) std::cout << "wrong here" << std::endl;
			if (temp == 0) std::cout << "wrong here" << std::endl;
			for (int k = 0; k < 2; k++)
			{
				if (envelope[cindex].size() == 12) {
					face[0] = p_face[cutp[i]][0];
					face[1] = p_face[cutp[i]][1];
					face[2] = p_face[cutp[i]][2];
				}
				if (envelope[cindex].size() == 8) {
					face[0] = c_face[cutp[i]][0];
					face[1] = c_face[cutp[i]][1];
					face[2] = c_face[cutp[i]][2];
				}
				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					seg00[k], seg01[k], seg02[k],
					seg10[k], seg11[k], seg12[k],
					envelope[cindex][face[0]][0], envelope[cindex][face[0]][1], envelope[cindex][face[0]][2],
					envelope[cindex][face[1]][0], envelope[cindex][face[1]][1], envelope[cindex][face[1]][2],
					envelope[cindex][face[2]][0], envelope[cindex][face[2]][1], envelope[cindex][face[2]][2],
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
				if (precom == false)
				{
					cut[cutp[i]] = true;
					break;
				}
				for (int j = 0; j < cutp.size(); j++)
				{
					if (i == j)
						continue;

					if (envelope[cindex].size() == 12) {
						face[0] = p_face[cutp[j]][0];
						face[1] = p_face[cutp[j]][1];
						face[2] = p_face[cutp[j]][2];
					}
					if (envelope[cindex].size() == 8) {
						face[0] = c_face[cutp[j]][0];
						face[1] = c_face[cutp[j]][1];
						face[2] = c_face[cutp[j]][2];
					}
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							seg00[k], seg01[k], seg02[k],
							envelope[cindex][face[0]][0], envelope[cindex][face[0]][1], envelope[cindex][face[0]][2],
							envelope[cindex][face[1]][0], envelope[cindex][face[1]][1], envelope[cindex][face[1]][2],
							envelope[cindex][face[2]][0], envelope[cindex][face[2]][1], envelope[cindex][face[2]][2]);

					if (ori == 1)
						break;
				}
				if (ori != 1)
				{
					cut[cutp[i]] = true;
					break;
				}
			}
	
			ori = 0;
		}

		if (cutp.size() <= 2)
		{
			for (int i = 0; i < number; i++)
			{
				if (cut[i] == true)
					cid.emplace_back(i);
			}
			return true;
		}
		// triangle-facet-facet intersection
		Scalar n1, n2, n3, max3, max4, max6, max7;
		int face1[3], face2[3];
		for (int i = 0; i < cutp.size(); i++)
		{
			for (int j = i + 1; j < cutp.size(); j++)
			{
				if (cut[cutp[i]] == true && cut[cutp[j]] == true)
					continue;

				int id,id0;
				if (envelope[cindex].size() == 12) {
					id= cutp[i] * 8 + cutp[j];
					id0 = prism_map[id][0];
					if (id0 == -1) continue;
					face1[0] = p_face[cutp[i]][0];
					face1[1] = p_face[cutp[i]][1];
					face1[2] = p_face[cutp[i]][2];
					face2[0] = p_face[cutp[j]][0];
					face2[1] = p_face[cutp[j]][1];
					face2[2] = p_face[cutp[j]][2];
					
				}
				if (envelope[cindex].size() == 8) {
					id = cutp[i] * 6 + cutp[j];
					id0 = cubic_map[id][0];
					if (id0 == -1) continue;
					face1[0] = c_face[cutp[i]][0];
					face1[1] = c_face[cutp[i]][1];
					face1[2] = c_face[cutp[i]][2];
					face2[0] = c_face[cutp[j]][0];
					face2[1] = c_face[cutp[j]][1];
					face2[2] = c_face[cutp[j]][2];
				}

				

				int inter = is_3_triangle_cut_float_fast(
					tri0, tri1, tri2,
					envelope[cindex][face1[0]],
					envelope[cindex][face1[1]],
					envelope[cindex][face1[2]],
					envelope[cindex][face2[0]],
					envelope[cindex][face2[1]],
					envelope[cindex][face2[2]]);
				if (inter == 2)
				{ //we dont know if point exist or if inside of triangle
					cut[cutp[i]] == true;
					cut[cutp[j]] == true;
					continue;
				}
				if (inter == 0)
					continue; // sure not inside

				bool pre =
					orient3D_TPI_prefilter(
						tri0[0], tri0[1], tri0[2],
						tri1[0], tri1[1], tri1[2],
						tri2[0], tri2[1], tri2[2],
						envelope[cindex][face1[0]][0], envelope[cindex][face1[0]][1], envelope[cindex][face1[0]][2],
						envelope[cindex][face1[1]][0], envelope[cindex][face1[1]][1], envelope[cindex][face1[1]][2],
						envelope[cindex][face1[2]][0], envelope[cindex][face1[2]][1], envelope[cindex][face1[2]][2],
						envelope[cindex][face2[0]][0], envelope[cindex][face2[0]][1], envelope[cindex][face2[0]][2],
						envelope[cindex][face2[1]][0], envelope[cindex][face2[1]][1], envelope[cindex][face2[1]][2],
						envelope[cindex][face2[2]][0], envelope[cindex][face2[2]][1], envelope[cindex][face2[2]][2],
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

				for (int k = 0; k < cutp.size(); k++)
				{

					if (k == i || k == j)
						continue;
					if (envelope[cindex].size() == 12) {
						
						face[0] = p_face[cutp[k]][0];
						face[1] = p_face[cutp[k]][1];
						face[2] = p_face[cutp[k]][2];	
					}
					if (envelope[cindex].size() == 8) {

						face[0] = c_face[cutp[k]][0];
						face[1] = c_face[cutp[k]][1];
						face[2] = c_face[cutp[k]][2];
					}
					ori =
						orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							envelope[cindex][face[0]][0], envelope[cindex][face[0]][1], envelope[cindex][face[0]][2],
							envelope[cindex][face[1]][0], envelope[cindex][face[1]][1], envelope[cindex][face[1]][2],
							envelope[cindex][face[2]][0], envelope[cindex][face[2]][1], envelope[cindex][face[2]][2]);

					if (ori == 1)
						break;
				}

				if (ori != 1)
				{
					cut[cutp[i]] = true;
					cut[cutp[j]] = true;
				}
			}
		}

		for (int i = 0; i < number; i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::is_seg_cut_cube(const int &cindex,
		const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const
	{
		bool cut[6];
		for (int i = 0; i < 6; i++)
		{
			cut[i] = false;
		}
		int o1[6], o2[6], ori = 0;
		std::vector<int> cutp;

		for (int i = 0; i < 6; i++)
		{

			o1[i] = Predicates::orient_3d(seg0, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);
			o2[i] = Predicates::orient_3d(seg1, envcubic[cindex][c_face[i][0]], envcubic[cindex][c_face[i][1]], envcubic[cindex][c_face[i][2]]);

			if (o1[i] + o2[i] >= 1)
			{
				return false;
			}

			if (o1[i] == 0 && o2[i] == 0)
			{
				return false;
			}

			if (o1[i] * o2[i] == -1)
				cutp.emplace_back(i);
		}
		if (cutp.size() == 0)
		{
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++)
		{

			bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
				seg0[0], seg0[1], seg0[2],
				seg1[0], seg1[1], seg1[2],
				envcubic[cindex][c_face[cutp[i]][0]][0], envcubic[cindex][c_face[cutp[i]][0]][1], envcubic[cindex][c_face[cutp[i]][0]][2],
				envcubic[cindex][c_face[cutp[i]][1]][0], envcubic[cindex][c_face[cutp[i]][1]][1], envcubic[cindex][c_face[cutp[i]][1]][2],
				envcubic[cindex][c_face[cutp[i]][2]][0], envcubic[cindex][c_face[cutp[i]][2]][1], envcubic[cindex][c_face[cutp[i]][2]][2],
				a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
			if (precom == false)
			{
				cut[cutp[i]] = true;
				continue;
			}
			for (int j = 0; j < cutp.size(); j++)
			{
				if (i == j)
					continue;
				ori =
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						seg0[0], seg0[1], seg0[2],
						envcubic[cindex][c_face[cutp[j]][0]][0], envcubic[cindex][c_face[cutp[j]][0]][1], envcubic[cindex][c_face[cutp[j]][0]][2],
						envcubic[cindex][c_face[cutp[j]][1]][0], envcubic[cindex][c_face[cutp[j]][1]][1], envcubic[cindex][c_face[cutp[j]][1]][2],
						envcubic[cindex][c_face[cutp[j]][2]][0], envcubic[cindex][c_face[cutp[j]][2]][1], envcubic[cindex][c_face[cutp[j]][2]][2]);

				if (ori == 1)
					break;
			}
			if (ori != 1)
			{
				cut[cutp[i]] = true;
			}
		}

		for (int i = 0; i < 6; i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::is_seg_cut_polyhedra(const int &cindex,
		const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const
	{

		int number;
		if (envelope[cindex].size() == 12) {
			number = 8;
		}
		if (envelope[cindex].size() == 8) {
			number = 6;
		}
		std::vector<bool> cut;
		cut.resize(number);
		for (int i = 0; i < number; i++)
		{
			cut[i] = false;
		}
		std::vector<int> o1, o2;
		o1.resize(number);
		o2.resize(number);
		int ori = 0;
		std::vector<int> cutp;
		int face[3];
		for (int i = 0; i < number; i++)
		{
			if (envelope[cindex].size() == 12) {
				face[0] = p_face[i][0];
				face[1] = p_face[i][1];
				face[2] = p_face[i][2];
			}
			if (envelope[cindex].size() == 8) {
				face[0] = c_face[i][0];
				face[1] = c_face[i][1];
				face[2] = c_face[i][2];
			}
			o1[i] = Predicates::orient_3d(seg0, envelope[cindex][face[0]], envelope[cindex][face[1]], envelope[cindex][face[2]]);
			o2[i] = Predicates::orient_3d(seg1, envelope[cindex][face[0]], envelope[cindex][face[1]], envelope[cindex][face[2]]);

			if (o1[i] + o2[i] >= 1)
			{
				return false;
			}

			if (o1[i] == 0 && o2[i] == 0)
			{
				return false;
			}

			if (o1[i] * o2[i] == -1)
				cutp.emplace_back(i);
		}
		if (cutp.size() == 0)
		{
			return false;
		}

		Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
		for (int i = 0; i < cutp.size(); i++)
		{
			if (envelope[cindex].size() == 12) {
				face[0] = p_face[cutp[i]][0];
				face[1] = p_face[cutp[i]][1];
				face[2] = p_face[cutp[i]][2];
			}
			if (envelope[cindex].size() == 8) {
				face[0] = c_face[cutp[i]][0];
				face[1] = c_face[cutp[i]][1];
				face[2] = c_face[cutp[i]][2];
			}
			bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
				seg0[0], seg0[1], seg0[2],
				seg1[0], seg1[1], seg1[2],
				envelope[cindex][face[0]][0], envelope[cindex][face[0]][1], envelope[cindex][face[0]][2],
				envelope[cindex][face[1]][0], envelope[cindex][face[1]][1], envelope[cindex][face[1]][2],
				envelope[cindex][face[2]][0], envelope[cindex][face[2]][1], envelope[cindex][face[2]][2],
				a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
			if (precom == false)
			{
				cut[cutp[i]] = true;
				continue;
			}
			for (int j = 0; j < cutp.size(); j++)
			{
				if (i == j)
					continue;
				if (envelope[cindex].size() == 12) {
					face[0] = p_face[cutp[j]][0];
					face[1] = p_face[cutp[j]][1];
					face[2] = p_face[cutp[j]][2];
				}
				if (envelope[cindex].size() == 8) {
					face[0] = c_face[cutp[j]][0];
					face[1] = c_face[cutp[j]][1];
					face[2] = c_face[cutp[j]][2];
				}
				ori =
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						seg0[0], seg0[1], seg0[2],
						envelope[cindex][face[0]][0], envelope[cindex][face[0]][1], envelope[cindex][face[0]][2],
						envelope[cindex][face[1]][0], envelope[cindex][face[1]][1], envelope[cindex][face[1]][2],
						envelope[cindex][face[2]][0], envelope[cindex][face[2]][1], envelope[cindex][face[2]][2]);

				if (ori == 1)
					break;
			}
			if (ori != 1)
			{
				cut[cutp[i]] = true;
			}
		}

		for (int i = 0; i < number; i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const
	{

		int ori,number,face[3];
		
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump)
				continue;
			if (envelope[prismindex[i]].size() == 12) {
				number = 8;
			}
			if (envelope[prismindex[i]].size() == 8) {
				number = 6;
			}
			for (int j = 0; j < number; j++)
			{
				if (envelope[prismindex[i]].size() == 12) {
					face[0] = p_face[j][0];
					face[1] = p_face[j][1];
					face[2] = p_face[j][2];
				}
				if (envelope[prismindex[i]].size() == 8) {
					face[0] = c_face[j][0];
					face[1] = c_face[j][1];
					face[2] = c_face[j][2];
				}
				ori = Predicates::orient_3d(envelope[prismindex[i]][face[0]], envelope[prismindex[i]][face[1]], envelope[prismindex[i]][face[2]], point);
				if (ori == -1 || ori == 0)
				{
					break;
				}
				if (j == number - 1)
				{

					return false;
				}
			}
			
		}

		return true;
	}

	int FastEnvelope::is_triangle_degenerated(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2)
	{
		const auto to_2d = [](const Vector3 &p, int t)
		{
			return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
		};

		Vector3 a = triangle0 - triangle1, b = triangle0 - triangle2;
		Vector3 normal = a.cross(b);
		Scalar nbr = normal.norm();

		if (nbr > SCALAR_ZERO)
		{
			return NOT_DEGENERATED;
		}
		int ori;
		std::array<Vector2, 3> p;
		for (int j = 0; j < 3; j++)
		{

			p[0] = to_2d(triangle0, j);
			p[1] = to_2d(triangle1, j);
			p[2] = to_2d(triangle2, j);

			ori = Predicates::orient_2d(p[0], p[1], p[2]);
			if (ori != 0)
			{
				return NERLY_DEGENERATED;
			}
		}

		if (triangle0[0] != triangle1[0] || triangle0[1] != triangle1[1] || triangle0[2] != triangle1[2])
		{
			return DEGENERATED_SEGMENT;
		}
		if (triangle0[0] != triangle2[0] || triangle0[1] != triangle2[1] || triangle0[2] != triangle2[2])
		{
			return DEGENERATED_SEGMENT;
		}
		return DEGENERATED_POINT;
	}

	void FastEnvelope::halfspace_generation(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, std::vector<std::array<Vector3, 12>> &envprism, std::vector<std::array<Vector3, 8>> &envbox, const Scalar &epsilon)
	{
		const auto dot_sign = [](const Vector3 &a, const Vector3 &b)
		{
			Scalar t = a.dot(b);
			if (t > SCALAR_ZERO)
				return 1;
			if (t < -1 * SCALAR_ZERO)
				return -1;
			Multiprecision a0(a[0]), a1(a[1]), a2(a[2]), b0(b[0]), b1(b[1]), b2(b[2]), dot;
			dot = a0 * b0 + a1 * b1 + a2 * b2;
			if (dot.get_sign() > 0)
				return 1;
			if (dot.get_sign() < 0)
				return -1;
			return 0;
		};

		static const fastEnvelope::Vector3 origin = fastEnvelope::Vector3(0, 0, 0);

		envprism.reserve(m_faces.size());
		Vector3 AB, AC, BC, normal, vector1, ABn;
		std::array<Vector3, 6> polygon;
		std::array<Vector3, 12> polygonoff;
		std::array<Vector3, 8> box;
		static const std::array<Vector3, 8> boxorder = {
			{
				{1, 1, 1},
				{-1, 1, 1},
				{-1, -1, 1},
				{1, -1, 1},
				{1, 1, -1},
				{-1, 1, -1},
				{-1, -1, -1},
				{1, -1, -1},
			} };

		Scalar tolerance = epsilon / sqrt(3);
		Scalar de;

		for (int i = 0; i < m_faces.size(); i++)
		{
			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			de = is_triangle_degenerated(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);

			if (de == DEGENERATED_POINT)
			{
				logger().debug("Envelope Triangle Degeneration- Point");
				for (int j = 0; j < 8; j++)
				{
					box[j] = m_ver[m_faces[i][0]] + boxorder[j] * tolerance;
				}
				envbox.emplace_back(box);
				continue;
			}
			if (de == DEGENERATED_SEGMENT)
			{
				logger().debug("Envelope Triangle Degeneration- Segment");
				Scalar length1 = AB.norm(), length2 = AC.norm(), length3 = BC.norm();
				if (length1 >= length2 && length1 >= length3)
				{
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], tolerance, box);
					envbox.emplace_back(box);
				}
				if (length2 >= length1 && length2 >= length3)
				{
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][2]], tolerance, box);
					envbox.emplace_back(box);
				}
				if (length3 >= length1 && length3 >= length2)
				{
					seg_cube(m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tolerance, box);
					envbox.emplace_back(box);
				}
				continue;
			}
			if (de == NERLY_DEGENERATED)
			{
				logger().debug("Envelope Triangle Degeneration- Nearly");
				//std::cout << std::setprecision(17) <<m_ver[m_faces[i][0]][0] << " " << m_ver[m_faces[i][0]][1] << " " << m_ver[m_faces[i][0]][2] << std::endl;
				//std::cout << std::setprecision(17) <<m_ver[m_faces[i][1]][0] << " " << m_ver[m_faces[i][1]][1] << " " << m_ver[m_faces[i][1]][2] << std::endl;
				//std::cout << std::setprecision(17) <<m_ver[m_faces[i][2]][0] << " " << m_ver[m_faces[i][2]][1] << " " << m_ver[m_faces[i][2]][2] << std::endl;

				normal = accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][0]], m_ver[m_faces[i][2]]);
				//std::cout << "pass1" << std::endl;
				vector1 = accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], origin, normal);
				//std::cout << "pass2" << std::endl;
			}
			else
			{
				normal = AB.cross(AC).normalized();
				vector1 = AB.cross(normal).normalized();
			}

			ABn = AB.normalized();
			polygon[0] = m_ver[m_faces[i][0]] + (vector1 - ABn) * tolerance;
			polygon[1] = m_ver[m_faces[i][1]] + (vector1 + ABn) * tolerance;
			if (dot_sign(AB, BC) < 0)
			{
				polygon[2] = m_ver[m_faces[i][1]] + (-vector1 + ABn) * tolerance;
				polygon[3] = m_ver[m_faces[i][2]] + (-vector1 + ABn) * tolerance;
				polygon[4] = m_ver[m_faces[i][2]] + (-vector1 - ABn) * tolerance;
				if (dot_sign(AB, AC) < 0)
				{
					polygon[5] = m_ver[m_faces[i][2]] + (vector1 - ABn) * tolerance;
				}
				else
				{
					polygon[5] = m_ver[m_faces[i][0]] + (-vector1 - ABn) * tolerance;
				}
			}
			else
			{
				polygon[2] = m_ver[m_faces[i][2]] + (vector1 + ABn) * tolerance;
				polygon[3] = m_ver[m_faces[i][2]] + (-vector1 + ABn) * tolerance;
				polygon[4] = m_ver[m_faces[i][2]] + (-vector1 - ABn) * tolerance;
				polygon[5] = m_ver[m_faces[i][0]] + (-vector1 - ABn) * tolerance;
			}

			for (int j = 0; j < 6; j++)
			{
				polygonoff[j] = polygon[j] + normal * tolerance;
			}
			for (int j = 6; j < 12; j++)
			{
				polygonoff[j] = polygon[j - 6] - normal * tolerance;
			}
			envprism.emplace_back(polygonoff);
		}
	}

	void FastEnvelope::seg_cube(const Vector3 &p1, const Vector3 &p2, const Scalar &width, std::array<Vector3, 8> &envbox)
	{
		Vector3 v1, v2, v = p2 - p1; //p12
		if (v[0] != 0)
		{
			v1 = Vector3((0 - v[1] - v[2]) / v[0], 1, 1);
		}
		else
		{
			if (v[1] != 0)
			{
				v1 = Vector3(1, (0 - v[2]) / v[1], 1);
			}
			else
			{
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
		envbox[3] = p2 + width * (v + v1 - v2); //right hand out direction
		envbox[4] = p1 + width * (-v + v1 + v2);
		envbox[5] = p1 + width * (-v - v1 + v2);
		envbox[6] = p1 + width * (-v - v1 - v2);
		envbox[7] = p1 + width * (-v + v1 - v2); //right hand in direction
	}

	Vector3 FastEnvelope::accurate_normal_vector(const Vector3 &p0, const Vector3 &p1,
		const Vector3 &q0, const Vector3 &q1)
	{

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
		if (ssum == 0)
		{
			logger().debug("divided by zero in accuratexxx");
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

		x = x / length;
		y = y / length;
		z = z / length;

		Scalar fx = x.to_double(), fy = y.to_double(), fz = z.to_double();
		return Vector3(fx, fy, fz);
	}





} // namespace fastEnvelope

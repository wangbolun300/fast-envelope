#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Predicates.hpp>
#include <fastenvelope/Logger.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_reorder.h>

#include <igl/Timer.h>
#include<igl/write_triangle_mesh.h>
#include <fstream>
#include <istream>

int dbg1=0, dbg2=0,dbg3=0,dbg4=0,dbgout1=0,dbgout2=0,dbgout3=0,dbgout4=0,dbgout5=0;
int ct1=0, ct2 = 0, ct3 = 0, ct4 = 0, ct5 = 0;
double time1 = 0, time2 = 0, time3 = 0, time4 = 0, time5 = 0, timesearch = 0, timecheck = 0,
time6 = 0, time7 = 0, time8 = 0, time9 = 0, time10 = 0, time11 = 0, time12 = 0, time13 = 0;
static const std::array<std::array<int, 2>, 3> triseg = {
	{{{0, 1}}, {{0, 2}}, {{1, 2}}}
};
std::vector<int> tempid;

namespace fastEnvelope
{
	namespace {
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

		//delete me in the future
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
	void FastEnvelope::printnumber() {
		//std::cout << "where checking is out " << dbgout1 << " " << dbgout2 << " " << dbgout3 << " " << dbgout4 << " " << dbgout5 << std::endl;
		std::cout << "time for search, "<< timesearch << std::endl;
		std::cout << "time for check, " << timecheck << std::endl;
		std::cout << "time for vertices, " << time1 << std::endl;
		std::cout << "time for degeneration, " << time2 << std::endl;
		std::cout << "time for intersection, " << time3 << std::endl;
		std::cout << "time for lpi, " << time4 << std::endl;
		std::cout << "time for local tree, " << time5 << std::endl;
		std::cout << "time for tpi all, " << time6 << std::endl;
		std::cout << "time for tpi post, " << time7 << std::endl;
		std::cout << "time for tpi part1 multi post, " << time12 << std::endl;
		std::cout << "time for tpi number, " << ct1 << std::endl;
		std::cout << "time for tpi part1, " << time8 << std::endl;
		std::cout << "time for tpi part2, " << time9 << std::endl;
		std::cout << "time for tpi multi post time, " << time10 << std::endl;
		std::cout << "time for tpi multi post number, " << ct2 << std::endl;
		std::cout << "time for see if point on prism, " << time11 << std::endl;
		std::cout << "time for local search, " << time13 << std::endl;
	}
	void FastEnvelope::reset_time() {
		time1 = 0; time2 = 0; time3 = 0; time4 = 0; time5 = 0; timesearch = 0; timecheck = 0;
		time6 = 0; time7 = 0; time8 = 0; time9 = 0; time10 = 0;
		time11 = 0; time12 = 0; time13 = 0;
		ct1 = 0; ct2 = 0;
		std::cout << "\n" << std::endl;
	}

	FastEnvelope::FastEnvelope(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, const Scalar eps)
	{
		igl::Timer timer;
		
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
		std::vector<Vector3> ver_new;
		std::vector<Vector3i> faces_new;
		timer.start();
		GEO::Mesh M;
		ver_new.clear();
		faces_new.clear();

		to_geogram_mesh(m_ver, m_faces, M);
		GEO::mesh_reorder(M, GEO::MESH_ORDER_MORTON);
		from_geogram_mesh(M, ver_new, faces_new);
		timer.stop();
		logger().info("Resorting mesh time {}s", timer.getElapsedTimeInSec());
		/////////////////just for debug


		/////////////////
		timer.start();
		
		halfspace_init(ver_new, faces_new, halfspace, cornerlist, epsilon);
		
		
		timer.stop();
		logger().info("halfspace and bounding box initialize time {}s", timer.getElapsedTimeInSec());


		timer.start();
		tree.init_envelope(cornerlist);
		timer.stop();
		logger().info("Tree init time {}s", timer.getElapsedTimeInSec());

		//initializing types
		initFPU();

		logger().debug("halfspace size {}", halfspace.size());
	}

	bool FastEnvelope::is_outside(const std::array<Vector3, 3> &triangle) const
	{
		
		igl::Timer timer;
		timer.start();
		std::vector<unsigned int> querylist;
		tree.facet_in_envelope(triangle[0], triangle[1], triangle[2], querylist);
		timesearch += timer.getElapsedTimeInSec();
		timer.start();
		const auto res = FastEnvelopeTestImplicit(triangle, querylist);
		timecheck += timer.getElapsedTimeInSec();

		return res;
	}

	/*void FastEnvelope::print_prisms(const std::array<Vector3, 3> &triangle, const std::string &path) const
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
	}*/

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
		igl::Timer timer;

		int jump1, jump2;
		std::cout << "prism size " << prismindex.size() << std::endl;
		
		std::vector<DATA_TPI> tpilist; 
		std::vector<unsigned int> filted_intersection; filted_intersection.reserve(prismindex.size() / 3);
		std::vector<std::vector<int>>intersect_face; intersect_face.reserve(prismindex.size() / 3);
		bool out, cut;

		int inter, inter1, record1, record2,

			tti; //triangle-triangle intersection

		jump1 = -1;

		

		timer.start();
		for (int i = 0; i < 3; i++)

		{

			out = point_out_prism(triangle[i], prismindex, jump1);

			if (out)
			{
				time1 += timer.getElapsedTimeInSec();
				return true;
			}
		}
		time1 += timer.getElapsedTimeInSec();

		if (prismindex.size() == 1)
			return false;

		////////////////////degeneration fix
		
		timer.start();
		int degeneration = is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

		if (degeneration == DEGENERATED_POINT)
		{ //case 1 degenerate to a point
			return false;
		}

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
						
						Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
						bool precom = orient3D_LPI_prefilter( //
							triangle[triseg[we][0]][0], triangle[triseg[we][0]][1], triangle[triseg[we][0]][2],
							triangle[triseg[we][1]][0], triangle[triseg[we][1]][1], triangle[triseg[we][1]][2],
							halfspace[prismindex[i]][cid[j]][0][0], halfspace[prismindex[i]][cid[j]][0][1], halfspace[prismindex[i]][cid[j]][0][2],
							halfspace[prismindex[i]][cid[j]][1][0], halfspace[prismindex[i]][cid[j]][1][1], halfspace[prismindex[i]][cid[j]][1][2],
							halfspace[prismindex[i]][cid[j]][2][0], halfspace[prismindex[i]][cid[j]][2][1], halfspace[prismindex[i]][cid[j]][2][2],
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
						if (precom == true)
						{
							inter = Implicit_Seg_Facet_interpoint_Out_Prism_double(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
								triangle[triseg[we][0]], triangle[triseg[we][1]],
								halfspace[prismindex[i]][cid[j]][0], halfspace[prismindex[i]][cid[j]][1], halfspace[prismindex[i]][cid[j]][2],
								prismindex, jump1);
							if (inter == 1)
							{
								time2 += timer.getElapsedTimeInSec();
								return true;
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
					time2 += timer.getElapsedTimeInSec();
					return 1;
			}
			time2 += timer.getElapsedTimeInSec();
			return 0;
		}
		//
		////////////////////////////////degeneration fix over
		time2 += timer.getElapsedTimeInSec();
		// logger().info("Bit time? time {}s", timer.getElapsedTimeInSec())

		// timer_bigpart.start();
		timer.start();
		std::vector<int> cidl; cidl.reserve(8);
		for (int i = 0; i < prismindex.size(); i++) {
			tti = is_triangle_cut_envelope_polyhedra(prismindex[i],
				triangle[0], triangle[1], triangle[2], cidl);
			if (tti == 2) {
				time3 += timer.getElapsedTimeInSec();
				return false;
			}
			if (tti == 1 && cidl.size() > 0) {
				 
				filted_intersection.emplace_back(prismindex[i]);
				intersect_face.emplace_back(cidl);
				dbg3++;
			}
		}
		if (filted_intersection.size() == 0) {
			time3 += timer.getElapsedTimeInSec();
			return false;//inside
		}
		std::vector<std::vector<bool>>face_flags;
		face_flags.resize(intersect_face.size());
		for (int i = 0; i < intersect_face.size(); i++) {
			face_flags[i].resize(intersect_face[i].size());
			for (int j = 0; j < intersect_face[i].size(); j++) {
				face_flags[i][j] = false;
			}
		}
		time3 += timer.getElapsedTimeInSec();

		timer.start();
		std::vector<int> ast, astf;
		ast.reserve(50); astf.reserve(50);

		int check_id;
		for (int i = 0; i < filted_intersection.size(); i++)
		{
			jump1 = filted_intersection[i];

			for (int j = 0; j < intersect_face[i].size(); j++) {
				if (face_flags[i][j]) continue;
				for (int k = 0; k < 3; k++) {
					tti = seg_cut_plane(triangle[triseg[k][0]], triangle[triseg[k][1]],
						halfspace[filted_intersection[i]][intersect_face[i][j]][0], halfspace[filted_intersection[i]][intersect_face[i][j]][1], halfspace[filted_intersection[i]][intersect_face[i][j]][2] );

					if (tti != CUT_FACE) continue;

					Scalar a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;
					bool precom = orient3D_LPI_prefilter( //
						triangle[triseg[k][0]][0], triangle[triseg[k][0]][1], triangle[triseg[k][0]][2],
						triangle[triseg[k][1]][0], triangle[triseg[k][1]][1], triangle[triseg[k][1]][2],
						halfspace[filted_intersection[i]][intersect_face[i][j]][0][0], halfspace[filted_intersection[i]][intersect_face[i][j]][0][1], halfspace[filted_intersection[i]][intersect_face[i][j]][0][2],
						halfspace[filted_intersection[i]][intersect_face[i][j]][1][0], halfspace[filted_intersection[i]][intersect_face[i][j]][1][1], halfspace[filted_intersection[i]][intersect_face[i][j]][1][2],
						halfspace[filted_intersection[i]][intersect_face[i][j]][2][0], halfspace[filted_intersection[i]][intersect_face[i][j]][2][1], halfspace[filted_intersection[i]][intersect_face[i][j]][2][2],
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
					if (precom)
					{

						inter = Implicit_Seg_Facet_interpoint_Out_Prism_double_return_id(a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							triangle[triseg[k][0]], triangle[triseg[k][1]],
							halfspace[filted_intersection[i]][intersect_face[i][j]][0], halfspace[filted_intersection[i]][intersect_face[i][j]][1], halfspace[filted_intersection[i]][intersect_face[i][j]][2],
							filted_intersection, jump1,check_id);
						if (inter == 1)
						{
							dbgout1++;
							time4 += timer.getElapsedTimeInSec();
							return true;
						}
						if (inter == 0) {
							for (int h = 0; h < 3; h++) {
								if (k == h) continue;
								tti = seg_cut_plane(triangle[triseg[h][0]], triangle[triseg[h][1]],
									halfspace[filted_intersection[i]][intersect_face[i][j]][0], halfspace[filted_intersection[i]][intersect_face[i][j]][1], halfspace[filted_intersection[i]][intersect_face[i][j]][2]);
								if (tti != CUT_FACE) continue;
								int check = Implicit_Seg_Facet_interpoint_Out_Prism_check_id(triangle[triseg[h][0]], triangle[triseg[h][1]],
									halfspace[filted_intersection[i]][intersect_face[i][j]][0], halfspace[filted_intersection[i]][intersect_face[i][j]][1], halfspace[filted_intersection[i]][intersect_face[i][j]][2], check_id);
								if (check == 0) face_flags[i][j] = true;
							}
						}
					}
					else
					{
						datalpi.segid = k;
						datalpi.prismid = filted_intersection[i];
						datalpi.facetid = intersect_face[i][j];
						datalpi.jump1 = jump1;
						lpi_list.emplace_back(datalpi);
						ast.emplace_back(i);
						astf.emplace_back(j);
					}
				}

			}

		}
		assert(ast.size() == lpi_list.size());
		


		for (int i = 0; i < lpi_list.size(); i++)
		{
			if (face_flags[ast[i]][astf[i]]) continue;
			inter = Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision_return_id(lpi_list[i], triangle, filted_intersection,check_id);
			if (inter == 1) {
				dbgout2++;
				time4 += timer.getElapsedTimeInSec();
				return true;
			}
				
			if (inter == 0) {
				for (int h = 0; h < 3; h++) {
					if (lpi_list[i].segid == h) continue;
					tti = seg_cut_plane(triangle[triseg[h][0]], triangle[triseg[h][1]],
						halfspace[lpi_list[i].prismid][lpi_list[i].facetid][0], 
						halfspace[lpi_list[i].prismid][lpi_list[i].facetid][1], 
						halfspace[lpi_list[i].prismid][lpi_list[i].facetid][2]);
					if (tti != CUT_FACE) continue;
					int check = Implicit_Seg_Facet_interpoint_Out_Prism_check_id(
						triangle[triseg[h][0]], triangle[triseg[h][1]],
						halfspace[lpi_list[i].prismid][lpi_list[i].facetid][0],
						halfspace[lpi_list[i].prismid][lpi_list[i].facetid][1],
						halfspace[lpi_list[i].prismid][lpi_list[i].facetid][2], check_id);
					if (check == 0) face_flags[ast[i]][astf[i]] = true;
				}
			}
		} 


		std::vector<unsigned int> filted_intersection_new; filted_intersection_new.reserve(filted_intersection.size());
		std::vector<std::vector<int>>intersect_face_new; intersect_face_new.reserve(intersect_face.size());
		std::vector<int> tempvector,oldtonew;
		oldtonew.resize(intersect_face.size());
		int tot;
		for (int i = 0; i < filted_intersection.size(); i++) {
			tot = 0;
			for (int j = 0; j < intersect_face[i].size(); j++) {
				if (face_flags[i][j] == false) tot++;
			}
			if (tot > 0) {
				filted_intersection_new.emplace_back(filted_intersection[i]);
				tempvector.clear();
				for (int j = 0; j < intersect_face[i].size(); j++) {
					
					if (face_flags[i][j] == false) tempvector.emplace_back(intersect_face[i][j]);
				}
				intersect_face_new.emplace_back(tempvector);
				oldtonew[i] = intersect_face_new.size() - 1;
			}
			else {
				oldtonew[i] = -1;
			}
		}
		time4 += timer.getElapsedTimeInSec();
		if (filted_intersection_new.size() == 0) return false;//inside
		///////////////////////////////////////////////////tpp
		igl::Timer timerm,timer1;
		timerm.start();
		AABB localtree;
		std::vector<std::array<Vector3, 2>>localcorners;
		localcorners.resize(filted_intersection.size());

		for (int i = 0; i < filted_intersection.size(); i++) {
			localcorners[i] = cornerlist[filted_intersection[i]];//caution: this must not use the new list
		}
		
		localtree.init_envelope(localcorners);
		time5 += timerm.getElapsedTimeInSec();
		//std::wcout << "time for local tree " << timerm.getElapsedTimeInSec() << std::endl;

		

		tpilist.clear();
		tpilist.reserve(filted_intersection_new.size()/50);
		
		int id, id0 = 0;
		std::vector<unsigned int> localist;
		std::vector<unsigned int> potential_prism;
		std::vector<std::vector<int>> potential_facets;
		timerm.start();
		
		timer.start();
		for (int i = 0; i < filted_intersection_new.size(); i++)
		{
			jump1 = filted_intersection_new[i];
			timer1.start();
			localtree.bbd_finding_in_envelope(cornerlist[filted_intersection_new[i]][0], cornerlist[filted_intersection_new[i]][1], localist);
			time13 += timer1.getElapsedTimeInSec();
			potential_prism.clear();
			potential_facets.clear();
			potential_prism.resize(localist.size());
			potential_facets.resize(localist.size());
			for (int k = 0; k < localist.size(); k++) {
				potential_prism[k] = filted_intersection[localist[k]];
				potential_facets[k] = intersect_face[localist[k]];
			}
			for (int k = 0; k < intersect_face_new[i].size(); k++) {
				for (int j = 0; j < localist.size(); j++) {
					if (oldtonew[localist[j]] == -1) continue;
					jump2= filted_intersection[localist[j]];
					if (jump1 >= jump2) continue;
					for (int h = 0; h < intersect_face_new[oldtonew[localist[j]]].size(); h++) {
						Scalar d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7;
						bool multiflag;
						
						bool pre = orient3D_TPI_prefilter(triangle[0][0], triangle[0][1], triangle[0][2],
							triangle[1][0], triangle[1][1], triangle[1][2], triangle[2][0], triangle[2][1], triangle[2][2],
							halfspace[jump1][intersect_face_new[i][k]][0][0],
							halfspace[jump1][intersect_face_new[i][k]][0][1],
							halfspace[jump1][intersect_face_new[i][k]][0][2],

							halfspace[jump1][intersect_face_new[i][k]][1][0],
							halfspace[jump1][intersect_face_new[i][k]][1][1],
							halfspace[jump1][intersect_face_new[i][k]][1][2],

							halfspace[jump1][intersect_face_new[i][k]][2][0],
							halfspace[jump1][intersect_face_new[i][k]][2][1],
							halfspace[jump1][intersect_face_new[i][k]][2][2],

							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][0][0],
							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][0][1],
							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][0][2],

							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][1][0],
							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][1][1],
							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][1][2],
											 
							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][2][0],
							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][2][1],
							halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][2][2],

							d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7);

						if (pre)
						{
							TPI_exact_suppvars s;
							
							cut = is_3_triangle_cut_double(d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7, triangle,
								halfspace[jump1][intersect_face_new[i][k]][0],
								halfspace[jump1][intersect_face_new[i][k]][1],
								halfspace[jump1][intersect_face_new[i][k]][2],
								
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][0],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][1],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][2],
								multiflag, s);
							
							if (!cut)
								continue;
							timer1.start();
							cut = is_tpp_on_polyhedra_double(triangle,
								halfspace[jump1][intersect_face_new[i][k]][0],
								halfspace[jump1][intersect_face_new[i][k]][1],
								halfspace[jump1][intersect_face_new[i][k]][2],

								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][0],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][1],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][2], jump1, intersect_face_new[i][k]);
							time11 += timer1.getElapsedTimeInSec();
							if (!cut) continue;
							timer1.start();
							cut = is_tpp_on_polyhedra_double(triangle,
								halfspace[jump1][intersect_face_new[i][k]][0],
								halfspace[jump1][intersect_face_new[i][k]][1],
								halfspace[jump1][intersect_face_new[i][k]][2],

								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][0],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][1],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][2], jump2, intersect_face_new[oldtonew[localist[j]]][h]);
							time11 += timer1.getElapsedTimeInSec();
							if (!cut) continue;

					

							dbg1++;
							inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double_with_face_order( //TODO takes most of time
								d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7, triangle,
								halfspace[jump1][intersect_face_new[i][k]][0],
								halfspace[jump1][intersect_face_new[i][k]][1],
								halfspace[jump1][intersect_face_new[i][k]][2],

								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][0],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][1],
								halfspace[jump2][intersect_face_new[oldtonew[localist[j]]][h]][2],
								potential_prism,potential_facets, jump1, jump2, multiflag, s);
							
							if (inter == 1)
							{
								dbgout3++;
								time8 += timer.getElapsedTimeInSec();
								time6 += timerm.getElapsedTimeInSec();
								return true;
							}
							
						}
						else
						{
							
							tpilist.emplace_back();
							auto &datatpi = tpilist.back();
							datatpi.prismid1 = jump1;
							datatpi.facetid1 = intersect_face_new[i][k];
							datatpi.prismid2 = jump2;
							datatpi.facetid2 = intersect_face_new[oldtonew[localist[j]]][h];
							datatpi.jump1 = jump1;
							datatpi.jump2 = jump2;
							
						}
					}
				}
			}
		}
		time8 += timer.getElapsedTimeInSec();
		igl::Timer timer2;
		timer.start();
		for (int i = 0; i < tpilist.size(); i++)
		{
			TPI_exact_suppvars s;
			bool premulti = orient3D_TPI_pre_exact(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],

				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][0][0],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][0][1],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][0][2],

				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][1][0],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][1][1],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][1][2],

				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][2][0],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][2][1],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][2][2],

				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][0][0],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][0][1],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][0][2],

				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][1][0],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][1][1],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][1][2],

				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][2][0],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][2][1],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][2][2],
				s);
			if (premulti == false) continue;
			cut = is_3_triangle_cut_pure_multiprecision(triangle, s);
			if (cut == false) continue;
			timer1.start();
			cut = is_tpp_on_polyhedra_double(triangle,
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][0],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][1],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][2],

				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][0],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][1],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][2], tpilist[i].prismid1, tpilist[i].facetid1);
			time11 += timer1.getElapsedTimeInSec();
			if (cut == false) continue;
			timer1.start();
			cut = is_tpp_on_polyhedra_double(triangle,
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][0],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][1],
				halfspace[tpilist[i].prismid1][tpilist[i].facetid1][2],

				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][0],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][1],
				halfspace[tpilist[i].prismid2][tpilist[i].facetid2][2], tpilist[i].prismid2, tpilist[i].facetid2);
			time11 += timer1.getElapsedTimeInSec();
			if (cut == false) continue;
			timer2.start();
			localtree.bbd_finding_in_envelope(cornerlist[tpilist[i].prismid1][0], cornerlist[tpilist[i].prismid1][1], localist);
			time13 += timer2.getElapsedTimeInSec();
			potential_prism.clear();
			potential_prism.resize(localist.size());
			for (int k = 0; k < localist.size(); k++) {
				potential_prism[k] = filted_intersection[localist[k]];
			}
			
			inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(tpilist[i], triangle, potential_prism,s); //is_3_intersection is already in it
			if (inter == 1) {
				dbgout4++;
				time9 += timer.getElapsedTimeInSec();
				time6 += timerm.getElapsedTimeInSec();
				return true;
			}
				
		}
		time9 += timer.getElapsedTimeInSec();
		time6 += timerm.getElapsedTimeInSec();
		std::cout << "how many times of tpp calculation " << dbg1 << std::endl;
		/*std::cout << "time for tpp " << timerm.getElapsedTimeInSec() << std::endl;
		std::cout << "how many prisms intersection " << dbg3 << std::endl;
		
		std::cout << "how many times of directly multi tpp calculation " << dbg2 << std::endl;
		*/

		//Eigen::MatrixXd V(ver_new.size(), 3);
		//for (int i = 0; i < ver_new.size(); ++i)
		//	V.row(i) = ver_new[i];
		//Eigen::MatrixXi F(tempid.size(), 3);

		//for (int i = 0; i < tempid.size(); ++i)
		//	F.row(i) = faces_new[tempid[i]];

		////igl::write_triangle_mesh("patch.stl", V, F);
		//Eigen::MatrixXd V1(3, 3);
		//Eigen::MatrixXi F1(1, 3);
		//V1(0,0) = triangle[0][0];
		//V1(0, 1) = triangle[0][1];
		//V1(0, 2) = triangle[0][2];
		//V1(1, 0) = triangle[1][0];
		//V1(1, 1) = triangle[1][1];
		//V1(1, 2) = triangle[1][2];
		//V1(2, 0) = triangle[2][0];
		//V1(2, 1) = triangle[2][1];
		//V1(2, 2) = triangle[2][2];

		//
		//F1 << 0, 1, 2;
		//igl::write_triangle_mesh("tri.stl", V1, F1);

		return false;
	}


	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const
	{
		int tot, ori;
		
		// timer.start();
		LPI_exact_suppvars s;
		bool premulti = orient3D_LPI_pre_exact(
			triangle[triseg[datalpi.segid][0]][0], 
			triangle[triseg[datalpi.segid][0]][1],
			triangle[triseg[datalpi.segid][0]][2],

			triangle[triseg[datalpi.segid][1]][0],
			triangle[triseg[datalpi.segid][1]][1],
			triangle[triseg[datalpi.segid][1]][2],

			halfspace[datalpi.prismid][datalpi.facetid][0][0],
			halfspace[datalpi.prismid][datalpi.facetid][0][1],
			halfspace[datalpi.prismid][datalpi.facetid][0][2],

			halfspace[datalpi.prismid][datalpi.facetid][1][0],
			halfspace[datalpi.prismid][datalpi.facetid][1][1],
			halfspace[datalpi.prismid][datalpi.facetid][1][2],

			halfspace[datalpi.prismid][datalpi.facetid][2][0],
			halfspace[datalpi.prismid][datalpi.facetid][2][1],
			halfspace[datalpi.prismid][datalpi.facetid][2][2],
			s);

		// time_multi += timer.getElapsedTimeInSec();

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == datalpi.jump1)
			{
				continue;
			}
			tot = 0;
			
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
				
				// timer.start();
				ori = orient3D_LPI_post_exact(s,
					triangle[triseg[datalpi.segid][0]][0],
					triangle[triseg[datalpi.segid][0]][1],
					triangle[triseg[datalpi.segid][0]][2],

					halfspace[prismindex[i]][j][0][0],
					halfspace[prismindex[i]][j][0][1],
					halfspace[prismindex[i]][j][0][2],

					halfspace[prismindex[i]][j][1][0],
					halfspace[prismindex[i]][j][1][1],
					halfspace[prismindex[i]][j][1][2],

					halfspace[prismindex[i]][j][2][0],
					halfspace[prismindex[i]][j][2][1],
					halfspace[prismindex[i]][j][2][2]);
				if (ori == 1 || ori == 0)
				{
					break;
				}

				if (ori == -1)
				{
					tot++;
				}
			}
			if (tot == prismindex.size())
			{
				return IN_PRISM;
			}
		}
			
			
		return OUT_PRISM;
	}
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision_return_id(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle, 
		const std::vector<unsigned int> &prismindex,int& id) const
	{
		int tot, ori;

		// timer.start();
		LPI_exact_suppvars s;
		bool premulti = orient3D_LPI_pre_exact(
			triangle[triseg[datalpi.segid][0]][0],
			triangle[triseg[datalpi.segid][0]][1],
			triangle[triseg[datalpi.segid][0]][2],

			triangle[triseg[datalpi.segid][1]][0],
			triangle[triseg[datalpi.segid][1]][1],
			triangle[triseg[datalpi.segid][1]][2],

			halfspace[datalpi.prismid][datalpi.facetid][0][0],
			halfspace[datalpi.prismid][datalpi.facetid][0][1],
			halfspace[datalpi.prismid][datalpi.facetid][0][2],

			halfspace[datalpi.prismid][datalpi.facetid][1][0],
			halfspace[datalpi.prismid][datalpi.facetid][1][1],
			halfspace[datalpi.prismid][datalpi.facetid][1][2],

			halfspace[datalpi.prismid][datalpi.facetid][2][0],
			halfspace[datalpi.prismid][datalpi.facetid][2][1],
			halfspace[datalpi.prismid][datalpi.facetid][2][2],
			s);

		// time_multi += timer.getElapsedTimeInSec();

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == datalpi.jump1)
			{
				continue;
			}
			tot = 0;

			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

				// timer.start();
				ori = orient3D_LPI_post_exact(s,
					triangle[triseg[datalpi.segid][0]][0],
					triangle[triseg[datalpi.segid][0]][1],
					triangle[triseg[datalpi.segid][0]][2],

					halfspace[prismindex[i]][j][0][0],
					halfspace[prismindex[i]][j][0][1],
					halfspace[prismindex[i]][j][0][2],

					halfspace[prismindex[i]][j][1][0],
					halfspace[prismindex[i]][j][1][1],
					halfspace[prismindex[i]][j][1][2],

					halfspace[prismindex[i]][j][2][0],
					halfspace[prismindex[i]][j][2][1],
					halfspace[prismindex[i]][j][2][2]);
				if (ori == 1 || ori == 0)
				{
					break;
				}

				if (ori == -1)
				{
					tot++;
				}
			}
			if (tot == prismindex.size())
			{
				id = prismindex[i];
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
		
		int tot;
		int ori, ori1;
		static INDEX index;
		static std::vector<INDEX> recompute;

		recompute.clear();
		int face[3];
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump) continue;

			index.FACES.clear();
			tot = 0;
			
			
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
				
				ori =
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						segpoint0[0], segpoint0[1], segpoint0[2],
						halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
						halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
						halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);

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
			if (tot == halfspace[prismindex[i]].size())
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
				
			// timer.start();
			LPI_exact_suppvars s;
			bool premulti = orient3D_LPI_pre_exact(
				segpoint0[0],segpoint0[1],segpoint0[2],
				segpoint1[0],segpoint1[1],segpoint1[2],
				triangle0[0],triangle0[1],triangle0[2],
				triangle1[0],triangle1[1],triangle1[2],
				triangle2[0],triangle2[1],triangle2[2],
				s);
			// time_multi += timer.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = recompute[k].Pi;
				
				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];
			
					ori = orient3D_LPI_post_exact(s, segpoint0[0], segpoint0[1], segpoint0[2],
						halfspace[in1][in2][0][0], halfspace[in1][in2][0][1], halfspace[in1][in2][0][2],
						halfspace[in1][in2][1][0], halfspace[in1][in2][1][1], halfspace[in1][in2][1][2],
						halfspace[in1][in2][2][0], halfspace[in1][in2][2][1], halfspace[in1][in2][2][2]);
					// time_multi += timer.getElapsedTimeInSec();
					if (ori == 1 || ori == 0)
						break;
					
				}
				if (ori == -1) return IN_PRISM;

			}
		}

		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_double_return_id(
		const Scalar &a11, const Scalar &a12, const Scalar &a13, const Scalar &d, const Scalar &fa11,
		const Scalar &fa12, const Scalar &fa13, const Scalar &max1, const Scalar &max2, const Scalar &max5,
		const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
		const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const
	{

		int tot;
		int ori, ori1;
		static INDEX index;
		static std::vector<INDEX> recompute;

		recompute.clear();
		int face[3];
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump) continue;

			index.FACES.clear();
			tot = 0;


			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

				ori =
					orient3D_LPI_postfilter(
						a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
						segpoint0[0], segpoint0[1], segpoint0[2],
						halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
						halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
						halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);

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
			if (tot == halfspace[prismindex[i]].size())
			{
				id = prismindex[i];
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

			// timer.start();
			LPI_exact_suppvars s;
			bool premulti = orient3D_LPI_pre_exact(
				segpoint0[0], segpoint0[1], segpoint0[2],
				segpoint1[0], segpoint1[1], segpoint1[2],
				triangle0[0], triangle0[1], triangle0[2],
				triangle1[0], triangle1[1], triangle1[2],
				triangle2[0], triangle2[1], triangle2[2],
				s);
			// time_multi += timer.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = recompute[k].Pi;

				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];

					ori = orient3D_LPI_post_exact(s, segpoint0[0], segpoint0[1], segpoint0[2],
						halfspace[in1][in2][0][0], halfspace[in1][in2][0][1], halfspace[in1][in2][0][2],
						halfspace[in1][in2][1][0], halfspace[in1][in2][1][1], halfspace[in1][in2][1][2],
						halfspace[in1][in2][2][0], halfspace[in1][in2][2][1], halfspace[in1][in2][2][2]);
					// time_multi += timer.getElapsedTimeInSec();
					if (ori == 1 || ori == 0)
						break;

				}
				
				if (ori == -1) {
					id = in1;
					return IN_PRISM;
				}
			}
		}

		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_check_id(
		const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
		const Vector3 &triangle1, const Vector3 &triangle2, const int &id) const
	{

		Scalar a11, a12, a13, d, fa11,
			fa12, fa13, max1, max2, max5;
		int tot;
		int ori;
		static INDEX index;
		static std::vector<INDEX> recompute;

		recompute.clear();

		index.FACES.clear();
		tot = 0;

		bool precom = orient3D_LPI_prefilter( //
			segpoint0[0], segpoint0[1], segpoint0[2],
			segpoint1[0], segpoint1[1], segpoint1[2],
			triangle0[0], triangle0[1], triangle0[2],
			triangle1[0], triangle1[1], triangle1[2],
			triangle2[0], triangle2[1], triangle2[2],
			a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5);
		if (precom == false) {
			LPI_exact_suppvars s;
			bool premulti = orient3D_LPI_pre_exact(
				segpoint0[0], segpoint0[1], segpoint0[2],
				segpoint1[0], segpoint1[1], segpoint1[2],
				triangle0[0], triangle0[1], triangle0[2],
				triangle1[0], triangle1[1], triangle1[2],
				triangle2[0], triangle2[1], triangle2[2],
				s);

			for (int j = 0; j < halfspace[id].size(); j++) {
				ori = orient3D_LPI_post_exact(s, segpoint0[0], segpoint0[1], segpoint0[2],
					halfspace[id][j][0][0], halfspace[id][j][0][1], halfspace[id][j][0][2],
					halfspace[id][j][1][0], halfspace[id][j][1][1], halfspace[id][j][1][2],
					halfspace[id][j][2][0], halfspace[id][j][2][1], halfspace[id][j][2][2]);
				if (ori == 1 || ori == 0) {
					break;
				}

			}
			if (ori == -1) {

				return IN_PRISM;
			}
			else {
				return OUT_PRISM;
			}

		}

		for (int j = 0; j < halfspace[id].size(); j++) {

			ori =
				orient3D_LPI_postfilter(
					a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
					segpoint0[0], segpoint0[1], segpoint0[2],
					halfspace[id][j][0][0], halfspace[id][j][0][1], halfspace[id][j][0][2],
					halfspace[id][j][1][0], halfspace[id][j][1][1], halfspace[id][j][1][2],
					halfspace[id][j][2][0], halfspace[id][j][2][1], halfspace[id][j][2][2]);

			if (ori == 1)
			{
				return OUT_PRISM;
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
		if (tot == halfspace[id].size())
		{

			return IN_PRISM;
		}


		assert(!index.FACES.empty());
		index.Pi = id;
		recompute.emplace_back(index);

		if (!recompute.empty())
		{

			// timer.start();
			LPI_exact_suppvars s;
			bool premulti = orient3D_LPI_pre_exact(
				segpoint0[0], segpoint0[1], segpoint0[2],
				segpoint1[0], segpoint1[1], segpoint1[2],
				triangle0[0], triangle0[1], triangle0[2],
				triangle1[0], triangle1[1], triangle1[2],
				triangle2[0], triangle2[1], triangle2[2],
				s);
			// time_multi += timer.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = recompute[k].Pi;

				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];

					ori = orient3D_LPI_post_exact(s, segpoint0[0], segpoint0[1], segpoint0[2],
						halfspace[in1][in2][0][0], halfspace[in1][in2][0][1], halfspace[in1][in2][0][2],
						halfspace[in1][in2][1][0], halfspace[in1][in2][1][1], halfspace[in1][in2][1][2],
						halfspace[in1][in2][2][0], halfspace[in1][in2][2][1], halfspace[in1][in2][2][2]);
					// time_multi += timer.getElapsedTimeInSec();
					if (ori == 1 || ori == 0)
						break;

				}
	
				if (ori == -1) return IN_PRISM;

			}
		}

		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision(const DATA_TPI &datatpi, const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex, TPI_exact_suppvars& s) const
	{
		igl::Timer timer;
		int tot, ori;
		int jump1 = datatpi.jump1, jump2 = datatpi.jump2;

		
		
		dbg1++;
		dbg2++;
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump1|| prismindex[i] == jump2)	continue;
			tot = 0;
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
			
				timer.start();
				ori = orient3D_TPI_post_exact(s, 
					halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
					halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
					halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);
				time10 += timer.getElapsedTimeInSec();
				ct2++;
				if (ori == 1 || ori == 0)
				{
					break;
				}

				if (ori == -1)
				{
					tot++;
				}
			}
			if (tot == halfspace[prismindex[i]].size())
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
		
		int ori;
		int tot;
		igl::Timer timer;
		static INDEX index;
		static std::vector<INDEX> recompute;
		recompute.clear();

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump1 || prismindex[i] == jump2)
				continue;
			
			index.FACES.clear();
			tot = 0;
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
				timer.start();
				ct1 += 1;
				ori =
					orient3D_TPI_postfilter(
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
						halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
						halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
						halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);
				// timetpp1 += timer_a.getElapsedTimeInSec();
				time7 += timer.getElapsedTimeInSec();
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
			if (tot == halfspace[prismindex[i]].size())
			{

				return IN_PRISM;
			}

			if (ori != 1)
			{
				index.Pi = prismindex[i];
				recompute.emplace_back(index);

			}
		}
		timer.start();
		if (recompute.size() > 0)
		{

			if (!multiflag)
			{
				// timer.start();
				bool premulti = orient3D_TPI_pre_exact(
					triangle[0][0],triangle[0][1],triangle[0][2],
					triangle[1][0],triangle[1][1],triangle[1][2],
					triangle[2][0],triangle[2][1],triangle[2][2],
								  
					facet10[0], facet10[1], facet10[2],
					facet11[0], facet11[1], facet11[2],
					facet12[0], facet12[1], facet12[2],

					facet20[0], facet20[1], facet20[2],
					facet21[0], facet21[1], facet21[2],
					facet22[0], facet22[1], facet22[2],
					s);
				// time_multi += timer.getElapsedTimeInSec();
			}

			// timetpp2 += timer_a.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = recompute[k].Pi;

				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];
					
					ori = orient3D_TPI_post_exact(s,
						halfspace[in1][in2][0][0], halfspace[in1][in2][0][1], halfspace[in1][in2][0][2],
						halfspace[in1][in2][1][0], halfspace[in1][in2][1][1], halfspace[in1][in2][1][2],
						halfspace[in1][in2][2][0], halfspace[in1][in2][2][1], halfspace[in1][in2][2][2]
						);

					if (ori == 1 || ori == 0)	break;
				}
				if (ori == -1) {
					time12+=timer.getElapsedTimeInSec();
					return IN_PRISM;
				}
			}
		}
		time12 += timer.getElapsedTimeInSec();
		return OUT_PRISM;
	}

	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_double_with_face_order(
		const Scalar &d, const Scalar &n1, const Scalar &n2, const Scalar &n3,
		const Scalar &max1, const Scalar &max2, const Scalar &max3, const Scalar &max4, const Scalar &max5, const Scalar &max6, const Scalar &max7,
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
		const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>>intersect_face, const int &jump1, const int &jump2, const bool &multiflag,
		TPI_exact_suppvars &s) const
	{

		int ori;
		int tot;
		igl::Timer timer;
		static INDEX index;
		static std::vector<INDEX> recompute;
		recompute.clear();
		bool flag = false;

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump1 || prismindex[i] == jump2)
				continue;
			
			index.FACES.clear();
			tot = 0;
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
				flag = false;
				for (int k = 0; k < intersect_face[i].size(); k++) {
					if (j == intersect_face[i][k]) {
						flag = true;
						break;
					}
				}
				if (flag == false) continue;
				timer.start();
				ct1 += 1;
				ori =
					orient3D_TPI_postfilter(
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
						halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
						halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
						halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);
				// timetpp1 += timer_a.getElapsedTimeInSec();
				time7 += timer.getElapsedTimeInSec();
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
			if (ori == 1) continue;
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
				flag = false;
				for (int k = 0; k < intersect_face[i].size(); k++) {
					if (j == intersect_face[i][k]) {
						flag = true;
						break;
					}
				}
				if (flag == true) continue;
				timer.start();
				ct1 += 1;
				ori =
					orient3D_TPI_postfilter(
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
						halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
						halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
						halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);
				// timetpp1 += timer_a.getElapsedTimeInSec();
				time7 += timer.getElapsedTimeInSec();
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
			if (ori == 1) continue;
			if (tot == halfspace[prismindex[i]].size())
			{

				return IN_PRISM;
			}

			if (ori != 1)
			{
				index.Pi = prismindex[i];
				recompute.emplace_back(index);

			}
		}
		timer.start();
		if (recompute.size() > 0)
		{

			if (!multiflag)
			{
				// timer.start();
				bool premulti = orient3D_TPI_pre_exact(
					triangle[0][0], triangle[0][1], triangle[0][2],
					triangle[1][0], triangle[1][1], triangle[1][2],
					triangle[2][0], triangle[2][1], triangle[2][2],

					facet10[0], facet10[1], facet10[2],
					facet11[0], facet11[1], facet11[2],
					facet12[0], facet12[1], facet12[2],

					facet20[0], facet20[1], facet20[2],
					facet21[0], facet21[1], facet21[2],
					facet22[0], facet22[1], facet22[2],
					s);
				// time_multi += timer.getElapsedTimeInSec();
			}

			// timetpp2 += timer_a.getElapsedTimeInSec();

			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = recompute[k].Pi;

				for (int j = 0; j < recompute[k].FACES.size(); j++) {
					int in2 = recompute[k].FACES[j];

					ori = orient3D_TPI_post_exact(s,
						halfspace[in1][in2][0][0], halfspace[in1][in2][0][1], halfspace[in1][in2][0][2],
						halfspace[in1][in2][1][0], halfspace[in1][in2][1][1], halfspace[in1][in2][1][2],
						halfspace[in1][in2][2][0], halfspace[in1][in2][2][1], halfspace[in1][in2][2][2]
					);

					if (ori == 1 || ori == 0)	break;
				}
				if (ori == -1) {
					time12 += timer.getElapsedTimeInSec();
					return IN_PRISM;
				}
			}
		}
		time12 += timer.getElapsedTimeInSec();
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




	int FastEnvelope::is_triangle_cut_envelope_polyhedra(const int &cindex,
		const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, std::vector<int> &cid) const
	{
		
		cid.clear();
		std::vector<bool> cut;
		cut.resize(halfspace[cindex].size());
		for (int i = 0; i < halfspace[cindex].size(); i++)
		{
			cut[i] = false;
		}
		std::vector<int> o1, o2, o3, cutp;
		o1.resize(halfspace[cindex].size());
		o2.resize(halfspace[cindex].size());
		o3.resize(halfspace[cindex].size());
		int  ori = 0, ct = 0;
		

		for (int i = 0; i < halfspace[cindex].size(); i++)
		{
			
			o1[i] = Predicates::orient_3d(tri0, halfspace[cindex][i][0], halfspace[cindex][i][1], halfspace[cindex][i][2]);
			o2[i] = Predicates::orient_3d(tri1, halfspace[cindex][i][0], halfspace[cindex][i][1], halfspace[cindex][i][2]);
			o3[i] = Predicates::orient_3d(tri2, halfspace[cindex][i][0], halfspace[cindex][i][1], halfspace[cindex][i][2]);
			if (o1[i] + o2[i] + o3[i] >= 3)
			{
				return 0;
			}
			if (o1[i] == -1) ct++;
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1)
			{
				return 0;
			}
			if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0)
			{
				return 0;
			}
			if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0)
			{
				return 0;
			}
			if (o1[i] == 0 && o2[i] == 0 && o3[i] == 0)
			{
				return 0;
			}

			if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1)
				cutp.emplace_back(i);
		}
		if (cutp.size() == 0 && ct == halfspace[cindex].size()) return 2;//means totally inside of this prism
		if (cutp.size() == 0)
		{
			return 0;
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
			
			for (int k = 0; k < 2; k++)
			{
				
				bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
					seg00[k], seg01[k], seg02[k],
					seg10[k], seg11[k], seg12[k],
					halfspace[cindex][cutp[i]][0][0], halfspace[cindex][cutp[i]][0][1], halfspace[cindex][cutp[i]][0][2],
					halfspace[cindex][cutp[i]][1][0], halfspace[cindex][cutp[i]][1][1], halfspace[cindex][cutp[i]][1][2],
					halfspace[cindex][cutp[i]][2][0], halfspace[cindex][cutp[i]][2][1], halfspace[cindex][cutp[i]][2][2],
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
					ori =
						orient3D_LPI_postfilter(
							a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5,
							seg00[k], seg01[k], seg02[k],
							halfspace[cindex][cutp[j]][0][0], halfspace[cindex][cutp[j]][0][1], halfspace[cindex][cutp[j]][0][2],
							halfspace[cindex][cutp[j]][1][0], halfspace[cindex][cutp[j]][1][1], halfspace[cindex][cutp[j]][1][2],
							halfspace[cindex][cutp[j]][2][0], halfspace[cindex][cutp[j]][2][1], halfspace[cindex][cutp[j]][2][2]);

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
			for (int i = 0; i < halfspace[cindex].size(); i++)
			{
				if (cut[i] == true)
					cid.emplace_back(i);
			}
			return 1;
		}
		// triangle-facet-facet intersection
		Scalar n1, n2, n3, max3, max4, max6, max7;

		for (int i = 0; i < cutp.size(); i++)
		{
			for (int j = i + 1; j < cutp.size(); j++)
			{
				if (cut[cutp[i]] == true && cut[cutp[j]] == true)
					continue;

				bool neib = is_two_facets_neighbouring(cindex, cutp[i], cutp[j]);
				if (neib == false) continue;
				

				int inter = is_3_triangle_cut_float_fast(
					tri0, tri1, tri2,
					halfspace[cindex][cutp[i]][0],
					halfspace[cindex][cutp[i]][1],
					halfspace[cindex][cutp[i]][2],
					halfspace[cindex][cutp[j]][0],
					halfspace[cindex][cutp[j]][1],
					halfspace[cindex][cutp[j]][2]);
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
						halfspace[cindex][cutp[i]][0][0], halfspace[cindex][cutp[i]][0][1], halfspace[cindex][cutp[i]][0][2],
						halfspace[cindex][cutp[i]][1][0], halfspace[cindex][cutp[i]][1][1], halfspace[cindex][cutp[i]][1][2],
						halfspace[cindex][cutp[i]][2][0], halfspace[cindex][cutp[i]][2][1], halfspace[cindex][cutp[i]][2][2],
						halfspace[cindex][cutp[j]][0][0], halfspace[cindex][cutp[j]][0][1], halfspace[cindex][cutp[j]][0][2],
						halfspace[cindex][cutp[j]][1][0], halfspace[cindex][cutp[j]][1][1], halfspace[cindex][cutp[j]][1][2],
						halfspace[cindex][cutp[j]][2][0], halfspace[cindex][cutp[j]][2][1], halfspace[cindex][cutp[j]][2][2],
						d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);

				for (int k = 0; k < cutp.size(); k++)
				{

					if (k == i || k == j)
						continue;
					
					ori =
						orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
							halfspace[cindex][cutp[k]][0][0], halfspace[cindex][cutp[k]][0][1], halfspace[cindex][cutp[k]][0][2],
							halfspace[cindex][cutp[k]][1][0], halfspace[cindex][cutp[k]][1][1], halfspace[cindex][cutp[k]][1][2],
							halfspace[cindex][cutp[k]][2][0], halfspace[cindex][cutp[k]][2][1], halfspace[cindex][cutp[k]][2][2]);

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

		for (int i = 0; i < halfspace[cindex].size(); i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return 1;
	}
	bool FastEnvelope::is_tpp_on_polyhedra_double(
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
		const int &prismid, const int &faceid)const {
		int ori;
		Scalar  d, n1, n2, n3,
			max1, max2, max3, max4, max5,
			max6, max7;
		TPI_exact_suppvars s;
		bool premulti = false;
		bool pre =
			orient3D_TPI_prefilter(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],
				facet10[0], facet10[1], facet10[2],
				facet11[0], facet11[1], facet11[2],
				facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2],
				facet21[0], facet21[1], facet21[2],
				facet22[0], facet22[1], facet22[2],
				d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);
		if (pre == true) {
			for (int i = 0; i < halfspace[prismid].size(); i++) {
				/*bool neib = is_two_facets_neighbouring(prismid, i, faceid);
				if (neib == false) continue;*/
				if (i == faceid) continue;
				ori =
					orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
						halfspace[prismid][i][0][0], halfspace[prismid][i][0][1], halfspace[prismid][i][0][2],
						halfspace[prismid][i][1][0], halfspace[prismid][i][1][1], halfspace[prismid][i][1][2],
						halfspace[prismid][i][2][0], halfspace[prismid][i][2][1], halfspace[prismid][i][2][2]);
				
				if (ori == 0) {
					if (premulti == false) {
						premulti = orient3D_TPI_pre_exact(
							triangle[0][0], triangle[0][1], triangle[0][2],
							triangle[1][0], triangle[1][1], triangle[1][2],
							triangle[2][0], triangle[2][1], triangle[2][2],

							facet10[0], facet10[1], facet10[2],
							facet11[0], facet11[1], facet11[2],
							facet12[0], facet12[1], facet12[2],

							facet20[0], facet20[1], facet20[2],
							facet21[0], facet21[1], facet21[2],
							facet22[0], facet22[1], facet22[2],
							s);
						if (premulti == false) return false;
					}
					ori = orient3D_TPI_post_exact(s,
						halfspace[prismid][i][0][0], halfspace[prismid][i][0][1], halfspace[prismid][i][0][2],
						halfspace[prismid][i][1][0], halfspace[prismid][i][1][1], halfspace[prismid][i][1][2],
						halfspace[prismid][i][2][0], halfspace[prismid][i][2][1], halfspace[prismid][i][2][2]);
				}
				//if (ori == 0) std::cout << "f id " << halfspace[prismid].size() <<" "<< faceid << " " << i << std::endl;
				if (ori == 1) return false;

			}
			
		}
		else {
			premulti = orient3D_TPI_pre_exact(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],

				facet10[0], facet10[1], facet10[2],
				facet11[0], facet11[1], facet11[2],
				facet12[0], facet12[1], facet12[2],

				facet20[0], facet20[1], facet20[2],
				facet21[0], facet21[1], facet21[2],
				facet22[0], facet22[1], facet22[2],
				s);
			if (premulti == false) return false;
			for (int i = 0; i < halfspace[prismid].size(); i++) {
				/*bool neib = is_two_facets_neighbouring(prismid, i, faceid);
				if (neib == false) continue;*/
				if (i == faceid) continue;
				ori = orient3D_TPI_post_exact(s,
					halfspace[prismid][i][0][0], halfspace[prismid][i][0][1], halfspace[prismid][i][0][2],
					halfspace[prismid][i][1][0], halfspace[prismid][i][1][1], halfspace[prismid][i][1][2],
					halfspace[prismid][i][2][0], halfspace[prismid][i][2][1], halfspace[prismid][i][2][2]);
				//if (ori == 0) std::cout << "f id " << faceid << " " << i << std::endl;
				if (ori == 1) return false;

			}
		}
		return true;
	}
	bool FastEnvelope::is_tpp_inside_bbd(
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
		const Vector3& min, const Vector3 &max)const {
		int ori;
		Scalar  d, n1, n2, n3,
			max1, max2, max3, max4, max5,
			max6, max7;
		TPI_exact_suppvars s;
		std::array<Vector3, 8> ver;
		ver[0] = Vector3(min[0], min[1], min[2]);
		ver[1] = Vector3(max[0], min[1], min[2]);
		ver[2] = Vector3(max[0], max[1], min[2]);
		ver[3] = Vector3(min[0], max[1], min[2]);
		ver[4] = Vector3(min[0], min[1], max[2]);
		ver[5] = Vector3(max[0], min[1], max[2]);
		ver[6] = Vector3(max[0], max[1], max[2]);
		ver[7] = Vector3(min[0], max[1], max[2]);

		std::array<std::array<Vector3, 3>, 6> box;
		box[0][0] = ver[0]; box[0][1] = ver[3]; box[0][2] = ver[2];
		box[1][0] = ver[4]; box[1][1] = ver[5]; box[1][2] = ver[6];
		box[2][0] = ver[5]; box[2][1] = ver[4]; box[2][2] = ver[0];
		box[3][0] = ver[5]; box[3][1] = ver[1]; box[3][2] = ver[6];
		box[4][0] = ver[7]; box[4][1] = ver[6]; box[4][2] = ver[2];
		box[5][0] = ver[4]; box[5][1] = ver[7]; box[5][2] = ver[3];
		bool premulti = false;
		bool pre =
			orient3D_TPI_prefilter(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],
				facet10[0], facet10[1], facet10[2],
				facet11[0], facet11[1], facet11[2],
				facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2],
				facet21[0], facet21[1], facet21[2],
				facet22[0], facet22[1], facet22[2],
				d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7);
		if (pre == true) {
			for (int i = 0; i < 6; i++) {
				
				ori =
					orient3D_TPI_postfilter(d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7,
						box[i][0][0], box[i][0][1], box[i][0][2],
						box[i][1][0], box[i][1][1], box[i][1][2],
						box[i][2][0], box[i][2][1], box[i][2][2]);

				if (ori == 0) {
					if (premulti == false) {
						premulti = orient3D_TPI_pre_exact(
							triangle[0][0], triangle[0][1], triangle[0][2],
							triangle[1][0], triangle[1][1], triangle[1][2],
							triangle[2][0], triangle[2][1], triangle[2][2],

							facet10[0], facet10[1], facet10[2],
							facet11[0], facet11[1], facet11[2],
							facet12[0], facet12[1], facet12[2],

							facet20[0], facet20[1], facet20[2],
							facet21[0], facet21[1], facet21[2],
							facet22[0], facet22[1], facet22[2],
							s);
						if (premulti == false) return false;
					}
					ori = orient3D_TPI_post_exact(s,
						box[i][0][0], box[i][0][1], box[i][0][2],
						box[i][1][0], box[i][1][1], box[i][1][2],
						box[i][2][0], box[i][2][1], box[i][2][2]);
				}
				if (ori == 1 || ori == 0) return false;

			}

		}
		else {
			premulti = orient3D_TPI_pre_exact(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],

				facet10[0], facet10[1], facet10[2],
				facet11[0], facet11[1], facet11[2],
				facet12[0], facet12[1], facet12[2],

				facet20[0], facet20[1], facet20[2],
				facet21[0], facet21[1], facet21[2],
				facet22[0], facet22[1], facet22[2],
				s);
			if (premulti == false) return false;
			for (int i = 0; i < 6; i++) {

				ori = orient3D_TPI_post_exact(s,
					box[i][0][0], box[i][0][1], box[i][0][2],
					box[i][1][0], box[i][1][1], box[i][1][2],
					box[i][2][0], box[i][2][1], box[i][2][2]);

				if (ori == 1 || ori == 0) return false;

			}
		}
		return true;
	}
	bool FastEnvelope::is_seg_cut_polyhedra(const int &cindex,
		const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const
	{

		
		
		std::vector<bool> cut;
		cut.resize(halfspace[cindex].size());
		for (int i = 0; i < halfspace[cindex].size(); i++)
		{
			cut[i] = false;
		}
		std::vector<int> o1, o2;
		o1.resize(halfspace[cindex].size());
		o2.resize(halfspace[cindex].size());
		int ori = 0;
		std::vector<int> cutp;
	
		for (int i = 0; i < halfspace[cindex].size(); i++)
		{
			
			
			o1[i] = Predicates::orient_3d(seg0, halfspace[cindex][i][0], halfspace[cindex][i][1], halfspace[cindex][i][2]);
			o2[i] = Predicates::orient_3d(seg1, halfspace[cindex][i][0], halfspace[cindex][i][1], halfspace[cindex][i][2]);

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
				halfspace[cindex][cutp[i]][0][0], halfspace[cindex][cutp[i]][0][1], halfspace[cindex][cutp[i]][0][2],
				halfspace[cindex][cutp[i]][1][0], halfspace[cindex][cutp[i]][1][1], halfspace[cindex][cutp[i]][1][2],
				halfspace[cindex][cutp[i]][2][0], halfspace[cindex][cutp[i]][2][1], halfspace[cindex][cutp[i]][2][2],
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
						halfspace[cindex][cutp[j]][0][0], halfspace[cindex][cutp[j]][0][1], halfspace[cindex][cutp[j]][0][2],
						halfspace[cindex][cutp[j]][1][0], halfspace[cindex][cutp[j]][1][1], halfspace[cindex][cutp[j]][1][2],
						halfspace[cindex][cutp[j]][2][0], halfspace[cindex][cutp[j]][2][1], halfspace[cindex][cutp[j]][2][2]);

				if (ori == 1)
					break;
			}
			if (ori != 1)
			{
				cut[cutp[i]] = true;
			}
		}

		for (int i = 0; i < halfspace[cindex].size(); i++)
		{
			if (cut[i] == true)
				cid.emplace_back(i);
		}

		return true;
	}
	bool FastEnvelope::point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const
	{

		int ori;
		
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump)
				continue;
			
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++)
			{
				
				ori = Predicates::orient_3d(halfspace[prismindex[i]][j][0], halfspace[prismindex[i]][j][1], halfspace[prismindex[i]][j][2], point);
				if (ori == -1 || ori == 0)
				{
					break;
				}
				if (j == halfspace[prismindex[i]].size() - 1)
				{

					return false;
				}
			}
			
		}

		return true;
	}
	bool FastEnvelope::point_out_prism_return_id(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump,int &id) const
	{

		int ori;

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump)
				continue;

			for (int j = 0; j < halfspace[prismindex[i]].size(); j++)
			{

				ori = Predicates::orient_3d(halfspace[prismindex[i]][j][0], halfspace[prismindex[i]][j][1], halfspace[prismindex[i]][j][2], point);
				if (ori == -1 || ori == 0)
				{
					break;
				}
				if (j == halfspace[prismindex[i]].size() - 1)
				{
					id = prismindex[i];
					return false;
				}
			}

		}

		return true;
	}
	bool FastEnvelope::point_out_prism_return_id_list(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump, std::vector<int> &idlist) const
	{
		idlist.clear();

		int ori;

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump)
				continue;

			for (int j = 0; j < halfspace[prismindex[i]].size(); j++)
			{

				ori = Predicates::orient_3d(halfspace[prismindex[i]][j][0], halfspace[prismindex[i]][j][1], halfspace[prismindex[i]][j][2], point);
				if (ori == -1 || ori == 0)
				{
					break;
				}
				if (j == halfspace[prismindex[i]].size() - 1)
				{
					idlist.emplace_back(prismindex[i]);
					
				}
			}

		}
		if (idlist.size() > 0) return false;
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

	void FastEnvelope::halfspace_init(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, 
		std::vector<std::vector<std::array<Vector3, 3>>>& halfspace, std::vector<std::array<Vector3, 2>>& cornerlist, const Scalar &epsilon) {
		
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
		const auto get_bb_corners_12 = [](const std::array<Vector3, 12> &vertices) {//TODO why use this one
			std::array<Vector3, 2> corners;
			corners[0] = vertices[0];
			corners[1] = vertices[0];

			for (size_t j = 0; j < 12; j++) {
				for (int i = 0; i < 3; i++) {
					corners[0][i] = std::min(corners[0][i], vertices[j][i]);
					corners[1][i] = std::max(corners[1][i], vertices[j][i]);
				}
			}

			//const Scalar dis = (max - min).minCoeff() * 3;//TODO  change to 1e-5 or sth
			const Scalar dis = 1e-4;
			for (int j = 0; j < 3; j++) {
				corners[0][j] -= dis;
				corners[1][j] += dis;
			}
			return corners;

		};
		const auto get_bb_corners_8 = [](const std::array<Vector3, 8> &vertices) {//TODO why use this one
			std::array<Vector3, 2> corners;
			corners[0] = vertices[0];
			corners[1] = vertices[0];

			for (size_t j = 0; j < 8; j++) {
				for (int i = 0; i < 3; i++) {
					corners[0][i] = std::min(corners[0][i], vertices[j][i]);
					corners[1][i] = std::max(corners[1][i], vertices[j][i]);
				}
			}

			//const Scalar dis = (max - min).minCoeff() * 3;//TODO  change to 1e-5 or sth
			const Scalar dis = 1e-4;
			for (int j = 0; j < 3; j++) {
				corners[0][j] -= dis;
				corners[1][j] += dis;
			}
			return corners;
		};
		static const fastEnvelope::Vector3 origin = fastEnvelope::Vector3(0, 0, 0);

		halfspace.resize(m_faces.size());
		cornerlist.resize(m_faces.size());
		Vector3 AB, AC, BC, normal, vector1, ABn, min, max;
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
		static const int p_face[8][3] = { {0, 1, 3}, {7, 6, 9}, {1, 0, 7}, {2, 1, 7}, {3, 2, 8}, {3, 9, 10}, {5, 4, 11}, {0, 5, 6} }; //prism triangle index. all with orientation.
		static const int c_face[6][3] = { {0, 1, 2}, {4, 7, 6}, {0, 3, 4}, {1, 0, 4}, {1, 5, 2}, {2, 6, 3} };
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
				halfspace[i].resize(6);
				for (int j = 0; j < 6; j++) {
					halfspace[i][j][0] = box[c_face[j][0]];
					halfspace[i][j][1] = box[c_face[j][1]];
					halfspace[i][j][2] = box[c_face[j][2]];
				}
				cornerlist[i] = (get_bb_corners_8(box));

				continue;
			}
			if (de == DEGENERATED_SEGMENT)
			{
				logger().debug("Envelope Triangle Degeneration- Segment");
				Scalar length1 = AB.norm(), length2 = AC.norm(), length3 = BC.norm();
				if (length1 >= length2 && length1 >= length3)
				{
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], tolerance, box);
					
				}
				if (length2 >= length1 && length2 >= length3)
				{
					seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][2]], tolerance, box);
					
				}
				if (length3 >= length1 && length3 >= length2)
				{
					seg_cube(m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tolerance, box);
				}
				halfspace[i].resize(6);
				for (int j = 0; j < 6; j++) {
					halfspace[i][j][0] = box[c_face[j][0]];
					halfspace[i][j][1] = box[c_face[j][1]];
					halfspace[i][j][2] = box[c_face[j][2]];
				}
				cornerlist[i] = (get_bb_corners_8(box));
				continue;
			}
			if (de == NERLY_DEGENERATED)
			{
				logger().debug("Envelope Triangle Degeneration- Nearly");
				
				normal = accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][0]], m_ver[m_faces[i][2]]);

				vector1 = accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], origin, normal);

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
			halfspace[i].resize(8);
			for (int j = 0; j < 8; j++) {
				halfspace[i][j][0] = polygonoff[p_face[j][0]];
				halfspace[i][j][1] = polygonoff[p_face[j][1]];
				halfspace[i][j][2] = polygonoff[p_face[j][2]];
			}
			cornerlist[i] = (get_bb_corners_12(polygonoff));
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

	Vector3 FastEnvelope::accurate_normal_vector(const Vector3 &p0, const Vector3 &p1,//TODO use marco's code
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

	bool FastEnvelope::is_two_facets_neighbouring(const int & pid, const int &i, const int &j)const {
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
		if (halfspace[pid].size() == 8) {
			int id = i * 8 + j;
			if (prism_map[id][0] == -1) return false;
		}
		if (halfspace[pid].size() == 6) {
			int id = i * 6 + j;
			if (cubic_map[id][0] == -1) return false;
		}
		return true;
	}



} // namespace fastEnvelope

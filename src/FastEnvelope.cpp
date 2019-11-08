#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Predicates.hpp>
#include <fastenvelope/Logger.hpp>
#include<fastenvelope/Morton.h>


#include <igl/Timer.h>
#include<igl/write_triangle_mesh.h>
#include <fstream>


int dbg1 = 0, dbg2 = 0, dbg3 = 0, dbg4 = 0, dbgout1 = 0, dbgout2 = 0, dbgout3 = 0, dbgout4 = 0, dbgout5 = 0, howmany = 0;
int ct1 = 0, ct2 = 0, ct3 = 0, ct4 = 0, ct5 = 0, count_int1 = 0, count_int2 = 0, count_int3 = 0, count_int4 = 0;
double time1 = 0, time2 = 0, time3 = 0, time4 = 0, time5 = 0, timesearch = 0, timecheck = 0,
time6 = 0, time7 = 0, time8 = 0, time9 = 0, time10 = 0, time11 = 0, time12 = 0, time13 = 0, time14 = 0, time15 = 0, time16 = 0, time17 = 0,
time18 = 0;
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
		std::cout << "where checking is out " << dbgout1 << " " << dbgout2 << " " << dbgout3 << " " << dbgout4 << " " << dbgout5 << std::endl;
		std::cout << "time for search, " << timesearch << std::endl;
		std::cout << "time for check, " << timecheck << std::endl;
		std::cout << "time for vertices, " << time1 << std::endl;
		std::cout << "time for degeneration, " << time2 << std::endl;
		std::cout << "time for intersection, " << time3 << std::endl;
		std::cout << "time for lpi, " << time4 << std::endl;
		std::cout << "time for local tree, " << time5 << std::endl;
		std::cout << "time for tpi all, " << time6 << std::endl;
		std::cout << "time for tpi part one check list, " << time14 << std::endl;
		std::cout << "time for tpi part two check list, " << time15 << std::endl;
		std::cout << "time for tpi part one check neighbours, " << time16 << std::endl;
		std::cout << "time for tpi part two check neighbours, " << time17 << std::endl;
		std::cout << "time for tpi post, " << time7 << std::endl;
		std::cout << "time for tpi part1 multi post, " << time12 << std::endl;
		std::cout << "time for tpi number, " << ct1 << std::endl;
		//std::cout << "time for tpi part1, " << time8 << std::endl;
		//std::cout << "time for tpi part2, " << time9 << std::endl;
		std::cout << "time for tpi multi post time, " << time10 << std::endl;
		std::cout << "time for tpi multi post number, " << ct2 << std::endl;
		std::cout << "time for see if point on prism, " << time11 << std::endl;
		std::cout << "time for local search, " << time13 << std::endl;
		std::cout << "time for 3 triangle intersection, " << time18 << std::endl;
		std::cout << "count intersections " << count_int1<<" "<<count_int2<< " "<<count_int3<<" "<< count_int4 << std::endl;
	}
	void FastEnvelope::reset_time() {
		time1 = 0; time2 = 0; time3 = 0; time4 = 0; time5 = 0; timesearch = 0; timecheck = 0;
		time6 = 0; time7 = 0; time8 = 0; time9 = 0; time10 = 0;
		time11 = 0; time12 = 0; time13 = 0;
		ct1 = 0; ct2 = 0;
		time14 = 0; time15 = 0; time16 = 0; time17 = 0; time18 = 0;
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
		/*std::vector<Vector3> ver_new;
		std::vector<Vector3i> faces_new;*/
		timer.start();
		
		//std::vector<Vector3> ver_new;
		std::vector<Vector3i> faces_new;
		//----
		/*ver_new.clear();
		faces_new.clear();
		GEO::Mesh M;
		to_geogram_mesh(m_ver, m_faces, M);
		GEO::mesh_reorder(M, GEO::MESH_ORDER_MORTON);
		from_geogram_mesh(M, ver_new, faces_new);
*/
		//----
		/*ver_new = m_ver;
		faces_new = m_faces;*/
		//----
		//ver_new = m_ver;
		resorting(m_ver, m_faces, faces_new);//resort the facets order
		//----
		timer.stop();
		logger().info("Resorting mesh time {}s", timer.getElapsedTimeInSec());

		timer.start();

		halfspace_init(m_ver, faces_new, halfspace, cornerlist, epsilon);


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
	void FastEnvelope::resorting(const std::vector<Vector3> &V, const std::vector<Vector3i> &F, std::vector<Vector3i> &fnew) {
		std::vector<std::array<int,3>> ct;
		struct sortstruct {
			int order;
			Resorting::MortonCode64 morton;
		};
		std::vector<sortstruct> list;
		const int multi = 1000;
		ct.resize(F.size());
		list.resize(F.size());
		
		for (int i = 0; i < F.size(); i++) {
			ct[i][0] = int(((V[F[i][0]] + V[F[i][1]] + V[F[i][2]])*multi)[0]);
			ct[i][1] = int(((V[F[i][0]] + V[F[i][1]] + V[F[i][2]])*multi)[1]);
			ct[i][2] = int(((V[F[i][0]] + V[F[i][1]] + V[F[i][2]])*multi)[2]);
			list[i].morton = Resorting::MortonCode64(ct[i][0], ct[i][1], ct[i][2]);
			list[i].order = i;
		}
		const auto morton_compare = [](const sortstruct &a, const sortstruct &b)
		{
			return (a.morton < b.morton);
		};
		std::sort(list.begin(), list.end(), morton_compare);

		fnew.resize(F.size());
		for (int i = 0; i < F.size(); i++) {
			fnew[i] = F[list[i].order];
		}


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
		//bool res=debugcode(triangle, querylist);
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
		//std::cout << "prism size " << prismindex.size() << std::endl;


		std::vector<unsigned int> filted_intersection; filted_intersection.reserve(prismindex.size() / 3);
		std::vector<std::vector<int>>intersect_face; intersect_face.reserve(prismindex.size() / 3);
		bool out, cut;

		int inter, inter1, record1, record2,

			tti; //triangle-triangle intersection

		jump1 = -1;

		int check_id;

		timer.start();
		for (int i = 0; i < 3; i++) {
			out = point_out_prism_return_id(triangle[i], prismindex, jump1, check_id);

			if (out) {
				dbgout1++;

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
			std::vector<unsigned int > queue, idlist;
			queue.emplace_back(check_id);//queue contains the id in prismindex
			idlist.emplace_back(prismindex[check_id]);
			std::vector<unsigned int > list;
			std::vector<int> cidl; cidl.reserve(8);
			for (int i = 0; i < queue.size(); i++) {

				jump1 = prismindex[queue[i]];

				for (int k = 0; k < 3; k++) {
					cut = is_seg_cut_polyhedra(jump1, triangle[triseg[k][0]], triangle[triseg[k][1]], cidl);
					if (cut&&cidl.size() == 0)continue;
					if (!cut) continue;
					for (int j = 0; j < cidl.size(); j++) {
						/*tti = seg_cut_plane(triangle[triseg[k][0]], triangle[triseg[k][1]],
							halfspace[prismindex[queue[i]]][cidl[j]][0], halfspace[prismindex[queue[i]]][cidl[j]][1], halfspace[prismindex[queue[i]]][cidl[j]][2]);
						if (tti != CUT_FACE) continue;*/
						inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(triangle[triseg[k][0]], triangle[triseg[k][1]],
							halfspace[prismindex[queue[i]]][cidl[j]][0], halfspace[prismindex[queue[i]]][cidl[j]][1], halfspace[prismindex[queue[i]]][cidl[j]][2],
							idlist, jump1, check_id);

						if (inter == 1)
						{
							inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(triangle[triseg[k][0]], triangle[triseg[k][1]],
								halfspace[prismindex[queue[i]]][cidl[j]][0], halfspace[prismindex[queue[i]]][cidl[j]][1], halfspace[prismindex[queue[i]]][cidl[j]][2],
								prismindex, jump1, check_id);

							if (inter == 1) {
								time2 += timer.getElapsedTimeInSec();
								dbgout4++;
								return true;
							}
							if (inter == 0) {
								queue.emplace_back(check_id);
								idlist.emplace_back(prismindex[check_id]);
							}
						}
					}
				}
			}

			return false;
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
			else if (tti == 1 && cidl.size() > 0) {
				//if(prismindex[i] == check_id) 
				filted_intersection.emplace_back(prismindex[i]);
				intersect_face.emplace_back(cidl);

				dbg3++;
			}
		}
		time3 += timer.getElapsedTimeInSec();

		if (filted_intersection.size() == 0) {
			time3 += timer.getElapsedTimeInSec();
			return false;//inside
		}


		//tree
		timer.start();
		AABB localtree;
		std::vector<std::array<Vector3, 2>> localcorner;
		localcorner.resize(filted_intersection.size());

		for (int i = 0; i < filted_intersection.size(); i++) {
			localcorner[i] = cornerlist[filted_intersection[i]];
		}

		localtree.init_envelope(localcorner);
		time5 += timer.getElapsedTimeInSec();


		std::vector<unsigned int > queue, idlist;
		queue.emplace_back(0);//queue contains the id in filted_intersection
		idlist.emplace_back(filted_intersection[queue[0]]);

		std::vector<std::vector<int>> neighbour_facets, idlistorder;

		idlistorder.emplace_back(intersect_face[queue[0]]);

		time3 += timer.getElapsedTimeInSec();






		std::vector<unsigned int> neighbours;//local id

		std::vector<unsigned int > list;
		timer.start();
		for (int i = 0; i < queue.size(); i++) {

			jump1 = filted_intersection[queue[i]];

			for (int k = 0; k < 3; k++) {
				for (int j = 0; j < intersect_face[queue[i]].size(); j++) {
					tti = seg_cut_plane(triangle[triseg[k][0]], triangle[triseg[k][1]],
						halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][0], halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][1], halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][2]);
					if (tti != CUT_FACE) continue;
					inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(triangle[triseg[k][0]], triangle[triseg[k][1]],
						halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][0], halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][1], halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][2],
						idlist, idlistorder, jump1, check_id);

					if (inter == 1)
					{
						inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(triangle[triseg[k][0]], triangle[triseg[k][1]],
							halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][0], halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][1], halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][2],
							filted_intersection, intersect_face, jump1, check_id);

						if (inter == 1) {
							dbgout2++;
							time4 += timer.getElapsedTimeInSec();
							return true;
						}
						if (inter == 0) {
							idlistorder.emplace_back(intersect_face[check_id]);
							queue.emplace_back(check_id);
							idlist.emplace_back(filted_intersection[check_id]);
						}
					}
				}
			}
		}


		time4 += timer.getElapsedTimeInSec();
		//std::cout << "the size of queue " << queue.size()<<" "<<idlist.size()  << std::endl;
		/*obb tobb;
		tobb = obb::build_triangle_obb_matrixs(triangle[0], triangle[1], triangle[2]);*/
	
		//tpi part
		igl::Timer timer1;
		timer.start();
		DATA_TPI datatpi;
		for (int i = 1; i < queue.size(); i++)
		{
			jump1 = filted_intersection[queue[i]];
			timer1.start();
			localtree.bbd_finding_in_envelope(cornerlist[jump1][0], cornerlist[jump1][1], list);
			neighbours.clear();
			neighbours.resize(list.size());
			neighbour_facets.resize(list.size());
			for (int j = 0; j < list.size(); j++) {
				neighbours[j] = filted_intersection[list[j]];
				neighbour_facets[j] = intersect_face[list[j]];
			}


			time13 += timer1.getElapsedTimeInSec();
			for (int j = 0; j < i; j++) {
				jump2 = filted_intersection[queue[j]];
				if (!box_box_intersection(cornerlist[jump1][0], cornerlist[jump1][1], cornerlist[jump2][0], cornerlist[jump2][1]))
					continue;
				for (int k = 0; k < intersect_face[queue[i]].size(); k++) {
					for (int h = 0; h < intersect_face[queue[j]].size(); h++) {
						//box on faces?
						count_int1++;
						/*cut = tobb.intersects(obblist[jump1][intersect_face[queue[i]][k]], obblist[jump2][intersect_face[queue[j]][h]]);
						if (cut == false) continue;*/
						timer1.start();
						cut = is_3_triangle_cut(triangle,
							halfspace[jump1][intersect_face[queue[i]][k]][0],
							halfspace[jump1][intersect_face[queue[i]][k]][1],
							halfspace[jump1][intersect_face[queue[i]][k]][2],

							halfspace[jump2][intersect_face[queue[j]][h]][0],
							halfspace[jump2][intersect_face[queue[j]][h]][1],
							halfspace[jump2][intersect_face[queue[j]][h]][2]);
						
						time18 += timer1.getElapsedTimeInSec();
						if (!cut) continue;
						
						/*cut = obblist[jump1][intersect_face[queue[i]][k]].intersects(obblist[jump2][intersect_face[queue[j]][h]]);
						if (cut) count_int3++;*/
						timer1.start();
						cut = is_tpp_on_polyhedra(triangle,
							halfspace[jump1][intersect_face[queue[i]][k]][0],
							halfspace[jump1][intersect_face[queue[i]][k]][1],
							halfspace[jump1][intersect_face[queue[i]][k]][2],

							halfspace[jump2][intersect_face[queue[j]][h]][0],
							halfspace[jump2][intersect_face[queue[j]][h]][1],
							halfspace[jump2][intersect_face[queue[j]][h]][2], jump1, intersect_face[queue[i]][k]);
						time11 += timer1.getElapsedTimeInSec();
						if (!cut) continue;

						timer1.start();
						cut = is_tpp_on_polyhedra(triangle,
							halfspace[jump1][intersect_face[queue[i]][k]][0],
							halfspace[jump1][intersect_face[queue[i]][k]][1],
							halfspace[jump1][intersect_face[queue[i]][k]][2],

							halfspace[jump2][intersect_face[queue[j]][h]][0],
							halfspace[jump2][intersect_face[queue[j]][h]][1],
							halfspace[jump2][intersect_face[queue[j]][h]][2], jump2, intersect_face[queue[j]][h]);
						time11 += timer1.getElapsedTimeInSec();
						if (!cut) continue;
						count_int2++;
						/*cut = tobb.intersects(obblist[jump1][intersect_face[queue[i]][k]], obblist[jump2][intersect_face[queue[j]][h]]);
						if (cut == false) std::cout << "here wrong for obb intersection" << std::endl;
						if (cut == true) count_int4++;*/
						timer1.start();
						inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(triangle,
							halfspace[jump1][intersect_face[queue[i]][k]][0],
							halfspace[jump1][intersect_face[queue[i]][k]][1],
							halfspace[jump1][intersect_face[queue[i]][k]][2],

							halfspace[jump2][intersect_face[queue[j]][h]][0],
							halfspace[jump2][intersect_face[queue[j]][h]][1],
							halfspace[jump2][intersect_face[queue[j]][h]][2],
							idlist, idlistorder, jump1, jump2, check_id);

						time14 += timer1.getElapsedTimeInSec();
						if (inter == 1) {
							timer1.start();
							inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(triangle,
								halfspace[jump1][intersect_face[queue[i]][k]][0],
								halfspace[jump1][intersect_face[queue[i]][k]][1],
								halfspace[jump1][intersect_face[queue[i]][k]][2],

								halfspace[jump2][intersect_face[queue[j]][h]][0],
								halfspace[jump2][intersect_face[queue[j]][h]][1],
								halfspace[jump2][intersect_face[queue[j]][h]][2],
								neighbours, neighbour_facets, jump1, jump2, check_id);

							time16 += timer1.getElapsedTimeInSec();
							if (inter == 1) {
								dbgout3++;
								time6 += timer.getElapsedTimeInSec();
								return true;
							}
							if (inter == 0) {
								idlistorder.emplace_back(intersect_face[list[check_id]]);
								queue.emplace_back(list[check_id]);
								idlist.emplace_back(filted_intersection[list[check_id]]);
							}
						}
					}
				}
			}
		}
		//std::cout << "the size of queue " << queue.size() << " " << idlist.size() << std::endl;
		time6 += timer.getElapsedTimeInSec();
		return false;


	}
	bool FastEnvelope::debugcode(const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex) const {
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
			if (tot == halfspace[prismindex[i]].size())
			{
				return IN_PRISM;
			}
		}


		return OUT_PRISM;
	}
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision_return_id(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle,
		const std::vector<unsigned int> &prismindex, int& id) const
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

	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_pure_multiprecision_return_local_id(const DATA_LPI &datalpi, const std::array<Vector3, 3> &triangle,
		const std::vector<unsigned int> &prismindex, int& id) const
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
			if (tot == halfspace[prismindex[i]].size())
			{
				id = i;
				return IN_PRISM;
			}
		}


		return OUT_PRISM;
	}


	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_double_return_id(
		LPI_filtered_suppvars& sf,
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
						sf,
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

	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
		const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
		const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const {
		LPI_filtered_suppvars sf;
		bool precom = orient3D_LPI_prefilter( //
			segpoint0[0], segpoint0[1], segpoint0[2],
			segpoint1[0], segpoint1[1], segpoint1[2],
			triangle0[0], triangle0[1], triangle0[2],
			triangle1[0], triangle1[1], triangle1[2],
			triangle2[0], triangle2[1], triangle2[2],
			sf);
		if (precom == false) {
			int tot, ori;


			LPI_exact_suppvars s;
			bool premulti = orient3D_LPI_pre_exact(
				segpoint0[0], segpoint0[1], segpoint0[2],
				segpoint1[0], segpoint1[1], segpoint1[2],

				triangle0[0], triangle0[1], triangle0[2],
				triangle1[0], triangle1[1], triangle1[2],
				triangle2[0], triangle2[1], triangle2[2], s);

			// time_multi += timer.getElapsedTimeInSec();
			if (premulti == false) return 2;
			for (int i = 0; i < prismindex.size(); i++)
			{
				if (prismindex[i] == jump)
				{
					continue;
				}
				tot = 0;

				for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

					// timer.start();
					ori = orient3D_LPI_post_exact(s,
						segpoint0[0], segpoint0[1], segpoint0[2],

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
					id = i;
					return IN_PRISM;
				}
			}


			return OUT_PRISM;
		}
		int tot;
		int ori;
		static INDEX index;
		static std::vector<INDEX> recompute;

		recompute.clear();
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump) continue;

			index.FACES.clear();
			tot = 0;


			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

				ori =
					orient3D_LPI_postfilter(
						sf,
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
				id = i;
				return IN_PRISM;
			}

			if (ori != 1)
			{
				assert(!index.FACES.empty());
				index.Pi = i;//local id
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
				int in1 = prismindex[recompute[k].Pi];

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
					id = recompute[k].Pi;
					return IN_PRISM;
				}
			}
		}

		return OUT_PRISM;
	}
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
		const Vector3 &segpoint0, const Vector3 &segpoint1, const Vector3 &triangle0,
		const Vector3 &triangle1, const Vector3 &triangle2, const std::vector<unsigned int> &prismindex,
		const std::vector<std::vector<int>>& intersect_face, const int &jump, int &id) const {
		LPI_filtered_suppvars sl;
		bool precom = orient3D_LPI_prefilter( //
			segpoint0[0], segpoint0[1], segpoint0[2],
			segpoint1[0], segpoint1[1], segpoint1[2],
			triangle0[0], triangle0[1], triangle0[2],
			triangle1[0], triangle1[1], triangle1[2],
			triangle2[0], triangle2[1], triangle2[2],
			sl);
		if (precom == false) {
			int tot, ori, fid;

			LPI_exact_suppvars s;
			bool premulti = orient3D_LPI_pre_exact(
				segpoint0[0], segpoint0[1], segpoint0[2],
				segpoint1[0], segpoint1[1], segpoint1[2],

				triangle0[0], triangle0[1], triangle0[2],
				triangle1[0], triangle1[1], triangle1[2],
				triangle2[0], triangle2[1], triangle2[2], s);

			// time_multi += timer.getElapsedTimeInSec();
			if (premulti == false) return 2;
			for (int i = 0; i < prismindex.size(); i++)
			{
				if (prismindex[i] == jump)
				{
					continue;
				}
				tot = 0; fid = 0;

				for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

					// timer.start();
					if (intersect_face[i][fid] == j)
					{
						if (fid + 1 < intersect_face[i].size()) fid++;
					}
					else continue;
					ori = orient3D_LPI_post_exact(s,
						segpoint0[0], segpoint0[1], segpoint0[2],

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
				if (ori == 1 || ori == 0) continue;
				fid = 0;
				for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

					// timer.start();
					if (intersect_face[i][fid] == j)
					{
						if (fid + 1 < intersect_face[i].size()) fid++;
						continue;
					}

					ori = orient3D_LPI_post_exact(s,
						segpoint0[0], segpoint0[1], segpoint0[2],

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
				if (ori == 1 || ori == 0) continue;
				if (tot == prismindex.size())
				{
					id = i;
					return IN_PRISM;
				}
			}


			return OUT_PRISM;
		}
		int tot, fid;
		int ori;
		static INDEX index;
		static std::vector<INDEX> recompute;

		recompute.clear();
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump) continue;

			index.FACES.clear();
			tot = 0; fid = 0;


			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
				if (intersect_face[i][fid] == j)
				{
					if (fid + 1 < intersect_face[i].size()) fid++;
				}
				else continue;
				ori =
					orient3D_LPI_postfilter(
						sl,
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
			if (ori == 1) continue;
			fid = 0;
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
				if (intersect_face[i][fid] == j)
				{
					if (fid + 1 < intersect_face[i].size()) fid++;
					continue;
				}

				ori =
					orient3D_LPI_postfilter(
						sl,
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
			if (ori == 1) continue;
			if (tot == halfspace[prismindex[i]].size())
			{
				id = i;
				return IN_PRISM;
			}

			if (ori != 1)
			{
				assert(!index.FACES.empty());
				index.Pi = i;//local id
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
				int in1 = prismindex[recompute[k].Pi];

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
					id = recompute[k].Pi;
					return IN_PRISM;
				}
			}
		}

		return OUT_PRISM;
	}
	int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_double_return_local_id(
		LPI_filtered_suppvars& sf,
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
						sf,
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
				id = i;
				return IN_PRISM;
			}

			if (ori != 1)
			{
				assert(!index.FACES.empty());
				index.Pi = i;//local id
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
				int in1 = prismindex[recompute[k].Pi];

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
					id = recompute[k].Pi;
					return IN_PRISM;
				}
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
			if (prismindex[i] == jump1 || prismindex[i] == jump2)	continue;
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
	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_pure_multiprecision_return_local_id(const DATA_TPI &datatpi,
		const std::array<Vector3, 3> &triangle, const std::vector<unsigned int> &prismindex, TPI_exact_suppvars& s, int &id) const
	{
		igl::Timer timer;
		int tot, ori;
		int jump1 = datatpi.jump1, jump2 = datatpi.jump2;



		dbg1++;
		dbg2++;
		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump1 || prismindex[i] == jump2)	continue;
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
				id = i;
				return IN_PRISM;
			}

		}
		return OUT_PRISM;
	}






	int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
		const std::vector<unsigned int> &prismindex, const std::vector<std::vector<int>>&intersect_face, const int &jump1, const int &jump2,
		int &id) const {
		igl::Timer timer;
		TPI_exact_suppvars s;
		TPI_filtered_suppvars st;
		Scalar d, n1d, n2d, n3d, max1, max2, max3, max4, max5, max6, max7;
		int tot, ori, fid;
		bool pre = orient3D_TPI_prefilter(triangle[0][0], triangle[0][1], triangle[0][2],
			triangle[1][0], triangle[1][1], triangle[1][2], triangle[2][0], triangle[2][1], triangle[2][2],
			facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
			facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
			st);

		if (pre == false) {
			bool premulti = orient3D_TPI_pre_exact(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2], triangle[2][0], triangle[2][1], triangle[2][2],
				facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
				s);
			if (premulti == false) return 2;
			for (int i = 0; i < prismindex.size(); i++)
			{

				if (prismindex[i] == jump1 || prismindex[i] == jump2)	continue;
				if (!box_box_intersection(cornerlist[prismindex[i]][0], cornerlist[prismindex[i]][1], cornerlist[jump1][0], cornerlist[jump1][1])) continue;
				if (!box_box_intersection(cornerlist[prismindex[i]][0], cornerlist[prismindex[i]][1], cornerlist[jump2][0], cornerlist[jump2][1])) continue;

				tot = 0;
				fid = 0;
				for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
					if (intersect_face[i][fid] == j)
					{
						if (fid + 1 < intersect_face[i].size()) fid++;
					}
					else continue;
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
				if (ori == 1 || ori == 0) continue;
				fid = 0;
				for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
					if (intersect_face[i][fid] == j)
					{
						if (fid + 1 < intersect_face[i].size()) fid++;
						continue;
					}

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
				if (ori == 1 || ori == 0) continue;
				if (tot == halfspace[prismindex[i]].size())
				{
					id = i;
					return IN_PRISM;
				}

			}

			return OUT_PRISM;
		}

		static INDEX index;
		static std::vector<INDEX> recompute;
		recompute.clear();


		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump1 || prismindex[i] == jump2)
				continue;
			if (!box_box_intersection(cornerlist[prismindex[i]][0], cornerlist[prismindex[i]][1], cornerlist[jump1][0], cornerlist[jump1][1])) continue;
			if (!box_box_intersection(cornerlist[prismindex[i]][0], cornerlist[prismindex[i]][1], cornerlist[jump2][0], cornerlist[jump2][1])) continue;

			index.FACES.clear();
			tot = 0;
			fid = 0;
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {


				if (intersect_face[i][fid] == j)
				{
					if (fid + 1 < intersect_face[i].size()) fid++;
				}
				else continue;
				timer.start();
				ct1 += 1;
				ori =
					orient3D_TPI_postfilter(
						st,
						halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
						halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
						halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);
				// timetpp1 += timer_a.getElapsedTimeInSec();
				time7 += timer.getElapsedTimeInSec();
				if (ori == 1)
				{
					//std::cout << "SHOULD HAPPEN" << std::endl;
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
			fid = 0;
			for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

				if (intersect_face[i][fid] == j)
				{
					if (fid + 1 < intersect_face[i].size()) fid++;
					continue;
				}

				timer.start();
				ct1 += 1;
				ori =
					orient3D_TPI_postfilter(
						st,
						halfspace[prismindex[i]][j][0][0], halfspace[prismindex[i]][j][0][1], halfspace[prismindex[i]][j][0][2],
						halfspace[prismindex[i]][j][1][0], halfspace[prismindex[i]][j][1][1], halfspace[prismindex[i]][j][1][2],
						halfspace[prismindex[i]][j][2][0], halfspace[prismindex[i]][j][2][1], halfspace[prismindex[i]][j][2][2]);
				// timetpp1 += timer_a.getElapsedTimeInSec();
				time7 += timer.getElapsedTimeInSec();
				if (ori == 1)
				{
					//std::cout << "should not happen" << std::endl;
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
				id = i;
				return IN_PRISM;
			}

			if (ori != 1)
			{
				index.Pi = i;
				recompute.emplace_back(index);

			}
		}
		timer.start();
		if (recompute.size() > 0)
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
			if (premulti == false) return 2;
			for (int k = 0; k < recompute.size(); k++)
			{
				int in1 = prismindex[recompute[k].Pi];

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
					id = recompute[k].Pi;
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


	bool FastEnvelope::is_3_triangle_cut(
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22)
	{
		TPI_filtered_suppvars st;
		TPI_exact_suppvars s;
		static Scalar
			t00, t01, t02,
			t10, t11, t12,
			t20, t21, t22,

			f100, f101, f102,
			f110, f111, f112,
			f120, f121, f122,

			f200, f201, f202,
			f210, f211, f212,
			f220, f221, f222;
		bool pre =
			orient3D_TPI_prefilter(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],
				facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
				st);
		if (pre == false) {
			bool premulti = orient3D_TPI_pre_exact(
				triangle[0][0], triangle[0][1], triangle[0][2],
				triangle[1][0], triangle[1][1], triangle[1][2],
				triangle[2][0], triangle[2][1], triangle[2][2],
				facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
				s);
			if (premulti == false) return false;
			return is_3_triangle_cut_pure_multiprecision(triangle, s);
		}

		Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];
		if (Predicates::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0)
		{
			logger().debug("Degeneration happens");
			//move this guy in constructor and use fixed seed
			n = { {Vector3(rand(), rand(), rand())} };
		}


		bool premulti = false;
		int o1 = orient3D_TPI_postfilter(st, n[0], n[1], n[2],
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

		int o2 = orient3D_TPI_postfilter(st, n[0], n[1], n[2],
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

		int o3 = orient3D_TPI_postfilter(st, n[0], n[1], n[2],
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
		TPI_filtered_suppvars st;
		bool pre =
			orient3D_TPI_prefilter(
				tri0[0], tri0[1], tri0[2],
				tri1[0], tri1[1], tri1[2],
				tri2[0], tri2[1], tri2[2],
				facet10[0], facet10[1], facet10[2], facet11[0], facet11[1], facet11[2], facet12[0], facet12[1], facet12[2],
				facet20[0], facet20[1], facet20[2], facet21[0], facet21[1], facet21[2], facet22[0], facet22[1], facet22[2],
				st);

		if (pre == false)
			return 2; // means we dont know
		int o1 = orient3D_TPI_postfilter(st, n[0], n[1], n[2],
			tri0[0], tri0[1], tri0[2],
			tri1[0], tri1[1], tri1[2]);
		int o2 = orient3D_TPI_postfilter(st, n[0], n[1], n[2],
			tri1[0], tri1[1], tri1[2],
			tri2[0], tri2[1], tri2[2]);
		if (o1 * o2 == -1)
			return 0;
		int o3 = orient3D_TPI_postfilter(st, n[0], n[1], n[2],
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

		LPI_filtered_suppvars sl;
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
					sl);
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
							sl,
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
		TPI_filtered_suppvars st;

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
					cut[cutp[i]] = true;
					cut[cutp[j]] = true;
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
						st);

				for (int k = 0; k < cutp.size(); k++)
				{

					if (k == i || k == j)
						continue;

					ori =
						orient3D_TPI_postfilter(st,
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
	bool FastEnvelope::is_tpp_on_polyhedra(
		const std::array<Vector3, 3> &triangle,
		const Vector3 &facet10, const Vector3 &facet11, const Vector3 &facet12, const Vector3 &facet20, const Vector3 &facet21, const Vector3 &facet22,
		const int &prismid, const int &faceid)const {
		int ori;
		TPI_filtered_suppvars st;
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
				st);
		if (pre == true) {
			for (int i = 0; i < halfspace[prismid].size(); i++) {
				/*bool neib = is_two_facets_neighbouring(prismid, i, faceid);
				if (neib == false) continue;*/
				if (i == faceid) continue;
				ori =
					orient3D_TPI_postfilter(st,
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

	bool FastEnvelope::is_seg_cut_polyhedra(const int &cindex,
		const Vector3 &seg0, const Vector3 &seg1, std::vector<int> &cid) const
	{


		cid.clear();
		std::vector<bool> cut;
		cut.resize(halfspace[cindex].size());
		for (int i = 0; i < halfspace[cindex].size(); i++)
		{
			cut[i] = false;
		}
		std::vector<int> o1, o2;
		o1.resize(halfspace[cindex].size());
		o2.resize(halfspace[cindex].size());
		int ori = 0, count = 0;
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
			if (o1[i] == -1) count++;
		}
		if (cutp.size() == 0 && count == halfspace[cindex].size()) return true;
		if (cutp.size() == 0)
		{
			return false;
		}

		LPI_filtered_suppvars sf;
		for (int i = 0; i < cutp.size(); i++)
		{

			bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
				seg0[0], seg0[1], seg0[2],
				seg1[0], seg1[1], seg1[2],
				halfspace[cindex][cutp[i]][0][0], halfspace[cindex][cutp[i]][0][1], halfspace[cindex][cutp[i]][0][2],
				halfspace[cindex][cutp[i]][1][0], halfspace[cindex][cutp[i]][1][1], halfspace[cindex][cutp[i]][1][2],
				halfspace[cindex][cutp[i]][2][0], halfspace[cindex][cutp[i]][2][1], halfspace[cindex][cutp[i]][2][2],
				sf);
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
						sf,
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
	bool FastEnvelope::point_out_prism_return_id(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump, int &id) const
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
					id = i;
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
			Rational a0(a[0]), a1(a[1]), a2(a[2]), b0(b[0]), b1(b[1]), b2(b[2]), dot;
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
			const Scalar dis = 1e-6;
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
			const Scalar dis = 1e-6;
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
		static const std::array<std::vector<int>, 8> p_facepoint = {
		{{0,1,2,3,4,5},
	{8,7,6,11,10,9},
	{7,1,0,6},
	{2,1,7,8},
	{3,2,8,9},
	{4,3,9,10},
	{4,10,11,5},
	{6,0,5,11}}
		};

		static const std::array<std::array<int, 4>, 6> c_facepoint = {
				{
					{0,1,2,3},
			{4,7,6,5},
			{4,0,3,7},
			{1,0,4,5},
			{2,1,5,6},
			{3,2,6,7}
				}
		};

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

				/*envelope_vertices.resize(8);
				for (int j = 0; j < 8; j++) {
					envelope_vertices[j][0] = box[j][0];
					envelope_vertices[j][1] = box[j][1];
					envelope_vertices[j][2] = box[j][2];
				}
				obblist[i] = obb::build_obb_matrixs(envelope_vertices, p_face, c_face, p_facepoint, c_facepoint);*/
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

				/*envelope_vertices.resize(8);
				for (int j = 0; j < 8; j++) {
					envelope_vertices[j][0] = box[j][0];
					envelope_vertices[j][1] = box[j][1];
					envelope_vertices[j][2] = box[j][2];
				}
				obblist[i] = obb::build_obb_matrixs(envelope_vertices, p_face, c_face, p_facepoint, c_facepoint);*/
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
			/*envelope_vertices.resize(12);
			for (int j = 0; j < 12; j++) {
				envelope_vertices[j][0] = polygonoff[j][0];
				envelope_vertices[j][1] = polygonoff[j][1];
				envelope_vertices[j][2] = polygonoff[j][2];
			}
			obblist[i] = obb::build_obb_matrixs(envelope_vertices, p_face, c_face, p_facepoint, c_facepoint);*/
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
		
		logger().debug("use accurate normal vector solver for nearly degenerated triangles");
		
		Rational p00r(p0[0]), p01r(p0[1]), p02r(p0[2]),
			p10r(p1[0]), p11r(p1[1]), p12r(p1[2]),
			q00r(q0[0]), q01r(q0[1]), q02r(q0[2]),
			q10r(q1[0]), q11r(q1[1]), q12r(q1[2]);
		Rational axr(p10r - p00r), ayr(p11r - p01r), azr(p12r - p02r),
			bxr(q10r - q00r), byr(q11r - q01r), bzr(q12r - q02r);
		Rational xr = ayr * bzr - azr * byr;
		Rational yr = azr * bxr - axr * bzr;
		Rational zr = axr * byr - ayr * bxr;//get the direction (x,y,z), now normalize
		int xsign, ysign, zsign;
		xsign = xr.get_sign();
		ysign = yr.get_sign();
		zsign = zr.get_sign();
		Rational ssumr = xr * xr + yr * yr + zr * zr;
		xr = xr * xr / ssumr;
		yr = yr * yr / ssumr;
		zr = zr * zr / ssumr;

		Scalar x, y, z;
		x = sqrt(xr.to_double())*xsign;
		y = sqrt(yr.to_double())*ysign;
		z = sqrt(zr.to_double())*zsign;
		return Vector3(x, y, z);
		
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
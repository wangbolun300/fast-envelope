#include<fastenvelope/common_algorithms.h>
#include<fastenvelope/Predicates.hpp>
#include<iostream>
#include<array>
#include<fastenvelope/Morton.h>
//#include<fastenvelope/ //logger.hpp>
#include <indirectpredicates/ip_filtered_ex.h>
namespace fastEnvelope {
	namespace algorithms {
		int seg_cut_plane(const Vector3 &seg0, const Vector3 &seg1, const Vector3 &t0, const Vector3 &t1, const Vector3 &t2)
		{
			int o1, o2;
			o1 = Predicates::orient_3d(seg0, t0, t1, t2);
			o2 = Predicates::orient_3d(seg1, t0, t1, t2);
			int op = o1 * o2;
			if (op >= 0)
			{
				return FE_CUT_COPLANAR; //in fact, coplanar and not cut this plane
			}
			return FE_CUT_FACE;
		}
		void get_tri_corners(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2, Vector3 &mint, Vector3 &maxt)
		{
			mint[0] = std::min(std::min(triangle0[0], triangle1[0]), triangle2[0]);
			mint[1] = std::min(std::min(triangle0[1], triangle1[1]), triangle2[1]);
			mint[2] = std::min(std::min(triangle0[2], triangle1[2]), triangle2[2]);
			maxt[0] = std::max(std::max(triangle0[0], triangle1[0]), triangle2[0]);
			maxt[1] = std::max(std::max(triangle0[1], triangle1[1]), triangle2[1]);
			maxt[2] = std::max(std::max(triangle0[2], triangle1[2]), triangle2[2]);
		}
		bool box_box_intersection(const Vector3 &min1, const Vector3 &max1, const Vector3 &min2, const Vector3 &max2)//TDOO
		{
			if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
				return 0;
			if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
				return 0;
			return 1;
		}
		Vector2 to_2d(const Vector3 &p, int t)

		{
			return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
		}


		int is_triangle_degenerated(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2)
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

		
		Vector3 accurate_normal_vector(const Vector3 &p0, const Vector3 &p1,
			const Vector3 &p2)
		{
			Vector3 v;
			triangle_normal_exact(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],v[0],v[1],v[2]);
			return v;
		}
		Vector3 accurate_cross_product_direction(const Vector3 &p0, const Vector3 &p1,
			const Vector3 &q0, const Vector3 &q1)
		{
			Vector3 v;
			cross_product_normalized_exact(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], q0[0], q0[1], q0[2], q1[0], q1[1], q1[2],
				v[0], v[1], v[2]);
			return v;
		}


		void resorting(const std::vector<Vector3> &V, const std::vector<Vector3i> &F, std::vector<Vector3i> &fnew) {
			std::vector<std::array<int, 3>> ct;
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

		void resorting(const std::vector<Vector3> &V, const std::vector<Vector3i> &F, std::vector<Vector3i> &fnew, std::vector<int>& new2old){
			std::vector<std::array<int, 3>> ct;
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
			new2old.resize(F.size());

			for (int i = 0; i < F.size(); i++) {
				new2old[i] = list[i].order;
			}

		}


		void seg_cube(const Vector3 &p1, const Vector3 &p2, const Scalar &width, std::array<Vector3, 8> &envbox)
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
	}
	
	
	
	
	// this function provide an algorithm build halfspaces for a list of triangles. each prism has 7-8 facets
	void algorithms::halfspace_generation(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
		std::vector<std::array<Vector3, 2>>& cornerlist, const Scalar &epsilon) {
		Scalar tolerance = epsilon / sqrt(3);// the envelope thickness, to be conservative
		Vector3 AB, AC, BC, normal;
		int de;
		std::array<Vector3, 3> plane;
		std::array<Vector3, 8> box;
		Vector3 tmin, tmax;
		Vector3 bbox_offset;
		bbox_offset[0] = 1;
		bbox_offset[1] = 1;
		bbox_offset[2] = 1;
		bbox_offset = bbox_offset * tolerance*sqrt(3)*(1 + 1e-6);

		const auto get_triangle_obtuse_angle= [](const Vector3& p0, const Vector3& p1, const Vector3& p2) {
			const auto dot_sign = [](const Vector3 &a, const Vector3 &b)
			{
				Scalar t = a.dot(b);
				if (t > SCALAR_ZERO)
					return 1;
				if (t < -1 * SCALAR_ZERO)
					return -1;

				return dot_product_sign(a[0], a[1], a[2], b[0], b[1], b[2]);
			};
			
			if (dot_sign(p1 - p0, p1 - p2) < 0) return 1;
			if (dot_sign(p2 - p1, p2 - p0) < 0) return 2;
			if (dot_sign(p0 - p1, p0 - p2) < 0) return 0;
			return -1;

		};
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
		const auto get_corner_plane=[](const Vector3& p0, const Vector3& midp, const Vector3 &normal, const Scalar& distance,
			Vector3& plane0, Vector3& plane1, Vector3& plane2, const bool use_exact) {
			Scalar distance_small = distance * 1;// to be conservative to reduce numerical error, can set the Scalar as 0.999
			Vector3 direction = (p0 - midp).normalized();
			plane0 = p0 + direction * distance_small;
			plane1 = plane0 + normal;
			Vector3 origin = Vector3(0, 0, 0);
			Vector3 axis;
			if (use_exact) axis = accurate_cross_product_direction(midp, p0, origin, normal);
			else axis = direction.cross(normal).normalized();
			plane2 = plane0 + axis;
		};
		bool use_accurate_cross = false;
		Vector3 origin = Vector3(0, 0, 0);

		static const int c_face[6][3] = { {0, 1, 2}, {4, 7, 6}, {0, 3, 4}, {1, 0, 4}, {1, 5, 2}, {2, 6, 3} };

		halfspace.resize(m_faces.size());
		cornerlist.resize(m_faces.size());
		for (int i = 0; i < m_faces.size(); i++)
		{
			algorithms::get_tri_corners(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tmin, tmax);
			cornerlist[i][0] = tmin - bbox_offset;
			cornerlist[i][1] = tmax + bbox_offset;

			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			de = algorithms::is_triangle_degenerated(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);

			if (de == DEGENERATED_POINT)
			{
				 //logger().debug("Envelope Triangle Degeneration- Point");
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
				//cornerlist[i] = (get_bb_corners_8(box));

				
				continue;
			}
			if (de == DEGENERATED_SEGMENT)
			{
				 //logger().debug("Envelope Triangle Degeneration- Segment");
				Scalar length1 = AB.dot(AB), length2 = AC.dot(AC), length3 = BC.dot(BC);
				if (length1 >= length2 && length1 >= length3)
				{
					algorithms::seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], tolerance, box);

				}
				if (length2 >= length1 && length2 >= length3)
				{
					algorithms::seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][2]], tolerance, box);

				}
				if (length3 >= length1 && length3 >= length2)
				{
					algorithms::seg_cube(m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tolerance, box);
				}
				halfspace[i].resize(6);
				for (int j = 0; j < 6; j++) {
					halfspace[i][j][0] = box[c_face[j][0]];
					halfspace[i][j][1] = box[c_face[j][1]];
					halfspace[i][j][2] = box[c_face[j][2]];
				}
				//cornerlist[i] = (get_bb_corners_8(box));

			
				continue;
			}
			if (de == NERLY_DEGENERATED)
			{
				//logger().debug("Envelope Triangle Degeneration- Nearly");
				use_accurate_cross = true;

				normal = algorithms::accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);

			}
			else
			{
				normal = AB.cross(AC).normalized();
			}
			halfspace[i].reserve(8);
			Vector3 normaldist = normal * tolerance;
			Vector3 edgedire, edgenormaldist;
			plane[0] = m_ver[m_faces[i][0]] + normaldist;
			plane[1] = m_ver[m_faces[i][1]] + normaldist;
			plane[2] = m_ver[m_faces[i][2]] + normaldist;
			halfspace[i].emplace_back(plane);// number 0

			plane[0] = m_ver[m_faces[i][0]] - normaldist;
			plane[1] = m_ver[m_faces[i][2]] - normaldist;
			plane[2] = m_ver[m_faces[i][1]] - normaldist;// order: 0, 2, 1
			halfspace[i].emplace_back(plane);// number 1
			
			int obtuse = get_triangle_obtuse_angle(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);


			edgedire = AB.normalized();
			if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(origin, edgedire, origin, normal)*tolerance;
			else edgenormaldist = edgedire.cross(normal).normalized()*tolerance;
			plane[0] = m_ver[m_faces[i][0]] + edgenormaldist;
			plane[1] = m_ver[m_faces[i][1]] + edgenormaldist;
			plane[2] = plane[0] + normal;
			halfspace[i].emplace_back(plane);// number 2

			
			if (obtuse != 1) {
				get_corner_plane(m_ver[m_faces[i][1]], (m_ver[m_faces[i][0]] + m_ver[m_faces[i][2]]) / 2, normal,
					tolerance, plane[0], plane[1], plane[2],use_accurate_cross);
				halfspace[i].emplace_back(plane);// number 3;

			}

			edgedire = BC.normalized();
			if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(origin, edgedire, origin, normal)*tolerance;
			else edgenormaldist = edgedire.cross(normal).normalized()*tolerance;

			plane[0] = m_ver[m_faces[i][1]] + edgenormaldist;
			plane[1] = m_ver[m_faces[i][2]] + edgenormaldist;
			plane[2] = plane[0] + normal;
			halfspace[i].emplace_back(plane);// number 4

			if (obtuse != 2) {
				get_corner_plane(m_ver[m_faces[i][2]], (m_ver[m_faces[i][0]] + m_ver[m_faces[i][1]]) / 2, normal,
					tolerance, plane[0], plane[1], plane[2],use_accurate_cross);
				halfspace[i].emplace_back(plane);// number 5;

			}

			edgedire = -AC.normalized();
			if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(origin, edgedire, origin, normal)*tolerance;
			else edgenormaldist = edgedire.cross(normal).normalized()*tolerance;
			
			plane[0] = m_ver[m_faces[i][2]] + edgenormaldist;
			plane[1] = m_ver[m_faces[i][0]] + edgenormaldist;
			plane[2] = plane[1] + normal;
			halfspace[i].emplace_back(plane);// number 6 

			if (obtuse != 0) {
				get_corner_plane(m_ver[m_faces[i][0]], (m_ver[m_faces[i][1]] + m_ver[m_faces[i][2]]) / 2, normal,
					tolerance, plane[0], plane[1], plane[2],use_accurate_cross);
				halfspace[i].emplace_back(plane);// number 7;

			}
			//std::cout << "envelope face nbr " << halfspace[i].size() << std::endl;
		
		}

	}

	// use user defined epsilon list to initialize adaptive envelope
	void algorithms::halfspace_generation(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces, std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
		std::vector<std::array<Vector3, 2>>& cornerlist, const std::vector<Scalar>& epsilons) {
		std::vector<Scalar> tolerance; tolerance.resize(epsilons.size());
		for (int i = 0; i < epsilons.size(); i++) {
			tolerance[i] = epsilons[i] / sqrt(3);// the envelope thickness, to be conservative
		}
		
		
		Vector3 AB, AC, BC, normal;
		int de;
		std::array<Vector3, 3> plane;
		std::array<Vector3, 8> box;
		Vector3 tmin, tmax;
		std::vector<Vector3> bbox_offset; bbox_offset.resize(epsilons.size());
		for (int i = 0; i < epsilons.size(); i++) {
			bbox_offset[i][0] = 1;
			bbox_offset[i][1] = 1;
			bbox_offset[i][2] = 1;
			bbox_offset[i] = bbox_offset[i] * tolerance[i]*sqrt(3)*(1 + 1e-6);
		}

		

		const auto get_triangle_obtuse_angle = [](const Vector3& p0, const Vector3& p1, const Vector3& p2) {
			const auto dot_sign = [](const Vector3 &a, const Vector3 &b)
			{
				Scalar t = a.dot(b);
				if (t > SCALAR_ZERO)
					return 1;
				if (t < -1 * SCALAR_ZERO)
					return -1;

				return dot_product_sign(a[0], a[1], a[2], b[0], b[1], b[2]);
			};

			if (dot_sign(p1 - p0, p1 - p2) < 0) return 1;
			if (dot_sign(p2 - p1, p2 - p0) < 0) return 2;
			if (dot_sign(p0 - p1, p0 - p2) < 0) return 0;
			return -1;

		};
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
		const auto get_corner_plane = [](const Vector3& p0, const Vector3& midp, const Vector3 &normal, const Scalar& distance,
			Vector3& plane0, Vector3& plane1, Vector3& plane2, const bool use_exact) {
			Scalar distance_small = distance * 1;// to be conservative to reduce numerical error, can set the Scalar as 0.999
			Vector3 direction = (p0 - midp).normalized();
			plane0 = p0 + direction * distance_small;
			plane1 = plane0 + normal;
			Vector3 origin = Vector3(0, 0, 0);
			Vector3 axis;
			if (use_exact) axis = accurate_cross_product_direction(midp, p0, origin, normal);
			else axis = direction.cross(normal).normalized();
			plane2 = plane0 + axis;
		};
		bool use_accurate_cross = false;
		Vector3 origin = Vector3(0, 0, 0);

		static const int c_face[6][3] = { {0, 1, 2}, {4, 7, 6}, {0, 3, 4}, {1, 0, 4}, {1, 5, 2}, {2, 6, 3} };

		halfspace.resize(m_faces.size());
		cornerlist.resize(m_faces.size());
		for (int i = 0; i < m_faces.size(); i++)
		{
			algorithms::get_tri_corners(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tmin, tmax);
			cornerlist[i][0] = tmin - bbox_offset[i];
			cornerlist[i][1] = tmax + bbox_offset[i];

			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			de = algorithms::is_triangle_degenerated(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);

			if (de == DEGENERATED_POINT)
			{
				//logger().debug("Envelope Triangle Degeneration- Point");
				for (int j = 0; j < 8; j++)
				{
					box[j] = m_ver[m_faces[i][0]] + boxorder[j] * tolerance[i];
				}
				halfspace[i].resize(6);
				for (int j = 0; j < 6; j++) {
					halfspace[i][j][0] = box[c_face[j][0]];
					halfspace[i][j][1] = box[c_face[j][1]];
					halfspace[i][j][2] = box[c_face[j][2]];
				}
				//cornerlist[i] = (get_bb_corners_8(box));


				continue;
			}
			if (de == DEGENERATED_SEGMENT)
			{
				//logger().debug("Envelope Triangle Degeneration- Segment");
				Scalar length1 = AB.dot(AB), length2 = AC.dot(AC), length3 = BC.dot(BC);
				if (length1 >= length2 && length1 >= length3)
				{
					algorithms::seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], tolerance[i], box);

				}
				if (length2 >= length1 && length2 >= length3)
				{
					algorithms::seg_cube(m_ver[m_faces[i][0]], m_ver[m_faces[i][2]], tolerance[i], box);

				}
				if (length3 >= length1 && length3 >= length2)
				{
					algorithms::seg_cube(m_ver[m_faces[i][1]], m_ver[m_faces[i][2]], tolerance[i], box);
				}
				halfspace[i].resize(6);
				for (int j = 0; j < 6; j++) {
					halfspace[i][j][0] = box[c_face[j][0]];
					halfspace[i][j][1] = box[c_face[j][1]];
					halfspace[i][j][2] = box[c_face[j][2]];
				}
				//cornerlist[i] = (get_bb_corners_8(box));


				continue;
			}
			if (de == NERLY_DEGENERATED)
			{
				//logger().debug("Envelope Triangle Degeneration- Nearly");

				normal = algorithms::accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);
				use_accurate_cross = true;
			}
			else
			{
				normal = AB.cross(AC).normalized();
			}
			halfspace[i].reserve(8);
			Vector3 normaldist = normal * tolerance[i];
			Vector3 edgedire, edgenormaldist;
			plane[0] = m_ver[m_faces[i][0]] + normaldist;
			plane[1] = m_ver[m_faces[i][1]] + normaldist;
			plane[2] = m_ver[m_faces[i][2]] + normaldist;
			halfspace[i].emplace_back(plane);// number 0

			plane[0] = m_ver[m_faces[i][0]] - normaldist;
			plane[1] = m_ver[m_faces[i][2]] - normaldist;
			plane[2] = m_ver[m_faces[i][1]] - normaldist;// order: 0, 2, 1
			halfspace[i].emplace_back(plane);// number 1

			int obtuse = get_triangle_obtuse_angle(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);


			edgedire = AB.normalized();
			if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(origin, edgedire, origin, normal)*tolerance[i];
			else edgenormaldist = edgedire.cross(normal).normalized()*tolerance[i];
			plane[0] = m_ver[m_faces[i][0]] + edgenormaldist;
			plane[1] = m_ver[m_faces[i][1]] + edgenormaldist;
			plane[2] = plane[0] + normal;
			halfspace[i].emplace_back(plane);// number 2


			if (obtuse != 1) {
				get_corner_plane(m_ver[m_faces[i][1]], (m_ver[m_faces[i][0]] + m_ver[m_faces[i][2]]) / 2, normal,
					tolerance[i], plane[0], plane[1], plane[2],use_accurate_cross);
				halfspace[i].emplace_back(plane);// number 3;

			}

			edgedire = BC.normalized();
			if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(origin, edgedire, origin, normal)*tolerance[i];
			else edgenormaldist = edgedire.cross(normal).normalized()*tolerance[i];
			plane[0] = m_ver[m_faces[i][1]] + edgenormaldist;
			plane[1] = m_ver[m_faces[i][2]] + edgenormaldist;
			plane[2] = plane[0] + normal;
			halfspace[i].emplace_back(plane);// number 4

			if (obtuse != 2) {
				get_corner_plane(m_ver[m_faces[i][2]], (m_ver[m_faces[i][0]] + m_ver[m_faces[i][1]]) / 2, normal,
					tolerance[i], plane[0], plane[1], plane[2],use_accurate_cross);
				halfspace[i].emplace_back(plane);// number 5;

			}

			edgedire = -AC.normalized();
			if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(origin, edgedire, origin, normal)*tolerance[i];
			else edgenormaldist = edgedire.cross(normal).normalized()*tolerance[i];
			plane[0] = m_ver[m_faces[i][2]] + edgenormaldist;
			plane[1] = m_ver[m_faces[i][0]] + edgenormaldist;
			plane[2] = plane[1] + normal;
			halfspace[i].emplace_back(plane);// number 6 

			if (obtuse != 0) {
				get_corner_plane(m_ver[m_faces[i][0]], (m_ver[m_faces[i][1]] + m_ver[m_faces[i][2]]) / 2, normal,
					tolerance[i], plane[0], plane[1], plane[2],use_accurate_cross);
				halfspace[i].emplace_back(plane);// number 7;

			}
			//std::cout << "envelope face nbr " << halfspace[i].size() << std::endl;

		}

	}

}
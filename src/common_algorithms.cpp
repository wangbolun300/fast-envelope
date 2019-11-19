#include<fastenvelope/common_algorithms.h>
#include<fastenvelope/Predicates.hpp>
#include<iostream>
#include<array>
#include<fastenvelope/Morton.h>
#include<fastenvelope/Logger.hpp>
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

		//Vector3 accurate_normal_vector(const Vector3 &p0, const Vector3 &p1,
		//	const Vector3 &q0, const Vector3 &q1)
		//{

		//	Rational p00r(p0[0]), p01r(p0[1]), p02r(p0[2]),
		//		p10r(p1[0]), p11r(p1[1]), p12r(p1[2]),
		//		q00r(q0[0]), q01r(q0[1]), q02r(q0[2]),
		//		q10r(q1[0]), q11r(q1[1]), q12r(q1[2]);
		//	Rational axr(p10r - p00r), ayr(p11r - p01r), azr(p12r - p02r),
		//		bxr(q10r - q00r), byr(q11r - q01r), bzr(q12r - q02r);
		//	Rational xr = ayr * bzr - azr * byr;
		//	Rational yr = azr * bxr - axr * bzr;
		//	Rational zr = axr * byr - ayr * bxr;//get the direction (x,y,z), now normalize
		//	int xsign, ysign, zsign;
		//	xsign = xr.get_sign();
		//	ysign = yr.get_sign();
		//	zsign = zr.get_sign();
		//	Rational ssumr = xr * xr + yr * yr + zr * zr;
		//	xr = xr * xr / ssumr;
		//	yr = yr * yr / ssumr;
		//	zr = zr * zr / ssumr;

		//	Scalar x, y, z;
		//	x = sqrt(xr.to_double())*xsign;
		//	y = sqrt(yr.to_double())*ysign;
		//	z = sqrt(zr.to_double())*zsign;
		//	return Vector3(x, y, z);
		//}
		Vector3 accurate_normal_vector(const Vector3 &p0, const Vector3 &p1,
			const Vector3 &p2)
		{
			Vector3 v;
			triangle_normal_exact(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],v[0],v[1],v[2]);
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
	
	
	void algorithms::halfspace_init(const std::vector<Vector3> &m_ver, const std::vector<Vector3i> &m_faces,
		std::vector<std::vector<std::array<Vector3, 3>>>& halfspace, std::vector<std::array<Vector3, 2>>& cornerlist, const Scalar &epsilon) {

		const auto dot_sign = [](const Vector3 &a, const Vector3 &b)
		{
			Scalar t = a.dot(b);
			if (t > SCALAR_ZERO)
				return 1;
			if (t < -1 * SCALAR_ZERO)
				return -1;
			std::cout << "need accurate dot sign" << std::endl;
			return dot_product_sign(a[0], a[1], a[2], b[0], b[1], b[2]);
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
		//	static const std::array<std::vector<int>, 8> p_facepoint = {
		//	{{0,1,2,3,4,5},
		//{8,7,6,11,10,9},
		//{7,1,0,6},
		//{2,1,7,8},
		//{3,2,8,9},
		//{4,3,9,10},
		//{4,10,11,5},
		//{6,0,5,11}}
		//	};

		//	static const std::array<std::array<int, 4>, 6> c_facepoint = {
		//			{
		//				{0,1,2,3},
		//		{4,7,6,5},
		//		{4,0,3,7},
		//		{1,0,4,5},
		//		{2,1,5,6},
		//		{3,2,6,7}
		//			}
		//	};

		Scalar tolerance = epsilon / sqrt(3);
		Scalar de;

		for (int i = 0; i < m_faces.size(); i++)
		{
			AB = m_ver[m_faces[i][1]] - m_ver[m_faces[i][0]];
			AC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][0]];
			BC = m_ver[m_faces[i][2]] - m_ver[m_faces[i][1]];
			de = algorithms::is_triangle_degenerated(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);

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
				cornerlist[i] = (get_bb_corners_8(box));


				continue;
			}
			if (de == NERLY_DEGENERATED)
			{
				logger().debug("Envelope Triangle Degeneration- Nearly");

				normal = algorithms::accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], m_ver[m_faces[i][2]]);

				vector1 = algorithms::accurate_normal_vector(m_ver[m_faces[i][0]], m_ver[m_faces[i][1]], normal + m_ver[m_faces[i][0]]);

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

	
	

}
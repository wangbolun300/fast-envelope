#pragma once

#include <fastenvelope/Types.hpp>
#include<fastenvelope/Predicates.hpp>
#include <vector>
#include <array>

namespace fastEnvelope {
	class AABB {
	private:
		static const int NOT_INTERSECTED = 2;
		static const int INTERSECTED = 1;
		static const int OUT_PRISM = 1;
		static const int IN_PRISM = 0;
		static const int CUT_COPLANAR = 4;
		static const int CUT_EMPTY = -1;
		static const int CUT_FACE = 3;

		static const int NOT_DEGENERATED = 0;
		static const int NERLY_DEGENERATED = 1;
		static const int DEGENERATED_SEGMENT = 2;
		static const int DEGENERATED_POINT = 3;
		static Vector2 to_2d(const Vector3 &p, int t)
		{
			return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
		}
		std::vector<std::array<Vector3, 2>> boxlist;
		size_t n_corners = -1;

		void init_envelope_boxes_recursive(
			const std::vector<std::array<Vector3, 2>> &cornerlist,
			int node_index,
			int b, int e);

		void facet_in_envelope_recursive(
			const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2,
			std::vector<unsigned int> &list,
			int n, int b, int e) const;
		static int is_triangle_degenerated(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2)
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
		static void get_tri_corners(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2, Vector3 &mint, Vector3 &maxt)
		{
			mint[0] = std::min(std::min(triangle0[0], triangle1[0]), triangle2[0]);
			mint[1] = std::min(std::min(triangle0[1], triangle1[1]), triangle2[1]);
			mint[2] = std::min(std::min(triangle0[2], triangle1[2]), triangle2[2]);
			maxt[0] = std::max(std::max(triangle0[0], triangle1[0]), triangle2[0]);
			maxt[1] = std::max(std::max(triangle0[1], triangle1[1]), triangle2[1]);
			maxt[2] = std::max(std::max(triangle0[2], triangle1[2]), triangle2[2]);
		}
		void bbd_searching_recursive(
			const Vector3 &bbd0, const Vector3 &bbd1,
			std::vector<unsigned int> &list,
			int n, int b, int e) const;

		static int envelope_max_node_index(int node_index, int b, int e);

		bool is_triangle_cut_bounding_box(const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, int index) const;
		bool is_bbd_cut_bounding_box(const Vector3 &bbd0, const Vector3 &bbd1, int index) const;
	public:
		void init_envelope(const std::vector<std::array<Vector3, 2>> &cornerlist);

		inline void facet_in_envelope(
			const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2,
			std::vector<unsigned int> &list) const
		{
			assert(n_corners >= 0);
			int de = is_triangle_degenerated(triangle0, triangle1, triangle2);
			if (de == DEGENERATED_SEGMENT || DEGENERATED_POINT) {
				Vector3 tmin, tmax;
				get_tri_corners(triangle0, triangle1, triangle2, tmin, tmax);
				bbd_finding_in_envelope(tmin, tmax, list);
			}
			else
				facet_in_envelope_recursive(triangle0, triangle1, triangle2, list, 1, 0, n_corners);
		}
		inline void bbd_finding_in_envelope(
			const Vector3 &bbd0, const Vector3 &bbd1,
			std::vector<unsigned int> &list) const
		{
			list.clear();
			assert(n_corners >= 0);
			bbd_searching_recursive(bbd0,bbd1, list, 1, 0, n_corners);
		}
	};
}
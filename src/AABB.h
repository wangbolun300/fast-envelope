#pragma once

#include <fastenvelope/Types.hpp>

#include <vector>
#include <array>

namespace fastEnvelope {
	class AABB {
	private:
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
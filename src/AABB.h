#pragma once

#include <fastenvelope/Types.hpp>
#include<fastenvelope/Predicates.hpp>
#include <vector>
#include <array>
#include<fastenvelope/common_algorithms.h>

namespace fastEnvelope {
	class AABB {
	private:
		
		
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
			int de = algorithms::is_triangle_degenerated(triangle0, triangle1, triangle2);
			if (de == DEGENERATED_SEGMENT || DEGENERATED_POINT) {
				Vector3 tmin, tmax;
				algorithms::get_tri_corners(triangle0, triangle1, triangle2, tmin, tmax);
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
#pragma once


#include <fastenvelope/Types.hpp>

#include <AABB.h>

#include <vector>
#include <array>
#include <cassert>
#include <memory>

namespace fastEnvelope {
	class AABB {
	private:
		bool use_aabbcc_ = false;

		std::vector<std::array<Vector3, 2>> boxlist;
		size_t n_corners = -1;

		std::unique_ptr<aabb::Tree> aabb;

		void init_envelope_boxes_recursive(
			const std::vector<std::array<Vector3, 2>> &cornerlist,
			int node_index,
			int b, int e);

		void facet_in_envelope_recursive(
			const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2,
			std::vector<unsigned int> &list,
			int n, int b, int e) const;

		static int envelope_max_node_index(int node_index, int b, int e);

	public:
		void init_envelope(const std::vector<std::array<Vector3, 2>> &cornerlist, bool use_aabbcc);

		inline void facet_in_envelope(
			const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2,
			std::vector<unsigned int> &list) const
		{
			if (use_aabbcc_)
			{
				assert(aabb);
				std::vector<double> minv = {
					std::min(std::min(triangle0[0], triangle1[0]), triangle2[0]),
					std::min(std::min(triangle0[1], triangle1[1]), triangle2[1]),
					std::min(std::min(triangle0[2], triangle1[2]), triangle2[2])};
				std::vector<double> maxv = {
					std::max(std::max(triangle0[0], triangle1[0]), triangle2[0]),
					std::max(std::max(triangle0[1], triangle1[1]), triangle2[1]),
					std::max(std::max(triangle0[2], triangle1[2]), triangle2[2])};
				aabb::AABB box(minv, maxv);

				list = aabb->query(box);
			}
			else
			{
				assert(n_corners >= 0);
				facet_in_envelope_recursive(triangle0, triangle1, triangle2, list, 1, 0, n_corners);
			}
		}

	};
}
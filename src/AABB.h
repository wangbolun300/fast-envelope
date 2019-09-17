#pragma once


#include <fastenvelope/Types.hpp>

#include <vector>
#include <array>
#include <cassert>

namespace fastEnvelope {
	class AABB {
	private:
		std::vector<std::array<fastEnvelope::Vector3, 2>> boxlist;
		size_t n_corners = -1;
	
		void init_envelope_boxes_recursive(
			const std::vector<std::array<fastEnvelope::Vector3, 2>> &cornerlist,
			int node_index,
			int b, int e);


		void facet_in_envelope_recursive(
			const fastEnvelope::Vector3& triangle0, const fastEnvelope::Vector3& triangle1, const fastEnvelope::Vector3& triangle2,
			std::vector<int>& list,
			int n, int b, int e
		) const;

		static int envelope_max_node_index(int node_index, int b, int e);

	public:
		void init_envelope(const std::vector<std::array<fastEnvelope::Vector3, 2>> &cornerlist);

		inline void facet_in_envelope(
			const fastEnvelope::Vector3& triangle0, const fastEnvelope::Vector3& triangle1, const fastEnvelope::Vector3& triangle2,
			std::vector<int>& list
		) const {
			assert(n_corners >= 0);
			facet_in_envelope_recursive(triangle0, triangle1, triangle2, list, 1, 0, n_corners);
		}

	};
}
#include <fastenvelope/AABB.h>
#include <fastenvelope/FastEnvelope.h>

namespace fastEnvelope {
	void AABB::init_envelope_boxes_recursive(
		const std::vector<std::array<fastEnvelope::Vector3, 2>> &cornerlist,
		int node_index,
		int b, int e)
	{
		assert(b != e);
		assert(node_index < boxlist.size());

		
		if (b + 1 == e) {
			boxlist[node_index] = cornerlist[b];
			return;
		}
		int m = b + (e - b) / 2;
		int childl = 2 * node_index;
		int childr = 2 * node_index + 1;

		assert(childl < boxlist.size());
		assert(childr < boxlist.size());

		init_envelope_boxes_recursive(cornerlist, childl, b, m);
		init_envelope_boxes_recursive(cornerlist, childr, m, e);

		assert(childl < boxlist.size());
		assert(childr < boxlist.size());
		for (int c = 0; c < 3; ++c) {
			boxlist[node_index][0][c] = std::min(boxlist[childl][0][c], boxlist[childr][0][c]);
			boxlist[node_index][1][c] = std::max(boxlist[childl][1][c], boxlist[childr][1][c]);
		}
	}


	void AABB::facet_in_envelope_recursive(
		const fastEnvelope::Vector3& triangle0, const fastEnvelope::Vector3& triangle1, const fastEnvelope::Vector3& triangle2,
		std::vector<int>& list,
		int n, int b, int e
	) const {
		assert(e != b);
		
		assert(n < boxlist.size());
		bool cut = fastEnvelope::FastEnvelope::is_triangle_cut_bounding_box(triangle0, triangle1, triangle2, boxlist[n][0], boxlist[n][1]);

		if (cut == false) return;

		// Leaf case
		if (e == b + 1) {
			list.emplace_back(b);
			return;
		}

		int m = b + (e - b) / 2;
		int childl = 2 * n;
		int childr = 2 * n + 1;

		//assert(childl < boxlist.size());
		//assert(childr < boxlist.size());



		// Traverse the "nearest" child first, so that it has more chances
		// to prune the traversal of the other child.
		facet_in_envelope_recursive(
			triangle0, triangle1, triangle2, list,
			childl, b, m
		);
		facet_in_envelope_recursive(
			triangle0, triangle1, triangle2, list,
			childr, m, e
		);
	}

	int AABB::envelope_max_node_index(int node_index, int b, int e) {
		assert(e > b);
		if (b + 1 == e) {
			return node_index;
		}
		int m = b + (e - b) / 2;
		int childl = 2 * node_index;
		int childr = 2 * node_index + 1;
		return std::max(
			envelope_max_node_index(childl, b, m),
			envelope_max_node_index(childr, m, e)
		);
	}

	void AABB::init_envelope(const std::vector<std::array<fastEnvelope::Vector3, 2>> &cornerlist)
	{
		n_corners = cornerlist.size();

		boxlist.resize(
			envelope_max_node_index(
				1, 0, n_corners
			) + 1 // <-- this is because size == max_index + 1 !!!
		);


		init_envelope_boxes_recursive(cornerlist, 1, 0, n_corners);
	}
}
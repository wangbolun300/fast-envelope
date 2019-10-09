#include <fastenvelope/AABB.h>
#include <fastenvelope/FastEnvelope.h>

#include <cassert>

namespace fastEnvelope {
	namespace {
		bool box_box_intersection(const Vector3 &min1, const Vector3 &max1, const Vector3 &min2, const Vector3 &max2)
		{
			if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
				return 0;
			if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
				return 0;
			return 1;
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
	}


	void AABB::init_envelope_boxes_recursive(
		const std::vector<std::array<Vector3, 2>> &cornerlist,
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
		const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2,
		std::vector<unsigned int> &list,
		int n, int b, int e) const
	{
		assert(e != b);

		assert(n < boxlist.size());
		bool cut = is_triangle_cut_bounding_box(triangle0, triangle1, triangle2, n);

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

	void AABB::init_envelope(const std::vector<std::array<Vector3, 2>> &cornerlist)
	{
		n_corners = cornerlist.size();

		boxlist.resize(
			envelope_max_node_index(
				1, 0, n_corners) +
			1 // <-- this is because size == max_index + 1 !!!
		);

		init_envelope_boxes_recursive(cornerlist, 1, 0, n_corners);
	}

	
#include<fastenvelope/Predicates.hpp>
	bool AABB::is_triangle_cut_bounding_box(
		const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, int index) const
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];
		Vector3 tmin, tmax;
		get_tri_corners(tri0, tri1, tri2, tmin, tmax);
		bool cut= box_box_intersection(tmin, tmax, bmin, bmax);
		if (cut) {
			
			std::array<Vector2,3> t;
			std::array<Vector2,4> m;
			int o0, o1, o2, o3, ori;
			for (int i = 0; i < 3; i++) {
				t[0] = to_2d(tri0, i);
				t[1] = to_2d(tri1, i);
				t[2] = to_2d(tri2, i);

				m[0] = to_2d(bmin, i);
				m[1] = to_2d(bmax, i);
				m[2][0] = m[0][0]; m[2][1] = m[1][1];
				m[3][0] = m[1][0]; m[3][1] = m[0][1];
		
				for (int j = 0; j < 3; j++) {
					o0 = fastEnvelope::Predicates::orient_2d(m[0], t[j % 3], t[(j + 1) % 3]);
					o1 = fastEnvelope::Predicates::orient_2d(m[1], t[j % 3], t[(j + 1) % 3]);
					o2 = fastEnvelope::Predicates::orient_2d(m[2], t[j % 3], t[(j + 1) % 3]);
					o3 = fastEnvelope::Predicates::orient_2d(m[3], t[j % 3], t[(j + 1) % 3]);
					ori = fastEnvelope::Predicates::orient_2d(t[(j + 2) % 3], t[j % 3], t[(j + 1) % 3]);
					if (ori*o0 <= 0 && ori*o1 <= 0 && ori*o2 <= 0 && ori*o3 <= 0) return false;
				}
			}
		}
		else {
			return false;
		}
		return true;
	}
}
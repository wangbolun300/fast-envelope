#include <fastenvelope/AABB.h>
#include <fastenvelope/FastEnvelope.h>
#include <cassert>

namespace fastEnvelope {
	


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
	
	void AABB::bbd_searching_recursive(
		const Vector3 &bbd0, const Vector3 &bbd1,
		std::vector<unsigned int> &list,
		int n, int b, int e) const
	{
		assert(e != b);

		assert(n < boxlist.size());
		bool cut = is_bbd_cut_bounding_box(bbd0, bbd1, n);

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
		bbd_searching_recursive(
			bbd0, bbd1, list,
			childl, b, m
		);
		bbd_searching_recursive(
			bbd0, bbd1, list,
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
	




	bool AABB::is_triangle_cut_bounding_box(
		const Vector3 &tri0, const Vector3 &tri1, const Vector3 &tri2, int index) const
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];
		Vector3 tmin, tmax;
		
		algorithms::get_tri_corners(tri0, tri1, tri2, tmin, tmax);
		bool cut=algorithms:: box_box_intersection(tmin, tmax, bmin, bmax);
		if (cut == false) return false;
		
		if (cut) {
			
			std::array<Vector2,3> tri;
			std::array<Vector2,4> mp;
			int o0, o1, o2, o3, ori;
			for (int i = 0; i < 3; i++) {
				tri[0] = algorithms::to_2d(tri0, i);
				tri[1] = algorithms::to_2d(tri1, i);
				tri[2] = algorithms::to_2d(tri2, i);

				mp[0] = algorithms::to_2d(bmin, i);
				mp[1] = algorithms::to_2d(bmax, i);
				mp[2][0] = mp[0][0]; mp[2][1] = mp[1][1];
				mp[3][0] = mp[1][0]; mp[3][1] = mp[0][1];
		
				for (int j = 0; j < 3; j++) {
					o0 = fastEnvelope::Predicates::orient_2d(mp[0], tri[j % 3], tri[(j + 1) % 3]);
					o1 = fastEnvelope::Predicates::orient_2d(mp[1], tri[j % 3], tri[(j + 1) % 3]);
					o2 = fastEnvelope::Predicates::orient_2d(mp[2], tri[j % 3], tri[(j + 1) % 3]);
					o3 = fastEnvelope::Predicates::orient_2d(mp[3], tri[j % 3], tri[(j + 1) % 3]);
					ori = fastEnvelope::Predicates::orient_2d(tri[(j + 2) % 3], tri[j % 3], tri[(j + 1) % 3]);
					if (ori == 0) continue;
					if (ori*o0 <= 0 && ori*o1 <= 0 && ori*o2 <= 0 && ori*o3 <= 0) return false;
				}
			}
		}

		return cut;
	}
	bool AABB::is_bbd_cut_bounding_box(
		const Vector3 &bbd0, const Vector3 &bbd1, int index) const
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];

		
		return algorithms::box_box_intersection(bbd0, bbd1, bmin, bmax);
	}
}

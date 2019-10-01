#pragma once

#include <fastenvelope/mesh_AABB.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <memory>
#include<fastenvelope/Types.hpp>
#include<array>
namespace floatTetWild {
	using namespace fastEnvelope;
	void sample_triangle(const std::array<Vector3, 3>& vs, std::vector<GEO::vec3>& ps, Scalar sampling_dist) {
		Scalar sqrt3_2 = std::sqrt(3) / 2;

		std::array<Scalar, 3> ls;
		for (int i = 0; i < 3; i++) {
			ls[i] = (vs[i] - vs[(i + 1) % 3]).squaredNorm();
		}
		auto min_max = std::minmax_element(ls.begin(), ls.end());
		int min_i = min_max.first - ls.begin();
		int max_i = min_max.second - ls.begin();
		Scalar N = sqrt(ls[max_i]) / sampling_dist;
		if (N <= 1) {
			for (int i = 0; i < 3; i++)
				ps.push_back(GEO::vec3(vs[i][0], vs[i][1], vs[i][2]));
			return;
		}
		if (N == int(N))
			N -= 1;

		GEO::vec3 v0(vs[max_i][0], vs[max_i][1], vs[max_i][2]);
		GEO::vec3 v1(vs[(max_i + 1) % 3][0], vs[(max_i + 1) % 3][1], vs[(max_i + 1) % 3][2]);
		GEO::vec3 v2(vs[(max_i + 2) % 3][0], vs[(max_i + 2) % 3][1], vs[(max_i + 2) % 3][2]);

		GEO::vec3 n_v0v1 = GEO::normalize(v1 - v0);
		for (int n = 0; n <= N; n++) {
			ps.push_back(v0 + n_v0v1 * sampling_dist * n);
		}
		ps.push_back(v1);

		Scalar h = GEO::distance(GEO::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
		int M = h / (sqrt3_2 * sampling_dist);
		if (M < 1) {
			ps.push_back(v2);
			return;
		}

		GEO::vec3 n_v0v2 = GEO::normalize(v2 - v0);
		GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
		Scalar tan_v0, tan_v1, sin_v0, sin_v1;
		sin_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / (GEO::distance(v0, v2) * GEO::distance(v0, v1));
		tan_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / GEO::dot((v2 - v0), (v1 - v0));
		tan_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / GEO::dot((v2 - v1), (v0 - v1));
		sin_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / (GEO::distance(v1, v2) * GEO::distance(v0, v1));

		for (int m = 1; m <= M; m++) {
			int n = sqrt3_2 / tan_v0 * m + 0.5;
			int n1 = sqrt3_2 / tan_v0 * m;
			if (m % 2 == 0 && n == n1) {
				n += 1;
			}
			GEO::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
			GEO::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
			if (GEO::distance(v0_m, v1_m) <= sampling_dist)
				break;

			Scalar delta_d = ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * sampling_dist;
			GEO::vec3 v = v0_m + delta_d * n_v0v1;
			int N1 = GEO::distance(v, v1_m) / sampling_dist;
			//        ps.push_back(v0_m);
			for (int i = 0; i <= N1; i++) {
				ps.push_back(v + i * n_v0v1 * sampling_dist);
			}
			//        ps.push_back(v1_m);
		}
		ps.push_back(v2);

		//sample edges
		N = sqrt(ls[(max_i + 1) % 3]) / sampling_dist;
		if (N > 1) {
			if (N == int(N))
				N -= 1;
			GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
			for (int n = 1; n <= N; n++) {
				ps.push_back(v1 + n_v1v2 * sampling_dist * n);
			}
		}

		N = sqrt(ls[(max_i + 2) % 3]) / sampling_dist;
		if (N > 1) {
			if (N == int(N))
				N -= 1;
			GEO::vec3 n_v2v0 = GEO::normalize(v0 - v2);
			for (int n = 1; n <= N; n++) {
				ps.push_back(v2 + n_v2v0 * sampling_dist * n);
			}
		}
	}





	class AABBWrapper {
	private:
		GEO::Mesh b_mesh;
		GEO::Mesh tmp_b_mesh;
		const GEO::Mesh &sf_mesh;

		std::shared_ptr<GEO::MeshFacetsAABBWithEps> b_tree;
		std::shared_ptr<GEO::MeshFacetsAABBWithEps> tmp_b_tree;
		GEO::MeshFacetsAABBWithEps sf_tree;

		void init_b_mesh(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces);

	public:
		AABBWrapper(const GEO::Mesh &sf_mesh) : sf_mesh(sf_mesh), sf_tree(sf_mesh) {}

		inline Scalar get_sf_diag() const { return GEO::bbox_diagonal(sf_mesh); }

		void init_b_mesh_and_tree(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces) {
			init_b_mesh(input_vertices, input_faces);
			b_tree = std::make_shared<GEO::MeshFacetsAABBWithEps>(b_mesh);
		}


		void init_tmp_b_mesh_and_tree(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
			const std::vector<std::array<int, 2>>& b_edges);

		inline Scalar project_to_sf(Vector3 &p) const {
			GEO::vec3 geo_p(p[0], p[1], p[2]);
			GEO::vec3 nearest_p;
			double sq_dist = std::numeric_limits<double>::max(); //??
			sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
			p[0] = nearest_p[0];
			p[1] = nearest_p[1];
			p[2] = nearest_p[2];

			return sq_dist;
		}

		inline Scalar project_to_sf(const GEO::vec3 geo_p) const {
			GEO::vec3 nearest_p;
			double sq_dist = std::numeric_limits<double>::max(); //??
			sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
			return sq_dist;
		}

		inline Scalar project_to_b(Vector3 &p) const {
			GEO::vec3 geo_p(p[0], p[1], p[2]);
			GEO::vec3 nearest_p;
			double sq_dist = std::numeric_limits<double>::max(); //?
			b_tree->nearest_facet(geo_p, nearest_p, sq_dist);
			p[0] = nearest_p[0];
			p[1] = nearest_p[1];
			p[2] = nearest_p[2];

			return sq_dist;
		}

		inline Scalar project_to_tmp_b(Vector3 &p) const {
			GEO::vec3 geo_p(p[0], p[1], p[2]);
			GEO::vec3 nearest_p;
			double sq_dist = std::numeric_limits<double>::max(); //?
			tmp_b_tree->nearest_facet(geo_p, nearest_p, sq_dist);
			p[0] = nearest_p[0];
			p[1] = nearest_p[1];
			p[2] = nearest_p[2];

			return sq_dist;
		}

		inline Scalar project_to_b(const GEO::vec3 geo_p) const {
			GEO::vec3 nearest_p;
			double sq_dist = std::numeric_limits<double>::max(); //??
			b_tree->nearest_facet(geo_p, nearest_p, sq_dist);
			return sq_dist;
		}

		inline bool is_out_sf_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
			GEO::index_t prev_facet = GEO::NO_FACET) const {
			GEO::vec3 nearest_point;
			double sq_dist = std::numeric_limits<double>::max();

			for (const GEO::vec3 &current_point : ps) {
				if (prev_facet != GEO::NO_FACET) {
					get_point_facet_nearest_point(sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					sf_tree.facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					return true;
				}
			}

			return false;
		}

		inline Scalar dist_sf_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
			GEO::index_t prev_facet = GEO::NO_FACET) const {
			GEO::vec3 nearest_point;
			double sq_dist = std::numeric_limits<double>::max();

			for (const GEO::vec3 &current_point : ps) {
				if (prev_facet != GEO::NO_FACET) {
					get_point_facet_nearest_point(sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					sf_tree.facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					return sq_dist;
				}
			}

			return 0;
		}

		inline bool is_out_b_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
			GEO::index_t prev_facet = GEO::NO_FACET) const {
			GEO::vec3 nearest_point;
			double sq_dist = std::numeric_limits<double>::max();

			for (const GEO::vec3 &current_point : ps) {
				if (prev_facet != GEO::NO_FACET) {
					get_point_facet_nearest_point(b_mesh, current_point, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					b_tree->facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					return true;
				}
			}

			return false;
		}

		inline bool is_out_tmp_b_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
			GEO::index_t prev_facet = GEO::NO_FACET) const {
			GEO::vec3 nearest_point;
			double sq_dist = std::numeric_limits<double>::max();

			for (const GEO::vec3 &current_point : ps) {
				if (prev_facet != GEO::NO_FACET) {
					get_point_facet_nearest_point(tmp_b_mesh, current_point, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					tmp_b_tree->facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
				}
				if (Scalar(sq_dist) > eps_2) {
					return true;
				}
			}

			return false;
		}

		inline bool is_out_sf_envelope(const Vector3& p, const Scalar eps_2, GEO::index_t& prev_facet) const {
			GEO::vec3 nearest_p;
			double sq_dist;
			GEO::vec3 geo_p(p[0], p[1], p[2]);
			prev_facet = sf_tree.facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

			if (Scalar(sq_dist) > eps_2)
				return true;
			return false;
		}

		inline bool is_out_b_envelope(const Vector3& p, const Scalar eps_2, GEO::index_t& prev_facet) const {
			GEO::vec3 nearest_p;
			double sq_dist;
			GEO::vec3 geo_p(p[0], p[1], p[2]);
			prev_facet = b_tree->facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

			if (Scalar(sq_dist) > eps_2)
				return true;
			return false;
		}

		inline bool is_out_tmp_b_envelope(const Vector3& p, const Scalar eps_2, GEO::index_t& prev_facet) const {
			GEO::vec3 nearest_p;
			double sq_dist;
			GEO::vec3 geo_p(p[0], p[1], p[2]);
			prev_facet = tmp_b_tree->facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

			if (Scalar(sq_dist) > eps_2)
				return true;
			return false;
		}
	};
}

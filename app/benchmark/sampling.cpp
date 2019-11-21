#include<fastenvelope/sampling.h>
#include<fastenvelope/common_algorithms.h>
#include <fastenvelope/Rational.hpp>
#include<fastenvelope/Logger.hpp>
namespace fastEnvelope {
	sampling::sampling(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps) {
		Vector3 min, max;
		min = m_ver.front();
		max = m_ver.front();

		for (size_t j = 0; j < m_ver.size(); j++)
		{
			for (int i = 0; i < 3; i++)
			{
				min(i) = std::min(min(i), m_ver[j](i));
				max(i) = std::max(max(i), m_ver[j](i));
			}
		}



		const Scalar bbd = (max - min).norm();
		const Scalar epsilon = bbd * eps; //eps*bounding box diagnal


		std::vector<Vector3i> faces_new;


		algorithms::resorting(m_ver, m_faces, faces_new);//resort the facets order
		//----




		algorithms::halfspace_init(m_ver, faces_new, halfspace, cornerlist, epsilon);



		tree.init(cornerlist);


		//initializing types
	}

	
	bool sampling::sample_triangle_outside(const std::array<Vector3, 3> &triangle, const int &pieces) const
	{
		const auto triangle_sample_segment=[](const std::array<Vector3, 3> &triangle, Vector3 &ps, const int &pieces, const int &nbr)
		{
			int t = pieces - 1;
			if (triangle[1] - triangle[0] == Vector3(0, 0, 0))
			{

				ps = (triangle[0] + (triangle[2] - triangle[0]) * nbr / t);

				return;
			}
			if (triangle[2] - triangle[0] == Vector3(0, 0, 0))
			{

				ps = (triangle[0] + (triangle[1] - triangle[0]) * nbr / t);

				return;
			}
			if (triangle[2] - triangle[1] == Vector3(0, 0, 0))
			{

				ps = (triangle[0] + (triangle[1] - triangle[0]) * nbr / t);

				return;
			}

			Scalar d1 = (triangle[1] - triangle[0]).norm(), d2 = (triangle[2] - triangle[0]).norm(), d3 = (triangle[1] - triangle[2]).norm();
			if (d1 >= d2 && d1 >= d3)
			{
				ps = (triangle[0] + (triangle[1] - triangle[0]) * nbr / t);

				return;
			}
			if (d2 >= d1 && d2 >= d3)
			{
				ps = (triangle[0] + (triangle[2] - triangle[0]) * nbr / t);

				return;
			}
			if (d3 >= d1 && d3 >= d2)
			{
				ps = (triangle[1] + (triangle[2] - triangle[1]) * nbr / t);

				return;
			}
		};

		const auto triangle_sample_normal = [](const std::array<Vector3, 3> &triangle, Vector3 &ps, const int &pieces, const int &nbr1, const int &nbr2)
		{
			int l1s = pieces - 1; //
			Vector3 p1 = triangle[0] + (triangle[1] - triangle[0]) * nbr1 / l1s, d = (triangle[2] - triangle[1]) / l1s;
			ps = p1 + d * nbr2;
		};

		const auto triangle_sample_normal_rational = [](const std::array<Vector3, 3> &triangle, Rational &ps0, Rational &ps1, Rational &ps2, const int &pieces, const int &nbr1, const int &nbr2)
		{
			int l1s = pieces - 1; //
			Rational t00(triangle[0][0]), t01(triangle[0][1]), t02(triangle[0][2]), t10(triangle[1][0]), t11(triangle[1][1]),
				t12(triangle[1][2]), t20(triangle[2][0]), t21(triangle[2][1]), t22(triangle[2][2]), nbr1r(nbr1), nbr2r(nbr2), l1sr(l1s);

			Rational p0 = t00 + (t10 - t00) * nbr1r / l1sr, d0 = (t20 - t10) / l1sr;
			Rational p1 = t01 + (t11 - t01) * nbr1r / l1sr, d1 = (t21 - t11) / l1sr;
			Rational p2 = t02 + (t12 - t02) * nbr1r / l1sr, d2 = (t22 - t12) / l1sr;
			ps0 = p0 + d0 * nbr2;
			ps1 = p1 + d1 * nbr2;
			ps2 = p2 + d2 * nbr2;
		};

		bool out;
		Vector3 point;
		std::vector<unsigned int> querylist;
		tree.triangle_find_bbox(triangle[0], triangle[1], triangle[2], querylist);
		int jump = -1;
		if (querylist.size() == 0)
			return 1;

		int deg = algorithms::is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

		if (deg == DEGENERATED_POINT)
		{
			out = point_out_prism(triangle[0], querylist, jump);
			if (out == true)
			{

				return 1;
			}
			return 0;
		}
		if (deg == DEGENERATED_SEGMENT)
		{
			for (int i = 0; i < pieces; i++)
			{
				triangle_sample_segment(triangle, point, pieces, i);
				out = point_out_prism(point, querylist, jump);
				if (out == true)
				{

					return 1;
				}
			}
			return 0;
		}

		for (int i = 0; i < pieces; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				triangle_sample_normal(triangle, point, pieces, i, j);
				out = point_out_prism(point, querylist, jump);
				if (out == true)
				{

					return 1;
				}
			}
		}

		return 0;
	}
	bool sampling::point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const
	{

		int ori;

		for (int i = 0; i < prismindex.size(); i++)
		{
			if (prismindex[i] == jump)
				continue;

			for (int j = 0; j < halfspace[prismindex[i]].size(); j++)
			{

				ori = Predicates::orient_3d(halfspace[prismindex[i]][j][0], halfspace[prismindex[i]][j][1], halfspace[prismindex[i]][j][2], point);
				if (ori == -1 || ori == 0)
				{
					break;
				}
				if (j == halfspace[prismindex[i]].size() - 1)
				{

					return false;
				}
			}

		}

		return true;
	}
}

	
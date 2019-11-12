#include<fastenvelope/common_algorithms.h>
#include<fastenvelope/Predicates.hpp>
#include<fastenvelope/Rational.hpp>
#include<array>
namespace fastEnvelope {
	namespace algorithms {
		int seg_cut_plane(const Vector3 &seg0, const Vector3 &seg1, const Vector3 &t0, const Vector3 &t1, const Vector3 &t2)
		{
			int o1, o2;
			o1 = Predicates::orient_3d(seg0, t0, t1, t2);
			o2 = Predicates::orient_3d(seg1, t0, t1, t2);
			int op = o1 * o2;
			if (op >= 0)
			{
				return FE_CUT_COPLANAR; //in fact, coplanar and not cut this plane
			}
			return FE_CUT_FACE;
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
		bool box_box_intersection(const Vector3 &min1, const Vector3 &max1, const Vector3 &min2, const Vector3 &max2)//TDOO
		{
			if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
				return 0;
			if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
				return 0;
			return 1;
		}
		Vector2 to_2d(const Vector3 &p, int t)

		{
			return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
		}

	}
	int algorithms::is_triangle_degenerated(const Vector3 &triangle0, const Vector3 &triangle1, const Vector3 &triangle2)
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

	Vector3 algorithms::accurate_normal_vector(const Vector3 &p0, const Vector3 &p1,
		const Vector3 &q0, const Vector3 &q1)
	{

		Rational p00r(p0[0]), p01r(p0[1]), p02r(p0[2]),
			p10r(p1[0]), p11r(p1[1]), p12r(p1[2]),
			q00r(q0[0]), q01r(q0[1]), q02r(q0[2]),
			q10r(q1[0]), q11r(q1[1]), q12r(q1[2]);
		Rational axr(p10r - p00r), ayr(p11r - p01r), azr(p12r - p02r),
			bxr(q10r - q00r), byr(q11r - q01r), bzr(q12r - q02r);
		Rational xr = ayr * bzr - azr * byr;
		Rational yr = azr * bxr - axr * bzr;
		Rational zr = axr * byr - ayr * bxr;//get the direction (x,y,z), now normalize
		int xsign, ysign, zsign;
		xsign = xr.get_sign();
		ysign = yr.get_sign();
		zsign = zr.get_sign();
		Rational ssumr = xr * xr + yr * yr + zr * zr;
		xr = xr * xr / ssumr;
		yr = yr * yr / ssumr;
		zr = zr * zr / ssumr;

		Scalar x, y, z;
		x = sqrt(xr.to_double())*xsign;
		y = sqrt(yr.to_double())*ysign;
		z = sqrt(zr.to_double())*zsign;
		return Vector3(x, y, z);

	}
}
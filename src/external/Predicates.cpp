#include <fastenvelope/Predicates.hpp>

#ifdef ENVELOPE_WITH_GEO
#ifdef ENVELOPE_WITH_GEO_PSM
#include <Predicates_psm.h>
#else
#include <geogram/delaunay/delaunay_3d.h>
#endif
#else
#include <igl/predicates/predicates.h>
#endif

namespace fastEnvelope
{
	namespace
	{
		void init()
		{
#ifndef ENVELOPE_WITH_GEO
			igl::predicates::exactinit();
#endif
		}
	}

	const int Predicates::ORI_POSITIVE;
	const int Predicates::ORI_ZERO;
	const int Predicates::ORI_NEGATIVE;
	const int Predicates::ORI_UNKNOWN;

	int Predicates::orient_3d(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4)
	{
		init();
#ifdef ENVELOPE_WITH_GEO
		const int result = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
		const int result = (int)igl::predicates::orient3d(p1, p2, p3, p4);
#endif

		//		if (result > SCALAR_ZERO)
		//			return ORI_POSITIVE;
		//		else if (result < -SCALAR_ZERO)
		//			return ORI_NEGATIVE;
		//		else
		//			return ORI_ZERO;

		if (result > 0)
			return ORI_POSITIVE;
		else if (result < 0)
			return ORI_NEGATIVE;
		else
			return ORI_ZERO;
	}

	// int Predicates::orient_3d_tolerance(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p) {
	// 	const Scalar result = orient3d(p1.data(), p2.data(), p3.data(), p.data());
	// 	if (result == 0)
	// 		return ORI_ZERO;

	// 	Vector3 n = ((p2 - p3).cross(p1 - p3)).normalized();
	// 	Scalar d = std::abs(n.dot(p - p1));
	// 	if (d <= SCALAR_ZERO)
	// 		return Predicates::ORI_ZERO;

	// 	if (result > 0)
	// 		return ORI_POSITIVE;
	// 	else
	// 		return ORI_NEGATIVE;
	// }

	Scalar Predicates::orient_3d_volume(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4)
	{
		init();
#ifdef ENVELOPE_WITH_GEO
		const int ori = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
		const int ori = (int)igl::predicates::orient3d(p1, p2, p3, p4);
#endif
		if (ori <= 0)
			return ori;
		else
			return (p1 - p4).dot((p2 - p4).cross(p3 - p4)) / 6;
	}

	int Predicates::orient_2d(const Vector2 &p1, const Vector2 &p2, const Vector2 &p3)
	{
		init();
#ifdef ENVELOPE_WITH_GEO
		const int result = -GEO::PCK::orient_2d(p1.data(), p2.data(), p3.data());
#else
		const int result = (int)igl::predicates::orient2d(p1, p2, p3);
#endif
		if (result > 0)
			return ORI_POSITIVE;
		else if (result < 0)
			return ORI_NEGATIVE;
		else
			return ORI_ZERO;
	}

} // namespace fastEnvelope

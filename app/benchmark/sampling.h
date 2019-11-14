#pragma once
#include<fastenvelope/Types.hpp>
#include<vector>
#include<array>
#include<fastenvelope/AABB.h>

namespace fastEnvelope {
	class sampling {
	public:
		sampling(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps);
		bool sample_triangle_outside(const std::array<Vector3, 3> &triangle, const int &pieces) const;
	private:
		std::vector<std::array<Vector3, 2>> cornerlist;
		std::vector<std::vector<std::array<Vector3, 3>>> halfspace;
		AABB tree;


		
		bool point_out_prism(const Vector3 &point, const std::vector<unsigned int> &prismindex, const int &jump) const;
		
	};


}
#pragma once

#include <fastenvelope/Types.hpp>

#include <geogram/mesh/mesh_geometry.h>

#include <vector>

namespace fastEnvelope
{
	class MeshIO
	{
	public:
		static bool load_mesh(const std::string &path, std::vector<Vector3> &points, std::vector<Vector3i> &faces, GEO::Mesh& input);
	};
}

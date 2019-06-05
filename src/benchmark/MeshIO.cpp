#include "MeshIO.hpp"

#include <fastenvelope/Logger.hpp>


#include <igl/Timer.h>

#include <igl/write_triangle_mesh.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>


#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_reorder.h>


#include <numeric>


namespace fastEnvelope
{
	bool MeshIO::load_mesh(const std::string &path, std::vector<Vector3> &points, std::vector<Vector3i> &faces, GEO::Mesh& input)
	{
		logger().debug("Loading mesh at {}...", path);
		igl::Timer timer; timer.start();

		input.clear(false,false);
		std::vector<int> flags;

		const bool ok = GEO::mesh_load(path, input);

		if(!ok)
			return false;

		bool is_valid = (flags.size() == input.facets.nb());
		if(is_valid)
		{
			assert(flags.size() == input.facets.nb());
			GEO::Attribute<int> bflags(input.facets.attributes(), "bbflags");
			for (int index = 0; index < (int) input.facets.nb(); ++index) {
				bflags[index] = flags[index];
			}
		}

		if(!input.facets.are_simplices()) {
			mesh_repair(
				input,
				GEO::MeshRepairMode(GEO::MESH_REPAIR_TRIANGULATE | GEO::MESH_REPAIR_QUIET)
				);
		}

		// #ifdef FLOAT_TETWILD_USE_FLOAT
		// 		input.vertices.set_single_precision();
		// #else
		// 		input.vertices.set_double_precision();
		// #endif

		GEO::mesh_reorder(input, GEO::MESH_ORDER_MORTON);

		if(is_valid)
		{
			flags.clear();
			flags.resize(input.facets.nb());
			GEO::Attribute<int> bflags(input.facets.attributes(), "bbflags");
			for (int index = 0; index < (int) input.facets.nb(); ++index) {
				flags[index] = bflags[index];
			}
		}

		points.resize(input.vertices.nb());
		for(size_t i=0; i<points.size(); i++)
			points[i]<<(input.vertices.point(i))[0], (input.vertices.point(i))[1], (input.vertices.point(i))[2];

		faces.resize(input.facets.nb());
		for(size_t i=0; i<faces.size(); i++)
			faces[i]<<input.facets.vertex(i, 0), input.facets.vertex(i, 1), input.facets.vertex(i, 2);


		return ok;
	}
}
set(SOURCES
	predicates.c
	Predicates.cpp
	Predicates.hpp

	mesh_AABB.h
	mesh_AABB.cpp


	#Exception.h

	triangle_triangle_intersection.cpp

)

prepend_current_path(SOURCES)
fast_envelope_copy_headers(${SOURCES})
fast_envelope_set_source_group(${SOURCES})
target_sources(FastEnvelope PRIVATE ${SOURCES})


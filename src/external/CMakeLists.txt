set(SOURCES
	Predicates.cpp
	Predicates.hpp

	triangle_triangle_intersection.cpp
)

prepend_current_path(SOURCES)
fast_envelope_copy_headers(${SOURCES})
fast_envelope_set_source_group(${SOURCES})
target_sources(FastEnvelope PRIVATE ${SOURCES})


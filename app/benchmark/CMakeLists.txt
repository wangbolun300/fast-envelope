set(SOURCES
	MeshIO.hpp
	MeshIO.cpp

	#EnvelopeTest.cpp
	#EnvelopeTest.h

	mesh_AABB.h
	mesh_AABB.cpp

	AABBWrapper.cpp
	AABBWrapper.h
)

prepend_current_path(SOURCES)
fast_envelope_copy_headers(${SOURCES})
fast_envelope_set_source_group(${SOURCES})
target_sources(FastEnvelope_bin PRIVATE ${SOURCES})


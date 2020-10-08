set(APP_SOURCES
	MeshIO.hpp
	MeshIO.cpp

	#EnvelopeTest.cpp
	#EnvelopeTest.h

	mesh_AABB.h
	mesh_AABB.cpp

	AABBWrapper.cpp
	AABBWrapper.h
csv_reader.h
csv_reader.cpp

getRSS.hpp
sampling.h
sampling.cpp
)

prepend_current_path(APP_SOURCES)
fast_envelope_copy_headers(${APP_SOURCES})
fast_envelope_set_source_group(${APP_SOURCES})
target_sources(FastEnvelope_bin PRIVATE ${APP_SOURCES})


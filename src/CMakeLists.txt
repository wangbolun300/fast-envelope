set(SOURCES
        EnvelopeTest.cpp
        EnvelopeTest.h
        # FastEnvelope.cpp
        # FastEnvelope.h
        FastEnvelopeV2.cpp
        FastEnvelopeV2.h
        MeshIO.hpp
        MeshIO.cpp
        intersections.h
        intersections.cpp
        Logger.hpp
		Parameters.h
        Types.hpp
)

prepend_current_path(SOURCES)
fast_envelope_copy_headers(${SOURCES})
fast_envelope_set_source_group(${SOURCES})
target_sources(FastEnvelope PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(external)
add_subdirectory(benchmark)

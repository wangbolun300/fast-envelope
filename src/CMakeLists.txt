set(SOURCES
        FastEnvelope.cpp
        FastEnvelope.h
        MeshIO.hpp
        MeshIO.cpp
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
add_subdirectory(Precision)
add_subdirectory(benchmark)
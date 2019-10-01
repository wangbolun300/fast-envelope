set(SOURCES
        AABB.h
        AABB.cpp
        FastEnvelope.cpp
        FastEnvelope.h
        Logger.hpp

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
add_subdirectory(types)
add_subdirectory(indirectPredicates)
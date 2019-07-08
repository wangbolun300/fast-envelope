set(SOURCES
        FastEnvelope.cpp
        FastEnvelope.h
        Logger.hpp
        Parameters.h
        Rational.hpp
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
add_subdirectory(indirectPredicates)
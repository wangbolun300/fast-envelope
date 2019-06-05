set(SOURCES
	fprecision.h
	intervalprecision.h
	iprecision.h
	precisioncore.cpp
)

prepend_current_path(SOURCES)
fast_envelope_copy_headers(${SOURCES})
fast_envelope_set_source_group(${SOURCES})
add_library(Precision ${SOURCES})



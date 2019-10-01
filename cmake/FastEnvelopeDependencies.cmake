# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

### Configuration
set(FAST_ENVELOPE_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(FAST_ENVELOPE_EXTERNAL ${THIRD_PARTY_DIR})

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(FastEnvelopeDownloadExternal)

################################################################################
# Required libraries
################################################################################

# Sanitizers
if(FAST_ENVELOPE_WITH_SANITIZERS)
    fast_envelope_download_sanitizers()
    find_package(Sanitizers)
endif()

# CL11
if(NOT TARGET CLI11::CLI11)
    fast_envelope_download_cli11()
    add_subdirectory(${FAST_ENVELOPE_EXTERNAL}/cli11)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
	fast_envelope_download_spdlog()
	add_subdirectory(${FAST_ENVELOPE_EXTERNAL}/spdlog)
endif()



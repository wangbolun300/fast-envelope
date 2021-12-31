################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FAST_ENVELOPE_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FAST_ENVELOPE_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(fast_envelope_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${FAST_ENVELOPE_EXTERNAL}/${name}
        DOWNLOAD_DIR ${FAST_ENVELOPE_EXTERNAL}/.cache/${name}
        QUIET
        ${FAST_ENVELOPE_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## libigl
function(fast_envelope_download_libigl)
    fast_envelope_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        fc7f4b99415e0b0e3c4109d8e55b7ff6a5fc8b72
    )
endfunction()


## CLI11
function(fast_envelope_download_cli11)
    fast_envelope_download_project(cli11
        URL     https://github.com/CLIUtils/CLI11/archive/v1.6.1.tar.gz
        URL_MD5 48ef97262adb0b47a2f0a7edbda6e2aa
    )
endfunction()

## tbb
function(fast_envelope_download_tbb)
    fast_envelope_download_project(tbb
        GIT_REPOSITORY https://github.com/wjakob/tbb.git
        GIT_TAG        08b4341a1893a72656467e96137f1f99d0112547
    )
endfunction()


## Sanitizers
function(fast_envelope_download_sanitizers)
    fast_envelope_download_project(sanitizers-cmake
        GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
        GIT_TAG        6947cff3a9c9305eb9c16135dd81da3feb4bf87f
    )
endfunction()

## spdlog
function(fast_envelope_download_spdlog)
    fast_envelope_download_project(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        188cff7d6567b80c6b99bc15899fef9637a8fe52
    )
endfunction()

## Geogram LGPL
function(fast_envelope_download_geogram)
    fast_envelope_download_project(geogram
        GIT_REPOSITORY  https://github.com/polyfem/geogram.git
        GIT_TAG        516c151c244d9019a9076a1a468d52a0f6dd195d
    )
endfunction()


# ## Geogram predicates LGPL
# function(fast_envelope_download_geogram_predicates)
#     fast_envelope_download_project(geogram_predicates
#         URL     https://gforge.inria.fr/frs/download.php/file/38199/Predicates_psm_1.7.2.tar.gz
#         URL_MD5 4635a4e19538514ee2d6183af791e280
#     )
# endfunction()

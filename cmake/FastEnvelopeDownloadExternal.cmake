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
        GIT_TAG        7f7551af1613f60c7e89bef3934baf0da6298b06
    )
endfunction()


## Json
function(fast_envelope_download_json)
    fast_envelope_download_project(json
        GIT_REPOSITORY https://github.com/jdumas/json
        GIT_TAG        0901d33bf6e7dfe6f70fd9d142c8f5c6695c6c5b
    )
endfunction()

## Catch2
function(fast_envelope_download_catch2)
    fast_envelope_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.4.2
    )
endfunction()

## Multi precision
function(fast_envelope_download_arbitrary_precision)
    fast_envelope_download_project(arbitrary_precision
        GIT_REPOSITORY      https://github.com/teseoch/arbitrary_precision.git
        GIT_TAG             f7b744ad07e131ef65cfde10f086b6f6682a0205
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
        GIT_REPOSITORY https://github.com/alicevision/geogram.git
        GIT_TAG        v1.6.8
    )
endfunction()


## aabbcc
function(fast_envelope_download_aabbcc)
    fast_envelope_download_project(aabbcc
            GIT_REPOSITORY https://github.com/lohedges/aabbcc.git
            GIT_TAG        0c85e61362d384d70c71946826bfed0fb24a74ba
            )
endfunction()

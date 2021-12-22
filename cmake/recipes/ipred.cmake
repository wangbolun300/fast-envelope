if(TARGET ipred)
    return()
endif()

message(STATUS "Third-party: creating target 'ipred'")

include(FetchContent)
FetchContent_Declare(
    ipred
    GIT_REPOSITORY https://github.com/teseoch/Indirect_Predicates.git
    GIT_TAG f0528c03f392e25fc077ace893eb233e7ebc704d
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipred)


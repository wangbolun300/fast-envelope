if(TARGET ipred)
    return()
endif()

message(STATUS "Third-party: creating target 'ipred'")

include(FetchContent)
FetchContent_Declare(
    ipred
    GIT_REPOSITORY https://github.com/teseoch/Indirect_Predicates.git
    GIT_TAG a0bf2cd73383202b491a3d8d89f65a8efff232b1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipred)


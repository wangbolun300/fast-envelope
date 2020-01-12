################################################################################
cmake_minimum_required(VERSION 3.1)
project(FastEnvelope)
################################################################################

# cmake_policy(SET CMP0063 NEW)
# set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(FAST_ENVELOPE_TOPLEVEL_PROJECT ON)
else()
	set(FAST_ENVELOPE_TOPLEVEL_PROJECT OFF)
endif()

if(INPUT_THIRD_PARTY_DIR)
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_THIRD_PARTY_DIR}/)
else()
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# FastEnvelope utils
include(FastEnvelopeUtils)

# Extra warnings
include(Warnings)

# Use C++11/14
include(CXXFeatures)


# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################


################################################################################

# FastEnvelope options
option(FAST_ENVELOPE_ENABLE_TBB        			"Enable TBB"                    ON)
# option(FAST_ENVELOPE_USE_FLOAT          		"Use floats instead of double"  OFF)
option(FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES        "Use whole Geogram for predicates"        OFF)
option(FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES    "Use Geogram predicates only"             ON)
option(FAST_ENVELOPE_WITH_GMP						"Use gmp to have rational predicates"	  ON)

# Sanitizer options
option(SANITIZE_ADDRESS   "Sanitize Address"       OFF)
option(SANITIZE_MEMORY    "Sanitize Memory"        OFF)
option(SANITIZE_THREAD    "Sanitize Thread"        OFF)
option(SANITIZE_UNDEFINED "Sanitize Undefined"     OFF)

# Options for libigl modules
option(LIBIGL_USE_STATIC_LIBRARY          "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR            "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL                   "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO                 "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK                   "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE                 "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM                    "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB                 "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK                  "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL                 "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW            "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI      "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG                    "Use PNG"                      OFF)
option(LIBIGL_WITH_PYTHON                 "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN                 "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE               "Use Triangle"                 OFF)
option(LIBIGL_WITH_VIEWER                 "Use OpenGL viewer"            OFF)
option(LIBIGL_WITH_XML                    "Use XML"                      OFF)




# Sanitizers
if(FAST_ENVELOPE_WITH_SANITIZERS)
	list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/sanitizers-cmake/cmake)
endif()

# Setup dependencies
include(FastEnvelopeDependencies)

################################################################################
# FastEnvelope library
################################################################################

# add_library() can only be called without any source since CMake 3.11 ...
add_library(${PROJECT_NAME} src/common_algorithms.cpp)
target_compile_features(${PROJECT_NAME} PRIVATE ${CXX11_FEATURES})

# Public include directory for FastEnvelope
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})


# libigl
if(NOT TARGET igl::core)
	fast_envelope_download_libigl()
	find_package(LIBIGL REQUIRED)
endif()
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)

if(FAST_ENVELOPE_WITH_SANITIZERS)
	add_sanitizers(${PROJECT_NAME})
endif()


################################################################################
# Required libraries
################################################################################


# TBB
if(FAST_ENVELOPE_ENABLE_TBB AND NOT TARGET tbb::tbb)
	fast_envelope_download_tbb()
	set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
	set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
	set(TBB_NO_DATE ON CACHE BOOL " " FORCE)

	add_subdirectory(${THIRD_PARTY_DIR}/tbb tbb)
	set_target_properties(tbb_static PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${THIRD_PARTY_DIR}/tbb/include"
	)
	if(NOT MSVC)
		set_target_properties(tbb_static PROPERTIES
			COMPILE_FLAGS "-Wno-implicit-fallthrough -Wno-missing-field-initializers -Wno-unused-parameter -Wno-keyword-macro"
		)
		set_target_properties(tbb_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
	endif()
	target_compile_definitions(tbb_static PUBLIC -DUSE_TBB)
	add_library(tbb::tbb ALIAS tbb_static)

	set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")
endif()


# if(${FAST_ENVELOPE_USE_FLOAT})
# 	target_compile_definitions(${PROJECT_NAME} PUBLIC -DFAST_ENVELOPE_USE_FLOAT)
# endif()


#spdlog
target_link_libraries(${PROJECT_NAME} PUBLIC spdlog::spdlog)
target_link_libraries(${PROJECT_NAME} PUBLIC CLI11::CLI11)

if(FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES)
	# fast_envelope_download_geogram_predicates()

	add_library(GeoPredicates geogram_predicates/Predicates_psm.cpp)
	target_include_directories(GeoPredicates PUBLIC geogram_predicates/)
	target_compile_features(GeoPredicates PRIVATE ${CXX11_FEATURES})
	target_compile_definitions(GeoPredicates PRIVATE GEO_STATIC_LIBS)

	target_link_libraries(${PROJECT_NAME} PUBLIC GeoPredicates)
	target_compile_definitions(${PROJECT_NAME} PUBLIC ENVELOPE_WITH_GEO_PSM)

elseif(FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES)
	if(NOT TARGET geogram AND NOT TARGET geogram::geogram)
		fast_envelope_download_geogram()
		include(geogram)
	endif()
	if(TARGET geogram)
		target_link_libraries(${PROJECT_NAME} PUBLIC geogram)
	elseif(TARGET geogram::geogram)
		target_link_libraries(${PROJECT_NAME} PUBLIC geogram::geogram)
	else()
		MESSAGE(FATAL_ERROR "Unable to get geogram or geogram::geogram targets")
	endif()
	target_compile_definitions(${PROJECT_NAME} PUBLIC ENVELOPE_WITH_GEO)
endif()



# Predicates
set(INDIRECTPREDICATES_SOURCES
	indirectPredicates/ip_filtered.cpp
	indirectPredicates/ip_filtered.h

	indirectPredicates/ip_filtered_ex.h
)
add_library(IndirectPredicates ${INDIRECTPREDICATES_SOURCES})
indirect_predicates_copy_headers(${INDIRECTPREDICATES_SOURCES})
target_include_directories(IndirectPredicates PRIVATE indirectPredicates/)


if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	target_compile_options(IndirectPredicates PRIVATE "/fp:strict")
elseif(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	target_compile_options(IndirectPredicates PRIVATE "-frounding-math")
endif()
target_link_libraries(${PROJECT_NAME} PRIVATE IndirectPredicates)



################################################################################
# Subdirectories
################################################################################

# Sources
add_subdirectory(src)

if(FAST_ENVELOPE_WITH_GMP)
	find_package(GMP)

	IF(NOT ${GMP_FOUND})
			MESSAGE(FATAL_ERROR "Cannot find GMP")
	ENDIF()
	target_include_directories(${PROJECT_NAME} PUBLIC ${GMP_INCLUDE_DIRS})
	target_link_libraries(${PROJECT_NAME} PUBLIC ${GMP_LIBRARIES})
	target_compile_definitions(${PROJECT_NAME} PUBLIC ENVELOPE_WITH_GMP)

endif()


################################################################################
# FastEnvelope binary
################################################################################

# Main executable
if(FAST_ENVELOPE_TOPLEVEL_PROJECT)
	add_executable(${PROJECT_NAME}_bin app/main.cpp)
	add_subdirectory(app)
	target_compile_features(${PROJECT_NAME}_bin PUBLIC ${CXX14_FEATURES})
	

	if(NOT FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES)
		# Geogram
		fast_envelope_download_geogram()
		include(geogram)
		target_link_libraries(${PROJECT_NAME}_bin PUBLIC geogram)
	endif()

	target_link_libraries(${PROJECT_NAME}_bin
		PUBLIC
			${PROJECT_NAME}
			warnings::all
	)

	if(FAST_ENVELOPE_WITH_SANITIZERS)
		add_sanitizers(${PROJECT_NAME}_bin)
	endif()

	
	if(FAST_ENVELOPE_ENABLE_TBB)
		target_link_libraries(${PROJECT_NAME}_bin PUBLIC tbb::tbb)
	endif()

	if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)
	endif()
endif()


if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	foreach(config ${CMAKE_CONFIGURATION_TYPES})
		string(TOUPPER ${config} config)
		string(REPLACE /MT /MD CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MT /MD CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
	endforeach()
endif()

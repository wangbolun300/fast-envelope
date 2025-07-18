################################################################################
cmake_minimum_required(VERSION 3.15)
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
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes)
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
option(FAST_ENVELOPE_WITH_UNIT_TESTS   			"unit test project" 			ON)
option(FAST_ENVELOPE_ENABLE_TBB        			"Enable TBB"                    ON)
# option(FAST_ENVELOPE_USE_FLOAT          		"Use floats instead of double"  OFF)

option(FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES        "Use whole Geogram for predicates"        ON)
option(FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES    "Use Geogram predicates only"             OFF)
option(FAST_ENVELOPE_WITH_GMP						"Use gmp to have rational predicates"	  OFF)
option(FAST_ENVELOPE_WITH_TIMER   					"use timer in the library for debug" 	  OFF)

# Sanitizer options
option(SANITIZE_ADDRESS   "Sanitize Address"       OFF)
option(SANITIZE_MEMORY    "Sanitize Memory"        OFF)
option(SANITIZE_THREAD    "Sanitize Thread"        OFF)
option(SANITIZE_UNDEFINED "Sanitize Undefined"     OFF)


if(NOT ${FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES} AND NOT ${FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES})
	set(LIBIGL_WITH_PREDICATES ON)
else()
	set(LIBIGL_WITH_PREDICATES OFF)
endif()



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
if(FAST_ENVELOPE_TOPLEVEL_PROJECT)
	if(NOT TARGET igl::core)
		fast_envelope_download_libigl()
	endif()
	igl_include(core)
	target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
	if(LIBIGL_WITH_PREDICATES)
		igl_include(predicates)
		target_link_libraries(${PROJECT_NAME} PUBLIC igl::predicates)
	endif()
else()
	if(LIBIGL_WITH_PREDICATES)
		if(NOT TARGET igl::core)
			fast_envelope_download_libigl()
		endif()
		igl_include(core)
		igl_include(predicates)
		target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
		target_link_libraries(${PROJECT_NAME} PUBLIC igl::predicates)
		# target_include_directories(${PROJECT_NAME} PUBLIC ${THIRD_PARTY_DIR}/libigl/include/igl/predicates)
	else()
		# Assume eigen exists in the other lib
		target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
	endif()
endif()

if(FAST_ENVELOPE_WITH_SANITIZERS)
	add_sanitizers(${PROJECT_NAME})
endif()


################################################################################
# Required libraries
################################################################################


# TBB
if(FAST_ENVELOPE_ENABLE_TBB AND NOT TARGET TBB::tbb)
include(tbb)
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
	target_compile_definitions(${PROJECT_NAME} PUBLIC ENVELOPE_WITH_GEO_PSM ENVELOPE_WITH_GEO)

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
include("ipred")
target_link_libraries(${PROJECT_NAME} PUBLIC indirectPredicates)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message(WARNING "On CLang, there is no compiler flag which is required for ensuring the correctness of the algorithm.")
endif()
set(INDIRECTPREDICATES_HEADERS
	indirectPredicates/ip_filtered.h
	indirectPredicates/ip_filtered_ex.h
)
indirect_predicates_copy_headers(${INDIRECTPREDICATES_HEADERS})

if(FAST_ENVELOPE_WITH_TIMER)
	target_compile_definitions(${PROJECT_NAME} PRIVATE FAST_ENVELOPE_USE_TIMER)
endif()





################################################################################
# Subdirectories
################################################################################

# Sources
add_subdirectory(src)
target_sources(${PROJECT_NAME} PRIVATE indirectPredicates/ip_filtered.cpp)

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
	# add an application library
	#add_subdirectory(app/benchmark)
	set(APP_SOURCES
		app/benchmark/MeshIO.hpp
		app/benchmark/MeshIO.cpp
		app/benchmark/mesh_AABB.h
		app/benchmark/mesh_AABB.cpp
		app/benchmark/AABBWrapper.cpp
		app/benchmark/AABBWrapper.h
		app/benchmark/csv_reader.h
		app/benchmark/csv_reader.cpp
		app/benchmark/getRSS.hpp
		app/benchmark/sampling.h
		app/benchmark/sampling.cpp
	)
	fast_envelope_copy_headers(${APP_SOURCES})
	add_library(${PROJECT_NAME}_app ${APP_SOURCES})
	target_include_directories(${PROJECT_NAME}_app PUBLIC ${PROJECT_BINARY_DIR}/include)
	target_include_directories(${PROJECT_NAME}_app PUBLIC app)

	# link libraries to the app library
	if(NOT FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES)
		# Geogram
		fast_envelope_download_geogram()
		include(geogram)
		if(TARGET geogram)
			target_link_libraries(${PROJECT_NAME}_app PUBLIC geogram)
		elseif(TARGET geogram::geogram)
			target_link_libraries(${PROJECT_NAME}_app PUBLIC geogram::geogram)
		else()
			MESSAGE(FATAL_ERROR "Unable to get geogram or geogram::geogram targets")
		endif()
	endif()
	target_link_libraries(${PROJECT_NAME}_app
		PUBLIC
			${PROJECT_NAME}

	)

	# executable
	add_executable(${PROJECT_NAME}_bin app/main.cpp)
	target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_app)

	target_compile_features(${PROJECT_NAME}_bin PUBLIC ${CXX14_FEATURES})
	if(FAST_ENVELOPE_WITH_SANITIZERS)
		add_sanitizers(${PROJECT_NAME}_bin)
	endif()


	if(FAST_ENVELOPE_ENABLE_TBB)
		target_link_libraries(${PROJECT_NAME}_bin PUBLIC TBB::tbb)
	endif()
	if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)
	endif()

	# unit tests
	if(FAST_ENVELOPE_WITH_UNIT_TESTS)

		add_executable(${PROJECT_NAME}_test app/test.cpp)
		target_compile_definitions(${PROJECT_NAME}_test PUBLIC
    ENVELOPE_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/app/test_data/")
		target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME}_app)

		target_compile_features(${PROJECT_NAME}_test PUBLIC ${CXX14_FEATURES})
		if(FAST_ENVELOPE_WITH_SANITIZERS)
			add_sanitizers(${PROJECT_NAME}_test)
		endif()


		if(FAST_ENVELOPE_ENABLE_TBB)
			target_link_libraries(${PROJECT_NAME}_test PUBLIC TBB::tbb)
		endif()
		if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
			set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_test)
		endif()
	endif()
	###################################################################



endif()

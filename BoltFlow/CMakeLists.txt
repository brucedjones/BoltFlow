cmake_minimum_required (VERSION 2.6)
project (BoltFlow)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

add_executable(runbolt2d index.cpp)
add_executable(runbolt3d index.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package (CGNS REQUIRED)
include_directories(${CGNS_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CGNS_LIBRARIES})

find_package (HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_LIBRARIES})

find_package (ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ZLIB_LIBRARIES})

find_package (SZIP REQUIRED)
include_directories(${SZIP_INCLUDE_DIRS})
set(LIBS ${LIBS} ${SZIP_LIBRARIES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET_TARGET_PROPERTIES(runbolt2d PROPERTIES COMPILE_DEFINITIONS "D2Q9")
SET_TARGET_PROPERTIES(runbolt3d PROPERTIES COMPILE_DEFINITIONS "D3Q15")

target_link_libraries (runbolt2d ${LIBS})
target_link_libraries (runbolt3d ${LIBS})

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	SET(CMAKE_INSTALL_PREFIX "/usr/local/bin" CACHE PATH "BoltFlow install prefix" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

install(TARGETS runbolt2d runbolt3d RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
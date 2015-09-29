cmake_minimum_required (VERSION 2.6)
project (BoltFlow)

add_executable(runbolt2d index.cpp)
add_executable(runbolt3d index.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package (CGNS REQUIRED)
include_directories(${CGNS_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CGNS_LIBRARIES})

find_package (HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_LIBRARIES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET_TARGET_PROPERTIES(runbolt2d PROPERTIES COMPILE_DEFINITIONS "D2Q9")
SET_TARGET_PROPERTIES(runbolt3d PROPERTIES COMPILE_DEFINITIONS "D3Q15")

target_link_libraries (runbolt2d ${LIBS})
target_link_libraries (runbolt3d ${LIBS})

install(TARGETS runbolt2d runbolt3d RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
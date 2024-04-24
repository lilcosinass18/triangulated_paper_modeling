cmake_minimum_required(VERSION 3.25)
project(paper_modeling)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_TOOLCHAIN_FILE "/home/big/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

add_subdirectory(src)

include_directories("/usr/include/nlohmann")


find_package(CGAL REQUIRED)
find_package(OpenMesh REQUIRED)

set(Eigen3_DIR /usr/local/share/eigen3/cmake)
#set(libigl_DIR /usr/local/lib/cmake/igl)

find_package(Eigen3 REQUIRED)
#find_package(libigl CONFIG REQUIRED)

#include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${IGL_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})

link_directories("/usr/local/include/OpenMesh/Core")
link_directories("/usr/local/include/OpenMesh/Tools")

find_package(nlohmann_json REQUIRED)

#file(GLOB_RECURSE HEADER_FILES "${CMAKE_SOURCE_DIR}/home/big/libs/libigl/include/igl/*.h")
#file(GLOB_RECURSE HEADER_FILES "${CMAKE_SOURCE_DIR}/home/big/libs/libigl/build/_deps/cgal-src/include/CGAL/*.h")

add_executable(paper_modeling main.cpp)

target_link_libraries(paper_modeling paper nlohmann_json::nlohmann_json CGAL::CGAL
        OpenMeshCore OpenMeshTools)

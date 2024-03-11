cmake_minimum_required(VERSION 3.25)
project(paper_modeling)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(src)

include_directories("lib/eigen/")
include_directories("/usr/include/nlohmann")

find_package(nlohmann_json REQUIRED)

add_executable(paper_modeling main.cpp)

target_link_libraries(paper_modeling nlohmann_json::nlohmann_json)
target_link_libraries(paper_modeling paper)
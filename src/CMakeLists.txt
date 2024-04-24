cmake_minimum_required(VERSION 3.25)
project(paper_modeling)

set(CMAKE_CXX_STANDARD 20)

include_directories("/usr/include/nlohmann")

find_package(nlohmann_json REQUIRED)

add_library(paper geometry.hpp geometry.cpp paper.hpp paper.cpp utils.hpp utils.cpp
            permution.hpp permution.cpp)
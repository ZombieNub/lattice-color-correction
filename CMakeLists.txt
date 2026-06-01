cmake_minimum_required(VERSION 4.2)
project(ilpfun)

set(CMAKE_CXX_STANDARD 20)

add_executable(ilpfun src/main.cpp)

target_include_directories(ilpfun PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_directories(ilpfun PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib)

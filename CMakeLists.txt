cmake_minimum_required(VERSION 3.30)
project(ilpfun)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configure OR-Tools features before pulling the source
set(BUILD_DEPS ON CACHE INTERNAL "")
set(USE_HIGHS ON CACHE INTERNAL "")

# Pull OR-Tools from GitHub
include(FetchContent)
FetchContent_Declare(
        or-tools
        GIT_REPOSITORY https://github.com/google/or-tools.git
        GIT_TAG        main
)
FetchContent_MakeAvailable(or-tools)

add_executable(ilpfun src/main.cpp)

target_include_directories(ilpfun PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ilpfun PRIVATE ortools::ortools)
target_compile_options(ilpfun PRIVATE -O3)
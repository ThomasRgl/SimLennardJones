cmake_minimum_required(VERSION 3.14)

project(LJS LANGUAGES C)

# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)
# find_package(OpenMP REQUIRED)

add_subdirectory(deps/argtable)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

set(SOURCES src/main.c src/eig.c)


add_subdirectory(src)
add_subdirectory(apps)

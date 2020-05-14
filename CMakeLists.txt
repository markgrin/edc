cmake_minimum_required(VERSION 3.9) # Проверка версии CMake.

project(quad_sieve)

add_executable(quad_sieve src/main.cpp src/common.cpp src/ressol.cpp)
target_link_libraries(quad_sieve PRIVATE gmp gmpxx)
target_compile_features(quad_sieve PRIVATE cxx_std_17)


cmake_minimum_required(VERSION 3.27)
project(Numerical_Optimization)

set(CMAKE_CXX_STANDARD 17)

add_executable(Numerical_Optimization main.cpp)

add_executable(matrix_test test/matrix_test.cpp
        Matrix/Matrix.hpp)

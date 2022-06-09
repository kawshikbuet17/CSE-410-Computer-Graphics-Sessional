cmake_minimum_required(VERSION 3.22)
project(Previous-Online-Codes)

set(CMAKE_CXX_STANDARD 14)

add_executable(A1-Online A1-Online.cpp)
target_link_libraries(A1-Online -lOpenGL32 -lfreeGLUT -lglu32)

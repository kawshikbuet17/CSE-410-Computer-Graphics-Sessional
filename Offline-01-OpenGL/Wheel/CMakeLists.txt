cmake_minimum_required(VERSION 3.22)
project(Wheel)

set(CMAKE_CXX_STANDARD 14)

add_executable(Wheel main.cpp)
target_link_libraries(Wheel -lOpenGL32 -lfreeGLUT -lglu32)
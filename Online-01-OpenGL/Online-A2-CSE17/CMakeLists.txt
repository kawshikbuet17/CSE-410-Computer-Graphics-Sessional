cmake_minimum_required(VERSION 3.22)
project(Online_A2_CSE17)

set(CMAKE_CXX_STANDARD 14)

add_executable(Online_A2_CSE17 main.cpp)
target_link_libraries(Online_A2_CSE17 -lOpenGL32 -lfreeGLUT -lglu32)
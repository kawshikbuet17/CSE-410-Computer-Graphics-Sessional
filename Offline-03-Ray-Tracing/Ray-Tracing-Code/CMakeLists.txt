cmake_minimum_required(VERSION 3.22)
project(Ray_Tracing_Code)

set(CMAKE_CXX_STANDARD 14)

add_executable(Ray_Tracing_Code 1705043_Main.cpp 1705043_Header.h)
target_link_libraries(Ray_Tracing_Code -lOpenGL32 -lfreeGLUT -lglu32)

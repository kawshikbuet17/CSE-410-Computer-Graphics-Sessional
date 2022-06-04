cmake_minimum_required(VERSION 3.22)
project(FullyControllableCamera)

set(CMAKE_CXX_STANDARD 14)

add_executable(CameraSphereCube main.cpp)
target_link_libraries(CameraSphereCube -lOpenGL32 -lfreeGLUT -lglu32)
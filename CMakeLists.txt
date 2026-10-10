cmake_minimum_required(VERSION 3.16)
PROJECT(jpeg_encoder)
add_executable(svjpeg svjpeg.cpp)
install(TARGETS svjpeg DESTINATION ${PROJECT_SOURCE_DIR})

cmake_minimum_required(VERSION 3.5)
project(IM)

set(CMAKE_CXX_STANDARD 14)
include_directories ("${PROJECT_SOURCE_DIR}/sfmt")


add_executable(IM main.cc allocator.h utils.cc TimGraph.cc anyoption.cc sfmt/SFMT.c allocator.cc)

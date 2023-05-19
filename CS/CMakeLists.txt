cmake_minimum_required(VERSION 3.5)
project(CS)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(CS main.cpp generate_data.h read_data.h GENSM_slow.h RandomClock.h Simultaneous.h RandomClock_arbitrary.h RandomClockMatroid.h GENSM_slow_matroid.h BFM_NM.h)
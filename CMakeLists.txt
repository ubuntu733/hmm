cmake_minimum_required(VERSION 3.8)
project(hmm)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp model.cc model.h real.h vector.cc vector.h matrix.cc matrix.h)
add_executable(hmm ${SOURCE_FILES})
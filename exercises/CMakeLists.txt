project("workshop")
cmake_minimum_required(VERSION 2.8)

find_package(AMDIS REQUIRED)

add_executable("exercise1" src/exercise1.cc)
target_link_libraries("exercise1" AMDiS)

add_executable("exercise2" src/exercise2.cc)
target_link_libraries("exercise2" AMDiS)

add_executable("exercise3" src/exercise3.cc)
target_link_libraries("exercise3" AMDiS)

add_executable("exercise4" src/exercise4.cc)
target_link_libraries("exercise4" AMDiS)

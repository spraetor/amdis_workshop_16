project("workshop")
cmake_minimum_required(VERSION 3.1)

find_package(AMDIS REQUIRED)

add_executable("poisson1" src/poisson1.cc)
target_link_libraries("poisson1" AMDiS)

add_executable("poisson2" src/poisson2.cc)
target_link_libraries("poisson2" AMDiS)

add_executable("poisson3" src/poisson3.cc)
target_link_libraries("poisson3" AMDiS)

add_executable("poisson4" src/poisson4.cc)
target_link_libraries("poisson4" AMDiS)

add_executable("heat1" src/heat1.cc)
target_link_libraries("heat1" AMDiS)

add_executable("heat2" src/heat2.cc)
target_link_libraries("heat2" AMDiS)

add_executable("heat3" src/heat3.cc)
target_link_libraries("heat3" AMDiS)

add_executable("heat4" src/heat4.cc)
target_link_libraries("heat4" AMDiS)

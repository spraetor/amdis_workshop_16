project("workshop")
cmake_minimum_required(VERSION 2.8)

find_package(AMDIS REQUIRED)
if(AMDIS_FOUND)
	message(STATUS "AMDiS was found.")
	include(${AMDIS_USE_FILE})
	SET(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

add_executable("poisson1" src/poisson1.cc)
target_link_libraries("poisson1" ${BASIS_LIBS})

add_executable("poisson2" src/poisson2.cc)
target_link_libraries("poisson2" ${BASIS_LIBS})

add_executable("poisson3" src/poisson3.cc)
target_link_libraries("poisson3" ${BASIS_LIBS})

add_executable("poisson4" src/poisson4.cc)
target_link_libraries("poisson4" ${BASIS_LIBS})

add_executable("heat1" src/heat1.cc)
target_link_libraries("heat1" ${BASIS_LIBS})

add_executable("heat2" src/heat2.cc)
target_link_libraries("heat2" ${BASIS_LIBS})

add_executable("heat3" src/heat3.cc)
target_link_libraries("heat3" ${BASIS_LIBS})

add_executable("heat4" src/heat4.cc)
target_link_libraries("heat4" ${BASIS_LIBS})
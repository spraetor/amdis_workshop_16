project("workshop")
cmake_minimum_required(VERSION 2.8)

find_package(AMDIS REQUIRED)
if(AMDIS_FOUND)
	message(STATUS "AMDiS was found.")
	include(${AMDIS_USE_FILE})
	SET(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

add_executable("exercise1" src/exercise1.cc)
target_link_libraries("exercise1" ${BASIS_LIBS})

add_executable("exercise2" src/exercise2.cc)
target_link_libraries("exercise2" ${BASIS_LIBS})

add_executable("exercise2b" src/exercise2b.cc)
target_link_libraries("exercise2b" ${BASIS_LIBS})

add_executable("exercise3" src/exercise3.cc)
target_link_libraries("exercise3" ${BASIS_LIBS})

add_executable("exercise3b" src/exercise3b.cc)
target_link_libraries("exercise3b" ${BASIS_LIBS})

add_executable("exercise4" src/exercise4.cc)
target_link_libraries("exercise4" ${BASIS_LIBS})

add_executable("exercise4b" src/exercise4b.cc)
target_link_libraries("exercise4b" ${BASIS_LIBS})

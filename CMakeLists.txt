cmake_minimum_required(VERSION 3.0)

project(bandmc)

set(CMAKE_CXX_FLAGS "-std=c++11")

enable_language(Fortran)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

FILE(GLOB headers ${CMAKE_SOURCE_DIR}/include/*.hh)
FILE(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/fortran/*.f)
#FILE(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cc)

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})


ROOT_GENERATE_DICTIONARY(
    GenTreeDict                         # path to dictionary to generate
    ${CMAKE_SOURCE_DIR}/include/GenTree.hh         # list of classes to process
    LINKDEF ${CMAKE_SOURCE_DIR}/include/LinkDef.hh  # ROOT linkDef file
)

add_executable(bandmc bandmc.cc ${headers} ${sources} GenTreeDict.cxx)
target_link_libraries(bandmc ${ROOT_LIBRARIES})

set(DAT_FILES upper_momentum_fit.dat lower_momentum_fit.dat)

foreach(_datfile ${DAT_FILES})
	configure_file(
    	${PROJECT_SOURCE_DIR}/${_datfile}
    	${PROJECT_BINARY_DIR}/${_datfile}
    	COPYONLY
   	)	
endforeach()



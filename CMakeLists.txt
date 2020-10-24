cmake_minimum_required(VERSION 3.0)

project(bandmc)

set(CMAKE_CXX_FLAGS "-std=c++11 -ggdb -O0")

enable_language(Fortran)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

FILE(GLOB headers ${CMAKE_SOURCE_DIR}/include/*.hh)
FILE(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/fortran/*.f)

#SET(BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libbandhit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libclashit_tree.dylib)
FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/*.dylib)

message("BAND LIBS = ${BAND_LIBS}")

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} $ENV{BANDSOFT_TOOLS_DIR}/include)

ROOT_GENERATE_DICTIONARY(
    bandmcDict                         # path to dictionary to generate
    ${CMAKE_SOURCE_DIR}/include/GenTree.hh         # list of classes to process
    LINKDEF ${CMAKE_SOURCE_DIR}/include/LinkDef.hh  # ROOT linkDef file
)

add_executable(bandmc bandmc.cc ${headers} ${sources} bandmcDict.cxx)
target_link_libraries(bandmc ${ROOT_LIBRARIES} Foam ${BAND_LIBS})

add_executable(wimgen wimgen.cc ${headers} ${sources} bandmcDict.cxx)
target_link_libraries(wimgen ${ROOT_LIBRARIES} Foam ${BAND_LIBS})

set(DAT_FILES dat/upper_momentum_fit.dat dat/lower_momentum_fit.dat dat/a09.dsfs_lNNC dat/a09.sfs_lNNC dat/edep_dist.root dat/RC_graph.root)

foreach(_datfile ${DAT_FILES})
	configure_file(
    	${PROJECT_SOURCE_DIR}/${_datfile}
    	${PROJECT_BINARY_DIR}/${_datfile}
    	COPYONLY
   	)	
endforeach()



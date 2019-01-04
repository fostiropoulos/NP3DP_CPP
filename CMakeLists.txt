cmake_minimum_required (VERSION 3.5)

# Executable name: 
##############################################################
SET(EXEC OUTFile)
##############################################################

project(${EXEC})

#include
include_directories(/usr/local/include/eigen3)
include_directories(/usr/local/include)
include_directories(/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/include)
include_directories(/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data)

#link
link_directories(/usrl/ocal/build/libnabo.a gomp)
link_directories(/usr/local/lib)

# support for Libnabo  
SET(CMAKE_CXX_FLAGS "-std=c++11 -O3")

set(UTLIB src/utilities.cpp
	src/NPAM_utilities.cpp 
	src/file_rw.cpp 
	src/transformation_utilities.cpp)

##############################################################
set(MAIN_FILE src/main_file.cpp)
##############################################################

add_executable(${EXEC} ${MAIN_FILE} ${UTLIB})

##LINKED LIBRARIES
# link library for libnabo
target_link_libraries(${EXEC} /usr/local/build/libnabo.a gomp)
# link library for -lpthread (used for compiling the libigl)
target_link_libraries(${EXEC} pthread)


set_property(TARGET ${EXEC} PROPERTY CXX_STANDARD 14)


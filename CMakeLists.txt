cmake_minimum_required (VERSION 3.2)

# Executable name:
##############################################################
SET(EXEC tool_path_library_test)
SET(PROJ_NAME tool_path_library)

##############################################################

project(${PROJ_NAME})

#include
include_directories(/usr/local/include/eigen3)
include_directories(/usr/local/include)
include_directories(include)
include_directories(data)
include_directories(../gnuplot-iostream)


#link
link_directories(/usr/local/lib)

# support for Libnabo
SET(CMAKE_CXX_FLAGS "-std=c++11 -O3")

set(UTLIB src/utilities.cpp
	src/NPAM_utilities.cpp
	src/file_rw.cpp
	src/transformation_utilities.cpp
	src/plot_utilities.cpp)

##############################################################
set(MAIN_FILE src/main_file.cpp)
add_library(${PROJECT_NAME} STATIC ${UTLIB})
##############################################################
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})
add_executable(${EXEC} ${MAIN_FILE} ${UTLIB})

##LINKED LIBRARIES
# link library for -lpthread (used for compiling the libigl)
target_link_libraries(${EXEC} pthread)

# target libraries for using the GNUPLOT
target_link_libraries(${EXEC} /usr/lib/x86_64-linux-gnu/libboost_iostreams.a)
target_link_libraries(${EXEC} /usr/lib/x86_64-linux-gnu/libboost_system.a)
target_link_libraries(${EXEC} /usr/lib/x86_64-linux-gnu/libboost_filesystem.a)
target_link_libraries(${EXEC} /usr/lib/x86_64-linux-gnu/libutil.a)
target_link_libraries(${EXEC} /usr/lib/x86_64-linux-gnu/libboost_thread.a)



set_property(TARGET ${EXEC} PROPERTY CXX_STANDARD 14)

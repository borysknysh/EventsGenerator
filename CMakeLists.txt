cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")

project(Gamapola)
set(CMAKE_BUILD_TYPE "Debug")

ADD_DEFINITIONS( "-Wl,--no-as-needed" )
ADD_DEFINITIONS( "`root-config --cflags`" )
ADD_DEFINITIONS( "-O2 -W -Wall -Wextra -Wno-long-long" )
ADD_DEFINITIONS( "-fno-common" )
ADD_DEFINITIONS( "-fopenmp" )
ADD_DEFINITIONS( "-Wno-deprecated" )
ADD_DEFINITIONS( "-std=c++11" )

include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})

link_directories($ENV{ROOTSYS}/lib)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_executable(main main.C ${sources} ${headers} )
target_link_libraries(main  -lm -lstdc++ -pthread -I/usr/local/include/root -L/usr/lib/x86_64-linux-gnu/ -L/usr/aarch64-linux-gnu/lib -L/usr/local/lib/root -L/usr/local/lib/ -lTree -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMinuit2 -lMathMore -lThread -lm -rdynamic -L/usr/local/lib -ldl)

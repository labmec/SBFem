# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 2.8.0)

project (SBFem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#builds and installs PZ before any project
ExternalProject_Add(PZ_proj
  URL http://github.com/labmec/neopz
  PREFIX /home/labmec/neopz/externos/
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make -j24
  INSTALL_COMMAND make install
  PREFIX=/home/labmec/neopz/neopz-install
)
#ExternalProject_Get_Property(PZ_proj install_dir)

#add_library(luajit STATIC IMPORTED)
#set_property(TARGET luajit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libluajit-5.1.a)
add_dependencies(PZ PZ_proj)
#add_executable(myapp main.c)
#include_directories(${install_dir}/include/luajit-2.0)
#target_link_libraries(myapp luajit)


IF (WIN32)
   find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
   find_package( PZ REQUIRED )
ENDIF (WIN32)

include_directories(${PZ_INCLUDE_DIRS})

if(NOT ${PZ_USING_MKL})
   MESSAGE(FATAL_ERROR "Found a PZ library not configured with Intel MKL library. SBFem needs PZ configured with MKL library to work!")
endif()

IF(APPLE)
   # Mac OS X specific code
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
ENDIF(APPLE)


# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include)

if(WIN32)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
   #add a define saying it is a VC compiler
   set(VC "#Will add a define VC on pz_config.h")
   add_definitions(-DVC)
   #define _USE_MATH_DEFINES for example: M_PI constant definitions.
   add_definitions(-D_USE_MATH_DEFINES)
   #disabling VC warnings
   add_definitions(-D_SCL_SECURE_NO_WARNINGS)
   add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	
   #define use of pthread static lib.
   add_definitions(-DPTW32_STATIC_LIB)
   #enabling /bigobj
   add_definitions("/bigobj")
 
   FOREACH(FLAG_TYPE EXE MODULE SHARED)
       STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
       STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
       STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
       STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
       SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
       MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
   ENDFOREACH ()
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
endif()

add_subdirectory(SBFem2D)
add_subdirectory(SBFem3D)

#Enabling nightly tests
ENABLE_TESTING()
include(CTest)

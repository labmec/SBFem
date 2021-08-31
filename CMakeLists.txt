# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.14)

project (SBFem LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF(APPLE)
   # Mac OS X specific code
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
ENDIF(APPLE)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/)

check_pz_opt(PZ_USING_MKL _PZ_HAS_MKL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing") #this flag is for preventing errors when casting from double to float

option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)
if(USING_BOOST)
	IF (WIN32)
		set(Boost_USE_STATIC_LIBS   ON)
		set(Boost_USE_MULTITHREADED ON)
        set(CMAKE_PREFIX_PATH "C:/local/boost_1_69_0" ${CMAKE_PREFIX_PATH})
	ENDIF()
	find_package( Boost COMPONENTS unit_test_framework date_time atomic REQUIRED )
endif(USING_BOOST)


#enabling Boost library
option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)

if(USING_BOOST)
	IF (WIN32)
		set(Boost_USE_STATIC_LIBS   ON)
		set(Boost_USE_MULTITHREADED ON)
        set(CMAKE_PREFIX_PATH "C:/local/boost_1_72_0" ${CMAKE_PREFIX_PATH})
	ENDIF()
	find_package( Boost COMPONENTS unit_test_framework date_time atomic REQUIRED )
endif(USING_BOOST)

option(USING_BLAZE "Whether the blaze package will be linked with the library" OFF)
if(USING_BLAZE)
	find_package( blaze REQUIRED)
	if( blaze_FOUND )
	    add_library( blaze_target INTERFACE )
	    target_link_libraries( blaze_target INTERFACE blaze::blaze )
	endif()
endif()

add_subdirectory(SBFem2D)
add_subdirectory(SBFem3D)
add_subdirectory(BubbleFunctions)

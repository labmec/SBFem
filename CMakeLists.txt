# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 2.8.0)

project (SBFem)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

# Boost library dependencies
IF (WIN32)
        set(Boost_USE_STATIC_LIBS   ON)
        set(Boost_USE_MULTITHREADED ON)
        set(CMAKE_PREFIX_PATH "C:/local/boost_1_72_0" ${CMAKE_PREFIX_PATH})
ENDIF()

#find_package( Boost 1.71 COMPONENTS unit_test_framework date_time atomic REQUIRED )
#enabling Boost library
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

option(USING_PETSC "Whether the blaze package will be linked with the library" OFF)
if(USING_PETSC)
  set(PETSC_DIR "/home/karol/Documents/GitHub/petsc")
  FIND_PACKAGE( PkgConfig REQUIRED )
  set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/arch-linux-c-debug/lib/pkgconfig")
  message("${PKG_CONFIG_PATH}")

  pkg_search_module (PETSc REQUIRED PETSc)
  set(_PETSC_LIBRARIES)
  foreach (lib ${PETSC_LIBRARIES})
    find_library(_LIB ${lib} PATHS ${PETSC_LIBRARY_DIRS} NO_DEFAULT_PATH)
    list(APPEND _PETSC_LIBRARIES ${_LIB})
  endforeach()

  set(PETSC_LIBRARIES ${_PETSC_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${PETSC_LIBRARIES})
  
  find_package_handle_standard_args(PETSc PETSC_LIBRARY_DIRS)
    add_definitions(-DUSING_PETSC)
endif()

add_subdirectory(SBFem2D)
add_subdirectory(SBFem3D)
add_subdirectory(BubbleFunctions)

#Enabling nightly tests
ENABLE_TESTING()
include(CTest)

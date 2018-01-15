# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 2.8.0)

project (SBFem)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")

IF (NOT WIN32)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
ENDIF()

IF (WIN32)
	find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
	find_package( PZ REQUIRED )
ENDIF (WIN32)

include_directories(${PZ_DIR}/../pzlib/include)

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include /usr/lib/atlas-base/ /usr/include/atlas /usr/local/atlas/include /usr/local/atlas/lib )

#enabling pthread library
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/include )
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/lib  ${PROJECT_SOURCE_DIR}/../pthread/lib)
if(PTHREAD_INCLUDE-NOTFOUND)
	set (PTHREAD_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
endif()
if(PTHREAD_LIB-NOTFOUND)
	set (PTHREAD_LIB "" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${PTHREAD_INCLUDE})


        add_definitions(-DUSING_MKL)

		if (WIN32)
			set(MKL_ROOT "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows")
			find_path(MKL_INCLUDE NAMES mkl.h PATHS ${MKL_ROOT}/mkl/include)
			find_library(MKL_LIB_INTEL NAMES mkl_rt.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
			find_library(MKL_LIB_CORE NAMES mkl_core.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
			find_library(MKL_LIB_THREAD NAMES mkl_intel_thread.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
            find_library(COMPOSER_OMP NAMES libiomp5md.lib PATHS ${MKL_ROOT}/compiler/lib/intel64_win)
		else()
			find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include /softwares/intel/mkl/include)
			find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
					${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64)
			find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
					/opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
			find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
					${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
			find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
					/opt/intel/composer_xe/compiler/lib
					/opt/intel/composerxe/lib/intel64
					/opt/intel/lib /softwares/intel/lib/intel64)
			#find_path(OMP_INCLUDE NAMES omp.h PATHS ${SEARCH_DIRS} /opt/intel/composer_xe_2013_sp1.1.103/compiler/include)
			#include_directories(${OMP_INCLUDE})
		endif()

		if(MKL_INCLUDE-NOTFOUND)
                set (MKL_INCLUDE "" CACHE PATH "Directory where mkl.h can be found")
        else()
                include_directories(${MKL_INCLUDE})
        endif()
        if(MKL_LIB-NOTFOUND)
                set (MKL_LIB_INTEL "" CACHE PATH "Directory where the mkl library can be found")
        else()
                link_libraries(${MKL_LIB_INTEL})
                link_libraries(${MKL_LIB_CORE})
                link_libraries(${MKL_LIB_THREAD})
                link_libraries(${COMPOSER_OMP})
        endif()





if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
	#disabling STL warning for std::_Vector_alloc when LOG4CXX is enabled
	if (USING_LOG4CXX)
	    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /wd4251")
	    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /wd4251")
	endif()
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

endif()

add_subdirectory(SBFem2D)
add_subdirectory(SBFem3D)

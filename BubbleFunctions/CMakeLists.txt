file(GLOB headers *.h)
file(GLOB sources *.cpp)

set (MainParabolic main2DParabolicEq.cpp Common.cpp ${headers})
add_executable(SB2DBubblesParabolic ${MainParabolic})
target_link_libraries(SB2DBubblesParabolic pz)

set (MainRegular mainregularsolution.cpp Common.cpp ${headers})
add_executable(SB2DBubblesRegular ${MainRegular})
target_link_libraries(SB2DBubblesRegular pz)


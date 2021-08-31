file(GLOB headers *.h)
file(GLOB sources *.cpp)

set (MainParabolic main2DParabolicEq.cpp Common.cpp ${headers})
add_executable(SB2DBubblesParabolic ${MainParabolic})
target_link_libraries(SB2DBubblesParabolic PUBLIC NeoPZ::pz)

set (MainRegular mainregularsolution.cpp Common.cpp ${headers})
add_executable(SB2DBubblesRegular ${MainRegular})
target_link_libraries(SB2DBubblesRegular PUBLIC NeoPZ::pz)

set (MainRegular3d mainregular3dsolution.cpp Common3D.cpp ${headers})
add_executable(SB3DBubblesRegular ${MainRegular3d})
target_link_libraries(SB3DBubblesRegular PUBLIC NeoPZ::pz)

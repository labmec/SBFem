file(GLOB headers *.h)
file(GLOB sources *.cpp)

# Regular solutions for Laplace and Elasticity with SBFEM and FEM
set (MainRegular mainregularsolution.cpp Common.cpp ${headers})
add_executable(SBFem2DRegular ${MainRegular})
target_link_libraries(SBFem2DRegular pz)

# Flow simulation - parabolic PDE
set (MainParabolic main2DParabolicEq.cpp Common.cpp ${headers})
add_executable(SBFem2DParabolic ${MainParabolic})
target_link_libraries(SBFem2DParabolic pz)

# Steklov problem
set (MainCircle maincircle.cpp Common.cpp ${headers})
add_executable(SBFem2DCircle ${MainCircle})
target_link_libraries(SBFem2DCircle pz)

# Heterogeneous circle domain
set (MainHeterogeneous mainheterogeneous.cpp Common.cpp ${headers})
add_executable(SBFem2DHeterogeneous ${MainHeterogeneous})
target_link_libraries(SBFem2DHeterogeneous pz)

# SBFem2D with scaling center at the boundary of the domain and Dirichlet imposed
set(MainSqrt mainsqrt.cpp Common.cpp ${headers})
add_executable(SBFem2DSqrt ${MainSqrt})
target_link_libraries(SBFem2DSqrt pz)

# Single crack
set(MainStressFiles mainstressintensity.cpp Common.cpp ${headers})
add_executable(SBFem2DStressIntensity ${MainStressFiles})
target_link_libraries(SBFem2DStressIntensity pz)

# Target for the paper Error Estimates
set(MainSqrtSquare mainsqrtsquare.cpp Common.cpp ${headers})
add_executable(SBFem2DSqrtSquare ${MainSqrtSquare})
target_link_libraries(SBFem2DSqrtSquare pz)

set (MainCrack2DFiles mainEmbeddedFractures.cpp Common.cpp Common3D.cpp ${headers})
add_executable(SBFem2DEmbeddedFractures ${MainCrack2DFiles})
target_link_libraries(SBFem2DEmbeddedFractures pz)
configure_file(TwoOrthogonalCracks.txt TwoOrthogonalCracks.txt COPYONLY)

set(MainSBFemFemFiles mainsbfem_fem.cpp Common.cpp ${headers})
add_executable(SBFem2DFem ${MainSBFemFemFiles})
target_link_libraries(SBFem2DFem pz)

set(MainHexagon hexagon.cpp Common.cpp Common3D.cpp ${headers})
add_executable(SBFemHexagon ${MainHexagon})
target_link_libraries(SBFemHexagon pz)

set(MainPolygons mainpolygons.cpp Common.cpp ${headers})
add_executable(SBFemPolygons ${MainPolygons})
target_link_libraries(SBFemPolygons pz)

configure_file(polygon1.txt polygon1.txt COPYONLY)
configure_file(polygon2.txt polygon2.txt COPYONLY)
configure_file(polygon3.txt polygon3.txt COPYONLY)
configure_file(polygon4.txt polygon4.txt COPYONLY)
configure_file(polygon5.txt polygon5.txt COPYONLY)

# Deprecated
#set(MymainOLD Main.cpp MainJSon.cpp Common.cpp ${headers})
#add_executable(SBFem2D ${MymainOLD})
#target_link_libraries(SBFem2D pz)

#set(MainCheckerboard maincheckerboard.cpp Common.cpp MainJSon.cpp ${headers})
#add_executable(SBFem2DCheckerboard ${MainCheckerboard})
#target_link_libraries(SBFem2DCheckerboard pz)

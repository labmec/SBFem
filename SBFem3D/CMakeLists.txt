file(GLOB headers *.h)
file(GLOB sources *.cpp)


set (MainRegular mainregular3dsolution.cpp Common3D.cpp ${headers})
add_executable(SBFem3DRegular ${MainRegular})
target_link_libraries(SBFem3DRegular PUBLIC NeoPZ::pz)

set (MainRegularFEM mainregular3Dfem.cpp Common3D.cpp ${headers})
add_executable(SBFem3DRegularFEM ${MainRegularFEM})
target_link_libraries(SBFem3DRegularFEM PUBLIC NeoPZ::pz)

set (MainCook maincooksimulation.cpp Common3D.cpp ${headers})
add_executable(SBFem3DCook ${MainCook})
target_link_libraries(SBFem3DCook PUBLIC NeoPZ::pz)

set (MainSphere mainspheresimulation.cpp Common3D.cpp ${headers})
add_executable(SBFem3DSphere ${MainSphere})
target_link_libraries(SBFem3DSphere PUBLIC NeoPZ::pz)

set (MainDragon maindragonsimulation.cpp Common3D.cpp ${headers})
add_executable(SBFem3DDragon ${MainDragon})
target_link_libraries(SBFem3DDragon PUBLIC NeoPZ::pz)

set (MainShangai mainshangai.cpp Common3D.cpp ${headers})
add_executable(SBFem3DShangai ${MainShangai})
target_link_libraries(SBFem3DShangai PUBLIC NeoPZ::pz)

set (MainPolygons mainpolygons.cpp Common3D.cpp ${headers})
add_executable(SBFem3DPolygons ${MainPolygons})
target_link_libraries(SBFem3DPolygons PUBLIC NeoPZ::pz)

configure_file(dragon_remesh_sbfemesh_128.txt dragon_remesh_sbfemesh_128.txt COPYONLY)
configure_file(dragon_sbfemesh_256.txt dragon_sbfemesh_256.txt COPYONLY)

configure_file(CooksMembrane_poly_16_1_1.txt CooksMembrane_poly_16_1_1.txt COPYONLY)
configure_file(CooksMembrane_sbfemesh_32_2_1.txt CooksMembrane_sbfemesh_32_2_1.txt COPYONLY)
configure_file(CooksMembrane_sbfemesh_64_2_1.txt CooksMembrane_sbfemesh_64_2_1.txt COPYONLY)
configure_file(CooksMembrane_sbfemesh_128_4_1.txt CooksMembrane_sbfemesh_128_4_1.txt COPYONLY)

configure_file(spheres_10_50_sbfemesh_64_8_1.txt spheres_10_50_sbfemesh_64_8_1.txt COPYONLY)
configure_file(spheres_10_50_sbfemesh_32_8_1.txt spheres_10_50_sbfemesh_32_8_1.txt COPYONLY)
configure_file(spheres_10_50_sbfemesh_128_8_1.txt spheres_10_50_sbfemesh_128_8_1.txt COPYONLY)

configure_file(Shanghai_Oriental_Pearl_Building_sbfemesh_256.txt Shanghai_Oriental_Pearl_Building_sbfemesh_256.txt COPYONLY)

configure_file(n64-id1-3.msh n64-id1-3.msh COPYONLY)
configure_file(n8-id1-2.msh n8-id1-2.msh COPYONLY)
configure_file(n512-id1-3.msh n512-id1-3.msh COPYONLY)

file(GLOB headers *.h)

# First test of the Hdiv code
set (FirstTest mainregularsolution.cpp ${headers})
add_executable(SBFemHdivV0 ${FirstTest})
target_link_libraries(SBFemHdivV0 pz)
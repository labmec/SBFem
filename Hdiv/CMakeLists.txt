file(GLOB headers *.h)

set(Material hybridpoissoncollapsed.h hybridpoissoncollapsed.cpp)

# First test of the Hdiv code
set (Trial mainregular.cpp SBFemHdiv.cpp hybridpoissoncollapsed.cpp ${headers})
add_executable(SBFemHdivV0 ${Trial})
target_link_libraries(SBFemHdivV0 pz)

# First test of the Hdiv code
set (FirstTest mainregularsolution.cpp ${Material})
add_executable(SBFemHdivV1 ${FirstTest})
target_link_libraries(SBFemHdivV1 pz)
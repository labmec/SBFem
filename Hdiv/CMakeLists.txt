file(GLOB headers *.h)

set(Material hybridpoissoncollapsed.h hybridpoissoncollapsed.cpp)

# Deprecated
# set (Trial mainregular.cpp SBFemHdiv.cpp hybridpoissoncollapsed.cpp ${headers})
# add_executable(SBFemHdivV0 ${Trial})
# target_link_libraries(SBFemHdivV0 PUBLIC NeoPZ::pz)

# First test of the Hdiv code
set (FirstTest mainregularsolution.cpp ${Material})
add_executable(SBFemHdivV1 ${FirstTest})
target_link_libraries(SBFemHdivV1 PUBLIC NeoPZ::pz)

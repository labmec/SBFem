#ifndef COMMON3DHPP
#define COMMON3DHPP

#include "pzanalysis.h"
#include "pzcmesh.h"
#include "TPZAnalyticSolution.h"

#ifdef USING_BOOST
#include "boost/crc.hpp"


static void printvec(const std::string &name, TPZVec<boost::crc_32_type::value_type> &vec)
{
    std::ofstream out(name);
    int64_t nel = vec.size();
    for (int64_t el=0; el<nel; el++) {
        if(vec[el] != 0)
        {
            out << el << " " << vec[el] << std::endl;
        }
    }
}

#endif

#ifdef _AUTODIFF
extern TLaplaceExample1 ExactLaplace;
extern TElasticity3DAnalytic ExactElast;
#endif

//    This Solve Different analysis
void SolveSist3D(TPZAnalysis *an, TPZCompMesh *fCmesh, int numthreads);

/// insert material objects in the computational mesh
void InsertMaterialObjects3D(TPZCompMesh *cmesh, bool scalarproblem);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupSquareMesh3D(int nelx, int nrefskeleton, int porder, bool elasticityproblem);

enum MMATID {Ebc1 = -1,Enomat, Emat1, Emat2, Emat3, Emat4, Ebc2, Ebc3, Ebc4, Ebc5, Ebcpoint1, Ebcpoint2, Ebcpoint3, Ewrap, ESkeleton, EInterfaceMat1, EInterfaceMat2, EGroup};

#ifdef _AUTODIFF
/// Function defining the exact elasticity solution
inline void Elasticity_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    ExactElast.Solution(xv,val,deriv);
}
#endif

/// Read a UNSWSBFem file
TPZGeoMesh *ReadUNSWSBGeoFile(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices);

#endif

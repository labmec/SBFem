#ifndef COMMONHPP
#define COMMONHPP

#include "pzanalysis.h"
#include "pzcmesh.h"
#include "TPZAnalyticSolution.h"

#ifdef _AUTODIFF
extern TElasticity2DAnalytic ElastExact;
extern TLaplaceExampleTimeDependent TimeLaplaceExact;
extern TLaplaceExample1 LaplaceExact;
#endif

//    Setup the system of equations and invert
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads);

/// insert material objects in the computational mesh
void InsertMaterialObjects(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact);

// Geometry of a quadrilateral uniform mesh
TPZAutoPointer<TPZGeoMesh> SetupGeom(int nelx);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder, bool elasticityproblem, bool applyexact);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupSquareH1Mesh(int nelx, int porder, bool elasticityproblem, bool applyexact);

// Generate output files for geometry and comp mesh
void OutputGmshCmsh(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmesh);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupCrackedOneElement(int nrefskeleton, int porder, bool applyexact, bool elastic);

enum MMATID {Enomat, Emat1, Emat2, Emat3, Emat4, Ebc1, Ebc2, Ebc3, Ebc4, EBCPoint1, EBCPoint2, Ewrap, ESkeleton, EInterfaceMat1, EInterfaceMat2, EGroup};

#ifdef _AUTODIFF
/// Function defining the exact elasticity solution
inline void Elasticity_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    ElastExact.Solution(xv, val, deriv);
}

inline void Laplace_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    LaplaceExact.Solution(xv, val, deriv);
}

inline void TimeLaplace_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    TimeLaplaceExact.Solution(xv, val, deriv);
}
#endif

/// Read a JSon File and generate a computational mesh
TPZCompMesh *ReadJSonFile(const std::string &filename, int numrefskeleton, int pOrder, REAL contrast);

/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZCompMesh *cmesh);

void PostProcessing(TPZAnalysis &Analysis, const std::string &filename, bool scalarproblem, int numthreads, int POrder, int nelxcount, int irefskeleton);

void PrintEigval(TPZAnalysis Analysis, std::string &filename);

TPZGeoMesh *ReadUNSWQuadtreeMesh(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices);
#endif

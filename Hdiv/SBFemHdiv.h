#pragma once

#ifdef HAVE_CONFIG_H
#include "pz_config.h"
#endif

#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "pzcmesh.h"
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "TPZAnalyticSolution.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "mixedpoisson.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZVecL2.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "TPZCompElHDivCollapsed.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"
#include "hybridpoissoncollapsed.h"
#include "TPZCompElHDivSBFem.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZNullMaterial.h"
#include "pzelchdivbound2.h"
#include "pzelmat.h"

#ifdef _AUTODIFF
extern TLaplaceExample1 LaplaceExact;
#endif

extern bool gHybrid, gSbfemhdiv, gHdivcollapsed;

enum MMATID {Enull, Emat0, Emat1, Ebc1, Ebc2, Ebc3, Ebc4, ESkeleton, Eleftflux, Erightflux, Eleftpressure, Erightpressure, Eint};

inline void Laplace_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
#ifdef _AUTODIFF
    LaplaceExact.Solution(xv, val, deriv);
#endif
}

TPZGeoMesh * SetupSquareMesh(int nelx);
TPZGeoMesh * SetupCollapsedMesh(int nelx);

TPZCompMesh * flux(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder);
TPZCompMesh * fluxhdivsbfem(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder);

TPZCompMesh * pressure(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder);

TPZMultiphysicsCompMesh * multiphysics(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder);
TPZMultiphysicsCompMesh *  multiphysicscollapsed(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder);

void AddInterfaceElements(TPZCompMesh * cmesh);

void AssociateElements(TPZCompMesh *cmesh, TPZManVector<int64_t> &elementgroup);
void GroupandCondense(TPZCompMesh * cmesh);

void AdjustExtPressureConnectivity(TPZCompMesh * cmeshm, TPZCompMesh * cmeshf, TPZManVector<int64_t> &perm);

void ComputeMatrices(TPZFMatrix<REAL> &E0fMat, TPZFMatrix<REAL> &E1fMat, TPZFMatrix<REAL> &E2fMat, TPZMatrix<REAL> &KCond);

void ComputeEigenvalues(TPZFMatrix<STATE> &E0fMat, TPZFMatrix<STATE> &E1fMat, TPZFMatrix<STATE> &E2fMat);
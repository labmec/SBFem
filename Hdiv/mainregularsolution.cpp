#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "pzgeoelbc.h"

#include "TPZVTKGeoMesh.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZBuildSBFemHdiv.h"

#include "TPZBndCondT.h"
#include "TPZMaterial.h"
#include "TPZNullMaterial.h"
#include "hybridpoissoncollapsed.h"
#include "TPZLagrangeMultiplier.h"

#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

#include "TPZAnalyticSolution.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
TLaplaceExample1 LaplaceExact;
#endif

using namespace std;

enum EMat {Emat1, Egroup, Ebc1, Ebc2, Ebc3, Ebc4, ESkeleton};

// Geometric mesh
TPZGeoMesh * GMeshRegular(int nelx, int irefskeleton, TPZManVector<int64_t> &scalingcenter, TPZManVector<int64_t> &elpartition);

// Computational mesh
TPZCompMesh * cmeshpressure(TPZAutoPointer<TPZGeoMesh> & gmesh, int POrder);
TPZCompMesh * cmeshflux(TPZAutoPointer<TPZGeoMesh> & gmesh, int POrder);
TPZMultiphysicsCompMesh *  cmeshmultiphysics(TPZAutoPointer<TPZGeoMesh> & gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder);

TLaplaceExample1 LaplaceExact;

int main(int argc, char *argv[])
{
    // Initial data
    auto minnelxcount = 1, maxnelxcount = 5;
    auto minrefskeleton = 0, maxrefskeleton = 1;
    auto usesbfem = true; // false for FEM simulations

    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;

    int countstep = 1;
    for (int porder = 1; porder < 4; porder++)
    {
        for(int nelxcount = minnelxcount; nelxcount < maxnelxcount; nelxcount ++)
        {
            int nelx = 1 << (nelxcount-1);
            
            // Creating geometric mesh
            TPZManVector<int64_t> scalingcenter, elpartition;
            TPZAutoPointer<TPZGeoMesh> gmesh = GMeshRegular(nelx, 0, scalingcenter, elpartition);

            // Creating pressure mesh
            auto cmeshp = cmeshpressure(gmesh, porder);

            // Creating flux mesh
            auto cmeshf = cmeshflux(gmesh, porder);

            // Creating multiphysics mesh
            TPZMultiphysicsCompMesh * cmeshm = cmeshmultiphysics(gmesh, cmeshp, cmeshf, porder);
            std::ofstream mout("teste.txt");
            cmeshm->Print(mout);
            
            map<int,int> matmap;
            matmap[Egroup] = Emat1;

            TPZBuildSBFemHdiv build(gmesh, ESkeleton, matmap);
            build.StandardConfiguration();
            build.BuildMultiphysicsCompMesh(*cmeshm);

#ifdef PZDEBUG
            std::ofstream gout("GeometrySBFEM.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(cmeshm->Reference(), gout, true);
            std::ofstream sout("CmeshSBFEM.txt");
            cmeshm->Print(sout);
#endif
            std::cout << "Analysis...\n";
            auto neq = cmeshm->NEquations();
            std::cout << "neq = " << neq << std::endl;
            TPZLinearAnalysis an(cmeshm, false);
            TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(cmeshm);
            // strmat.SetNumThreads(4);
            an.SetStructuralMatrix(strmat);

            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            an.SetSolver(step);

            an.Run();

            // {
                ofstream ssout("globalmatrices.txt");
            //     auto ssolv = dynamic_cast<TPZStepSolver<STATE> *>(an.Solver());
            //     ssolv->Matrix()->Print("ekglob = ", sout, EMathematicaInput);
            //     // an.Rhs().Print("rhs = ", sout, EMathematicaInput);
                auto sol = cmeshm->Solution();
                sol.Print("sol", ssout, EMathematicaInput);
            // }

            {
                ofstream sout("postprocessing.vtk");
                TPZStack<std::string> vecnames,scalnames;
                scalnames.Push("Pressure");
                vecnames.Push("Flux");
                an.DefineGraphMesh(2, scalnames, vecnames, "postprocessing.vtk");
                int res = porder+1;
                if (res >5) {
                    res = 5;
                }
                an.PostProcess(res);
            }
            an.SetExact(LaplaceExact.ExactSolution());

            TPZManVector<REAL> ervec(6,0.);
            an.PostProcessError(ervec, false);
        }
    }
    cout << "Check:: Calculation finished successfully" << endl;
    return EXIT_SUCCESS;
}

TPZGeoMesh * GMeshRegular(int nelx, int irefskeleton, TPZManVector<int64_t> &scalingcenters, TPZManVector<int64_t> &elpartitions)
{
    // FEM GEO MESH:
    TPZGeoMesh * gmesh = new TPZGeoMesh();

    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);

    // Setting boundary elements
    {
        gengrid.Read(gmesh, Egroup);
        gengrid.SetBC(gmesh, 4, Ebc1);
        gengrid.SetBC(gmesh, 5, Ebc1);
        gengrid.SetBC(gmesh, 6, Ebc1);
        gengrid.SetBC(gmesh, 7, Ebc1);
    }
    gmesh->BuildConnectivity();

#ifdef PZDEBUG
    std::ofstream out("GeometryFEM.vtk");
    TPZVTKGeoMesh gvtk;
    gvtk.PrintGMeshVTK(gmesh, out, true);
#endif

    return gmesh;
}

TPZCompMesh * cmeshpressure(TPZAutoPointer<TPZGeoMesh> & gmesh, int POrder)
{        
    auto dim = 1; auto nstate = 1;

    auto cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    
    // Volumetric material
    auto mat = new TPZNullMaterial(Emat1, 1, nstate);
    cmesh->InsertMaterialObject(mat);

    // Boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        bcond->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, ESkeleton, 1, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    set<int> matid = {Emat1, Ebc1, ESkeleton};
    cmesh->AutoBuild(matid);

    for(auto newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZCompMesh * cmeshflux(TPZAutoPointer<TPZGeoMesh> & gmesh, int POrder)
{
    auto dim = 2; auto nstate = 1;
    auto cmeshcollapsed = new TPZCompMesh(gmesh);
    cmeshcollapsed->SetDefaultOrder(POrder);
    cmeshcollapsed->SetDimModel(dim);
    cmeshcollapsed->CleanUp();
    
    // Volumetric material
    auto mat = new TPZNullMaterial(Emat1, 1, nstate);
    cmeshcollapsed->InsertMaterialObject(mat);

    cmeshcollapsed->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshcollapsed->AutoBuild();
    
    return cmeshcollapsed;
}

TPZMultiphysicsCompMesh *  cmeshmultiphysics(TPZAutoPointer<TPZGeoMesh> & gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder)
{
    int dim = gmesh->Dimension();
    int nstate = 1;

    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    auto mat = new TPZHybridPoissonCollapsed(Emat1,2);
    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);

    // auto matskeleton = new TPZMixedDarcyFlow(ESkeleton,1);
    // mat->SetPermeabilityFunction(1.);
    // cmesh->InsertMaterialObject(matskeleton);

    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        bcond->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, ESkeleton, 1, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }

    cout << "Creating multiphysics mesh\n";
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshf;
    meshvector[1] = cmeshp;

    TPZManVector<int> active(2,1);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}
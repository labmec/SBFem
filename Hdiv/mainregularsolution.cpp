#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "pzgeoelbc.h"

#include "TPZVTKGeoMesh.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZBuildSBFemHdiv.h"

#include "pzbndcond.h"
#include "TPZMaterial.h"
#include "TPZNullMaterial.h"
#include "hybridpoissoncollapsed.h"
#include "TPZLagrangeMultiplier.h"

#include "pzanalysis.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

using namespace std;

enum EMat {Emat1, Egroup, Ebc1, Ebc2, Ebc3, Ebc4, ESkeleton};

// Geometric mesh
TPZAutoPointer<TPZGeoMesh> GMeshRegular(int nelx, int irefskeleton, TPZManVector<int64_t> &scalingcenter, TPZManVector<int64_t> &elpartition);

// Computational mesh
TPZCompMesh * cmeshpressure(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder);
TPZCompMesh * cmeshflux(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder);
TPZMultiphysicsCompMesh *  cmeshmultiphysics(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initial data
    auto minnelxcount = 1, maxnelxcount = 5;
    auto minrefskeleton = 0, maxrefskeleton = 4;
    auto usesbfem = true; // false for FEM simulations
    auto porder = 1;

    int countstep = 1;
    for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
    {
        for(int nelxcount = minnelxcount; nelxcount < maxnelxcount; nelxcount ++)
        {
            int nelx = 1 << (nelxcount-1);
            
            // Creating geometric mesh
            TPZManVector<int64_t> scalingcenter, elpartition;
            auto gmesh = GMeshRegular(nelx, irefskeleton, scalingcenter, elpartition);

            // Creating pressure mesh
            auto cmeshp = cmeshpressure(gmesh, porder);

            // Creating flux mesh
            auto cmeshf = cmeshflux(gmesh, porder);

            // Creating multiphysics mesh
            TPZCompMesh * cmeshm = cmeshmultiphysics(gmesh, cmeshp, cmeshf, porder);
            
            map<int,int> matmap;
            matmap[Egroup] = Emat1;

            TPZBuildSBFemHdiv build(gmesh, ESkeleton, matmap);
            build.StandardConfiguration();
            build.BuildMultiphysicsCompMesh(*cmeshm);

#ifdef PZDEBUG
            std::ofstream gout("GeometrySBFEM.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(cmeshm->Reference(), gout, true);
#endif
            // TPZAnalysis an(cmeshm);


        }
    }
    cout << "Check:: Calculation finished successfully" << endl;
    return EXIT_SUCCESS;
}

TPZAutoPointer<TPZGeoMesh> GMeshRegular(int nelx, int irefskeleton, TPZManVector<int64_t> &scalingcenters, TPZManVector<int64_t> &elpartitions)
{
    // FEM GEO MESH:
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh();

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
        gengrid.SetBC(gmesh, 5, Ebc2);
        gengrid.SetBC(gmesh, 6, Ebc3);
        gengrid.SetBC(gmesh, 7, Ebc4);
    }
    gmesh->BuildConnectivity();

#ifdef PZDEBUG
    std::ofstream out("GeometryFEM.vtk");
    TPZVTKGeoMesh gvtk;
    gvtk.PrintGMeshVTK(gmesh, out, true);
#endif

    return gmesh;
}

TPZCompMesh * cmeshpressure(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{        
    auto dim = 1; auto nstate = 1;

    auto cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    
    auto mat = new TPZNullMaterial(Egroup, dim, nstate);
    cmesh->InsertMaterialObject(mat);
    // TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    // {
    //     auto bcond = mat->CreateBC(mat, ESkeleton, 1, val1, val2);
    //     cmesh->InsertMaterialObject(bcond);
    // }
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();

    for(auto newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZCompMesh * cmeshflux(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = gmesh->Dimension(); auto nstate = 1;
    auto cmeshcollapsed = new TPZCompMesh(gmesh);
    cmeshcollapsed->SetDefaultOrder(POrder);
    cmeshcollapsed->SetDimModel(dim);
    cmeshcollapsed->CleanUp();

    auto mat = new TPZNullMaterial(Egroup, dim, nstate);
    cmeshcollapsed->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, ESkeleton, 1, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }

    cmeshcollapsed->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshcollapsed->AutoBuild(); 

#ifdef PZDEBUG
    std::ofstream cout("CMeshFlux.vtk");
    TPZVTKGeoMesh cvtk;
    cvtk.PrintCMeshVTK(cmeshcollapsed, cout, true);
#endif
    
    return cmeshcollapsed;
}

TPZMultiphysicsCompMesh *  cmeshmultiphysics(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder)
{
    int dim = gmesh->Dimension();
    int nstate = 1;

    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    auto mat = new TPZHybridPoissonCollapsed(Egroup,dim);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2); 
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2); 
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2); 
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2); 
        cmesh->InsertMaterialObject(bcond);
    }
    // {
    //     TPZMaterial * bcond = mat->CreateBC(mat, ESkeleton, 1, val1, val2); 
    //     cmesh->InsertMaterialObject(bcond);
    // }

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
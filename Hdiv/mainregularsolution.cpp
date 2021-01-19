#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "TPZBuildSBFemHdiv.h"

#include "TPZMatLaplacian.h"
#include "pzbndcond.h"
#include "TPZMaterial.h"

#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

using namespace std;

enum EMat {Emat1, Egroup, Ebc1, Ebc2, Ebc3, Ebc4, Eskeleton};

TPZAutoPointer<TPZGeoMesh> GMeshRegular(int nelx, int irefskeleton);
TPZCompMesh * SBFemHdivMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int porder);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initial data
    auto minnelxcount = 2, maxnelxcount = 5;
    auto minrefskeleton = 0, maxrefskeleton = 4;
    auto usesbfem = true; // false for FEM simulations
    auto porder = 1;

    int countstep = 1;
    for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
    {
        for(int nelxcount = minnelxcount; nelxcount < maxnelxcount; nelxcount ++)
        {
            int nelx = 1 << (nelxcount-1);
            
            auto gmesh = GMeshRegular(nelx, irefskeleton);

            auto cmesh = SBFemHdivMesh(gmesh, porder);
        }
    }
    cout << "Check:: Calculation finished successfully" << endl;
    return EXIT_SUCCESS;
}

TPZAutoPointer<TPZGeoMesh> GMeshRegular(int nelx, int nrefskeleton)
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

    std::ofstream out("GeometryFEM.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintGMeshVTK(gmesh, out, true);

    return gmesh;
}

TPZCompMesh * SBFemHdivMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int porder)
{
    // Defining configuration for SBFEM mesh
    std::map<int, int> matmap;
    matmap[Egroup] = Emat1;
    TPZBuildSBFemHdiv build(gmesh, Eskeleton, matmap);

    build.StandardConfiguration();

    // Defining computational mesh and material data
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);

    TPZMatLaplacian *material = new TPZMatLaplacian(Emat1);
    material->SetDimension(gmesh->Dimension());
    material->SetSymmetric();
    SBFem->InsertMaterialObject(material);

    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    TPZMaterial *BCond1 = material->CreateBC(material, Ebc1, 0, val1, val2);
    TPZMaterial *BCond2 = material->CreateBC(material, Ebc2, 0, val1, val2);
    TPZMaterial *BCond3 = material->CreateBC(material, Ebc3, 0, val1, val2);
    TPZMaterial *BCond4 = material->CreateBC(material, Ebc4, 0, val1, val2);
    TPZMaterial *BCSkeleton = material->CreateBC(material, Eskeleton, 0, val1, val2);
    SBFem->InsertMaterialObject(BCond1);
    SBFem->InsertMaterialObject(BCond2);
    SBFem->InsertMaterialObject(BCond3);
    SBFem->InsertMaterialObject(BCond4);
    SBFem->InsertMaterialObject(BCSkeleton);

    // Generating SBFEM mesh
    build.BuildComputationMeshHdiv(*SBFem);

    std::ofstream out("GeometrySBFEM.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintGMeshVTK(SBFem->Reference(), out, true);

    return SBFem;
}
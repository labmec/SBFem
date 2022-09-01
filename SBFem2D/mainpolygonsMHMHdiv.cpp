#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZMultiphysicsCompMesh.h"
#include "pzgeoelbc.h"

#include "Elasticity/TPZElasticity2D.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"

#include "TPZVTKGeoMesh.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"

#include "tpzintpoints.h"
#include "pzgeoelrefless.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "TPZRefPatternDataBase.h"
#include "TPZAnalyticSolution.h"


void AddBoundaryElements(TPZAutoPointer<TPZGeoMesh> & gmesh);

void CreateTriangleElements(TPZAutoPointer<TPZGeoMesh> gmesh, std::map<int,int> &matmap, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenter);

/// divide the triangles and boundary elements
void DivideTriangles(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition);

/// divide the triangles and boundary elements
void DivideSkeleton(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition, int matid);

/// generate a mesh with HDiv elements
TPZCompMesh *GenerateFluxMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder_vol, int porder_skel);

/// generate a mesh with L2 elements
TPZCompMesh *GeneratePressureMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, TPZVec<int64_t> &elpartition);

/// generate a mesh with constant elements
TPZCompMesh *GenerateConstantMesh(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition, int level);

/// adjust the element partition of the boundary elements
void AdjustElPartition(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition);

/// build and MHM multiphysics mesh
TPZMultiphysicsCompMesh *BuildMultiphysicsMesh(int POrder_vol, int POrder_skel, TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition);

void InsertMaterialObjects(TPZMultiphysicsCompMesh &cmesh);

void PutinSubstructures(TPZCompMesh &cmesh, TPZVec<int64_t> &elpartition);

void CondenseElements(TPZCompMesh &cmesh);

void SolveSistPolygons(TPZLinearAnalysis *an, TPZCompMesh *Cmesh, int numthreads);

int average_pressure_lagrange_level = 5;
int distribute_flux_lagrange_level = 4;
int main(int argc, char *argv[])
{
#if PZ_LOG
    TPZLogger::InitializePZLOG();
#endif // PZ_LOG
    // Initial data
    int minnelxcount = 1, maxnelxcount = 2;
    int minporder = 1, maxporder = 4;
    int maxrefskeleton = 4;
    int numthreads = 6;
    bool scalarproblem = true; // false for elasticity 2D problems
    bool usesbfem = true; // false for FEM simulations
    bool useexact = true;
    LaplaceExact.fExact = TLaplaceExample1::E2SinSin;
    int counter = 1;
    
    std::string filename;
    for ( int POrder = 1; POrder < 2; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < 4; irefskeleton++) {
            for (int nelx = 1; nelx <=1 ; nelx++) {
                if(nelx == 1) filename = "polygon1.txt";
                else if (nelx == 2) filename = "polygon2.txt";
                else if (nelx == 3) filename = "polygon3.txt";
                else if (nelx == 4) filename = "polygon4.txt";
                else if (nelx == 5) filename = "polygon5.txt";
                else DebugStop();
    #ifdef MACOSX
                filename = "../"+filename;
    #endif
                std::string vtkfilename;
                std::string vtkfilegeom;
                std::string vtkfilecmesh;
                std::string rootname;
                std::string boundaryname;
                {
                    int pos = filename.find(".txt");
                    std::string truncate = filename;
                    truncate.erase(pos);
                    rootname = truncate;
                    std::stringstream sout;
                    sout << truncate << "_p" << POrder << ".vtk";
                    vtkfilename = sout.str();
                    std::stringstream boundstr;
                    boundstr << truncate << "_boundary";
                    boundaryname = boundstr.str();
                    std::stringstream vtkgeom;
                    vtkgeom << truncate << "_"<< irefskeleton << "_geom.vtk";
                    vtkfilegeom = vtkgeom.str();
                    std::stringstream vtkcmesh;
                    vtkcmesh << truncate << "_cmesh.vtk";
                    vtkfilecmesh = vtkcmesh.str();
                }
                
                TPZVec<int64_t> elpartition;
                TPZVec<int64_t> scalingcenterindices;
                TPZAutoPointer<TPZGeoMesh> gmesh = ReadUNSWQuadtreeMesh(filename, elpartition, scalingcenterindices);
                
                AddBoundaryElements(gmesh);
                elpartition.Resize(gmesh->NElements(), -1);
                scalingcenterindices.Resize(gmesh->NElements(), -1);

                std::cout << "Building computational mesh\n";
                std::map<int,int> matmap;
                matmap[ESkeleton] = Emat1;
                int EPoly = 100;
                matmap[EPoly] = Emat2;
                
                CreateTriangleElements(gmesh, matmap, elpartition, scalingcenterindices);
                for (int i=0; i< irefskeleton+1; i++) {
                    DivideTriangles(gmesh, elpartition);
                }
                for (int i=0; i<irefskeleton; i++) {
                    DivideSkeleton(gmesh, elpartition, 100);
                }
                if(1) {
                    std::cout << "Plotting the geometric mesh\n";
                    std::ofstream outtxt("gmesh.txt");
                    gmesh->Print(outtxt);
                    std::ofstream out(vtkfilegeom);
                    TPZVTKGeoMesh vtk;
                    vtk.PrintGMeshVTK(gmesh.operator->(), out,elpartition);
                }
                int POrder_vol = POrder+1;
                int POrder_skel = POrder;
                auto cmesh = BuildMultiphysicsMesh(POrder_vol, POrder_skel, gmesh, elpartition);
                if(0)
                {
                    std::ofstream out("mphysics.txt");
                    cmesh->Print(out);
                }
                PutinSubstructures(*cmesh, elpartition);
                CondenseElements(*cmesh);
                TPZCompMesh *SBFem = cmesh;
                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
                Analysis->SetStep(counter++);

                SolveSistPolygons(Analysis, SBFem, numthreads);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;

                std::cout << "Post processing\n";
                std::string sout;
                sout.append("../PolygonsSolution_MHM");
                
                PostProcessing(*Analysis, sout, scalarproblem, numthreads, POrder, nelx, irefskeleton);
                TPZVec<TPZCompMesh *> meshvec = cmesh->MeshVector();
                delete SBFem;
                for (int i = 0; i<meshvec.size(); i++) {
                    delete meshvec[i];
                }
                delete Analysis;
            }
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void AddBoundaryElements(TPZAutoPointer<TPZGeoMesh> & gmesh)
{
    std::set<int64_t> setbottom,setright,settop,setleft;
    int64_t nnodes = gmesh->NNodes();
    int dim = gmesh->Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[in].GetCoordinates(xco);
        if (fabs(xco[1]+1.) < 1.e-3)
        {
            setbottom.insert(in);
        }
        if (fabs(xco[0]-1.) < 1.e-3)
        {
            setright.insert(in);
        }
        if (fabs(xco[1]-1.) < 1.e-3)
        {
            settop.insert(in);
        }
        if (fabs(xco[0]+1.) < 1.e-3)
        {
            setleft.insert(in);
        }
    }
    int64_t nelem = gmesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            int nsidenodes = gel->NSideNodes(is);
            int nfoundbottom = 0;
            int nfoundright = 0;
            int nfoundtop = 0;
            int nfoundleft = 0;
            for (int in=0; in<nsidenodes; in++) {
                int64_t nodeindex = gel->SideNodeIndex(is, in);
                if (setbottom.find(nodeindex) != setbottom.end()) {
                    nfoundbottom++;
                }
                if (setright.find(nodeindex) != setright.end()) {
                    nfoundright++;
                }
                if (settop.find(nodeindex) != settop.end()) {
                    nfoundtop++;
                }
                if (setleft.find(nodeindex) != setleft.end()) {
                    nfoundleft++;
                }
            }
            if (nfoundbottom == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc1);
            }
            if (nfoundright == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc2);
            }
            if (nfoundtop == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc3);
            }
            if (nfoundleft == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc4);
            }
            else
            {
                TPZGeoElSide gelside(gel,is);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if (neighbour == gelside) {
                    int EPoly = 100;
                    TPZGeoElBC(gelside,EPoly);
                }
            }
        }
    }
}


void SolveSistPolygons(TPZLinearAnalysis *an, TPZCompMesh *Cmesh, int numthreads)
{
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(Cmesh);
    strmat.SetNumThreads(numthreads);
    an->SetStructuralMatrix(strmat);
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    an->Assemble();
//    try {
//        an->Assemble();
//    } catch (...) {
//        exit(-1);
//    }
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
}

void CreateTriangleElements(TPZAutoPointer<TPZGeoMesh> gmesh, std::map<int,int> &matmap, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenter)
{
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t> extpartition(elpartition);
    extpartition.Expand(nel*3);
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int gelmat = gel->MaterialId();
        if(gel->Type() != EOned) DebugStop();
        auto partition = elpartition[el];
        if(partition >= 0)
        {
            if(matmap.find(gelmat) == matmap.end()) DebugStop();
            int64_t center = scalingcenter[partition];
            TPZManVector<int64_t,3> nodes(3);
            nodes[0] = gel->NodeIndex(0);
            nodes[1] = gel->NodeIndex(1);
            nodes[2] = center;
            int64_t index;
            gmesh->CreateGeoElement(ETriangle, nodes, matmap[gelmat], index);
            if(extpartition.size() < index+1) extpartition.resize(index+1);
            extpartition[index] = partition;
            // set the oned element partition as -1. It will serve as skeleton element
            extpartition[el] = -1;
        }
    }
    gmesh->BuildConnectivity();
    elpartition = extpartition;
}
void DivideTriangles(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition)
{
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t> extpartition(elpartition);
    extpartition.Expand(nel*3);
    static int once = 1;
    if(once) {
        gRefDBase.InitializeUniformRefPattern(EOned);
        gRefDBase.InitializeUniformRefPattern(ETriangle);
        once = 0;
    }
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int gelmat = gel->MaterialId();
        bool boundmat = false;
        if(gelmat == Ebc1 || gelmat == Ebc2 || gelmat == Ebc3 || gelmat == Ebc4) boundmat  = true;
        if(boundmat)
        {
            TPZGeoElSide gelside(gel);
            for(TPZGeoElSide neigh = gelside.Neighbour(); neigh != gelside; neigh++)
            {
                if(neigh.Element()->Dimension() == 2)
                {
                    elpartition[el] = elpartition[neigh.Element()->Index()];
                    break;
                }
            }
            auto skel = gelside.HasNeighbour(100);
            if(skel)
            {
                skel.Element()->RemoveConnectivities();
                delete skel.Element();
            }
        }
        if(gel->Type() != ETriangle && !boundmat) continue;
        auto partition = elpartition[el];
        if(partition < 0) DebugStop();
        TPZManVector<TPZGeoEl *> subels(4);
        gel->Divide(subels);
        int nsubels = subels.size();
        for (int sub=0; sub<nsubels; sub++) {
            int64_t index = subels[sub]->Index();
            if(extpartition.size() < index+1) extpartition.resize(index+1);
            extpartition[index] = partition;
        }
    }
    elpartition = extpartition;
}

void DivideSkeleton(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition, int matid)
{
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t> extpartition(elpartition);
    extpartition.Expand(nel*3);
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        int gelmat = gel->MaterialId();
        if(gelmat != matid) continue;
        auto partition = elpartition[el];
        if(partition >= 0) DebugStop();
        TPZManVector<TPZGeoEl *> subels(4);
        gel->Divide(subels);
        int nsubels = subels.size();
        for (int sub=0; sub<nsubels; sub++) {
            int64_t index = subels[sub]->Index();
            if(extpartition.size() < index+1) extpartition.resize(index+1);
            extpartition[index] = partition;
        }
    }
    elpartition = extpartition;
}

/// generate a mesh with HDiv elements
TPZCompMesh *GenerateFluxMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder_vol,int porder_skel)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    auto *mat = new TPZNullMaterial<>(Emat2,dim,1);
    cmesh->InsertMaterialObject(mat);
    mat = new TPZNullMaterial<>(Ebc1,dim-1,1);
    cmesh->InsertMaterialObject(mat);
    mat = new TPZNullMaterial<>(Ebc2,dim-1,1);
    cmesh->InsertMaterialObject(mat);
    mat = new TPZNullMaterial<>(Ebc3,dim-1,1);
    cmesh->InsertMaterialObject(mat);
    mat = new TPZNullMaterial<>(Ebc4,dim-1,1);
    cmesh->InsertMaterialObject(mat);
    mat = new TPZNullMaterial<>(100,dim-1,1);
    cmesh->InsertMaterialObject(mat);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    std::set<int> matvol = {Emat2,Ebc1,Ebc2,Ebc3,Ebc4};
    std::set<int> matskel = {100};
    cmesh->SetDefaultOrder(porder_vol);
    cmesh->AutoBuild(matvol);
    cmesh->SetDefaultOrder(porder_skel);
    cmesh->AutoBuild(matskel);
    cmesh->InitializeBlock();
    int64_t nc = cmesh->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        auto &c = cmesh->ConnectVec()[ic];
        c.SetLagrangeMultiplier(0);
    }
    // skeleton elements have lagrange multiplier 4
    int64_t nel = cmesh->NElements();
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() == Emat1)
        {
            cel->Connect(0).SetLagrangeMultiplier(4);
        }
    }
    gmesh->ResetReference();
    return cmesh;
}

/// generate a mesh with L2 elements
TPZCompMesh *GeneratePressureMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, TPZVec<int64_t> &elpartition)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    int dim = gmesh->Dimension();
    auto *mat = new TPZNullMaterial<>(Emat2,dim,1);
    cmesh->InsertMaterialObject(mat);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    int64_t ncon = cmesh->ConnectVec().NElements();
    for (int64_t ic = 0; ic<ncon; ic++) {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    
//    std::set<int64_t> partitionprocessed;
//    int64_t nel = cmesh->NElements();
//    for (int el = 0; el<nel; el++) {
//        TPZCompEl *cel = cmesh->Element(el);
//        TPZGeoEl *gel = cel->Reference();
//        int64_t partition = elpartition[gel->Index()];
//        if(partitionprocessed.find(partition) == partitionprocessed.end())
//        {
//            partitionprocessed.insert(partition);
//            cel->Connect(0).SetLagrangeMultiplier(5);
//        }
//    }
    gmesh->ResetReference();
    return cmesh;
}

/// generate a mesh with constant elements
TPZCompMesh *GenerateConstantMesh(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition, int level)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(0);
    int dim = gmesh->Dimension();
    auto *mat = new TPZNullMaterial<>(Emat2,dim,1);
    cmesh->InsertMaterialObject(mat);
    cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        auto index = gel->Index();
        cel->Connect(0).SetLagrangeMultiplier(level);
    }
    gmesh->ResetReference();
    return cmesh;

}

/// build and MHM multiphysics mesh
TPZMultiphysicsCompMesh *BuildMultiphysicsMesh(int POrder_vol, int POrder_skel, TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &elpartition)
{
    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    auto flux = GenerateFluxMesh(gmesh,POrder_vol,POrder_skel);
    auto pressure = GeneratePressureMesh(gmesh, POrder_vol,elpartition);
    auto distflux = GenerateConstantMesh(gmesh, elpartition,distribute_flux_lagrange_level);
    auto avpressure = GenerateConstantMesh(gmesh, elpartition,average_pressure_lagrange_level);
    if(0)
    {
        std::ofstream outflux("flux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(flux, outflux);
        std::ofstream outpressure("pressure.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(pressure, outpressure);
        std::ofstream outconstflux("constflux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(distflux, outconstflux);
        std::ofstream outconstpress("constpressure.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(avpressure, outconstpress);
    }
    if(0)
    {
        std::ofstream outflux("flux.txt");
        flux->Print(outflux);
        std::ofstream outpressure("pressure.txt");
        pressure->Print(outpressure);
        std::ofstream outconstflux("constflux.txt");
        distflux->Print(outconstflux);
        std::ofstream outconstpress("constpressure.txt");
        avpressure->Print(outconstpress);
    }
    TPZManVector<TPZCompMesh *> meshvec = {flux,pressure,distflux,avpressure};
    InsertMaterialObjects(*mphys);
    mphys->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    mphys->BuildMultiphysicsSpace(meshvec);
    return mphys;
}
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZBndCondT.h"

//int dim = gmesh->Dimension();
//auto *mat = new TPZNullMaterial<>(Emat2,dim,1);
//cmesh->InsertMaterialObject(mat);
//mat = new TPZNullMaterial<>(Ebc1,dim-1,1);
//cmesh->InsertMaterialObject(mat);
//mat = new TPZNullMaterial<>(Ebc2,dim-1,1);
//cmesh->InsertMaterialObject(mat);
//mat = new TPZNullMaterial<>(Ebc3,dim-1,1);
//cmesh->InsertMaterialObject(mat);
//mat = new TPZNullMaterial<>(Ebc4,dim-1,1);
//cmesh->InsertMaterialObject(mat);
//mat = new TPZNullMaterial<>(Emat1,dim-1,1);
//cmesh->InsertMaterialObject(mat);

void ForceOne(const TPZVec<REAL> &x, TPZVec<REAL> &force)
{
    force[0] = 1.;
}

void InsertMaterialObjects(TPZMultiphysicsCompMesh &cmesh)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    int dim = gmesh->Dimension();
    auto *darcy = new TPZMixedDarcyFlow(Emat2,dim);
    
    darcy->SetForcingFunction(LaplaceExact.ForceFunc(), 3);
    darcy->SetExactSol(LaplaceExact.ExactSolution(), 3);
    cmesh.InsertMaterialObject(darcy);
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZVec<STATE> val2(1,1.);
    auto *bnd1 = darcy->CreateBC(darcy, Ebc1, 0, val1, val2);
    auto *bnd1T = dynamic_cast<TPZBndCondT<STATE>*>(bnd1);
    bnd1T->SetForcingFunctionBC(LaplaceExact.ExactSolution(), 3);
    cmesh.InsertMaterialObject(bnd1);
    auto *bnd2 = darcy->CreateBC(darcy, Ebc2, 0, val1, val2);
    bnd2->SetForcingFunctionBC(LaplaceExact.ExactSolution(), 3);
    cmesh.InsertMaterialObject(bnd2);
    auto *bnd3 = darcy->CreateBC(darcy, Ebc3, 0, val1, val2);
    bnd3->SetForcingFunctionBC(LaplaceExact.ExactSolution(), 3);
    cmesh.InsertMaterialObject(bnd3);
    auto *bnd4 = darcy->CreateBC(darcy, Ebc4, 0, val1, val2);
    bnd4->SetForcingFunctionBC(LaplaceExact.ExactSolution(), 3);
    cmesh.InsertMaterialObject(bnd4);
    auto *nullmat = new TPZNullMaterialCS<>(100,dim,1);
    cmesh.InsertMaterialObject(nullmat);
}

#include "pzsubcmesh.h"

void PutinSubstructures(TPZCompMesh &cmesh, TPZVec<int64_t> &elpartition)
{
    // create subcompmeshes
    std::map<int64_t,TPZSubCompMesh *> submeshes;
    for(auto part : elpartition) {
        if(part == -1) continue;
        if(submeshes.find(part) == submeshes.end())
        {
            auto *sub = new TPZSubCompMesh(cmesh);
            submeshes[part] = sub;
        }
    }
    int64_t nel = cmesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        auto *cel = cmesh.Element(el);
        auto *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) continue;
        auto *gel = cel->Reference();
        auto index = gel->Index();
        auto part = elpartition[index];
        if(part >= 0) {
            auto *sub = submeshes[part];
            sub->TransferElement(&cmesh, el);
        }
    }
    cmesh.ComputeNodElCon();
    for(auto itsub : submeshes) {
        TPZCompEl *submesh = itsub.second;
        int64_t ncon = submesh->NConnects();
        for(int64_t ic = 0; ic<ncon; ic++) {
            auto &c = submesh->Connect(ic);
            if(c.LagrangeMultiplier() == average_pressure_lagrange_level) {
                c.IncrementElConnected();
                break;
            }
        }
    }
    for(auto itsub : submeshes) {
        TPZSubCompMesh *submesh = itsub.second;
        submesh->MakeAllInternal();
        submesh->CleanUpUnconnectedNodes();
    }
    cmesh.ComputeNodElCon();
    {
        int64_t ncon = cmesh.NConnects();
        for(int64_t ic = 0; ic<ncon; ic++) {
            auto &c = cmesh.ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency()) {
                c.RemoveDepend();
            }
        }
    }
    cmesh.CleanUpUnconnectedNodes();
    if(0)
    {
        std::ofstream out("cmesh_substruct.txt");
        cmesh.Print(out);
    }
}

#include "pzcondensedcompel.h"

void CondenseElements(TPZCompMesh &cmesh)
{
    cmesh.ComputeNodElCon();
    int64_t nel = cmesh.NElements();
    cmesh.ElementSolution().Resize(nel, 5);
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        auto *subcel = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcel) {
            CondenseElements(*subcel);
            subcel->CleanUpUnconnectedNodes();
            subcel->SetAnalysisSparse(2);
//            subcel->SetAnalysisSkyline();
        }
        else if(cel)
        {
            TPZGeoEl *gel = cel->Reference();
            if(gel && gel->Dimension() == 2) {
                // prevent the pressure connect from being condensed
                int nc = cel->NConnects();
                for (int ic = 0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    if(c.LagrangeMultiplier() == average_pressure_lagrange_level) {
                        c.IncrementElConnected();
                    }
                }
                auto *cond = new TPZCondensedCompEl(cel);
            }
        }
    }
    auto *sub = dynamic_cast<TPZSubCompMesh *>(&cmesh);
    if(0 && !sub)
    {
        std::ofstream out("CondensedCompMesh.txt");
        cmesh.Print(out);
    }
}



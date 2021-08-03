#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZBndCond.h"
#include "TPZVTKGeoMesh.h"
#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZGenGrid2D.h"
#include "pzgeoelbc.h"
#include "pzgnode.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExactLower;
TElasticity2DAnalytic ElastExactUpper;
TLaplaceExample1 LaplaceExactLower;
TLaplaceExample1 LaplaceExactUpper;
#endif

TPZCompMesh *SetupCrackedOneSquareElement(int nrefskeleton, int porder, bool applyexact, bool elastic);

TPZCompMesh * CreateCMesh(int nelx, int porder);

TPZCompMesh *BuildSBFem(TPZAutoPointer<TPZGeoMesh> gmesh, int nx, int porder);

void AddBoundaryElements(TPZGeoMesh *gmesh);

void AddSkeletonElements(TPZGeoMesh *gmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;
    bool hasexact = true;

    int numrefskeleton = 5;
    int maxporder = 7;
    int counter = 1;
    bool hrefinement = true;
    int numthreads = 1;
#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::ESquareRoot;
    LaplaceExactLower = LaplaceExact;
    LaplaceExactUpper = LaplaceExact;
    LaplaceExactLower.fExact = TLaplaceExample1::ESquareRootLower;
    LaplaceExactUpper.fExact = TLaplaceExample1::ESquareRootUpper;
#endif
    for ( int POrder = 3; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 2; irefskeleton < numrefskeleton; irefskeleton++)
        {
            bool elastic = !scalarproblem;
            TPZCompMesh *SBFem;
            if (hrefinement)
            {    
                SBFem = CreateCMesh(irefskeleton, POrder);
            }
            else
            {
                SBFem = SetupCrackedOneSquareElement(irefskeleton, POrder, hasexact, elastic);
            }
#ifdef _AUTODIFF
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(0);
                mat->SetForcingFunction(LaplaceExact.TensorFunction());
            }
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(0);
                mat->SetForcingFunction(LaplaceExact.TensorFunction());
            }
            // Dirichlet at the singularity
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(0);
                mat->SetForcingFunction(LaplaceExact.TensorFunction());
            }
#endif
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                std::ofstream gout("gmesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(SBFem->Reference(), gout,true);
                std::stringstream sout;
                SBFem->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            auto neq = SBFem->NEquations();
            std::cout << "neq = " << neq << std::endl;
            SolveSist(*Analysis, SBFem, numthreads);
            
            std::cout << "Post processing\n";
#ifdef _AUTODIFF
            Analysis->SetExact(Laplace_exact);
#endif
            
            TPZManVector<REAL> errors(3,0.);
            std::stringstream filename;
            filename << "SquareRootOneElement_NR_" << irefskeleton << "_P_" << POrder << ".vtk";
            TPZStack<std::string> vecnames,scalnames;
            // scalar
            scalnames.Push("State");
            vecnames.Push("Flux");
            Analysis->DefineGraphMesh(2, scalnames, vecnames, filename.str());
            Analysis->PostProcess(4);
            

            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            if(hasexact)
            {
            
                std::cout << "Compute errors\n";
                
                Analysis->PostProcessError(errors);
                std::stringstream sout;
                sout << "../CrackRestrainedShapeScalar.txt";
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(*  numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "Errmat[[" << irefskeleton+1 << "," << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            }
            
            delete Analysis;
            delete SBFem;   
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}
TPZCompMesh *SetupCrackedOneSquareElement(int nrefskeleton, int porder, bool applyexact, bool elastic)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    REAL coor[][3] = {
        {-1, 0},
        {-1, 1},
        {0, 1},
        {1, 1},
        {1, 0},
        {0, 0}};
    gmesh->NodeVec().Resize(6);
    for (int i = 0; i < 6; i++)
    {
        TPZManVector<REAL, 3> co(3, 0);
        co[0] = coor[i][0];
        co[1] = coor[i][1];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    {
        TPZManVector<int64_t, 4> nodeindices(4);
        int64_t index;
        nodeindices[0] = 0;
        nodeindices[1] = 1;
        nodeindices[2] = 5;
        nodeindices[3] = 5;
        gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
        nodeindices[0] = 1;
        nodeindices[1] = 2;
        nodeindices[2] = 5;
        nodeindices[3] = 5;
        gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
        nodeindices[0] = 2;
        nodeindices[1] = 3;
        nodeindices[2] = 5;
        nodeindices[3] = 5;
        gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
        nodeindices[0] = 3;
        nodeindices[1] = 4;
        nodeindices[2] = 5;
        nodeindices[3] = 5;
        gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
        nodeindices[0] = 0;
        nodeindices[1] = 5;
        gmesh->CreateGeoElement(EOned, nodeindices, EGroup+1, index);

        nodeindices[0] = 0;
        nodeindices[1] = 1;
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc1, index);
        nodeindices[0] = 1;
        nodeindices[1] = 2;
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc3, index);
        nodeindices[0] = 2;
        nodeindices[1] = 3;
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc3, index);
        nodeindices[0] = 3;
        nodeindices[1] = 4;
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc1, index);
        nodeindices[0] = 0;
        gmesh->CreateGeoElement(EPoint, nodeindices, Ebc2, index);
    }
    gmesh->BuildConnectivity();
    
    std::map<int, int> matidtranslation;
    matidtranslation[EGroup] = Emat1;
    matidtranslation[EGroup+1] = Emat2;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);

    TPZStack<int64_t> elids;
    elids.Push(0);
    elids.Push(1);
    elids.Push(2);
    elids.Push(3);
    elids.Push(4);
    build.AddPartition(elids,5);

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    InsertMaterialObjects(cmesh, !elastic, false);

    TPZMaterial *mat = cmesh->FindMaterial(Emat1);

    TPZMaterial *mat2 = mat->NewMaterial();
    mat2->SetId(Emat2);
    TPZMatPoisson<STATE> *mat2lapl = dynamic_cast<TPZMatPoisson<STATE> *>(mat2);
    mat2lapl->SetScaleFactor(1.e12);
    cmesh->InsertMaterialObject(mat2);

    std::set<int> volmatids,boundmatids;
    volmatids.insert(Emat1);
    volmatids.insert(Emat2);
    boundmatids.insert(Ebc1);
    boundmatids.insert(Ebc2);
    boundmatids.insert(Ebc3);
    boundmatids.insert(ESkeleton);
    build.DivideSkeleton(nrefskeleton);

    build.BuildComputationMesh(*cmesh,volmatids,boundmatids);
    
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        cmesh->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return cmesh;
}

TPZCompMesh * CreateCMesh(int nelx, int porder)
{
    TPZManVector<REAL,4> x0(3,0.),x1(3,0.);
    x0[0] = -1, x0[1] = 0;
    x1[0] = 1, x1[1] = 1;
    TPZManVector<int,4> nx(2,nelx);
    nx[0] = nelx*2;

    TPZGenGrid2D gengrid(nx,x0,x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    
    TPZGeoMesh * gmesh = new TPZGeoMesh();
    
    gengrid.Read(gmesh, Emat1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
    TPZManVector<REAL,3> start(3,0.), end(3,0.);
    start[0] = -1.;
    start[1] = 0.;
    end[0] = -0.5;
    end[1] = 0.;
    gengrid.SetBC(gmesh, start, end, Ebc1);
    start[0] = 0.5;
    start[1] = 0.;
    end[0] = 1;
    end[1] = 0.;
    gengrid.SetBC(gmesh, start, end, Ebc1);

    int64_t index = nelx-1;
    TPZGeoEl *gel1 = gmesh->Element(index);
    TPZGeoElBC(gel1,7,ESkeleton);
    TPZGeoElBC(gel1,6,ESkeleton);
    gel1->RemoveConnectivities();
    delete gel1;

    TPZGeoEl *gel2 = gmesh->Element(index+1);
    TPZGeoElBC(gel2,5,ESkeleton);
    TPZGeoElBC(gel2,6,ESkeleton);
    gel2->RemoveConnectivities();
    delete gel2;

    TPZManVector<int64_t,10> scalingcenters(1);
    scalingcenters[0] = nelx;
    // TPZManVector<int64_t, 4> nodeindices(4);
    // nodeindices[0] = nelx*2+nelx;
    // nodeindices[1] = nelx-1;
    // nodeindices[2] = nelx;
    // nodeindices[3] = nelx;
    // gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
    // nodeindices[0] = nelx*2+nelx+1;
    // nodeindices[1] = nelx*2+nelx;
    // nodeindices[2] = nelx;
    // nodeindices[3] = nelx;
    // gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
    // nodeindices[0] = nelx*2+nelx+2;
    // nodeindices[1] = nelx*2+nelx+1;
    // nodeindices[2] = nelx;
    // nodeindices[3] = nelx;
    // gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
    // nodeindices[0] = nelx+1;
    // nodeindices[1] = nelx*2+nelx+2;
    // nodeindices[2] = nelx;
    // nodeindices[3] = nelx;
    // gmesh->CreateGeoElement(EQuadrilateral, nodeindices, EGroup, index);
    // nodeindices[0] = nelx-1;
    // nodeindices[1] = nelx;
    // gmesh->CreateGeoElement(EOned, nodeindices, EGroup+1, index);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    InsertMaterialObjects(cmesh, true, true);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    {
        std::ofstream outc("CMesh.txt");
        cmesh->Print(outc);
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
    }

    std::map<int, int> matidtranslation;
    matidtranslation[EGroup] = Emat1;
    matidtranslation[EGroup+1] = Emat2;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);

    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t,10> elementgroup(nel,-1);
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel && (gel->MaterialId() == ESkeleton) ) {
            elementgroup[el] = 0;
        }
    }
    build.SetPartitions(elementgroup, scalingcenters);

    TPZMaterial *mat = cmesh->FindMaterial(Emat1);

    // TPZMaterial *mat2 = mat->NewMaterial();
    // mat2->SetId(Emat2);
    // TPZMatLaplacian *mat2lapl = dynamic_cast<TPZMatLaplacian *>(mat2);
    // mat2lapl->SetDimension(1);
    // mat2lapl->SetParameters(1.e9, 0);
    // cmesh->InsertMaterialObject(mat2);

    // std::set<int> volmatids,boundmatids;
    // volmatids.insert(Emat1);
    // volmatids.insert(Emat2);
    // boundmatids.insert(Ebc1);
    // boundmatids.insert(Ebc2);
    // boundmatids.insert(Ebc4);
    // boundmatids.insert(ESkeleton);
    build.DivideSkeleton(0);

    // build.BuildComputationMesh(*cmesh,volmatids,boundmatids);
    build.BuildComputationalMeshFromSkeleton(*cmesh);
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        cmesh->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(cmesh->Reference(), out, true);
        std::ofstream outcmesh("CMesh.vtk");
        TPZVTKGeoMesh vtkcmesh;
        vtkcmesh.PrintCMeshVTK(cmesh, outcmesh, true);
    }
    return cmesh;
}
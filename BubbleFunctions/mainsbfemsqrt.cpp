#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZBndCond.h"
#include "TPZVTKGeoMesh.h"
#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZGenGrid2D.h"
#include "pzgeoelbc.h"
#include "pzgnode.h"

TPZCompMesh *SetupCrackedOneSquareElement(int nrefskeleton, int porder, bool applyexact, bool elastic);

TPZCompMesh * CreateCMesh(int nelx, int porder);

TPZCompMesh *BuildSBFem(TPZAutoPointer<TPZGeoMesh> gmesh, int nx, int porder);

void AddBoundaryElements(TPZGeoMesh *gmesh);

void AddSkeletonElements(TPZGeoMesh *gmesh);

int main(int argc, char *argv[])
{
    bool scalarproblem = true;
    bool hasexact = true;

    int numrefskeleton = 5;
    int maxporder = 7;
    int counter = 1;
    bool hrefinement = true;
    
    LaplaceExact.fExact = TLaplaceExample1::ESquareRoot;

    int numthreads = 1;
    for ( int POrder = 6; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 1; irefskeleton < numrefskeleton; irefskeleton++)
        {
            bool elastic = !scalarproblem;
            TPZCompMesh *SBFem;
            if (hrefinement)
            {    
                SBFem = CreateCMesh(irefskeleton, POrder);
            }
            else
            {
                SBFem = SetupCrackedOneSquareElement(irefskeleton-1, POrder, hasexact, elastic);
            }
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            auto neq = SBFem->NEquations();
            std::cout << "neq = " << neq << std::endl;
            SolveSist(Analysis, SBFem);
            
            std::cout << "Post processing\n";

            // scalar
            if (0)
            {
                std::stringstream filename;
                filename << "SquareRootOneElement_NR_" << irefskeleton << "_P_" << POrder << ".vtk";
                TPZStack<std::string> vecnames,scalnames;
                scalnames.Push("Solution");
                vecnames.Push("Derivative");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, filename.str());
                Analysis->PostProcess(4);
            }
            
            {
                std::cout << "Compute errors\n";
                
                Analysis->SetExact(LaplaceExact.ExactSolution());
            
                TPZManVector<REAL> errors(3,0.);

                Analysis->SetThreadsForError(numthreads);
                Analysis->PostProcessError(errors, false);
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

    {
        int nstate = 1;
        TPZFMatrix<STATE> val1(nstate, nstate, 0.);
        const TPZManVector<double> val2(nstate, 0.);

        auto matloc = new TPZDarcyFlow(Emat1, 2);
        matloc->SetConstantPermeability(1.);
        cmesh->InsertMaterialObject(matloc);

        auto matloc2 = new TPZDarcyFlow(Emat2, 1);
        matloc2->SetConstantPermeability(1.e12);
        cmesh->InsertMaterialObject(matloc2);
        
        auto BCond1 = matloc->CreateBC(matloc, Ebc1, 0, val1, val2);
        auto BCond2 = matloc->CreateBC(matloc, Ebc2, 0, val1, val2);
        auto BCond3 = matloc->CreateBC(matloc, Ebc3, 0, val1, val2);
        
        BCond1->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        BCond2->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        BCond3->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
    }

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

TPZCompMesh * CreateCMesh(int nc, int porder)
{
    TPZManVector<REAL,4> x0(3,0.),x1(3,0.);
    x0[0] = -1, x0[1] = 0;
    x1[0] = 1, x1[1] = 1;
    int nelx = pow(2,nc);
    TPZManVector<int,4> nx(2,nelx);
    nx[0] = nelx*2;

    TPZGenGrid2D gengrid(nx,x0,x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh();
    
    gengrid.Read(gmesh, Emat3);
    gengrid.SetBC(gmesh, 5, Ebc3);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc3);
    TPZManVector<REAL,3> start(3,0.), end(3,0.);
    start[0] = -1.;
    start[1] = 0.;
    end[0] = -0.5;
    end[1] = 0.;
    gengrid.SetBC(gmesh, start, end, Ebc3);
    start[0] = 0.5;
    start[1] = 0.;
    end[0] = 1;
    end[1] = 0.;
    gengrid.SetBC(gmesh, start, end, Ebc3);
    gmesh->BuildConnectivity();

    int neltoremovex = nelx;

    int nelcollapsed = nelx + 2*nc;
    TPZManVector<int64_t, 4> nodescollapsed(4);
    TPZStack<int64_t> elids;

    int64_t index;
    for (int j = 0; j < nelx/2; j++)
    {
        for (int i = 0; i < neltoremovex; i++)
        {
            int index2 = nelx/2 + j*(nelx*2) + i;
            TPZGeoEl *gel1 = gmesh->Element(index2);
            if (j == nelx/2-1)
            {
                nodescollapsed[0] = gel1->SideNodeIndex(6,0);
                nodescollapsed[1] = gel1->SideNodeIndex(6,1);
                nodescollapsed[2] = nelx;
                nodescollapsed[3] = nelx;
                gmesh->CreateGeoElement(EQuadrilateral, nodescollapsed, EGroup, index);
                elids.Push(index);
                gmesh->CreateGeoElement(EOned, nodescollapsed, Ebc1, index);
            }
            if (i == 0)
            {
                nodescollapsed[0] = gel1->SideNodeIndex(7,0);
                nodescollapsed[1] = gel1->SideNodeIndex(7,1);
                nodescollapsed[2] = nelx;
                nodescollapsed[3] = nelx;
                gmesh->CreateGeoElement(EQuadrilateral, nodescollapsed, EGroup, index);
                elids.Push(index);
                gmesh->CreateGeoElement(EOned, nodescollapsed, Ebc1, index);
            }
            if (i == neltoremovex-1)
            {
                nodescollapsed[0] = gel1->SideNodeIndex(5,0);
                nodescollapsed[1] = gel1->SideNodeIndex(5,1);
                nodescollapsed[2] = nelx;
                nodescollapsed[3] = nelx;
                gmesh->CreateGeoElement(EQuadrilateral, nodescollapsed, EGroup, index);
                elids.Push(index);
                gmesh->CreateGeoElement(EOned, nodescollapsed, Ebc1, index);
            }
            gel1->RemoveConnectivities();
            delete gel1;
        }
    }

    TPZManVector<int64_t, 2> nodeindices(2);
    nodeindices[0] = nelx/2;
    nodeindices[1] = nelx;
    gmesh->CreateGeoElement(EOned, nodeindices, EGroup+1, index);
    elids.Push(index);

    nodeindices.Resize(1);
    nodeindices[0] = nelx/2;
    gmesh->CreateGeoElement(EPoint, nodeindices, Ebc2, index);

    gmesh->BuildConnectivity();
    {
        TPZVTKGeoMesh vtk;
        std::ofstream out("Geometry.vtk");
        vtk.PrintGMeshVTK(gmesh, out, true);
    }
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    {
        int nstate = 1;
        TPZFMatrix<STATE> val1(nstate, nstate, 0.);
        const TPZManVector<double> val2(nstate, 0.);

        auto matloc = new TPZDarcyFlow(Emat1, 2);
        matloc->SetConstantPermeability(1.);
        cmesh->InsertMaterialObject(matloc);
        
        auto BCond1 = matloc->CreateBC(matloc, Ebc1, 0, val1, val2);
        BCond1->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        cmesh->InsertMaterialObject(BCond1);

        auto matloc2 = new TPZDarcyFlow(Emat2, 1);
        matloc2->SetConstantPermeability(1.e12);
        cmesh->InsertMaterialObject(matloc2);
        
        auto BCond2 = matloc2->CreateBC(matloc2, Ebc2, 0, val1, val2);
        BCond2->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        cmesh->InsertMaterialObject(BCond2);

        auto matloc3 = new TPZDarcyFlow(Emat3, 2);
        matloc3->SetConstantPermeability(1.);
        cmesh->InsertMaterialObject(matloc3);
        
        auto BCond3 = matloc3->CreateBC(matloc3, Ebc3, 0, val1, val2);
        BCond3->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        cmesh->InsertMaterialObject(BCond3);
    }

    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuous();

    std::map<int, int> matidtranslation;
    matidtranslation[EGroup] = Emat1;
    matidtranslation[EGroup+1] = Emat2;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);

    build.AddPartition(elids,nelx);

    std::set<int> volmatids,boundmatids;
    volmatids.insert(Emat1);
    volmatids.insert(Emat2);
    boundmatids.insert(Ebc1);
    boundmatids.insert(Ebc2);
    build.DivideSkeleton(0);
    build.BuildComputationMesh(*cmesh,volmatids,boundmatids);

    // std::set<int> matfem = {Emat3,Ebc3};
    // {
    // }
    // cmesh->AutoBuild(matfem);

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
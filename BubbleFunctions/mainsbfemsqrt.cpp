#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "pzgeoelbc.h"
#include "TPZGenGrid2D.h"
#include "pzcheckgeom.h"
#include "TPZBndCondT.h"
#include "TPZBuildSBFem.h"

#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

TPZAutoPointer<TPZGeoMesh> CreateGMesh(int nelx);

TPZCompMesh *BuildSBFem(TPZAutoPointer<TPZGeoMesh> & gmesh, int nx, int porder);

void IntegrateDirect(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    ElastExact.fProblemType = TElasticity2DAnalytic::ESquareRoot;
    ElastExact.gE = 10;
    ElastExact.gPoisson = 0.3;
    ElastExact.fPlaneStress = 1;
    ElastExactLower = ElastExact;
    ElastExactUpper = ElastExact;
    ElastExactLower.fProblemType = TElasticity2DAnalytic::ESquareRootLower;
    ElastExactUpper.fProblemType = TElasticity2DAnalytic::ESquareRootUpper;

    int maxnelxcount = 4;
    int maxporder = 5;
    int counter = 1;
    int numthreads = 4;
    int nx = 4;

    TPZAutoPointer<TPZGeoMesh> gmesh = CreateGMesh(nx);
    if(0)
    {
        std::cout << "Plotting the geometric mesh\n";
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }

    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for(int nelxcount = 0; nelxcount < maxnelxcount; nelxcount += 1)
        {
            TPZAutoPointer<TPZGeoMesh> locgmesh = new TPZGeoMesh(gmesh);
            {
                TPZCheckGeom check(locgmesh.operator->());
                check.UniformRefine(nelxcount);
            }
            TPZCompMesh *SBFem = BuildSBFem(locgmesh, nx, POrder);

            if(1)
            {
                std::cout << "Plotting the geometric mesh\n";
                std::stringstream sout;
                sout << "SBFem_Fem_Geometry." << counter << ".vtk";
                std::ofstream out(sout.str());
                TPZVTKGeoMesh vtk;
                vtk.PrintGMeshVTK(SBFem->Reference(), out,true);
            }

            std::cout << "nelx = " << nelxcount << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis Analysis(SBFem,mustOptimizeBandwidth);
            Analysis.SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(&Analysis, SBFem);
            
            std::cout << "Post processing\n";
            
            Analysis.SetExact(ElastExact.ExactSolution());
            
            TPZManVector<REAL> errors(3,0.);
            
            int64_t neq = SBFem->Solution().Rows();
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                vecnames.Push("Displacement");
                vecnames.Push("Strain");
                scalnames.Push("SigmaX");
                scalnames.Push("SigmaY");
                scalnames.Push("TauXY");
                Analysis.DefineGraphMesh(2, scalnames, vecnames, "../EmbeddedSBFemElasticity2DSolution.vtk");
                Analysis.PostProcess(3);
            }

            std::cout << "Compute errors\n";
            
            Analysis.PostProcessError(errors,false);
            Analysis.SetThreadsForError(numthreads);
            
            std::stringstream sout;
            sout << "../EmbeddedSBFem";
            sout << "Elastic2D.txt";
            
            std::ofstream results(sout.str(),std::ios::app);
            results.precision(15);
            results << "(* nx " << nx  << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
            TPZFMatrix<double> errmat(1,6);
            for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
            errmat(0,3) = 1./(10*pow(2.,nelxcount));
            errmat(0,4) = neq;
            errmat(0,5) = POrder;
            std::stringstream varname;
            varname << "Errmat[[" << nelxcount+1 << "," << POrder << "]] = (1/1000000)*";
            errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            
            delete SBFem;
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

TPZAutoPointer<TPZGeoMesh> CreateGMesh(int nelx)
{
    // only even number of elements are allowed
    if(nelx%2 != 0)
    {
        DebugStop();
    }
    TPZManVector<REAL,4> x0(3,-1.),x1(3,1.);
    x0[0] = -1;
    x0[1] = -1;
    x1[0] = 1;
    x1[1] = 1;
    x0[2] = 0.;
    x1[2] = 0.;
    TPZManVector<int,4> nx(2,nelx);
    TPZGenGrid2D gengrid(nx,x0,x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    
    gengrid.Read(gmesh,Emat1);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc1);
    gengrid.SetBC(gmesh, 6, Ebc1);
    TPZManVector<REAL,3> start(3,0.), end(3,0.);
    start[0] = -1.;
    start[1] = 1.;
    end[0] = -1.;
    gengrid.SetBC(gmesh, start, end, Ebc2);
    start = end;
    end[1] = -1.;
    gengrid.SetBC(gmesh, start, end, Ebc3);
    int64_t nnodes = gmesh->NNodes();
    gmesh->NodeVec().Resize(nnodes+nelx/2);
    REAL delx = 2./nelx;
    for (int in=0; in < nelx/2; in++) {
        TPZManVector<REAL,3> xco(3,0.);
        xco[0] = -1.+in*delx;
        gmesh->NodeVec()[nnodes+in].Initialize(xco, gmesh);
    }
    int minel = nelx*nelx/2;
    int maxel = nelx*(nelx+1)/2;
    for (int64_t el = minel; el < maxel; el++) {
        int64_t firstnode = el-nelx*nelx/2+nnodes;
        TPZGeoEl *gel = gmesh->Element(el);
        gel->SetNodeIndex(0, firstnode);
        if(firstnode+1 < nnodes+nelx/2)
        {
            gel->SetNodeIndex(1, firstnode+1);
        }
        if(el == nelx*nelx/2)
        {
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neighbour(gelside.Neighbour());
            while (neighbour != gelside) {
                if (neighbour.Element()->MaterialId() == Ebc2) {
                    int sidenodelocindex = neighbour.SideNodeLocIndex(0);
                    neighbour.Element()->SetNodeIndex(sidenodelocindex, firstnode);
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    int64_t index = (nelx-1)*(nelx/2)-1;
    TPZGeoEl *gel1 = gmesh->Element(index);
    TPZGeoElBC(gel1,7,ESkeleton);
    TPZGeoElBC(gel1,4,ESkeleton);
    TPZGeoEl *gel2 = gmesh->Element(index+1);
    TPZGeoElBC(gel2,4,ESkeleton);
    TPZGeoElBC(gel2,5,ESkeleton);
    TPZGeoEl *gel3 = gmesh->Element(index+nelx);
    TPZGeoElBC(gel3,6,ESkeleton);
    TPZGeoElBC(gel3,7,ESkeleton);
    TPZGeoEl *gel4 = gmesh->Element(index+nelx+1);
    TPZGeoElBC(gel4,5,ESkeleton);
    TPZGeoElBC(gel4,6,ESkeleton);
    gel1->RemoveConnectivities();
    delete gel1;
    gel2->RemoveConnectivities();
    delete gel2;
    gel3->RemoveConnectivities();
    delete gel3;
    gel4->RemoveConnectivities();
    delete gel4;
    return gmesh;
}

TPZCompMesh *BuildSBFem(TPZAutoPointer<TPZGeoMesh> & gmesh, int nx, int porder)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    InsertMaterialObjects(cmesh, false, true);
    {
        TPZMaterial *mat = cmesh->FindMaterial(Ebc1);
        auto bndcond = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        bndcond->SetType(0);
        bndcond->SetForcingFunctionBC(ElastExact.ExactSolution());
    }
    {
        TPZMaterial *mat = cmesh->FindMaterial(Ebc2);
        auto bndcond = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        bndcond->SetType(0);
        bndcond->SetForcingFunctionBC(ElastExactUpper.ExactSolution());
    }
    {
        TPZMaterial *mat = cmesh->FindMaterial(Ebc3);
        auto bndcond = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        bndcond->SetType(0);
        bndcond->SetForcingFunctionBC(ElastExactLower.ExactSolution());
    }
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    std::map<int,int> matidtranslation;
    matidtranslation[ESkeleton] = Emat1;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
    
    TPZManVector<int64_t,10> scalingcenters(1);
    scalingcenters[0] = ((nx+1)*(nx+1)-1)/2;

    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t,10> elementgroup(nel,-1);
    for (int64_t el=0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel && gel->MaterialId() == ESkeleton) {
            elementgroup[el] = 0;
        }
    }

    build.SetPartitions(elementgroup, scalingcenters);
    build.BuildComputationalMeshFromSkeleton(*cmesh);

    return cmesh;
}

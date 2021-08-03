#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common3D.h"
#include "TPZBndCond.h"
#include "TPZBuildSBFem.h"
#include "TPZBuildSBFem.h"
#include "TPZVTKGeoMesh.h"

#include "Poisson/TPZMatPoisson.h"
#include "TPZBndCond.h"
#include "pzgeoelbc.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExactLower;
TElasticity2DAnalytic ElastExactUpper;
#endif

TPZCompMesh * Hexagon(int porder);
void InsertMaterialObjectsDFN(TPZCompMesh *cmesh);
void SolveSistHexagon(TPZAnalysis *an, TPZCompMesh *Cmesh, int numthreads);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = false;
    bool hasexact = true;

    int numrefskeleton = 1;
    int maxporder = 3;
    int counter = 1;

    int numthreads = 0;

#ifdef _AUTODIFF
    // ExactLaplace.fExact = TLaplaceExample1::EConst;
    // ExactElast.fProblemType = TElasticity2DAnalytic::ELoadedBeam;
#endif
    for ( int POrder = 2; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            bool elastic = !scalarproblem;

            std::string filename("TwoOrthogonalCracks.txt");
            std::string vtkfilename;
            std::string vtkfilegeom;
            std::string rootname;
            std::string boundaryname;

            {
                int pos = filename.find(".txt");
                std::string truncate = filename;
                truncate.erase(pos);
                rootname = truncate;
                std::stringstream sout;
                sout << truncate << "_t" << numthreads << "_p" << POrder << "_href" << irefskeleton << ".vtk";
                vtkfilename = sout.str();
                std::stringstream boundstr;
                boundstr << truncate << "_boundary";
                boundaryname = boundstr.str();
                std::stringstream vtkgeom;
                vtkgeom << truncate << "_geom.vtk";
                vtkfilegeom = vtkgeom.str();
            }

            TPZManVector<int64_t,1000> elpartitions;
            TPZVec<int64_t> scalingcenterindices;
            std::cout << "Reading " << filename << std::endl;
            TPZAutoPointer<TPZGeoMesh> gmesh =ReadUNSWSBGeoFile(filename, elpartitions, scalingcenterindices);
            gmesh->SetDimension(2);
            {
                std::cout << "Plotting the geometric mesh\n";
//                std::ofstream outg("GMesh3D.txt");
//                gmesh->Print(outg);
                std::ofstream out("Geometry3D.vtk");
                TPZVTKGeoMesh vtk;
                vtk.PrintGMeshVTK(gmesh, out,true);
            }
            int nel = gmesh->NElements();
            for (int i = 0; i < nel; i++)
            {
                TPZGeoEl * gel = gmesh->Element(i);
                TPZGeoElBC(gel,2,Ebc1);
            }
            
            elpartitions.Resize(gmesh->NElements(), -1);

            std::cout << "Building the computational mesh\n";
            std::map<int,int> matidtranslation;
            matidtranslation[ESkeleton] = Emat1;
            matidtranslation[Ebc1] = Ebc1;
            TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
            build.SetPartitions(elpartitions, scalingcenterindices);
            build.DivideSkeleton(irefskeleton);

            TPZCompMesh *fem = new TPZCompMesh(gmesh);
            fem->SetDefaultOrder(POrder);
            InsertMaterialObjectsDFN(fem);
            build.BuildComputationalMeshFromSkeleton(*fem);

            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(fem, mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << fem->NEquations() << std::endl;
            SolveSistHexagon(Analysis, fem, numthreads);
            
            int64_t neq = fem->Solution().Rows();
            
            std::cout << "Plotting shape functions\n";
            TPZFMatrix<STATE> sol0 = fem->Solution();
            for (int i=1; i<sol0.Rows() ;i++)
            {        
                TPZFMatrix<STATE> sol = fem->Solution();
                sol.Zero();
                sol(i,0) = 1;
                
                fem->LoadSolution(sol);
                Analysis->LoadSolution(sol);
                
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("State");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, "../ShapeFunctions.vtk");
                Analysis->PostProcess(3);
            }
            
            delete Analysis;
            delete fem;
            //                exit(-1);
        }
        //            exit(-1);
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

TPZCompMesh * Hexagon(int POrder)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    REAL coor[][3] = {
        {0,0},
        {1.,0},
        {1.5,0.866025},
        {1,1.73205},
        {0,1.73205},
        {-0.5,0.866025},
        {0.5,0.866025}
    };
    gmesh->NodeVec().Resize(7);
    for (int i=0; i<7; i++) {
        TPZManVector<REAL,3> co(3,0);
        co[0] = coor[i][0];
        co[1] = coor[i][1];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    {
        int64_t index;
        TPZManVector<int64_t,2> nodeindices(4);
        for (int i=0; i<5; i++) {
            nodeindices[0] = i;
            nodeindices[1] = i+1;
            nodeindices[2] = 6;
            nodeindices[3] = 6;
            gmesh->CreateGeoElement(EQuadrilateral, nodeindices, ESkeleton, index);
        }
        nodeindices[0] = 5;
        nodeindices[1] = 0;
        nodeindices[2] = 6;
        nodeindices[3] = 6;
        gmesh->CreateGeoElement(EQuadrilateral, nodeindices, ESkeleton, index);
    }
    gmesh->BuildConnectivity();
    int64_t nel = gmesh->NElements();
    for (int i=0; i<nel; i++) {
        gmesh->Element(i)->CreateBCGeoEl(4, Ebc1);
    }
    gmesh->BuildConnectivity();
            
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    // InsertMaterialObjects(cmesh, true, true);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    std::ofstream out("CMesh.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintCMeshVTK(cmesh, out,true);

    return cmesh;

    
}



/// insert material objects in the computational mesh
void InsertMaterialObjectsDFN(TPZCompMesh *cmesh)
{
    
    // Getting mesh dimension
    int matId1 = Emat1;
    int nstate = 1;
    
    TPZMatPoisson<STATE> *matloc = new TPZMatPoisson<STATE>(matId1,2);
    nstate = 1;
    cmesh->InsertMaterialObject(matloc);
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.);
    TPZManVector<STATE> val2(nstate,0.);
    {
        val2[0] = 1.;
        auto BCond1 = matloc->CreateBC(matloc, Ebc1, 0, val1, val2);
        cmesh->InsertMaterialObject(BCond1);
    }
    
    {
        val1.Zero();
        val2[0] = 0.;
        auto BCond2 = matloc->CreateBC(matloc, Ebc2, 0, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    {
        val1.Zero();
        val2[0] = 0.;
        auto BCond2 = matloc->CreateBC(matloc, Ebc3, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    {
        TPZMaterial *BCond2 = matloc->NewMaterial();
        TPZMatPoisson<STATE> *matlap = dynamic_cast<TPZMatPoisson<STATE> *>(BCond2);
        matlap->SetScaleFactor(0.04e5);
        matlap->SetDimension(1);
        matlap->SetId(Ebc4);
        cmesh->InsertMaterialObject(BCond2);
    }
    
    
    val2[0] = 0.0;
    auto BSkeleton = matloc->CreateBC(matloc, ESkeleton, 1, val1, val2);
    cmesh->InsertMaterialObject(BSkeleton);

    
    
}
void SolveSistHexagon(TPZAnalysis *an, TPZCompMesh *Cmesh, int numthreads)
{
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix<REAL> strmat(Cmesh);
#endif

    strmat.SetNumThreads(numthreads);
    an->SetStructuralMatrix(strmat);
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    an->Assemble();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif
    
    
}

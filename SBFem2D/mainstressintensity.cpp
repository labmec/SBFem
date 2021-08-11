#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZBndCondT.h"
#include "TPZVTKGeoMesh.h"
#include "TPZSBFemElementGroup.h"
#include "Elasticity/TPZElasticity2D.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

TElasticity2DAnalytic ElastExactLower;
TElasticity2DAnalytic ElastExactUpper;
TLaplaceExample1 LaplaceExactLower;
TLaplaceExample1 LaplaceExactUpper;

void IntegrateDirect(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = false; // always elastic
    bool hasexact = true;

    int numrefskeleton = 4;
    int maxporder = 4;
    int counter = 1;

    int numthreads = 1;
    ElastExact.fProblemType = TElasticity2DAnalytic::ESquareRoot;
    ElastExact.gE = 10;
    ElastExact.gPoisson = 0.3;
    ElastExact.fPlaneStress = 1;
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            bool elastic = !scalarproblem;
            TPZCompMesh *SBFem = SetupCrackedOneElement(irefskeleton, POrder, hasexact, elastic);
            std::ofstream out("Geometry.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(SBFem->Reference(), out, true);
            

            
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
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(*Analysis, SBFem, numthreads);
            
            std::cout << "Post processing\n";
            Analysis->SetExact(Laplace_exact);
            if (elastic)
            {
                Analysis->SetExact(Elasticity_exact);
            }
            
            TPZManVector<REAL> errors(3,0.);
            
            int64_t neq = SBFem->Solution().Rows();
            
            if(!scalarproblem)
            {
                std::stringstream filename;
                filename << "SquareRootOneElement_NR_" << irefskeleton << "_P_" << POrder << ".vtk";
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                vecnames.Push("Displacement");
                scalnames.Push("SigmaX");
                scalnames.Push("SigmaY");
                scalnames.Push("TauXY");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, filename.str());
                Analysis->PostProcess(3);
            }
            

            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            if(hasexact)
            {
            
                std::cout << "Compute errors\n";
                
                Analysis->PostProcessError(errors, false);
                std::stringstream sout;
                sout << "../CrackRestrainedShape";
                if (scalarproblem) {
                    sout << "Scalar.txt";
                }
                else
                    sout << "Elastic2D.txt";
                
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
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZBndCondT.h"
#include "TPZVTKGeoMesh.h"
#include "TPZSBFemElementGroup.h"
#include "Elasticity/TPZElasticity2D.h"

int main(int argc, char *argv[])
{
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
    ElastExactLower = ElastExact;
    ElastExactUpper = ElastExact;
    ElastExactLower.fProblemType = TElasticity2DAnalytic::ESquareRootLower;
    ElastExactUpper.fProblemType = TElasticity2DAnalytic::ESquareRootUpper;

    for ( int POrder = 1; POrder < maxporder; POrder ++)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            bool elastic = !scalarproblem;
            TPZCompMesh *SBFem = SetupCrackedOneElement(irefskeleton, POrder, hasexact, elastic);
            std::ofstream out("Geometry.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(SBFem->Reference(), out, true);
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem, mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(*Analysis, SBFem, numthreads);
            
            std::cout << "Post processing\n";
            
            int64_t neq = SBFem->Solution().Rows();
            
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
            
            {
                std::cout << "Compute errors\n";
                Analysis->SetExact(ElastExact.ExactSolution());
                // Analysis->SetThreadsForError(numthreads);

                TPZManVector<REAL> errors(3,0.);
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
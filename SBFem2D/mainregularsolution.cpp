#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

int main(int argc, char *argv[])
{
    // Initial data
    int minnelxcount = 1, maxnelxcount = 5;
    int minporder = 1, maxporder = 4;
    int numrefskeleton = 3;
    int numthreads = 32;
    bool scalarproblem = false; // false for elasticity 2D problems
    bool usesbfem = true; // false for FEM simulations
    if (usesbfem == false) 
    {
        numrefskeleton = 1;
    }
    bool useexact = true;
    
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
    ElastExact.fProblemType = TElasticity2DAnalytic::ELoadedBeam;
    ElastExact.gE = 10;
    ElastExact.gPoisson = 0.3;
    ElastExact.fPlaneStress = 1;

    int countstep = 1;
    for ( int POrder = minporder; POrder <= maxporder; POrder ++)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            for(int nelxcount = minnelxcount; nelxcount < maxnelxcount; nelxcount ++)
            {
                int nelx = 1 << (nelxcount-1);
                
                TPZCompMesh *SBFem;
                if(usesbfem)
                {
                    SBFem = SetupSquareMesh(nelx,irefskeleton,POrder, scalarproblem, useexact);
                }
                else
                {
                    SBFem = SetupSquareH1Mesh(nelx, POrder, scalarproblem, useexact);
                }
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                std::cout << "Entering Analysis \n";
                bool mustOptimizeBandwidth = true;
                TPZLinearAnalysis Analysis(SBFem,mustOptimizeBandwidth);
                Analysis.SetStep(countstep++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(Analysis, SBFem, numthreads);
                
                bool printcmeshwsol = false;
                if(printcmeshwsol)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }

                std::cout << "Post processing\n";
                std::string sout;
                if (scalarproblem)
                {
                    sout.append("../RegularSolution");
                } else
                {
                    sout.append("../RegularElasticity2DSolution");
                }
                if(usesbfem)
                {
                    sout.append("_SBFem");
                }
                else
                {
                    sout.append("_H1");
                }
                PostProcessing(Analysis, sout, scalarproblem, numthreads, POrder, nelxcount, irefskeleton);
                
                delete SBFem;
            }
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

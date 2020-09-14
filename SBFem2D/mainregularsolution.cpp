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
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif

#ifndef _AUTODIFF
    std::cout << "This program needs FAD to run \n";
    DebugStop();
#endif

    // Initial data
    int minnelxcount = 1, maxnelxcount = 2;
    int minporder = 5, maxporder = 6;
    int numrefskeleton = 5;
    int numthreads = 1;
    bool scalarproblem = true; // false for elasticity 2D problems
    bool usesbfem = true; // false for FEM simulations
    if (usesbfem == false) 
    {
        numrefskeleton = 1;
    }
    bool useexact = true;
#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
    ElastExact.fProblemType = TElasticity2DAnalytic::ELoadedBeam;
    ElastExact.gE = 10;
    ElastExact.gPoisson = 0.3;
    ElastExact.fPlaneStress = 1;
#endif

    int countstep = 1;
    for ( int POrder = minporder; POrder <= maxporder; POrder ++)
    {
        for (int irefskeleton = 4; irefskeleton < numrefskeleton; irefskeleton++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
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
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    SBFem->Print(sout);
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                std::cout << "Entering Analysis \n";
#ifdef USING_BOOST
                boost::posix_time::ptime t01 = boost::posix_time::microsec_clock::local_time();
#endif		
                bool mustOptimizeBandwidth = true;
                TPZAnalysis Analysis(SBFem,mustOptimizeBandwidth);
                Analysis.SetStep(countstep++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(Analysis, SBFem, numthreads);
#ifdef USING_BOOST
                boost::posix_time::ptime t02 = boost::posix_time::microsec_clock::local_time();
                std::cout << "Time for analysis " << t02-t01 << std::endl;
#endif
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

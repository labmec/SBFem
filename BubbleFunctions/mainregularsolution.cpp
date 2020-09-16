#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZSBFemElementGroup.h"
#include "pzcondensedcompel.h"

#include <chrono>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;

    int maxnelxcount = 6;
    int numrefskeleton = 0;
    int maxporder = 3;
    int counter = 1;
    bool usesbfem = true;
    if (usesbfem == false) {
        numrefskeleton = 1;
    }
#ifdef _AUTODIFF
    ElastExact.fProblemType = TElasticity2DAnalytic::ELoadedBeam;
    LaplaceExact.fExact = TLaplaceExample1::ECosCos;
#endif
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
            int irefskeleton = 0;
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;
            for(int nelxcount = 1; nelxcount < maxnelxcount; nelxcount += 1)
            {
                int nelx = 1 << (nelxcount-1);
                bool useexact = true;
                if(!scalarproblem)
                {
#ifdef _AUTODIFF
                    ElastExact.gE = 10;
                    ElastExact.gPoisson = 0.3;
                    ElastExact.fPlaneStress = 0;
#endif
                }
                
                TPZCompMesh *SBFem;
                if(usesbfem)
                {
                    SBFem = SetupSquareMesh(nelx,irefskeleton,POrder, scalarproblem,useexact);
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
                int64_t nel = SBFem->NElements();
                for (int64_t el = 0; el<nel; el++) {
                    TPZCompEl *cel = SBFem->Element(el);
                    if(!cel) continue;
                    TPZSBFemElementGroup *sbgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
                    if(!sbgr) continue;
                    TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel,true);
                }
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;

                std::cout << "Entering on Analysis \n";
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                Analysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                
                std::clock_t begin_analysis = clock();
                SolveSist(Analysis, SBFem);

                std::clock_t end_analysis = clock();
                double elapsed_time = double(end_analysis - begin_analysis)/CLOCKS_PER_SEC;
		        std::cout << "Time taken for solving: " << elapsed_time << std::endl;

                std::cout << "Post processing\n";
#ifdef _AUTODIFF
                if(scalarproblem)
                {
                    Analysis->SetExact(Laplace_exact);
                }
                else
                {
                    Analysis->SetExact(Elasticity_exact);
                }
#endif                
                int64_t neq = SBFem->Solution().Rows();
                
                if(0)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    scalnames.Push("State");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularSolution.vtk");
                    Analysis->PostProcess(3);
                }
                if(0)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    vecnames.Push("Displacement");
                    scalnames.Push("SigmaX");
                    scalnames.Push("SigmaY");
                    scalnames.Push("TauXY");
                    scalnames.Push("EpsX");
                    scalnames.Push("EpsY");
                    scalnames.Push("EpsXY");
                    std::stringstream sout;
                    sout << "../RegularElasticity2DSolution";
                    if(usesbfem)
                    {
                        sout << "_SBFem.vtk";
                    }
                    else
                    {
                        sout << "_H1.vtk";
                    }
                    Analysis->DefineGraphMesh(2, scalnames, vecnames,sout.str());
                    Analysis->PostProcess(3);
                }

                if(0)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    vecnames.Push("Scalar");
                    std::stringstream sout;
                    sout << "../RegularLaplace2DSolution";
                    if(usesbfem)
                    {
                        sout << "_SBFem.vtk";
                    }
                    else
                    {
                        sout << "_H1.vtk";
                    }
                    Analysis->DefineGraphMesh(2, scalnames, vecnames,sout.str());
                    Analysis->PostProcess(3);
                }
                
                std::cout << "Compute errors\n";
                
                TPZManVector<REAL,10> errors(3,0.);
                Analysis->SetThreadsForError(4);
                Analysis->PostProcessError(errors, false);
                
                std::stringstream sout;
                sout << "../RegularSolution";
                if (scalarproblem) {
                    sout << "Scalar";
                }
                else
                    sout << "Elastic2D";
                
                if(usesbfem)
                {
                    sout << "_SBFem.txt";
                }
                else
                {
                    sout << "_H1.txt";
                }
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(* nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "Errmat[[" << nelxcount << "]][[" << irefskeleton+1 << "]][[" << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                
                TPZFMatrix<REAL> sol = Analysis->Solution();
                
                bool plotshape = false;
                if(plotshape)
                {
                    TPZFMatrix<REAL> sol0 = sol;
                    for (int i=0; i<sol0.Rows() ;i++){
                        
                        TPZFNMatrix<3,REAL> sol = SBFem->Solution();
                        sol.Zero();
                        sol(i,0) = 1;
                        
                        SBFem->LoadSolution(sol);
                        Analysis->LoadSolution(sol);
                        
                        TPZStack<std::string> vecnames,scalnames;
                        // scalar
                        scalnames.Push("State");
                        Analysis->DefineGraphMesh(2, scalnames, vecnames, "../ShapeFunctions.vtk");
                        Analysis->PostProcess(3);
                    }
                }
                
                delete Analysis;
                delete SBFem;
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

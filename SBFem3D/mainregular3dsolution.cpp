#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common3D.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzaxestools.h"
#include "TPZMaterial.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void AnalyseSolution(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minnelx = 1, maxnelx = 5;
    int minrefskeleton = 0, maxrefskeleton = 1;
    int minporder = 1, maxporder = 2;
    int counter = 1;
	int numthreads = 8;
    bool scalarproblem = true;
    bool usesbfem = false;

    ExactLaplace.fExact = TLaplaceExample1::EHarmonic2;
    ExactElast.fProblemType = TElasticity3DAnalytic::ETestShearMoment;
    ExactElast.fE = 1.;
    ExactElast.fPoisson = 0.2;
    // TPZMaterial::gBigNumber = 1e18;
    
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            for(int nelxcount = minnelx; nelxcount < maxnelx; nelxcount++)
            {
                // TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;
                // TPZSBFemElementGroup::gPolynomialShapeFunctions = true;
                int nelx = 1 << (nelxcount-1);
                TPZCompMesh *SBFem = SetupSquareMesh3D(nelx, irefskeleton, POrder, scalarproblem, usesbfem);
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    SBFem->Print(sout);
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif 
                // int64_t nel = SBFem->NElements();
                // for (int64_t el = 0; el<nel; el++) {
                //     TPZCompEl *cel = SBFem->Element(el);
                //     if(!cel) continue;
                //     TPZSBFemElementGroup *sbgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
                //     if(!sbgr) continue;
                //     TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel,true);
                // }
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                Analysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(Analysis, SBFem, numthreads);
                
                std::cout << "Plotting\n";
                
                if (scalarproblem)
                {
                    Analysis->SetExact(Laplace_exact);
                } else
                {
                    Analysis->SetExact(Elasticity_exact);
                }
                
                int64_t neq = SBFem->Solution().Rows();
                
                std::string vtkfilename;
                if (!scalarproblem) {
                    vtkfilename = "../Elast3DSolution.vtk";
                }
                else
                {
                    vtkfilename = "../Scalar3DSolution.vtk";
                }
                
                if(0)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    if(!scalarproblem)
                    {
                        vecnames.Push("State");
                        scalnames.Push("StressX");
                        scalnames.Push("StressY");
                        scalnames.Push("StressZ");
                    } else
                    {
                        scalnames.Push("State");
                    }
                    Analysis->DefineGraphMesh(3, scalnames, vecnames, vtkfilename);
                    Analysis->PostProcess(1);
                }
                
                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }

                std::cout << "Post processing\n";

                TPZManVector<REAL> errors(3,0.);
                Analysis->SetThreadsForError(numthreads);
                Analysis->PostProcessError(errors,false);
                
                std::stringstream sout;
                if (!scalarproblem) {
                    
                    sout << "../Elast3DSolutionBeam.txt";
                }
                else
                {
                    sout << "../Scalar3DSolution.txt";
                }
                
                {
                    std::ofstream results(sout.str(),std::ios::app);
                    results.precision(15);
                    results << "(* nelx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq <<
                    " elast " << ExactElast.fE << " poisson "  << ExactElast.fPoisson << " *)" << std::endl;
                    TPZFMatrix<double> errmat(1,3);
                    for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                    std::stringstream varname;
                    varname << "Errmat[[" << nelxcount << "," << irefskeleton+1 << "," << POrder << "]] = (1/1000000)*";
                    errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                }
                
                std::cout << "Plotting shape functions\n";
                bool plotshape = false;
                TPZFNMatrix<3,REAL> sol = SBFem->Solution();
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
                        Analysis->DefineGraphMesh(3, scalnames, vecnames, "../ShapeFunctions.vtk");
                        Analysis->PostProcess(3);
                    }
                }
                
                delete Analysis;
                delete SBFem;
            }
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}
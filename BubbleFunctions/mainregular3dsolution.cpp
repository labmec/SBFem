#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common3D.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzaxestools.h"
#include "pzcondensedcompel.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void AnalyseSolution(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    int minnelx = 1;
    int maxnelx = 5;
    int minrefskeleton = 0;
    int maxrefskeleton = 1;
    int minporder = 1;
    int maxporder = 4;
    int counter = 1;
	int numthreads = 32;
    bool elast = true;

    ExactElast.fProblemType = TElasticity3DAnalytic::Etest2;
    ExactLaplace.fExact = TLaplaceExample1::ESinSin;
    ExactElast.fE = 1.;
    ExactElast.fPoisson = 0.2;
    
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            for(int nelxcount = minnelx; nelxcount < maxnelx; nelxcount++)
            {
                TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;

                int nelx =  (1 << nelxcount-1);
                TPZCompMesh *SBFem = SetupSquareMesh3D(nelx,irefskeleton,POrder, elast);
                int64_t nel = SBFem->NElements();
                // if (TPZSBFemElementGroup::gDefaultPolynomialOrder != 0)
                // {
                //     for (int64_t el = 0; el<nel; el++) {
                //         TPZCompEl *cel = SBFem->Element(el);
                //         if(!cel) continue;
                //         TPZSBFemElementGroup *sbgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
                //         if(!sbgr) continue;
                //         TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel,false);
                //         if (nelxcount == 1)
                //         {
                //             std::cout << "el = " << sbgr->Index() << "," << sbgr->EigenValues() << std::endl;
                //         }
                //     }
                // }

                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZLinearAnalysis Analysis(SBFem,mustOptimizeBandwidth);
                Analysis.SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(&Analysis, SBFem, numthreads);
                
                std::cout << "Plotting\n";
                
                if(elast)
                {
                    Analysis.SetExact(ExactElast.ExactSolution());
                }
                else
                {
                    Analysis.SetExact(ExactLaplace.ExactSolution());
                }
                
                int64_t neq = SBFem->Solution().Rows();
                
                std::string vtkfilename;
                if (elast) {
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
                    if(elast)
                    {
                        vecnames.Push("State");
                        scalnames.Push("StressX");
                        scalnames.Push("StressY");
                        scalnames.Push("StressZ");
                    } else
                    {
                        scalnames.Push("State");
                    }
                    Analysis.DefineGraphMesh(3, scalnames, vecnames, vtkfilename);
                    Analysis.PostProcess(1);
                }
                
                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }

                std::cout << "Post processing\n";

                TPZManVector<REAL> errors(3,0.);
                Analysis.SetThreadsForError(numthreads);
                Analysis.PostProcessError(errors);
                
                std::stringstream sout;
                if (elast) {
                    
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
                    varname << "Errmat[[" << nelxcount+1 << "," << irefskeleton+1 << "," << POrder << "]] = (1/1000000)*";
                    errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                }
                
                std::cout << "************** END OF SIMULATION **************\n\n" << std::endl;
                
                delete SBFem;
            }
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}






void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

int64_t SBFemGroup(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *grp = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(grp) return el;
    }
    return -1;
}
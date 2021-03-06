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

#ifdef _AUTODIFF
void AnalyseSolution(TPZCompMesh *cmesh);
#endif

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minnelx = 1;
    int maxnelx = 5;
    int minrefskeleton = 0;
    int maxrefskeleton = 1;
    int minporder = 1;
    int maxporder = 4;
    int counter = 1;
	int numthreads = 4;
    bool elast = false;

#ifdef _AUTODIFF
    ExactElast.fProblemType = TElasticity3DAnalytic::Etest2;
    ExactLaplace.fExact = TLaplaceExample1::ESinSin;
    ExactElast.fE = 1.;
    ExactElast.fPoisson = 0.2;
#endif
    
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            for(int nelxcount = minnelx; nelxcount < maxnelx; nelxcount++)
            {
                TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;

                int nelx =  (1 << nelxcount-1);
                TPZCompMesh *SBFem = SetupSquareMesh3D(nelx,irefskeleton,POrder, elast);
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    SBFem->Print(sout);
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                int64_t nel = SBFem->NElements();
                if (TPZSBFemElementGroup::gDefaultPolynomialOrder != 0)
                {
                    for (int64_t el = 0; el<nel; el++) {
                        TPZCompEl *cel = SBFem->Element(el);
                        if(!cel) continue;
                        TPZSBFemElementGroup *sbgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
                        if(!sbgr) continue;
                        TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel,false);
                        if (nelxcount == 1)
                        {
                            std::cout << "el = " << sbgr->Index() << "," << sbgr->EigenValues() << std::endl;
                        }
                    }
                }

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
                
#ifdef _AUTODIFF
                if(elast){
                    Analysis->SetExact(Elasticity_exact);
                }
                else{
                    Analysis->SetExact(Laplace_exact);
                }

                // ElasticAnalysis->SetExact(Singular_exact);
#endif

                
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
                Analysis->SetThreadsForError(4);
                Analysis->PostProcessError(errors);
                
                
#ifdef _AUTODIFF
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
                    varname << "Errmat[[" << nelxcount+1 << "]][[" << irefskeleton+1 << "]][[" << POrder << "]] = (1/1000000)*";
                    errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                }
#endif
                
                cout << "************** END OF SIMULATION **************\n\n" << endl;
                
                delete Analysis;
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

#ifdef _AUTODIFF
void AnalyseSolution(TPZCompMesh *cmesh)
{
    int64_t el = SBFemGroup(cmesh);
    TPZSBFemElementGroup *elgrp = dynamic_cast<TPZSBFemElementGroup *>(cmesh->Element(el));
    auto subels = elgrp->GetElGroup();
    int nsub = subels.size();
    for (int isub = 0; isub < nsub; isub++) {
        TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(subels[isub]);
        TPZGeoEl *gel = vol->Reference();
        std::cout << "\n\n\n*********** ELEMENT " << isub << " ****************\n\n\n\n";
        for (int i=-1; i<2; i+=2) {
            for (int j=-1; j<2; j+=2) {
                TPZManVector<REAL,3> x(3,-1.), xco(3);
                x[0] = i;
                x[1] = j;
                gel->X(x, xco);
                TPZManVector<STATE,3> solex(3);
                TPZFNMatrix<9,STATE> dsolex(3,3);
                ExactElast.Solution(xco, solex, dsolex);
                
                TPZSolVec sol;
                TPZGradSolVec dsolax;
                TPZFNMatrix<9,REAL> axes(3,3);
                TPZFNMatrix<9,STATE> dsol(3,3), diff(3,3);
                vol->ComputeSolution(x, sol, dsolax, axes);
//                static void Axes2XYZ(const TPZFMatrix<TVar> &dudaxes, TPZFMatrix<TVar> &dudx, const TPZFMatrix<REAL> &axesv, bool colMajor = true){
                TPZAxesTools<STATE>::Axes2XYZ(dsolax[0], dsol, axes);
                diff = dsol-dsolex;
                REAL err = Norm(diff);
                if (err > 1.e-8)
                {
                    std::cout << "xco = " << x << " sol fem " << sol[0] << " dsol fem " << dsol << " axes " << axes << std::endl;
                    std::cout << "xco = " << x << " sol exa " << solex <<  " dsol exa " << dsolex << std::endl;
                    std::cout << "diff " << diff << std::endl;
                }
            }
        }
    }
}
#endif

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"
#include "TPZSBFemElementGroup.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExactLower;
TElasticity2DAnalytic ElastExactUpper;
TLaplaceExample1 LaplaceExactLower;
TLaplaceExample1 LaplaceExactUpper;
#endif

void IntegrateDirect(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = false;
    bool hasexact = true;

    int numrefskeleton = 4;
    int maxporder = 4;
    int counter = 1;

    int numthreads = 8;
#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::ESquareRoot;
    ElastExact.fProblemType = TElasticity2DAnalytic::ESquareRoot;
    ElastExact.gE = 10;
    ElastExact.gPoisson = 0.3;
    ElastExact.fPlaneStress = 1;
#endif
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            bool elastic = !scalarproblem;
            TPZCompMesh *SBFem = SetupCrackedOneElement(irefskeleton, POrder, hasexact, elastic);
            std::ofstream out("Geometry.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(SBFem->Reference(), out, true);
#ifdef _AUTODIFF
            ElastExactLower = ElastExact;
            ElastExactUpper = ElastExact;
            ElastExactLower.fProblemType = TElasticity2DAnalytic::ESquareRootLower;
            ElastExactUpper.fProblemType = TElasticity2DAnalytic::ESquareRootUpper;
            LaplaceExactLower = LaplaceExact;
            LaplaceExactUpper = LaplaceExact;
            LaplaceExactLower.fExact = TLaplaceExample1::ESquareRootLower;
            LaplaceExactUpper.fExact = TLaplaceExample1::ESquareRootUpper;
            
            if (elastic)
            {
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Emat1);
                    mat->SetForcingFunction(ElastExactLower.ForcingFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Emat2);
                    mat->SetForcingFunction(ElastExact.ForcingFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Emat3);
                    mat->SetForcingFunction(ElastExactUpper.ForcingFunction());
                }
            }
            else
            {
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Emat1);
                    mat->SetForcingFunction(LaplaceExactLower.ForcingFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Emat2);
                    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Emat3);
                    mat->SetForcingFunction(LaplaceExactUpper.ForcingFunction());
                }
            }
            
            if (elastic)
            {
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
                    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                    bc->SetType(0);
                    mat->SetForcingFunction(ElastExactLower.TensorFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
                    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                    bc->SetType(0);
                    mat->SetForcingFunction(ElastExact.TensorFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
                    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                    bc->SetType(0);
                    mat->SetForcingFunction(ElastExactUpper.TensorFunction());
                }
            }
            else
            {
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
                    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                    bc->SetType(0);
                    mat->SetForcingFunction(LaplaceExactLower.TensorFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
                    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                    bc->SetType(0);
                    mat->SetForcingFunction(LaplaceExact.TensorFunction());
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
                    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                    bc->SetType(0);
                    mat->SetForcingFunction(LaplaceExactUpper.TensorFunction());
                }
            }
            
#endif
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
            TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(*Analysis, SBFem, numthreads);
            
            std::cout << "Post processing\n";
#ifdef _AUTODIFF
            Analysis->SetExact(Laplace_exact);
            if (elastic)
            {
                Analysis->SetExact(Elasticity_exact);
            }
#endif
            
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
                scalnames.Push("EpsX");
                scalnames.Push("EpsY");
                scalnames.Push("EpsXY");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, filename.str());
                Analysis->PostProcess(3);
            }
            else
            {
                std::stringstream filename;
                filename << "SquareRootOneElement_NR_" << irefskeleton << "_P_" << POrder << ".vtk";
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("State");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, filename.str());
                int res = POrder+1;
                if (res >5) {
                    res = 5;
                }
                Analysis->PostProcess(res);
            }
            

            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            if(hasexact)
            {
            
                std::cout << "Compute errors\n";
                
                Analysis->PostProcessError(errors);
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
            if(0)
            {
                std::cout << "Plotting shape functions\n";
                int numshape = 25;
                if (numshape > SBFem->NEquations()) {
                    numshape = SBFem->NEquations();
                }
                TPZVec<int64_t> eqindex(numshape);
                for (int i=0; i<numshape; i++) {
                    eqindex[i] = i;
                }
                Analysis->ShowShape("OneElementCracked.vtk", eqindex);
            }
            
            delete Analysis;
            delete SBFem;
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void IntegrateDirect(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    bool scalarproblem = true;

    int minnelxcount = 1;
    int maxnelxcount = 5;
    
    int minporder = 1;
    int maxporder = 3;
    
    int numrefskeleton = 1;
    
    bool usesbfem = true;
    if (usesbfem == false) {
        numrefskeleton = 1;
    }

    int numthreads = 0;
#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
    ElastExact.fProblemType = TElasticity2DAnalytic::ELoadedBeam;
#endif

    int counter = 1;
    for ( int POrder = minporder; POrder <= maxporder; POrder ++)
    {
        for (int irefskeleton = 0; irefskeleton <= numrefskeleton; irefskeleton++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            for(int nelxcount = minnelxcount; nelxcount < maxnelxcount; nelxcount ++)
            {
                int nelx = 1 << (nelxcount-1);
                bool useexact = true;
                if(!scalarproblem)
                {
#ifdef _AUTODIFF
                    ElastExact.gE = 10;
                    ElastExact.gPoisson = 0.3;
                    ElastExact.fPlaneStress = 1;
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
                if(0 && !scalarproblem)
                {
#ifdef _AUTODIFF
                    TPZManVector<REAL,3> x(3,0.);
                    TPZFNMatrix<4,STATE> tensor(2,2);
                    for(int i=-1; i<3; i+=2)
                    {
                        for (int j=-1; j<3; j+=2) {
                                x[0] = i;
                                x[1] = j;
                                ElastExact.Sigma(x, tensor);
                                std::cout << "x = " << x << " tensor " << tensor << std::endl;
                        }
                    }
#endif
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
                
                std::cout << "Entering on Analysis \n";
#ifdef USING_BOOST
        	boost::posix_time::ptime t01 = boost::posix_time::microsec_clock::local_time();
#endif		
 		    bool mustOptimizeBandwidth = true;
            TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(Analysis, SBFem, numthreads);
                
#ifdef USING_BOOST
	        boost::posix_time::ptime t02 = boost::posix_time::microsec_clock::local_time();
        	std::cout << "Time for analysis " << t02-t01 << std::endl;
#endif

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
                              
                if(scalarproblem)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    scalnames.Push("State");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularSolution.vtk");
                    Analysis->PostProcess(3);
                }
                else
                {
                    TPZStack<std::string> vecnames,scalnames;
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
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                std::cout << "Compute errors\n";
                
		        int64_t neq = SBFem->NEquations();

                TPZManVector<REAL,10> errors(3,0.);
                Analysis->SetThreadsForError(numthreads);
                Analysis->PostProcessError(errors);
                
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
                    Analysis->ShowShape("Heterogeneous.vtk", eqindex);
                }
                
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

#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

void IntegrateDirect(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr) {
            TPZStack<TPZCompEl *,5> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            TPZElementMatrix ekvol, efvol, ekgrp, efgrp;
            elgr->CalcStiff(ekgrp, efgrp);
            for (int iv=0; iv<nvol; iv++) {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                TPZElementMatrix ek,ef;
                elvol->CalcStiff(ek, ef);
                if (iv==0) {
                    ekvol = ek;
                    efvol = ef;
                }
                else
                {
                    ekvol.fMat += ek.fMat;
                    efvol.fMat += ef.fMat;
                }
            }
//            ekgrp.fMat.Print("EKGRP = ",std::cout,EMathematicaInput);
//            ekvol.fMat.Print("EKVOL = ",std::cout,EMathematicaInput);
            ekvol.fMat -= ekgrp.fMat;
            std::cout << "IntegrateDirect Norm of difference " << Norm(ekvol.fMat) << std::endl;
            break;
        }
    }

    
}


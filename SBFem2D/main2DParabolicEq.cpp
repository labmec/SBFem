#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZSBFemElementGroup.h"
#include "pzbndcond.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

struct locconfig
{
    int porder;
    int refskeleton;
    int nelxcount;
    int nelx;
    int neq;
    REAL delt;
    int postprocfreq;
    int64_t nsteps;
};

locconfig LocalConfig;

// Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZAnalysis an, REAL delt, int nsteps, int numthreads);

// Setting the boundary conditions
void SetBoundaryConditions(TPZCompMesh *SBFem);

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
    bool scalarproblem = true;
    int maxnelxcount = 5;
    int maxrefskeleton = 3;
    int maxporder = 4;
    int numthreads = 4;
    bool useexact = true;
    LocalConfig.delt = 1./20000.;
    LocalConfig.postprocfreq = 2000;
    LocalConfig.nsteps = 20001;
#ifdef _AUTODIFF
    TimeLaplaceExact.fProblemType = TLaplaceExampleTimeDependent::ECos;
    TimeLaplaceExact.fTime = 0.;
    TimeLaplaceExact.fDelt = LocalConfig.delt;
#endif

    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        LocalConfig.porder = POrder;
        for (int irefskeleton = 0; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            LocalConfig.refskeleton = irefskeleton;
            for(int nelxcount = 3; nelxcount < maxnelxcount; nelxcount += 1)
            {
                int nelx = 2 << (nelxcount-1);

                TPZCompMesh *SBFem = SetupSquareMesh(nelx, irefskeleton, POrder, scalarproblem, useexact);
                SetBoundaryConditions(SBFem);
                
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

                LocalConfig.nelxcount = nelxcount;
                LocalConfig.nelx = nelx;
                LocalConfig.neq = SBFem->NEquations();

                bool mustOptimizeBandwidth = true;
                TPZAnalysis Analysis(SBFem, mustOptimizeBandwidth);
                std::cout << "neq = " << LocalConfig.neq << std::endl;
#ifdef _AUTODIFF
                Analysis.SetExact(TimeLaplace_exact);
#endif
                SolveParabolicProblem(Analysis, LocalConfig.delt, LocalConfig.nsteps, numthreads);          
                
                delete SBFem;
            }
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

/// set the timestep of all SBFem Element groups
void SetSBFemTimestep(TPZCompMesh *CMesh, REAL delt)
{
    int64_t nel = CMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = CMesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!elgr) {
            continue;
        }
        if (delt > 0.) {
            elgr->SetComputeTimeDependent(delt);
        } else
        {
            elgr->SetComputeOnlyMassMatrix();
        }
    }
}


void PostProcess(TPZAnalysis *Analysis, int step)
{
    TPZManVector<REAL,10> errors;
    std::cout << "Compute errors\n";

    Analysis->PostProcessError(errors);
    
    std::stringstream sout;
    sout << "../ParabolicSolutionErrors.txt";

    std::ofstream results(sout.str(),std::ios::app);
    results.precision(15);
    results << "(* nx " << LocalConfig.nelx << " numrefskel " << LocalConfig.refskeleton << " " << " POrder " << LocalConfig.porder << " neq " << LocalConfig.neq << "*)" << std::endl;
    TPZFMatrix<double> errmat(1,6);
    for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
    errmat(0,3) = LocalConfig.neq;
    errmat(0,4) = LocalConfig.nelx;
#ifdef _AUTODIFF
    errmat(0,5) = TimeLaplaceExact.fTime;
#endif
    std::stringstream varname;
    varname << "Errmat[[" << step << "]][[" << LocalConfig.nelxcount << "]][[" << LocalConfig.refskeleton+1 << "]][[" << LocalConfig.porder << "]] = (1/1000000)*";
    errmat.Print(varname.str().c_str(),results,EMathematicaInput);

}

//    Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZAnalysis an, REAL delt, int nsteps, int numthreads)
{
    TPZCompMesh *Cmesh = an.Mesh();
    
#ifdef _AUTODIFF
    TimeLaplaceExact.fDelt = delt;
#endif
    SetSBFemTimestep(Cmesh, delt);
    
#ifndef USING_MKL
    TPZSkylineStructMatrix strmat(Cmesh);
#else
    TPZSymetricSpStructMatrix strmat(Cmesh);
#endif
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    
    std::set<int> matids;
    matids.insert(Emat1);
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    matids.insert(Ebc4);
    strmat.SetMaterialIds(matids);
    strmat.SetNumThreads(numthreads);
    an.SetStructuralMatrix(strmat);
    an.Assemble();

    TPZAutoPointer<TPZMatrix<STATE> > stiff = an.Solver().Matrix();
    TPZFMatrix<STATE> rhs = an.Rhs();
    
    matids.clear();
    matids.insert(Emat1);
    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);
    
    SetSBFemTimestep(Cmesh, 0.);
    
    an.Solver().ResetMatrix();
    an.Assemble();

    TPZAutoPointer<TPZMatrix<STATE> > mass = an.Solver().Matrix();

    an.Solver().ResetMatrix();
    an.Solver().SetMatrix(stiff);

    // project the initial solution
    if(an.GetStep() == 0)
    {
        TPZAnalysis an2(an.Mesh(),false);
#ifndef USING_MKL
        TPZSkylineStructMatrix strmat(Cmesh);
#else
        TPZSymetricSpStructMatrix strmat(Cmesh);
#endif
        strmat.SetNumThreads(numthreads);
        std::set<int> matids;
        matids.insert(ESkeleton);
        matids.insert(Ebc1);
        matids.insert(Ebc2);
        matids.insert(Ebc3);
        matids.insert(Ebc4);
        strmat.SetMaterialIds(matids);
        an2.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        an2.SetSolver(step);
        an2.Run();
        bool plotinitsol = true;
        if(plotinitsol)
        {
            TPZStack<std::string> vecnames,scalnames;
            scalnames.Push("State");
            an2.DefineGraphMesh(2, scalnames, vecnames, "InitialSolution.vtk");
            an2.PostProcess(2);
        }
        an.LoadSolution(an2.Solution());
        std::cout << "compmesh solution norm " << Norm(Cmesh->Solution()) << std::endl;
    }
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << std::endl;
#endif
    
    for (int istep = 0; istep < nsteps; istep++)
    {
        
        if(istep%LocalConfig.postprocfreq == 0)
        {
            std::cout << "\n";
            int postprocindex = istep/LocalConfig.postprocfreq + 1;
            PostProcess(&an, postprocindex);
        }
        TPZFMatrix<STATE> rhstimestep;
        mass->MultAdd(an.Solution(), rhs, rhstimestep, 1./delt, 1.);
        an.Rhs() = rhstimestep;
        an.Solve();
        if(istep%LocalConfig.postprocfreq == 0)
        {
            std::stringstream sout;
            sout << "../Parabolic_nelx" << LocalConfig.nelx << "_p" << LocalConfig.porder << "_refsk" <<
                LocalConfig.refskeleton << ".vtk";
            TPZStack<std::string> vecnames,scalnames;
            // scalar
            scalnames.Push("State");
            an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
            an.PostProcess(2);
        }
#ifdef _AUTODIFF
        TimeLaplaceExact.fTime += delt;
#endif
    }
}

void SetBoundaryConditions(TPZCompMesh * SBFem)
{
    {
        TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        bnd->SetType(0);
#ifdef _AUTODIFF
        bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
    }
    {
        TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        bnd->SetType(0);
#ifdef _AUTODIFF
        bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
    }
    {
        TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        bnd->SetType(0);
#ifdef _AUTODIFF
        bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
    }
    {
        TPZMaterial *mat = SBFem->FindMaterial(Ebc4);
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        bnd->SetType(0);
#ifdef _AUTODIFF
        bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
    }
    {
        TPZMaterial *mat = SBFem->FindMaterial(ESkeleton);
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        bnd->SetType(0);
#ifdef _AUTODIFF
        bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
    }
}
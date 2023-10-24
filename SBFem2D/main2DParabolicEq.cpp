#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZSBFemElementGroup.h"
#include "TPZBndCondT.h"

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

void SubstituteMaterialObjects(TPZCompMesh *cmesh);

void InitializeSolution(TPZCompMesh *cmesh);

//    Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZLinearAnalysis *an, REAL delt, int nsteps, int numthreads);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;

    int maxnelxcount = 5;
    int maxrefskeleton = 3;
    int maxporder = 4;
    int counter = 1;
    
    TimeLaplaceExact.fProblemType = TLaplaceExampleTimeDependent::ECos;
    for ( int POrder = 2; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            for(int nelxcount = 1; nelxcount < maxnelxcount; nelxcount += 1)
            {
                int nelx = 2 << (nelxcount-1);
                bool useexact = true;
                

                TPZCompMesh *SBFem = SetupSquareMesh(nelx,irefskeleton,POrder, scalarproblem,useexact);
                auto exactsol = [](const TPZVec<REAL>&x, TPZVec<STATE>&u,
                                    TPZFMatrix<STATE>&du){
                    TimeLaplaceExact.Exact()->Execute(x, u, du);
                };
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
                    auto bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
                    bnd->SetType(0);
                    bnd->SetForcingFunctionBC(exactsol);
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
                    auto bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
                    bnd->SetType(0);
                    bnd->SetForcingFunctionBC(exactsol);
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
                    auto bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
                    bnd->SetType(0);
                    bnd->SetForcingFunctionBC(exactsol);
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc4);
                    auto bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
                    bnd->SetType(0);
                    bnd->SetForcingFunctionBC(exactsol);
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(ESkeleton);
                    auto bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
                    bnd->SetType(0);
                    bnd->SetForcingFunctionBC(exactsol);
                }
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                LocalConfig.porder = POrder;
                LocalConfig.refskeleton = irefskeleton;
                LocalConfig.nelxcount = nelxcount;
                LocalConfig.nelx = nelx;
                LocalConfig.neq = SBFem->NEquations();
                LocalConfig.delt = 1./20000.;
                LocalConfig.postprocfreq = 2000;
                LocalConfig.nsteps = 20001;
                TimeLaplaceExact.fTime = 0.;
                TimeLaplaceExact.fDelt = LocalConfig.delt;

                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                auto Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
                std::cout << "neq = " << LocalConfig.neq << std::endl;
                int numthreads = 0;
                Analysis->SetExact(TimeLaplaceExact.ExactSolution());
                SolveParabolicProblem(Analysis, LocalConfig.delt, LocalConfig.nsteps, numthreads);
                
                
                delete Analysis;
                delete SBFem;
            }
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void SwitchComputationMode(TPZCompMesh *cmesh, TPZSBFemElementGroup::EComputationMode mode, REAL delt)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(!elgr) continue;
        switch (mode) {
            case TPZSBFemElementGroup::EStiff:
                elgr->SetComputeStiff();
                break;
            case TPZSBFemElementGroup::EMass:
                elgr->SetComputeTimeDependent(delt);
                break;
            case TPZSBFemElementGroup::EOnlyMass:
                elgr->SetComputeOnlyMassMatrix();
                break;
            default:
                DebugStop();
                break;
        }
    }
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

void PostProcess(TPZLinearAnalysis *Analysis, int step)
{
    TPZManVector<REAL,10> errors;
    std::cout << "Compute errors\n";

    Analysis->PostProcessError(errors, false);

    std::stringstream sout;
    sout << "../ParabolicSolutionErrors.txt";

    std::ofstream results(sout.str(),std::ios::app);
    results.precision(15);
    results << "(* nx " << LocalConfig.nelx << " numrefskel " << LocalConfig.refskeleton << " " << " POrder " << LocalConfig.porder << " neq " << LocalConfig.neq << "*)" << std::endl;
    TPZFMatrix<double> errmat(1,6);
    for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
    errmat(0,3) = LocalConfig.neq;
    errmat(0,4) = LocalConfig.nelx;
    errmat(0,5) = TimeLaplaceExact.fTime;
    std::stringstream varname;
    varname << "Errmat[[" << step << "," << LocalConfig.nelxcount << "," << LocalConfig.refskeleton+1 << "," << LocalConfig.porder << "]] = (1/1000000)*";
    errmat.Print(varname.str().c_str(),results,EMathematicaInput);

}

//    Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZLinearAnalysis *an, REAL delt, int nsteps, int numthreads)
{
    TPZCompMesh *Cmesh = an->Mesh();
    
    TimeLaplaceExact.fDelt = delt;
    
    SetSBFemTimestep(Cmesh, delt);
    TPZSkylineStructMatrix<STATE> strmat(Cmesh);
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    std::set<int> matids;
    matids.insert(Emat1);
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    matids.insert(Ebc4);
    strmat.SetMaterialIds(matids);
    strmat.SetNumThreads(numthreads);
    an->SetStructuralMatrix(strmat);
    an->Assemble();

    TPZSolver * solver = an->Solver();
    auto ssolver = dynamic_cast<TPZMatrixSolver<STATE> * >(solver);
    auto stiff = ssolver->Matrix();
    auto rhs = an->Rhs();
    
    matids.clear();
    matids.insert(Emat1);
    strmat.SetMaterialIds(matids);
    an->SetStructuralMatrix(strmat);
    
    SetSBFemTimestep(Cmesh, 0.);
    
    an->Solver()->ResetMatrix();
    an->Assemble();
    solver = an->Solver();
    ssolver = dynamic_cast<TPZMatrixSolver<STATE> * >(solver);
    auto mass = ssolver->Matrix();

    an->Solver()->ResetMatrix();
    
    
    ssolver->SetMatrix(stiff);

    // project the initial solution
    if(an->GetStep() == 0)
    {
        TPZLinearAnalysis an2(an->Mesh(),false);
        TPZSkylineStructMatrix<STATE> strmat(Cmesh);

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
        step.SetDirect(ELDLt);
        an2.SetSolver(step);
        an2.Run();
        if(1)
        {
            TPZStack<std::string> vecnames,scalnames;
            // scalar
            scalnames.Push("Solution");
            an2.DefineGraphMesh(2, scalnames, vecnames, "InitialSolution.vtk");
            an2.PostProcess(2);
        }
        TPZFMatrix<STATE> sol = an2.Solution();
        an->LoadSolution(sol);
        std::cout << "compmesh solution norm " << Norm(Cmesh->Solution()) << std::endl;
    }
    
    for (int istep = 0; istep < nsteps; istep++)
    {
        if(istep%LocalConfig.postprocfreq == 0)
        {
            std::cout << "\n";
            int postprocindex = istep/LocalConfig.postprocfreq + 1;
            an->SetThreadsForError(numthreads);
            PostProcess(an, postprocindex);
        }
        
        TPZFMatrix<STATE> rhstimestep;
        mass->MultAdd(an->Solution(), rhs, rhstimestep, 1./delt, 1.);
        an->Rhs() = rhstimestep;
        an->Solve();
        if(istep%LocalConfig.postprocfreq == 0)
        {
            std::stringstream sout;
            sout << "../Parabolic_nelx" << LocalConfig.nelx << "_p" << LocalConfig.porder << "_refsk" <<
                LocalConfig.refskeleton << ".vtk";
            TPZStack<std::string> vecnames,scalnames;
            // scalar
            scalnames.Push("Solution");
            an->DefineGraphMesh(2, scalnames, vecnames, sout.str());
            an->PostProcess(2);
        }
        TimeLaplaceExact.fTime += delt;
    }
}


#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include <ctime>

#include "TPZGenGrid2D.h"

#include "TPZBuildSBFem.h"

#include "TPZSBFemElementGroup.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZBndCond.h"

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

TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;

    int maxnelxcount = 4;
    int maxrefskeleton = 1;
    int maxporder = 4;
    int counter = 1;
    bool useexact = true;
    LocalConfig.delt = 1./10.;
    LocalConfig.postprocfreq = 5;
    LocalConfig.nsteps = 11;
    TimeLaplaceExact.fTime = 0.;
    TimeLaplaceExact.fDelt = LocalConfig.delt;
    
    TimeLaplaceExact.fProblemType = TLaplaceExampleTimeDependent::ESin;
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            for(int nelxcount = 1; nelxcount < maxnelxcount; nelxcount++)
            {
                int nelx = 2 << (nelxcount-1);
                
                TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;

                LocalConfig.porder = POrder;
                LocalConfig.refskeleton = irefskeleton;
                LocalConfig.nelxcount = nelxcount;
                LocalConfig.nelx = nelx;

                TPZCompMesh *SBFem = SetupSquareMesh(nelx,irefskeleton,POrder);
                LocalConfig.neq = SBFem->NEquations();
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;

                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
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

TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder)
{
    TPZManVector<REAL,4> x0(3,-1.),x1(3,1.);
    x0[0] = -1;
    x0[1] = -1;
    x1[0] = 1;
    x1[1] = 1;
    x0[2] = 0.;
    x1[2] = 0.;
    TPZManVector<int,4> nx(2,nelx);
    TPZGenGrid2D gengrid(nx,x0,x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    
    gengrid.Read(gmesh,EGroup);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc1);
    gengrid.SetBC(gmesh, 6, Ebc1);
    gengrid.SetBC(gmesh, 7, Ebc1);
    
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    build.StandardConfiguration();
    build.DivideSkeleton(nrefskeleton);
    
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);

    auto forcingfunction = [](const TPZVec<REAL>&x, TPZVec<STATE>&u){
        TimeLaplaceExact.ForcingFunction()->Execute(x, u);
    };
    auto exactsol = [](const TPZVec<REAL>&x, TPZVec<STATE>&u,
                            TPZFMatrix<STATE>&du){
        TimeLaplaceExact.Exact()->Execute(x, u, du);
    };
    
    TPZDarcyFlow *matloc = new TPZDarcyFlow(Emat1,SBFem->Dimension());
    matloc->SetForcingFunction(forcingfunction, porder);
    SBFem->InsertMaterialObject(matloc);

    int nstate = 1;
    TPZFMatrix<STATE> val1(nstate,1,0.);
    TPZManVector<STATE> val2(nstate,0.);
    auto BCond1 = matloc->CreateBC(matloc, Ebc1, 0, val1, val2);
    BCond1->SetForcingFunctionBC(exactsol);
    SBFem->InsertMaterialObject(BCond1);

    auto BSkeleton = matloc->CreateBC(matloc, ESkeleton, 1, val1, val2);
    SBFem->InsertMaterialObject(BSkeleton);

    build.BuildComputationMesh(*SBFem);
    return SBFem;
}

/// set the timestep of all SBFem Element groups
void SetSBFemTimestep(TPZCompMesh *CMesh, REAL delt)
{
    for (auto cel : CMesh->ElementVec())
    {
        if (!cel) continue;
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!elgr) continue;
        elgr->SetComputeFullBubbleStiff();
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
    varname << "Errmat[[" << step << "," << LocalConfig.nelxcount << "," << LocalConfig.porder << "]] = (1/1000000)*";
    errmat.Print(varname.str().c_str(),results,EMathematicaInput);

}

//    Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZLinearAnalysis *an, REAL delt, int nsteps, int numthreads)
{
    TPZCompMesh *Cmesh = an->Mesh();
    
    TimeLaplaceExact.fDelt = delt;
    TimeLaplaceExact.fTime = 0;
    
    SetSBFemTimestep(Cmesh, 0.);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(Cmesh);
#endif
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    strmat.SetNumThreads(numthreads);

    an->SetStructuralMatrix(strmat);
    an->Assemble();
    an->Solve();
    
    for (int istep = 0; istep < nsteps; istep++)
    {
        an->AssembleResidual();
        an->Solve();

        if(istep%LocalConfig.postprocfreq == 0)
        {
            std::stringstream sout;
            sout << "../Parabolic_nelx" << LocalConfig.nelx << "_p" << LocalConfig.porder << "_refsk" <<
                LocalConfig.refskeleton << ".vtk";
            TPZStack<std::string> vecnames,scalnames;
            scalnames.Push("Solution");
            an->DefineGraphMesh(2, scalnames, vecnames, sout.str());
            an->PostProcess(2);
            int postprocindex = istep/LocalConfig.postprocfreq + 1;
            PostProcess(an, postprocindex);
        }
        TimeLaplaceExact.fTime += delt;
    }
    std::cout << "*******\n";
}


#ifdef HAVE_CONFIG_H
#include "pz_config.h"
#endif

#include "SBFemHdiv.h"
#include "TPZHybridizeHDiv.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

using namespace std;

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif

//-----------------------------------------

    // Initial data
    int maxporder = 6;
    bool useexact = true;
    int numthreads = 8;

#ifdef PZDEBUG
    // Outputs
    bool printvtk = true;
    bool printcmesh = true;
    bool printstifmatrix = true;
#endif
//-----------------------------------------

#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
#endif
    
    for (int POrder = 1; POrder <= maxporder; POrder ++)
    {
        std::cout << "POrder = " << POrder << "\n";

        std::cout << "Building geometric mesh...\n";
        int nelx = 1;

        auto gmesh = SetupCollapsedMesh(nelx);
        TPZCheckGeom gcheck(gmesh);
        gcheck.CheckIds();
        gcheck.CheckUniqueId();

        std::cout << "Building computational mesh...\n";
        auto cmeshp = pressure(gmesh, POrder);
        auto cmeshf = fluxhdivsbfem(gmesh, POrder);
        auto cmeshm = multiphysicscollapsed(gmesh, cmeshp, cmeshf, POrder);

        // Adjusting mesh
        AddInterfaceElements(cmeshm);
        ofstream sout("cmeshmultiphysics0.txt");
        cmeshm->Print(sout);

        // TPZManVector<int64_t> perm(cmeshm->NConnects(),-1);
        // AdjustExtPressureConnectivity(cmeshm, cmeshf, perm);
        // cout << perm << endl;
        // cmeshm->Permute(perm);
        ofstream sout1("cmeshmultiphysics1.txt");
        cmeshm->Print(sout1);
        
        GroupandCondense(cmeshm);

        cmeshm->ComputeNodElCon();
        cmeshm->CleanUpUnconnectedNodes();

        ofstream sout2("cmeshmultiphysics2.txt");
        cmeshm->Print(sout2);
        
        std::cout << "Analysis...\n";
        std::cout << "neq = " << cmeshm->NEquations() << std::endl;
        
        bool optimizeBandwidth = false;
        TPZAnalysis an(cmeshm, optimizeBandwidth);
        an.Assemble();

        TPZAutoPointer<TPZMatrix<REAL>> K = an.Solver().Matrix();
        std::ofstream sout3("correctstiffness.txt");
        K->Print("ekC = ", sout3, EMathematicaInput);
        
        TPZFMatrix<STATE> E0, E1, E2;
        ComputeMatrices(E0, E1, E2, K);
        ComputeEigenvalues(E0, E1, E2);
        
    }
        
    std::cout << "Check:: Calculation finished successfully\n";
    return EXIT_SUCCESS;
}

// Start: elementos h1 descontinuo.
// Volume - pegar elementos vizinhos
// Criar uma geometria correta - derivado do SBFEMBuildHdiv - OK
// Malha do fluxo criada no main - elementos HdivCollapsed d-1
// SBFEMVolumeHdiv: derivado do SBFEMVolume, criar um método pra me dar a lista de elementos
// SBFEMElementGroupHdiv - conjunto de elementos multifísicos p/ um grupo de Volumes
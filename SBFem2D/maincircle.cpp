#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

//#include "pzgengrid.h"
#include "Elasticity/TPZElasticity2D.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZBndCond.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZLinearAnalysis.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzmultiphysicselement.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

#include "TPZSSpStructMatrix.h"

#include "TPZBuildSBFem.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "tpzgeoelrefpattern.h"

#include "JSON.hpp"
void rect_mesh();

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

// Create a one element mesh going from angle = 0 to angle
TPZCompMesh *SetupOneArc(int numrefskeleton, int porder, REAL angle);

void PostProcessing(TPZLinearAnalysis Analysis, int POrder, int irefskeleton);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initial data:
    int minrefskeleton = 2;
    int maxrefskeleton = 3;
    int minporder = 2;
    int maxporder = 9;
    int countstep = 1;
    int numthreads = 4;
    REAL angle = M_PI*6./4.;
    bool scalarproblem = true;
    int nelxcount = 1;
    // LaplaceExact.fExact = TLaplaceExample1::ESingularCircle;

    for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
    {
        for ( int POrder = minporder; POrder < maxporder; POrder ++)
        {
            TPZCompMesh *SBFem = SetupOneArc(irefskeleton, POrder, angle);
            {
                auto BCond2 = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc2));
                BCond2->SetType(1);
                BCond2->SetForcingFunctionBC(LaplaceExact.ExactSolution());
                auto BC1 = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc1));
                TPZManVector<REAL> val2(1,0.);
                BC1->SetVal2(val2);
            }
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Analysis
            std::cout << "Entering on Analysis\n";
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis Analysis(SBFem,mustOptimizeBandwidth);
            Analysis.SetStep(countstep++);
	        std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(Analysis, SBFem, numthreads);
            bool printcmeshwsol = false;
            if(printcmeshwsol)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            std::cout << "Post processing\n";
            std::string filename = "../SingularSolution";
            PostProcessing(Analysis, filename, scalarproblem, numthreads, POrder, nelxcount, irefskeleton);

            bool printeigval = false;
            if(printeigval)
            {
                PrintEigval(Analysis, filename);
            }

            delete SBFem;
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}


TPZCompMesh *SetupOneArc(int numrefskeleton, int porder, REAL angle)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(4);
    TPZManVector<REAL,3> co(3,0.);
    gmesh->NodeVec()[0].Initialize(co, gmesh);
    co[0] = 1.;
    gmesh->NodeVec()[1].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle);
    co[1] = sin(angle);
    gmesh->NodeVec()[2].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle/2.);
    co[1] = sin(angle/2.);
    gmesh->NodeVec()[3].Initialize(co, gmesh);
    co.Fill(0.);

    TPZManVector<int64_t,4> nodeindex(1,0);
    nodeindex[0] = 1;
    int64_t elementid = 1;
    gmesh->CreateGeoElement(EPoint, nodeindex, Ebc1, elementid);
    
    // Definition of Arc coordenates
    // Create Geometrical Arc #1
    nodeindex.Resize(3);
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 3;
    elementid = 1;
    TPZGeoEl *arc = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (nodeindex, Ebc2, gmesh,elementid);
    
    nodeindex.Resize(4);
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 0;
    nodeindex[3] = 0;
    elementid = 2;
    TPZGeoEl *gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodeindex, EGroup, gmesh,elementid);
    
    gmesh->BuildConnectivity();

    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    TPZManVector<int64_t,5> elids(1,gblend->Index());
    build.AddPartition(elids, 0);
    
    build.DivideSkeleton(numrefskeleton);

    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
    bool applyexact = false;
    InsertMaterialObjects(SBFem,problemtype,applyexact);
    
    build.BuildComputationMesh(*SBFem);
    
    bool outputcmshgmsh = false;
    if (outputcmshgmsh)
    {
        OutputGmshCmsh(gmesh, SBFem);
    }

    return SBFem;
}
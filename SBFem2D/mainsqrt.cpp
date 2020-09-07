#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "TPZVTKGeoMesh.h"

#include "pzbndcond.h"
#include "TPZMatLaplacian.h"

#include "pzfunction.h"
#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

TPZCompMesh *SetupOneArcWithRestraint(int numrefskeleton, int porder, REAL angle);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif

#ifndef _AUTODIFF
    std::cout << "This program needs FAD to run \n";
    DebugStop();
#endif

    int minrefskeleton = 1;
    int maxrefskeleton = 5;
    int minporder = 1;
    int maxporder = 9;
    int counter = 1;
    int numthreads = 4;
    bool scalarproblem = true;
    int nelxcount = 1;
#ifdef _AUTODIFF
    // LaplaceExact.fExact = TLaplaceExample1::ESingularCircle;
#endif

    for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
    {
        for ( int POrder = minporder; POrder < maxporder; POrder += 1)
        {

            REAL angle = M_PI;
            TPZCompMesh *SBFem = SetupOneArcWithRestraint(irefskeleton,POrder, angle);
            {
                TPZBndCond *BCond2 = dynamic_cast<TPZBndCond *>(SBFem->FindMaterial(Ebc2));
                BCond2->SetType(1);
#ifdef _AUTODIFF
                BCond2->SetForcingFunction(0,LaplaceExact.TensorFunction());
#endif
                TPZBndCond *BC1 = dynamic_cast<TPZBndCond *>(SBFem->FindMaterial(Ebc1));
                BC1->Val2()(0,0) = 0;
            }
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Analysis
            std::cout << "Entering on Analysis\n";
            bool mustOptimizeBandwidth = true;
            TPZAnalysis Analysis(SBFem,mustOptimizeBandwidth);
            Analysis.SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(Analysis, SBFem, numthreads);
                        
            std::cout << "Post processing\n";
#ifdef _AUTODIFF
            Analysis.SetExact(Laplace_exact);
#endif
            int64_t neq = SBFem->Solution().Rows();
            std::string filename = "../SqrtSolution";
            PostProcessing(Analysis, filename, scalarproblem, numthreads, POrder, nelxcount, irefskeleton);
            
            bool printeigvalfile = false;
            if(printeigvalfile)
            {
                PrintEigval(Analysis, filename);
            }
            
            std::cout << "Plotting shape functions\n";
            if(POrder == maxporder-1 && irefskeleton == 0)
            {
                int numshape = 25;
                if (numshape > SBFem->NEquations()) {
                    numshape = SBFem->NEquations();
                }
                TPZVec<int64_t> eqindex(numshape);
                for (int i=0; i<numshape; i++) {
                    eqindex[i] = i;
                }
                Analysis.SetStep(0);
                Analysis.ShowShape("Sqrt.vtk", eqindex);
            }
            
            delete SBFem;
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

TPZCompMesh *SetupOneArcWithRestraint(int numrefskeleton, int porder, REAL angle)
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
    
    nodeindex.Resize(2);
    nodeindex[0] = 1;
    nodeindex[1] = 0;
    TPZGeoEl *oned = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (nodeindex, EGroup+1, gmesh, elementid);
    
    gmesh->BuildConnectivity();
    
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    matmap[EGroup+1] = Emat2;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    TPZStack<int64_t> elids;
    elids.Push(gblend->Index());
    elids.Push(oned->Index());
    build.AddPartition(elids, 0);
    
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
	bool apply_exact = false;
    InsertMaterialObjects(SBFem,problemtype, apply_exact);
    
    TPZMaterial *mat1 = SBFem->FindMaterial(Emat1);
    TPZMaterial *mat2 = mat1->NewMaterial();
    TPZMatLaplacian *mat2lapl = dynamic_cast<TPZMatLaplacian *>(mat2);
    mat2->SetId(Emat2);
    mat2lapl->SetDimension(1);
    mat2lapl->SetParameters(1.e9, 0);
    SBFem->InsertMaterialObject(mat2);
    
    std::set<int> volmatids,boundmatids;
    volmatids.insert(Emat1);
    volmatids.insert(Emat2);
    boundmatids.insert(Ebc1);
    boundmatids.insert(Ebc2);
    boundmatids.insert(ESkeleton);
    build.DivideSkeleton(numrefskeleton);
    build.BuildComputationMesh(*SBFem,volmatids,boundmatids);
    
    OutputGmshCmsh(gmesh, SBFem);

    return SBFem;
    
}

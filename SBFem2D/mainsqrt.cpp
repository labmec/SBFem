#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "TPZVTKGeoMesh.h"

#include "TPZBndCond.h"
#include "DarcyFlow/TPZDarcyFlow.h"

#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


REAL mult[] = {1.,10./45.,9./45.,8./45.,7./45.,6./45.,5./45.};

auto SingularNeumann = [](const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL Lambda0 = 1./2.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if(theta < 0.) theta += 2.*M_PI;
    val[0] = 0;
    for (int i=0; i<1; i++) {
        REAL Lambda = Lambda0+i;
        val[0] += mult[i]*Lambda*pow(r,Lambda-1.)*sin(Lambda*theta);
    }
};

auto SingularExact = [](const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL Lambda0 = 1./2;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if (theta < 0.)
    {
        theta += 2.*M_PI;
    }
    
    val[0] = 0;
    deriv.Resize(2,1);
    deriv.Zero();
    for (int i=0; i<1; i++)
    {
        REAL Lambda = Lambda0+i;
        val[0] += mult[i]*pow(r,Lambda)*sin(Lambda*theta);
        deriv(0,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[0]*sin(Lambda*theta)-pow(r,Lambda-2)*(Lambda)*cos(Lambda*theta)*(x[1]));
        deriv(1,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[1]*sin(Lambda*theta)+pow(r,Lambda-2)*(Lambda)*cos(Lambda*theta)*(x[0]));
    }  
};

TPZCompMesh *SetupOneArcWithRestraint(int numrefskeleton, int porder, REAL angle);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 1;
    int maxrefskeleton = 5;
    int minporder = 1;
    int maxporder = 9;
    int counter = 1;
    int numthreads = 2;
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            REAL angle = M_PI;
            TPZCompMesh *SBFem = SetupOneArcWithRestraint(irefskeleton,POrder, angle);
            {
                auto BCond2 = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc2));
                BCond2->SetType(0);
                BCond2->SetForcingFunctionBC(SingularExact,POrder);
                auto BC1 = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc1));
                BCond2->SetType(0);
                TPZManVector<STATE> v2(1,0);
                BC1->SetVal2(v2);
            }

            TPZSBFemElementGroup *celgrp = 0;
            int64_t nel = SBFem->NElements();
            for (int64_t el=0; el<nel; el++) {
                TPZSBFemElementGroup *cel = dynamic_cast<TPZSBFemElementGroup *>(SBFem->Element(el));
                if(cel)
                {
                    celgrp = cel;
                    break;
                }
            }
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(*Analysis, SBFem, numthreads);
                        
            std::cout << "Post processing\n";
            Analysis->SetExact(SingularExact);
            
            
            int64_t neq = SBFem->Solution().Rows();
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("Solution");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, "../SqrtSolution.vtk");
                Analysis->PostProcess(3);
            }
            TPZManVector<REAL> errors(3,0.);
            Analysis->SetThreadsForError(numthreads);
            Analysis->PostProcessError(errors);
            
            std::stringstream sout;
            sout << "../SqrtSolution.txt";
            
            std::ofstream results(sout.str(),std::ios::app);
            results.precision(15);
            results << "(* numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
            TPZFMatrix<double> errmat(1,4);
            for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
            errmat(0,3) = neq;
            std::stringstream varname;
            varname << "Errmat[[" << irefskeleton+1 << "," << POrder << "]] = (1/1000000)*";
            errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            
            if(0)
            {
                std::multimap<REAL,REAL> eigmap;
                TPZManVector<double> eigval = celgrp->EigenvaluesReal();
                TPZFMatrix<double> coef = celgrp->CoeficientsReal();
                for (int i=0; i<eigval.size(); i++) {
                    eigmap.insert(std::pair<REAL,REAL>(eigval[i],coef(i,0)));
                }
                for (std::multimap<REAL, REAL>::reverse_iterator it = eigmap.rbegin(); it!=eigmap.rend(); it++) {
                    results << it->first << "|" << it->second << " ";
                }
            }
            
            delete Analysis;
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
    
    nodeindex.Resize(3);
    // Definition of Arc coordenates
    // Create Geometrical Arc #1
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
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
	bool apply_exact = false;
    InsertMaterialObjects(SBFem,problemtype, apply_exact);
    
    
    TPZMaterial *mat2 = SBFem->FindMaterial(Emat2);
    auto mat2lapl = dynamic_cast<TPZDarcyFlow *>(mat2);
    mat2lapl->SetDimension(1);
    mat2lapl->SetConstantPermeability(1.e12);
    
    std::set<int> volmatids,boundmatids;
    volmatids.insert(Emat1);
    volmatids.insert(Emat2);
    boundmatids.insert(Ebc1);
    boundmatids.insert(Ebc2);
    boundmatids.insert(ESkeleton);
    build.DivideSkeleton(numrefskeleton);

    build.BuildComputationMesh(*SBFem,volmatids,boundmatids);
    
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
    
}
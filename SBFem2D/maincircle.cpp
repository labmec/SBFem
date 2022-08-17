#ifdef HAVE_CONFIG_H
#include <config.h>
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

#include "TPZGenGrid2D.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "tpzautopointer.h"
#include "TPZBndCond.h"

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


REAL mult[] = {10./45.,9./45.,8./45.,7./45.,6./45.,5./45.};

auto SingularNeumann = [](const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL Lambda0 = 2./3.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if(theta < 0.) theta += 2.*M_PI;
    val[0] = 0;
    for (int i=0; i<6; i++) {
        REAL Lambda = Lambda0*(i+1);
        val[0] += mult[i]*Lambda*pow(r,Lambda-1.)*cos(Lambda*theta);
    }
};

TPZCompMesh *SetupOneArc(int numrefskeleton, int porder, REAL angle);

auto Singular_exact = [](const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL Lambda0 = 2./3.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if (theta<0.) {
        theta += 2.*M_PI;
    }
    
    val[0] = 0;
    deriv.Resize(2,1);
    deriv.Zero();
    for (int i=0; i<6; i++) {
        REAL Lambda = Lambda0*(i+1);
        val[0] += mult[i]*pow(r,Lambda)*cos(Lambda*theta);
        deriv(0,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[0]*cos(Lambda*theta)+pow(r,Lambda-2)*(Lambda)*sin(Lambda*theta)*(x[1]));
        deriv(1,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[1]*cos(Lambda*theta)-pow(r,Lambda-2)*(Lambda)*sin(Lambda*theta)*(x[0]));
    }   
};

int main(int argc, char *argv[])
{
    int minrefskeleton = 2;
    int maxrefskeleton = 3;
    int minporder = 2;
    int maxporder = 9;
    int counter = 1;
    int numthreads = 32;
    for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
    {
        for ( int POrder = minporder; POrder < maxporder; POrder += 1)
        {
            REAL angle = M_PI*6./4.;
            TPZCompMesh *SBFem = SetupOneArc(irefskeleton,POrder,angle);
            {
                auto BCond2 = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc2));
                BCond2->SetType(0);
                BCond2->SetForcingFunctionBC(Singular_exact,POrder);
                auto BC1 = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc1));
                BCond2->SetType(0);
                TPZManVector<STATE> v2(1,1);
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
            Analysis->SetExact(Singular_exact);
            
            int64_t neq = SBFem->Solution().Rows();
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("Solution");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, "../SingularSolution.vtk");
                int res = POrder+1;
                if (res >5) {
                    res = 5;
                }
                Analysis->PostProcess(res);
            }
            TPZManVector<REAL> errors(3,0.);
            Analysis->SetThreadsForError(numthreads);
            Analysis->PostProcessError(errors, false);
            
            std::stringstream sout;
            sout << "../SingularSolution.txt";
            
            std::ofstream results(sout.str(),std::ios::app);
            results.precision(15);
            results << "(* numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
            TPZFMatrix<double> errmat(1,3);
            for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
            std::stringstream varname;
            varname << "Errmat[[" << POrder << "," << irefskeleton+1 << "]] = (1/1000000)*";
            errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            
            if(0)
            {
                std::multimap<REAL,REAL> eigmap;
                TPZManVector<double> eigval = celgrp->EigenvaluesReal();
                TPZFMatrix<double> coef = celgrp->CoeficientsReal();
                for (int i=0; i<eigval.size(); i++) {
                    eigmap.insert(std::pair<double,double>(eigval[i],coef(i,0)));
                }
                for (std::multimap<double, double>::reverse_iterator it = eigmap.rbegin(); it!=eigmap.rend(); it++) {
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
    
    gmesh->BuildConnectivity();
    
    //    gmesh->Print(std::cout);
    TPZManVector<REAL,3> xi(1),x(3);
    for (REAL s=-1.; s<=1.; s+= 1./10.) {
        xi[0] = s;
        arc->X(xi, x);
        std::cout << "xi " << xi << " x " << x << std::endl;
    }
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    TPZManVector<int64_t,5> elids(1,gblend->Index());
    build.AddPartition(elids, 0);
    
    build.DivideSkeleton(numrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
    bool applyexact = false;
    InsertMaterialObjects(SBFem,problemtype,applyexact);
    
    
    build.BuildComputationMesh(*SBFem);
    
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

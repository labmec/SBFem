#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "tpzgeoelrefpattern.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZVTKGeoMesh.h"

#include "pzinterpolationspace.h"

#include "TPZBuildSBFem.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

#include "pzfunction.h"

#include "Poisson/TPZMatPoisson.h"


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void OutputFourtyFive(TPZCompMesh *cmesh, REAL radius);

auto DirichletTestProblem = [](const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE>& dval)
{
    dval.Zero();
    REAL theta = atan2(x[1],x[0]);
    if(theta < 0.) theta += 2.*M_PI;
    
    if(theta < M_PI/2.)
    {
        val[0] = 1.;
    }
    else if(theta < M_PI)
    {
        val[0] = 1.-(theta-M_PI/2.)/(M_PI/2);
    }
    else if(theta < 3.*M_PI/2.)
    {
        val[0] = 0.;
    }
    else
    {
        val[0] = (theta -3.*M_PI/2.)/(M_PI/2.);
    }
};

TPZCompMesh *TestHeterogeneous(int numquadrant,TPZVec<REAL> &contrast, REAL radius, int numref, int porder);

int main(int argc, char *argv[])
{
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 1;
    int numrefskeleton = 2;
    int minporder = 1;
    int maxporder = 9;
    int counter = 1;
    int numthreads  = 2;

    for (int irefskeleton = minrefskeleton; irefskeleton < numrefskeleton; irefskeleton++)
    {
        for ( int POrder = minporder; POrder < maxporder; POrder += 1)
        {
            
            int numquadrant = 4;
            REAL radius = 1.;
            TPZManVector<REAL> contrast(4,1.);
            contrast[0] = 100;
            contrast[2] = 100;

            TPZCompMesh *SBFem = TestHeterogeneous(numquadrant , contrast, radius, irefskeleton, POrder);

            TPZSBFemElementGroup *celgrp = 0;
            int64_t nel = SBFem->NElements();
            for (auto cel : SBFem->ElementVec()) {
                auto celsbfem = dynamic_cast<TPZSBFemElementGroup *>(cel);
                if(celsbfem)
                {
                    celgrp = celsbfem;
                    break;
                }
            }
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis Analysis(SBFem,mustOptimizeBandwidth);
            Analysis.SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            // SolveSist(Analysis, SBFem, numthreads);
            
            std::cout << "Post processing\n";
            TPZManVector<STATE> errors(3,0.);
            int64_t neq = SBFem->Solution().Rows();
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("Solution");
                Analysis.DefineGraphMesh(2, scalnames, vecnames, "../Heterogeneous.vtk");
                Analysis.PostProcess(3);
            }
            
            std::stringstream sout;
            sout << "../Heterogeneous.txt";
            
            std::ofstream results(sout.str(),std::ios::app);
            results.precision(15);
            // for circular domain with contrast
            results << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << " contrast " << contrast << std::endl;
            
            if(1)
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
            results << std::endl;
            results << celgrp->EigenValues() << std::endl;
            
            delete SBFem;
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

TPZCompMesh *TestHeterogeneous(int numquadrant,TPZVec<REAL> &contrast, REAL radius, int numref, int porder)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);

    int64_t nodind;
    TPZManVector<REAL,3> co(3,0.);
    nodind = gmesh->NodeVec().AllocateNewElement();
    int64_t centernode = nodind;
    gmesh->NodeVec()[nodind].Initialize(co, gmesh);

    co[0] = radius;
    nodind = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodind].Initialize(co, gmesh);
    TPZManVector<int64_t> scalecenters(2*numquadrant,-1);
    
    int64_t lastnode = nodind;
    int64_t firstnode = nodind;

    for (int quadrant=0; quadrant<numquadrant; quadrant++)
    {
        TPZManVector<int64_t,4> nodes(3,0);
        nodes[0] = lastnode;
        REAL angle = M_PI*(quadrant+1)/2;
        co[0] = radius*cos(angle);
        co[1] = radius*sin(angle);
        if (quadrant == 3)
        {
            nodind = firstnode;
        }
        else
        {
            nodind = gmesh->NodeVec().AllocateNewElement();
        }
        gmesh->NodeVec()[nodind].Initialize(co, gmesh);

        nodes[1] = nodind;
        lastnode = nodind;
        angle = M_PI*(quadrant+1)/2-M_PI/4.;
        co[0] = radius*cos(angle);
        co[1] = radius*sin(angle);
        nodind = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[nodind].Initialize(co, gmesh);

        nodes[2] = nodind;
        int64_t elementindex;
        TPZGeoEl *arc = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (nodes, Ebc1, gmesh,elementindex);
        nodes.resize(4);
        nodes[2] = centernode;
        nodes[3] = centernode;
        int matid = EGroup;
        matid = EGroup+quadrant;
        TPZGeoEl *gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodes, matid, gmesh,elementindex);
        scalecenters[elementindex] = centernode;
    }
    
    gmesh->BuildConnectivity();

    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    matmap[EGroup+1] = Emat2;
    matmap[EGroup+2] = Emat3;
    matmap[EGroup+3] = Emat4;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    build.Configure(scalecenters);
    build.DivideSkeleton(numref);
    
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
	bool applyexact = false;
    InsertMaterialObjects(SBFem,problemtype,applyexact);
    
    {
        TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
        auto BCond1 = dynamic_cast<TPZBndCondT<STATE> * >(mat);
        BCond1->SetForcingFunctionBC(DirichletTestProblem);
    }
    {   
        TPZMatPoisson<STATE> *mat1 = dynamic_cast<TPZMatPoisson<STATE> *> (SBFem->FindMaterial(Emat1));
        TPZMatPoisson<STATE> *mat2 = dynamic_cast<TPZMatPoisson<STATE> *> (SBFem->FindMaterial(Emat2));
        TPZMatPoisson<STATE> *mat3 = dynamic_cast<TPZMatPoisson<STATE> *> (mat1->NewMaterial());
        TPZMatPoisson<STATE> *mat4 = dynamic_cast<TPZMatPoisson<STATE> *> (mat1->NewMaterial());
        STATE K = 1;
        mat1->SetScaleFactor(K*contrast[0]);
        mat2->SetScaleFactor(K*contrast[1]);
        mat3->SetScaleFactor(K*contrast[2]);
        mat4->SetScaleFactor(K*contrast[3]);
        mat3->SetId(Emat3);
        mat4->SetId(Emat4);
        SBFem->InsertMaterialObject(mat3);
        SBFem->InsertMaterialObject(mat4);
    }
    
    build.BuildComputationMesh(*SBFem);
    
    {
        int64_t nel = SBFem->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = SBFem->Element(el);
            TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if (elgr) {
                TPZElementMatrixT<STATE> ek,ef;
                elgr->CalcStiff(ek, ef);
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (intel && intel->NConnects() ==3) {
                TPZGeoEl *ref = intel->Reference();
                TPZManVector<REAL,3> co(3);
                TPZManVector<STATE,3> val(1);
                TPZFMatrix<STATE> dval(3,1);
                ref->NodePtr(0)->GetCoordinates(co);
                DirichletTestProblem(co, val, dval);
                int64_t seqnum = intel->Connect(0).SequenceNumber();
                {
                    TPZFMatrix<STATE> block(1, 1,val[0]);
                    SBFem->Block().PutBlock(seqnum, 0, block);
                }
                ref->NodePtr(1)->GetCoordinates(co);
                DirichletTestProblem(co, val, dval);
                seqnum = intel->Connect(1).SequenceNumber();
                {
                    TPZFMatrix<STATE> block(1, 1,val[0]);
                    SBFem->Block().PutBlock(seqnum, 0, block);
                }
                SBFem->LoadSolution(SBFem->Solution());
            }
        }
    }
    if(1)
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    if(numref)
    {
        OutputFourtyFive(SBFem,radius);
    }
    return SBFem;
    
}

void OutputFourtyFive(TPZCompMesh *cmesh, REAL radius)
{
    TPZGeoMesh *gmesh = cmesh->Reference();
    cmesh->LoadReferences();
    int64_t nelg = gmesh->NElements();
    TPZSBFemVolume *elfound = NULL;
    for (int64_t el=0; el<nelg; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->NodeIndex(0) ==9) {
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neighbour = gelside.Neighbour();
            TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(gel->Reference());
            if (vol) {
                elfound = vol;
                break;
            }
            while (neighbour != gelside)
            {
                TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(neighbour.Element()->Reference());
                if (vol) {
                    elfound = vol;
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if(elfound) break;
        }
    }
    if(!elfound)
    {
        DebugStop();
    }
    else
    {
        TPZFMatrix<REAL> phi,coef,eigvalmatrix;
        phi = elfound->PhiReal();
        coef = elfound->CoeficientsReal();
        eigvalmatrix.Resize(1, phi.Cols());
        TPZManVector<REAL > eig = elfound->EigenvaluesReal();
        for (int64_t i=0; i< eig.size(); i++) {
            eigvalmatrix(0,i) = eig[i];
        }
        std::ofstream out("Diagonal.nb");
        phi.Print("phi = ",out,EMathematicaInput);
        coef.Print("coef = ",out,EMathematicaInput);
        eigvalmatrix.Print("eig = ",out,EMathematicaInput);
        out << "a = Sum[Chop[phi[[1]][[i]] coef[[i]][[1]] ksi^(-eig[[1]][[i]])], {i,1, Length[phi[[1]]]}]\n"
        << "Plot[a, {ksi, 0, 1}, PlotRange -> All]\n"
        << "Plot[a, {ksi, 0, 0.00001}, PlotRange -> All]\n";
    }
}

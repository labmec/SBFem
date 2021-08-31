#include "Common.h"

#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "Elasticity/TPZElasticity2D.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZBndCond.h"

#include "TPZGenGrid2D.h"
#include "TPZBuildSBFem.h"

#include "TPZVTKGeoMesh.h"

#include "JSON.hpp"
#include "TPZSBFemElementGroup.h"
#include "pzinterpolationspace.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"

TLaplaceExample1 LaplaceExact;
TElasticity2DAnalytic ElastExact;
TElasticity2DAnalytic ElastExactUpper;
TElasticity2DAnalytic ElastExactLower;
TLaplaceExampleTimeDependent TimeLaplaceExact;

void SolveSist(TPZLinearAnalysis *an, TPZCompMesh *Cmesh)
{
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(Cmesh);
#endif
    
    // strmat.SetNumThreads(4);
    an->SetStructuralMatrix(strmat);
    
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
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    an->Assemble();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif
    
    
}

void BodyLoads(const TPZVec<REAL> &x, TPZVec<REAL> &val)
{
    val[0] = 0.;
    val[1] = -x[1];
}

void HarmonicNeumannLeft(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = -M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void HarmonicNeumannRight(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void Harmonic_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    val[0] = exp(M_PI*xv[0])*sin(M_PI*xv[1]);
    deriv(0,0) = M_PI*val[0];
    deriv(1,0) = M_PI*exp(M_PI*xv[0])*cos(M_PI*xv[1]);
    
}


void InsertMaterialObjects(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact)
{
    
    // Getting mesh dimension
    int dim = 2;

    TPZMaterial *material;
    int nstate = 1;
    if (!scalarproblem)
    {
        TPZElasticity2D *matloc1 = new TPZElasticity2D(Emat1);
        TPZElasticity2D *matloc2 = new TPZElasticity2D(Emat2);
        material = matloc1;
        nstate = 2;
		if (applyexact)
        {
			matloc1->SetPlaneStrain();
            matloc1->SetElasticity(ElastExact.gE, ElastExact.gPoisson);
            matloc2->SetPlaneStrain();
            matloc2->SetElasticity(ElastExact.gE, ElastExact.gPoisson);
		}
        REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
        matloc1->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
        matloc2->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);

        auto forcingfunction = [](const TPZVec<REAL>&x, TPZVec<STATE>&u){
            ElastExact.ForcingFunction()->Execute(x, u);
        };
        auto exactsol = [](const TPZVec<REAL>&x, TPZVec<STATE>&u,
                            TPZFMatrix<STATE>&du){
            ElastExact.Exact()->Execute(x, u, du);
        };
        if(applyexact)
        {
            constexpr int porder{3};
            matloc1->SetForcingFunction(forcingfunction, porder);
            matloc2->SetForcingFunction(forcingfunction, porder);
        }
        
        TPZFMatrix<STATE> val1(nstate,1,0.);
        TPZManVector<STATE> val2(nstate,0.);
        auto BCond1 = matloc1->CreateBC(matloc1, Ebc1, 0, val1, val2);
        if (applyexact)
        {
            BCond1->SetForcingFunctionBC(exactsol);
        }
        auto BSkeleton = matloc1->CreateBC(matloc1, ESkeleton, 1, val1, val2);

        cmesh->InsertMaterialObject(matloc1);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BSkeleton);
    }
    else
    {
        TPZMatPoisson<STATE> *matloc = new TPZMatPoisson<STATE>(Emat1,cmesh->Dimension());
        constexpr int porder{3};
        
        nstate = 1;

        TPZFMatrix<STATE> val1(nstate,1,0.);
        TPZManVector<STATE> val2(nstate,0.);
        auto BCond1 = matloc->CreateBC(matloc, Ebc1, 0, val1, val2);
        if (applyexact)
        {
            auto forcingfunction = [](const TPZVec<REAL>&x, TPZVec<STATE>&u){
                LaplaceExact.ForcingFunction()->Execute(x, u);
            };
            auto exactsol = [](const TPZVec<REAL>&x, TPZVec<STATE>&u,
                                TPZFMatrix<STATE>&du){
                LaplaceExact.Exact()->Execute(x, u, du);
            };
            BCond1->SetForcingFunctionBC(exactsol);
            matloc->SetForcingFunction(forcingfunction, porder);
        }
        auto BSkeleton = matloc->CreateBC(matloc, ESkeleton, 1, val1, val2);

        cmesh->InsertMaterialObject(matloc);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BSkeleton);
    }    
}

TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder, bool scalarproblem, bool useexact)
{
    bool elasticityproblem = !scalarproblem;
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
    
    //        OneQuad(gmesh);
    gengrid.Read(gmesh,EGroup);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc1);
    gengrid.SetBC(gmesh, 6, Ebc1);
    gengrid.SetBC(gmesh, 7, Ebc1);
    {
        TPZManVector<int64_t,2> nodeindex(1);
		int64_t index;
        nodeindex[0] = 0;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
        nodeindex[0] = nelx;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
        gmesh->BuildConnectivity();
    }
    
    std::map<int,int> matmap;
    matmap[EGroup] = 1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    build.StandardConfiguration();
    build.DivideSkeleton(nrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 0;
    if (elasticityproblem) {
        problemtype = 0;
    }
    else
    {
        problemtype = 1;
    }
    InsertMaterialObjects(SBFem,!elasticityproblem, useexact);
    
    build.BuildComputationMesh(*SBFem);
    
    if(1)
    {
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
}

TPZCompMesh *SetupSquareH1Mesh(int nelx, int porder, bool scalarproblem, bool useexact)
{
    bool elasticityproblem = !scalarproblem;
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
    
    //        OneQuad(gmesh);
    gengrid.Read(gmesh,Emat1);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
    {
        TPZManVector<int64_t,2> nodeindex(1);
		int64_t index;
        nodeindex[0] = 0;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
        nodeindex[0] = nelx;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
        gmesh->BuildConnectivity();
    }
    
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 0;
    if (elasticityproblem) {
        problemtype = 0;
    }
    else
    {
        problemtype = 1;
    }
    InsertMaterialObjects(SBFem,!elasticityproblem, useexact);

    SBFem->SetAllCreateFunctionsContinuous();
    SBFem->AutoBuild();
    
    if(1)
    {
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
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

void ElGroupEquations(TPZSBFemElementGroup *elgr, TPZVec<int64_t> &equations)
{
    equations.Resize(0, 0);
    TPZCompMesh *cmesh = elgr->Mesh();
    int nc = elgr->NConnects();
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = elgr->Connect(ic);
        int blsize = c.NDof();
        int64_t eqsize = equations.size();
        equations.Resize(eqsize+blsize, 0);
        int64_t seqnum = c.SequenceNumber();
        for (int idf = 0; idf<blsize; idf++) {
            equations[eqsize+idf] = cmesh->Block().Position(seqnum)+idf;
        }
    }
}

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupCrackedOneElement(int nrefskeleton, int porder, bool applyexact, bool elastic)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    REAL coor[][3] = {
        {0,0},
        {-1,0},
        {-1,-1},
        {1,-1},
        {1,1},
        {-1,1},
        {-1,0}
    };
    gmesh->NodeVec().Resize(7);
    for (int i=0; i<7; i++) {
        TPZManVector<REAL,3> co(3,0);
        co[0] = coor[i][0];
        co[1] = coor[i][1];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    {
        
        TPZManVector<int64_t,2> nodeindices(2);
        nodeindices[0] = 1;
        nodeindices[1] = 2;
        int64_t index;
        TPZManVector<int64_t,2> nodeindex(1);

        nodeindex[0] = 2;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
        nodeindex[0] = 3;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Emat1, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc1, index);
        for (int i=1; i<4; i++) {
            nodeindices[0] = i+1;
            nodeindices[1] = i+2;
            gmesh->CreateGeoElement(EOned, nodeindices, Emat2, index);
            gmesh->CreateGeoElement(EOned, nodeindices, Ebc2, index);
        }
        nodeindices[0] = 5;
        nodeindices[1] = 6;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat3, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc3, index);


        nodeindices[0] = 2;
        nodeindices[1] = 3;
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc4, index);
    }
    gmesh->BuildConnectivity();
    std::map<int,int> matidtranslation;
    matidtranslation[Emat1] = Emat1;
    matidtranslation[Emat2] = Emat2;
    matidtranslation[Emat3] = Emat3;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
    TPZManVector<int64_t,10> scalingcenters(1);
    scalingcenters[0] = 0;
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t,10> elementgroup(nel,-1);
    for (int i=0; i<nel; i+=2) {
        elementgroup[i] = 0;
    }
    build.SetPartitions(elementgroup, scalingcenters);
    std::set<int> matids;
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    matids.insert(Ebc4);
    build.DivideSkeleton(nrefskeleton,matids);
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    InsertMaterialObjects(cmesh, !elastic, true);
    TPZMaterial *mat = cmesh->FindMaterial(Emat1);
    TPZMaterial *mat2 = mat->NewMaterial();
    mat2->SetId(Emat2);
    cmesh->InsertMaterialObject(mat2);
    TPZMaterial *mat3 = mat->NewMaterial();
    mat3->SetId(Emat3);
    cmesh->InsertMaterialObject(mat3);
    build.BuildComputationalMeshFromSkeleton(*cmesh);
    {
        int64_t nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->Dimension() != gmesh->Dimension()-1) {
                continue;
            }
            if (gel->MaterialId() == Emat1 || gel->MaterialId() == Emat2 || gel->MaterialId() == Emat3) {
                gel->SetMaterialId(ESkeleton);
            }
        }
    }
    return cmesh;
}

TPZCompMesh *SetupCrackedTwoElements(int nrefskeleton, int porder, bool applyexact, bool elastic)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    REAL coor[][3] = {
        {-1, 0}, //0
        { 1, 0}, //1
        { 0, 0}, //2
        { 0,-2}, //3
        { 2,-2}, //4
        { 2, 2}, //5
        { 0, 2}, //6
        {-2, 2}, //7
        {-2,-2},  //8
        { 0, 0}  //9
    };
    gmesh->NodeVec().Resize(10);
    for (int i=0; i<7; i++) {
        TPZManVector<REAL,3> co(3,0);
        co[0] = coor[i][0];
        co[1] = coor[i][1];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    {
        // Element 1
        TPZManVector<int64_t,2> nodeindices(2);
        nodeindices[0] = 2;
        nodeindices[1] = 3;
        int64_t index;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat1, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc1, index);
        
        for (int i=0; i<3; i++) {
            nodeindices[0] = i+3;
            nodeindices[1] = i+4;
            gmesh->CreateGeoElement(EOned, nodeindices, Emat2, index);
            gmesh->CreateGeoElement(EOned, nodeindices, Ebc2, index);
        }

        nodeindices[0] = 6;
        nodeindices[1] = 9;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat3, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc3, index);

        // Element 2
        nodeindices[0] = 9;
        nodeindices[1] = 6;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat3, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc3, index);

        nodeindices[0] = 8;
        nodeindices[1] = 3;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat2, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc2, index);

        nodeindices[0] = 7;
        nodeindices[1] = 8;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat2, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc2, index);

        nodeindices[0] = 6;
        nodeindices[1] = 7;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat2, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc2, index);

        nodeindices[0] = 3;
        nodeindices[1] = 2;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat1, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc1, index);
        // nodeindices[0] = 2;
        // nodeindices[1] = 3;
        // gmesh->CreateGeoElement(EOned, nodeindices, Ebc4, index);
    }
    gmesh->BuildConnectivity();
    std::map<int,int> matidtranslation;
    matidtranslation[Emat1] = Emat1;
    matidtranslation[Emat2] = Emat2;
    matidtranslation[Emat3] = Emat3;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
    TPZManVector<int64_t,10> scalingcenters(2);
    scalingcenters[0] = 0;
    scalingcenters[1] = 1;
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t,10> elementgroup(nel,-1);
    for (int i=0; i<nel/2; i+=2) {
        elementgroup[i] = 0;
    }
    for (int i=nel/2; i<nel; i+=2) {
        elementgroup[i] = 1;
    }
    build.SetPartitions(elementgroup, scalingcenters);
    std::set<int> matids;
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    build.DivideSkeleton(nrefskeleton,matids);
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    InsertMaterialObjects(cmesh, !elastic, true);
    TPZMaterial *mat = cmesh->FindMaterial(Emat1);
    TPZMaterial *mat2 = mat->NewMaterial();
    mat2->SetId(Emat2);
    cmesh->InsertMaterialObject(mat2);
    TPZMaterial *mat3 = mat->NewMaterial();
    mat3->SetId(Emat3);
    cmesh->InsertMaterialObject(mat3);
    build.BuildComputationalMeshFromSkeleton(*cmesh);
    {
        int64_t nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->Dimension() != gmesh->Dimension()-1) {
                continue;
            }
            if (gel->MaterialId() == Emat1 || gel->MaterialId() == Emat2 || gel->MaterialId() == Emat3) {
                gel->SetMaterialId(ESkeleton);
            }
        }
    }

    if(1)
    {
        std::cout << "Plotting the geometric mesh\n";
        std::ofstream out("cmesh.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintCMeshVTK(cmesh, out,true);
    }
    return cmesh;
}

TPZGeoMesh *ReadUNSWQuadtreeMesh(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices){
    int maxvol = -1;
    
    std::ifstream file(filename);

    std::map<std::set<int64_t> , int64_t> midnode;
    std::string buf;
    getline(file,buf);
    if(!file) DebugStop();

    int64_t nnodes, nvolumes;
    file >> nnodes >> nvolumes;
    
    elpartition.Resize(nvolumes*6, -1);
    scalingcenterindices.Resize(nvolumes,0);
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(nnodes);
    
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3,0.);
        for (int i=0; i<2; i++) {
            file >> xco[i];
        }
        gmesh->NodeVec()[in].Initialize(xco, *gmesh);
    }
#ifdef PZDEBUG
    std::set<int64_t> badvolumes;
#endif
    int64_t nsurfaces;
    file >> nsurfaces;
    for (int64_t iv=0; iv<nsurfaces; iv++) {
#ifdef PZDEBUG
        std::map<std::set<int64_t>,int64_t> nodepairs;
#endif
        int elnnodes;
        file >> elnnodes;
        TPZManVector<int64_t,10> nodes(elnnodes);
        for (int i=0; i<elnnodes; i++) {
            file >> nodes[i];
            nodes[i]--;
        }
        
        if (elnnodes == 1)
        {
            int64_t index;
            MElementType eltype = EPoint;
            gmesh->CreateGeoElement(eltype, nodes, EGroup, index);
            elpartition[index] = iv;
            
        }
        else if (elnnodes == 2)
        {
            int64_t index;
            MElementType eltype = EOned;
            gmesh->CreateGeoElement(eltype, nodes, EGroup, index);
            elpartition[index] = iv;

        }
        else if (elnnodes == 3 || elnnodes == 4)
        {
            int64_t index;
            MElementType eltype = EQuadrilateral;
            if (elnnodes == 3) {
                eltype = ETriangle;
            }
            for (int i=0; i<4; i++){
                TPZVec<int64_t> nodeside(2);
                nodeside[0] = nodes[i];
                if (i==3){
                    nodeside[1] = nodes[0];
                } else{
                    nodeside[1] = nodes[i+1];
                }
                gmesh->CreateGeoElement(EOned, nodeside, ESkeleton, index);
                elpartition[index] = iv;
            }

            TPZManVector<REAL,3> midxco(3,0.);
            std::set<int64_t>  elnodes;
            for (int i=0; i<elnnodes; i++) {
                TPZManVector<REAL,3> x(3);
                gmesh->NodeVec()[nodes[i]].GetCoordinates(x);
                for(int j=0; j<3; j++) midxco[j] += x[j]/elnnodes;
            }
            int64_t midindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[midindex].Initialize(midxco, *gmesh);

            midnode[elnodes] = midindex;
            scalingcenterindices[iv] = midindex;
        }
        else if(elnnodes > 4)
        {
            std::set<int64_t>  elnodes;
            TPZManVector<REAL,3> midxco(3,0.);
            for (int i=0; i<elnnodes; i++) {
                elnodes.insert(nodes[i]);
                TPZManVector<REAL,3> x(3);
                gmesh->NodeVec()[nodes[i]].GetCoordinates(x);
                for(int j=0; j<3; j++) midxco[j] += x[j]/elnnodes;
            }
            int64_t midindex = -1;
            if (midnode.find(elnodes) == midnode.end()) {
                midindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[midindex].Initialize(midxco, *gmesh);
                midnode[elnodes] = midindex;
                TPZManVector<int64_t,10> nodeindices(1,midindex);
            }
            else
            {
                midindex = midnode[elnodes];
            }
            for (int triangle = 0; triangle <elnnodes; triangle++) {
                TPZManVector<int64_t,3> nodeindices(3);
                for (int in=0; in<2; in++) {
                    nodeindices[in] = nodes[(triangle+in)%elnnodes];
                }
                nodeindices[2] = midindex;
                int64_t index;
                TPZVec<int64_t> nodeside(2);
                nodeside[0] = nodeindices[0];
                nodeside[1] = nodeindices[1];
                gmesh->CreateGeoElement(EOned, nodeside, ESkeleton, index);
                elpartition[index] = iv;
                scalingcenterindices[iv] = midindex;
            }
        }
        else
        {
            DebugStop();
        }
        if (elpartition.size() < gmesh->NElements()+100) {
            elpartition.Resize(elpartition.size()*2, -1);
        }
    }
    TPZManVector<int64_t> matidelpartition(nvolumes);
    for (int64_t in=0; in<nvolumes; in++) {
        int64_t matid;
        file >> matid;
        matidelpartition[in] = matid;
    }
    for (int64_t i = 0; i < gmesh->NElements(); i++)
    {
        if (elpartition[i] == -1){
            continue;
        }
        int64_t matindex = elpartition[i];
        gmesh->Element(i)->SetMaterialId(matidelpartition[matindex]);
        
    }
    
    {
        std::ofstream mirror("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, mirror);
    }
    elpartition.Resize(gmesh->NElements(), -1);
    std::cout << "Building element connectivity\n";
    gmesh->BuildConnectivity();

    return gmesh;
}

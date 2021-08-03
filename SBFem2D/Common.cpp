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

#include <chrono>

TElasticity2DAnalytic ElastExact;
TLaplaceExample1 LaplaceExact;
TLaplaceExampleTimeDependent TimeLaplaceExact;

auto forcingfunctionelast = [](const TPZVec<REAL>&x, TPZVec<STATE>&u){
  ElastExact.ForcingFunction()->Execute(x, u);
};

auto forcingfunctionlaplace = [](const TPZVec<REAL>&x, TPZVec<STATE>&u){
  LaplaceExact.ForcingFunction()->Execute(x, u);
};


void SolveSist(TPZLinearAnalysis &an, TPZCompMesh *Cmesh, int numthreads)
{
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(Cmesh);
#endif
    
    strmat.SetNumThreads(numthreads);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    an.Run();
}

void InsertMaterialObjects(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact)
{
    int nstate = 1;
    int dim = 2;
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    const TPZManVector<double> val2(nstate, 0.);
    
    if (scalarproblem)
    {
        TPZMatPoisson<STATE> *matloc = new TPZMatPoisson<STATE>(Emat1, dim);
        cmesh->InsertMaterialObject(matloc);
        //TPZMatPoisson<STATE> *matloc2 = new TPZMatPoisson<STATE>(Emat2, dim);
        // matloc2->SetDimension(2);
        // cmesh->InsertMaterialObject(matloc2);
        // material = matloc2;
        constexpr int pOrder{2};
        matloc->SetForcingFunction(forcingfunctionlaplace, pOrder);
        
        auto BCond1 = matloc->CreateBC(matloc, Ebc1, 0, val1, val2);
        auto BCond2 = matloc->CreateBC(matloc, Ebc2, 0, val1, val2);
        auto BCond3 = matloc->CreateBC(matloc, Ebc3, 0, val1, val2);
        auto BCond4 = matloc->CreateBC(matloc, Ebc4, 0, val1, val2);
        
        BCond1->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        BCond2->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        BCond3->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        BCond4->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        
        val1.Zero();
        auto BCond5 = matloc->CreateBC(matloc, EBCPoint1, 0, val1, val2);
        BCond5->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        auto BCond6 = matloc->CreateBC(matloc, EBCPoint2, 0, val1, val2);
        BCond6->SetForcingFunctionBC(LaplaceExact.ExactSolution());
        
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
        cmesh->InsertMaterialObject(BCond4);
        cmesh->InsertMaterialObject(BCond5);
        cmesh->InsertMaterialObject(BCond6);
        
        auto BSkeleton = matloc->CreateBC(matloc, ESkeleton, 1, val1, val2);
        cmesh->InsertMaterialObject(BSkeleton);
    }
    else
    {
        TPZElasticity2D *matloc1 = new TPZElasticity2D(Emat1);
        TPZElasticity2D *matloc2 = new TPZElasticity2D(Emat2);
        nstate = 2;
        int porder = 4;
        if (applyexact)
        {
            matloc1->SetPlaneStress();
            matloc1->SetElasticity(ElastExact.gE, ElastExact.gPoisson);
            matloc2->SetPlaneStress();
            matloc2->SetElasticity(ElastExact.gE, ElastExact.gPoisson);
            constexpr int pOrder{2};
            matloc1->SetForcingFunction(forcingfunctionelast, pOrder);
            matloc2->SetForcingFunction(forcingfunctionelast, pOrder);
        }
        
    	auto BCond1 = matloc1->CreateBC(matloc1, Ebc1, 0, val1, val2);
    	auto BCond2 = matloc1->CreateBC(matloc1, Ebc2, 0, val1, val2);
    	auto BCond3 = matloc1->CreateBC(matloc1, Ebc3, 0, val1, val2);
    	auto BCond4 = matloc1->CreateBC(matloc1, Ebc4, 0, val1, val2);
    	
        constexpr int pOrder{2};
    	BCond1->SetForcingFunctionBC(ElastExact.ExactSolution());
    	BCond2->SetForcingFunctionBC(ElastExact.ExactSolution());
    	BCond3->SetForcingFunctionBC(ElastExact.ExactSolution());
    	BCond4->SetForcingFunctionBC(ElastExact.ExactSolution());
    	
    	val1.Zero();
    	auto BCond5 = matloc1->CreateBC(matloc1, EBCPoint1, 0, val1, val2);
    	BCond5->SetForcingFunctionBC(ElastExact.ExactSolution());
    	auto BCond6 = matloc1->CreateBC(matloc1, EBCPoint2, 0, val1, val2);
    	BCond6->SetForcingFunctionBC(ElastExact.ExactSolution());
    	
    	cmesh->InsertMaterialObject(matloc1);
    	cmesh->InsertMaterialObject(BCond1);
    	cmesh->InsertMaterialObject(BCond2);
    	cmesh->InsertMaterialObject(BCond3);
    	cmesh->InsertMaterialObject(BCond4);
    	cmesh->InsertMaterialObject(BCond5);
    	cmesh->InsertMaterialObject(BCond6);
    	
        int EPoly = 100;
        auto BPoly = matloc1->CreateBC(matloc1, EPoly, 1, val1, val2);
        cmesh->InsertMaterialObject(BPoly);
            
        auto BSkeleton = matloc1->CreateBC(matloc1, ESkeleton, 1, val1, val2);
        cmesh->InsertMaterialObject(BSkeleton);
    }
}

TPZAutoPointer<TPZGeoMesh> SetupGeom(int nelx)
{
    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh * gmesh = new TPZGeoMesh;

    // Setting boundary elements
    {
        gengrid.Read(gmesh, EGroup);
        gengrid.SetBC(gmesh, 4, Ebc1);
        gengrid.SetBC(gmesh, 5, Ebc2);
        gengrid.SetBC(gmesh, 6, Ebc3);
        gengrid.SetBC(gmesh, 7, Ebc4);
        TPZManVector<int64_t, 2> nodeindex(1);
        int64_t index;
        nodeindex[0] = 0;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
        nodeindex[0] = nelx;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
    }
    gmesh->BuildConnectivity();

    return gmesh;
}

TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder, bool scalarproblem, bool useexact)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = SetupGeom(nelx);

    // Defining configuration for SBFEM mesh
    std::map<int, int> matmap;
    matmap[EGroup] = Emat1;
    TPZBuildSBFem build(gmesh, ESkeleton, matmap);
    // this method will create center node elements and skeleton elements
    build.StandardConfiguration();
    build.DivideSkeleton(nrefskeleton);

    // Defining computational mesh and material data
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    InsertMaterialObjects(SBFem, scalarproblem, useexact);
    
    // Generating SBFEM mesh
    build.BuildComputationMesh(*SBFem);

    bool outputcmshgmsh = true;
    if (outputcmshgmsh)
    {
        OutputGmshCmsh(gmesh, SBFem);
    }

    return SBFem;
}

void OutputGmshCmsh(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmesh)
{
    std::ofstream outc("CMesh.txt");
    cmesh->Print(outc);
    std::ofstream outg("GMesh.txt");
    gmesh->Print(outg);
    std::ofstream out("Geometry.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintGMeshVTK(gmesh, out, true);
}

TPZCompMesh *SetupSquareH1Mesh(int nelx, int porder, bool scalarproblem, bool useexact)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = SetupGeom(nelx);

    /// put sbfem pyramids into the element groups
    TPZCompMesh *FEM = new TPZCompMesh(gmesh);
    FEM->SetDefaultOrder(porder);
    InsertMaterialObjects(FEM, scalarproblem, useexact);
    FEM->SetAllCreateFunctionsContinuous();
    FEM->AutoBuild();

    bool outputcmshgmsh = false;
    if (outputcmshgmsh)
    {
        OutputGmshCmsh(gmesh, FEM);
    }
    return FEM;
}

TPZCompMesh *ReadJSonFile(const std::string &filename, int numrefskeleton, int pOrder, REAL contrast)
{
    // read in json file
    std::ifstream myfile(filename);
    nlohmann::json json;
    myfile >> json;

    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    // Part coor
    std::vector<std::vector<double>> coor = json["coor"]; // "coor" are in 2d vector
    int nnodesTotal = coor.size();
    gmesh->NodeVec().Resize(nnodesTotal);
    std::cout << "Coordinates ( " << nnodesTotal << " in total )" << std::endl;
    for (int i = 0; i < nnodesTotal; i++)
    {
        TPZManVector<REAL, 3> co(3, 0.);
        for (int j = 0; j < 2; j++)
        {
            co[j] = coor[i][j];
        }
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }

    // Part elem
    std::vector<nlohmann::json> elem = json["elem"]; // "elem" are in 1d vector with json object
    int nElem = elem.size();
    TPZVec<int64_t> scalecenter(nElem, -1);
    for (int64_t el = 0; el < nElem; el++)
    {

        // sc idx
        int sc = elem[el]["sc"]; //sc idx
        scalecenter[el] = sc;
        // node idx
        std::vector<int> nodes = elem[el]["nodes"]; // nodes list
        int nnodes = nodes.size();
        TPZManVector<int64_t> nodeindices(nnodes, -1);
        for (int k = 0; k < nnodes; k++)
        {
            nodeindices[k] = nodes[k];
        }
        // mat idx
        int mat = elem[el]["mat"];
        if (mat == Emat1)
        {
            mat = EGroup;
        }
        if (mat == Emat2)
        {
            mat = EGroup + 1;
        }
        int64_t index;
        switch (nnodes)
        {
        case 1:
            gmesh->CreateGeoElement(EPoint, nodeindices, mat, index);
            break;
        case 2:
            gmesh->CreateGeoElement(EOned, nodeindices, mat, index);
            break;
        case 3:
            gmesh->CreateGeoElement(ETriangle, nodeindices, mat, index);
            break;
        case 4:
            gmesh->CreateGeoElement(EQuadrilateral, nodeindices, mat, index);
            break;
        default:
            DebugStop();
            break;
        }
        gmesh->BuildConnectivity();
    }
    if (0)
    {
        std::map<int, TPZStack<int>> elementset;
        for (int64_t el = 0; el < gmesh->NElements(); el++)
        {
            if (scalecenter[el] == -1)
            {
                continue;
            }
            elementset[scalecenter[el]].Push(el);
        }
        int materialindex = 100;
        for (std::map<int, TPZStack<int>>::iterator it = elementset.begin(); it != elementset.end(); it++)
        {
            int64_t nel = it->second.NElements();
            int64_t nodeindex = it->first;
            TPZManVector<REAL, 3> xcenter(3);
            gmesh->NodeVec()[nodeindex].GetCoordinates(xcenter);
            for (int64_t el = 0; el < nel; el++)
            {
                int64_t elindex = it->second[el];
                TPZGeoEl *gel = gmesh->Element(elindex);
                TPZManVector<REAL, 3> xi(2), xco(3);
                gel->CenterPoint(gel->NSides() - 1, xi);
                gel->X(xi, xco);
                int64_t newnode = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newnode].Initialize(xco, gmesh);
                TPZManVector<int64_t, 3> cornerindexes(2);
                cornerindexes[0] = nodeindex;
                cornerindexes[1] = newnode;
                int64_t index;
                gmesh->CreateGeoElement(EOned, cornerindexes, materialindex, index);
            }
            materialindex++;
        }
    }
    gmesh->BuildConnectivity();
    scalecenter.Resize(gmesh->NElements(), -1);
    std::map<int, int> matmap;
    matmap[EGroup] = Emat1;
    matmap[EGroup + 1] = Emat2;
    TPZBuildSBFem build(gmesh, ESkeleton, matmap);
    build.Configure(scalecenter);
    build.DivideSkeleton(numrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups

    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(pOrder);

    // problemtype - 1 laplace equation
    int problemtype = 0;

    bool applyexact = false;
    InsertMaterialObjects(SBFem, problemtype, applyexact);

    {
        auto mat = SBFem->FindMaterial(Emat1);
        REAL elast, poisson, lambda, G;
        {
            TPZElasticity2D *matelas = dynamic_cast<TPZElasticity2D *>(mat);
            elast = matelas->E();
            poisson = matelas->Nu();
            lambda = matelas->GetLambda(elast, poisson);
            G = matelas->GetMU(elast, poisson);
        }

        auto mat2 = mat->NewMaterial();
        mat2->SetId(Emat2);
        TPZElasticity2D *matelas2 = dynamic_cast<TPZElasticity2D *>(mat2);
        REAL elast2 = elast * contrast;
        matelas2->SetElasticity(elast2, poisson);
        SBFem->InsertMaterialObject(mat2);
        TPZFNMatrix<4, STATE> val1(2, 2, 0.);
        const TPZVec<double> val2(1, 0.);
        // zero neumann at the bottom
        auto bnd = matelas2->CreateBC(matelas2, -1, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        // zero neumann at the top
        bnd = matelas2->CreateBC(matelas2, -3, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val2[0] = 1.;
        bnd = matelas2->CreateBC(matelas2, -2, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val2[0] = 0;
        val1(1, 1) = 1.;
        val1(0, 0) = 1.;
        // remove rigid body modes
        bnd = matelas2->CreateBC(matelas2, -5, 2, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val1(0, 0) = 1.;
        val1(1, 1) = 0.;
        bnd = matelas2->CreateBC(matelas2, -6, 2, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val1.Zero();
        // traction to the left
        val2[0] = -1.;
        bnd = matelas2->CreateBC(matelas2, -4, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
    }

    build.BuildComputationMesh(*SBFem);
    if (0)
    {
        int64_t nel = SBFem->NElements();
        for (int64_t el = 0; el < nel; el++)
        {
            TPZCompEl *cel = SBFem->Element(el);
            TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if (elgr)
            {
                TPZElementMatrixT<STATE> ek, ef;
                elgr->CalcStiff(ek, ef);
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (intel && intel->NConnects() == 3)
            {
                TPZGeoEl *ref = intel->Reference();
                TPZManVector<REAL, 3> co(3), val(2, 0.);
                ref->NodePtr(0)->GetCoordinates(co);
                val[0] = co[0] * 0.01;
                int64_t seqnum = intel->Connect(0).SequenceNumber();
                SBFem->Block().at(seqnum, 0, 0, val[0]);
                ref->NodePtr(1)->GetCoordinates(co);
                val[0] = co[0] * 0.01;
                seqnum = intel->Connect(1).SequenceNumber();
                SBFem->Block().at(seqnum, 0, 0, val[0]);
            }
        }
        SBFem->LoadSolution(SBFem->Solution());
    }
    SBFem->LoadReferences();

    {
        std::ofstream out("JSonGeometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out, true);
        std::ofstream outg("JSonGeometry.txt");
        SBFem->Reference()->Print(outg);
        std::ofstream outc("JSonComp.txt");
        SBFem->Print(outc);
    }
    return SBFem;
}

void ElGroupEquations(TPZSBFemElementGroup *elgr, TPZVec<int64_t> &equations)
{
    equations.Resize(0, 0);
    TPZCompMesh *cmesh = elgr->Mesh();
    int nc = elgr->NConnects();
    for (int ic = 0; ic < nc; ic++)
    {
        TPZConnect &c = elgr->Connect(ic);
        int blsize = c.NDof();
        int64_t eqsize = equations.size();
        equations.Resize(eqsize + blsize, 0);
        int64_t seqnum = c.SequenceNumber();
        for (int idf = 0; idf < blsize; idf++)
        {
            equations[eqsize + idf] = cmesh->Block().Position(seqnum) + idf;
        }
    }
}
/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZSBFemVolume *celv)
{
    TPZGeoEl *gel = celv->Reference();
    int dim = gel->Dimension();
    int nstate = celv->Connect(0).NState();
    TPZCompMesh *cmesh = celv->Mesh();
    int volside = gel->NSides() - 1;
    TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cmesh->Element(celv->ElementGroupIndex()));
    TPZManVector<int64_t> globeq;
    ElGroupEquations(elgr, globeq);
    TPZIntPoints *intpoints = gel->CreateSideIntegrationRule(volside, 3);
    cmesh->Solution().Zero();
    for (int ip = 0; ip < intpoints->NPoints(); ip++)
    {
        TPZManVector<REAL, 3> xi(gel->Dimension(), 0.);
        TPZFNMatrix<32, REAL> phi, dphidxi;
        REAL weight;
        intpoints->Point(ip, xi, weight);
        celv->Shape(xi, phi, dphidxi);
        int64_t neq = globeq.size();
        for (int64_t eq = 0; eq < neq; eq++)
        {
            int64_t globindex = globeq[eq];
            cmesh->Solution().Zero();
            TPZFMatrix<STATE> solcmesh = (cmesh->Solution());
            solcmesh(globindex, 0) = 1.;
            cmesh->LoadSolution(solcmesh);

            TPZManVector<STATE> sol;
            celv->Solution(xi, 0, sol);
            REAL diffphi = 0.;
            for (int istate = 0; istate < nstate; istate++)
            {
                diffphi += (sol[istate] - phi(eq * nstate + istate)) * (sol[istate] - phi(eq * nstate + istate));
            }
            diffphi = sqrt(diffphi);
            if (diffphi > 1.e-8)
            {
                std::cout << "Wrong shape function diffphi = " << diffphi << "\n";
            }
        }
    }
    delete intpoints;
}

/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr)
        {
            TPZVec<TPZCompEl *> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            for (int iv = 0; iv < nvol; iv++)
            {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                VerifyShapeFunctionIntegrity(elvol);
            }
        }
    }
}

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupCrackedOneElement(int nrefskeleton, int porder, bool applyexact, bool elastic)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    REAL coor[][3] = {
        {0, 0},
        {-1, 0},
        {-1, -1},
        {1, -1},
        {1, 1},
        {-1, 1},
        {-1, 0}};
    gmesh->NodeVec().Resize(7);
    for (int i = 0; i < 7; i++)
    {
        TPZManVector<REAL, 3> co(3, 0);
        co[0] = coor[i][0];
        co[1] = coor[i][1];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    {
        TPZManVector<int64_t, 2> nodeindices(2);
        nodeindices[0] = 1;
        nodeindices[1] = 2;
        int64_t index;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat1, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc1, index);
        for (int i = 1; i < 4; i++)
        {
            nodeindices[0] = i + 1;
            nodeindices[1] = i + 2;
            gmesh->CreateGeoElement(EOned, nodeindices, Emat2, index);
            gmesh->CreateGeoElement(EOned, nodeindices, Ebc2, index);
        }
        nodeindices[0] = 5;
        nodeindices[1] = 6;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat3, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc3, index);
        nodeindices[0] = 6;
        nodeindices[1] = 0;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat4, index);

        nodeindices[0] = 6;
        gmesh->CreateGeoElement(EPoint, nodeindices, Ebc4, index);

        nodeindices[0] = 0;
        nodeindices[1] = 1;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat4, index);

        nodeindices[0] = 1;
        gmesh->CreateGeoElement(EPoint, nodeindices, Ebc4, index);

        nodeindices[0] = 3;
        gmesh->CreateGeoElement(EPoint, nodeindices, EBCPoint1, index);

        nodeindices[0] = 4;
        gmesh->CreateGeoElement(EPoint, nodeindices, EBCPoint2, index);
    }
    gmesh->BuildConnectivity();
    std::map<int, int> matidtranslation;
    matidtranslation[Emat1] = Emat1;
    matidtranslation[Emat2] = Emat2;
    matidtranslation[Emat3] = Emat3;
    matidtranslation[Emat4] = Emat3;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
    TPZManVector<int64_t, 10> scalingcenters(1);
    scalingcenters[0] = 0;
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t, 10> elementgroup(nel, -1);
    for (int i = 0; i < nel-2; i += 2)
    {
        elementgroup[i] = 0;
    }
    build.SetPartitions(elementgroup, scalingcenters);
    std::set<int> matids;
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    build.DivideSkeleton(nrefskeleton);
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

    std::set<int> volmatids,boundmatids;
    volmatids.insert(Emat1);
    volmatids.insert(Emat2);
    volmatids.insert(Emat3);
    boundmatids.insert(Ebc1);
    boundmatids.insert(Ebc2);
    boundmatids.insert(Ebc3);
    boundmatids.insert(Ebc4);
    boundmatids.insert(ESkeleton);
    build.DivideSkeleton(nrefskeleton);

    build.BuildComputationMesh(*cmesh,volmatids,boundmatids);
    // build.BuildComputationalMeshFromSkeleton(*cmesh);
    {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->Dimension() != gmesh->Dimension() - 1)
            {
                continue;
            }
            if (gel->MaterialId() == Emat1 || gel->MaterialId() == Emat2 || gel->MaterialId() == Emat3)
            {
                gel->SetMaterialId(ESkeleton);
            }
        }
    }
    return cmesh;
}

using namespace std;
void PostProcessing(TPZLinearAnalysis & Analysis, const std::string &filename, bool scalarproblem, int numthreads, int POrder, int nelxcount, int irefskeleton)
{
    // Generating Paraview file
#ifdef _AUTODIFF
    if(scalarproblem)
    {
        Analysis.SetExact(Laplace_exact);
    }
    else
    {
        Analysis.SetExact(Elasticity_exact);
    }
#endif

    
    std::stringstream soutvtk("RegularSolution.vtk");
    if (0)
    {
        if(scalarproblem)
        {
            TPZStack<std::string> vecnames,scalnames;
            scalnames.Push("State");
            Analysis.DefineGraphMesh(2, scalnames, vecnames, soutvtk.str());
            int res = POrder+1;
            if (res >5) {
                res = 5;
            }
            Analysis.PostProcess(res);
        }
        else
        {
            TPZStack<std::string> vecnames,scalnames;
            vecnames.Push("Displacement");
            scalnames.Push("SigmaX");
            scalnames.Push("SigmaY");
            scalnames.Push("TauXY");
            scalnames.Push("EpsX");
            scalnames.Push("EpsY");
            scalnames.Push("EpsXY");
            Analysis.DefineGraphMesh(2, scalnames, vecnames,soutvtk.str());
            Analysis.PostProcess(3);
        }
    }

    // Computing errors
    auto start = chrono::steady_clock::now();
    std::cout << "Compute errors\n";
    int64_t neq = Analysis.Mesh()->NEquations();
    TPZManVector<REAL,10> errors(3,0.);
    // Analysis.SetThreadsForError(numthreads);
    Analysis.PostProcessError(errors);
    
    std::stringstream sout("RegularSolution.txt");
    
    std::ofstream results(sout.str(),std::ios::app);
    results.precision(15);
    int nelx = 1 << (nelxcount-1);
    results << "(* nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
    TPZFMatrix<double> errmat(1,3);
    for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
    std::stringstream varname;
    varname << "Errmat[[" << nelxcount << "]][[" << irefskeleton+1 << "]][[" << POrder << "]] = (1/1000000)*";
    errmat.Print(varname.str().c_str(),results,EMathematicaInput);
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time to compute error (miliseconds): " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << "\n";
}

void PrintEigval(TPZLinearAnalysis Analysis, std::string &filename)
{
    std::stringstream sout(filename);
    sout << ".txt";
    std::ofstream results(sout.str(),std::ios::app);
    results.precision(15);

    TPZCompMesh * SBFem  = Analysis.Mesh();
    TPZSBFemElementGroup *celgrp = 0;
    int64_t nel = Analysis.Mesh()->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZSBFemElementGroup *cel = dynamic_cast<TPZSBFemElementGroup *>(SBFem->Element(el));
        if(cel)
        {
            celgrp = cel;
            break;
        }
    }
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
TPZGeoMesh *ReadUNSWQuadtreeMesh(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices){
    int maxvol = -1;
    
    std::ifstream file(filename);

    std::map<set<int64_t> , int64_t> midnode;
    std::string buf;
    std::getline(file,buf);
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
        std::map<set<int64_t>,int64_t> nodepairs;
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
        matidelpartition[in] = 100;
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

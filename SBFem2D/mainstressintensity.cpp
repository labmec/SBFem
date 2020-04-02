#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
// #include "Common3D.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"
#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzaxestools.h"

#include "pzgeoelbc.h"
#include "pzbndcond.h"
#include "pzelast3d.h"
#include "TPZVTKGeoMesh.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"

#include "tpzintpoints.h"
#include "pzgeoelrefless.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

// #ifdef _AUTODIFF
// TElasticity2DAnalytic ElastExactLower;
// TElasticity2DAnalytic ElastExactUpper;
// #endif

void IntegrateDirect(TPZCompMesh *cmesh);

TPZGeoMesh *ReadUNSWSBGeoFile(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = false;
    bool hasexact = false;

    int numrefskeleton = 1;
    int maxporder = 8;
    int counter = 1;
    int numproblem = 0;
#ifdef _AUTODIFF
    ElastExact.fProblemType = TElasticity2DAnalytic::ESquareRoot;
#endif
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            int polynomialorder = 6;
            TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(BodyLoads, polynomialorder);
            TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
            TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;
#ifdef _AUTODIFF
            ElastExact.fE = 400000;
            ElastExact.fPoisson = 0.3;
            ElastExact.fPlaneStress = 0;
            ElastExact = ElastExact;
            ElastExactLower.fProblemType = TElasticity2DAnalytic::ESquareRootLower;
#endif
            bool elastic = !scalarproblem;
            TPZCompMesh *SBFem = SetupCrackedTwoElements(irefskeleton, POrder, hasexact, elastic);

            // std::string filename("../../../src/SBFem/SBFem2D/OneCrack.txt");
            // std::cout << "Reading " << filename << std::endl;
            // TPZManVector<int64_t,1000> elpartitions;
            // TPZVec<int64_t> scalingcenterindices;
            // TPZGeoMesh *gmesh =ReadUNSWSBGeoFile(filename, elpartitions, scalingcenterindices);
            if(1)
            {
                std::cout << "Plotting the geometric mesh\n";
                std::ofstream out("gmesh.vtk");
                TPZVTKGeoMesh vtk;
                vtk.PrintGMeshVTK(SBFem->Reference(), out,true);
            }
            // int64_t nel = gmesh.NElements();
            // for (int64_t el=0; el<nel; el++) {
            //     TPZGeoEl *gel = gmesh.Element(el);
            //     if(gel->Type() == EPoint)
            //     {
            //         gel->SetMaterialId(Ebc3);
            //         continue;
            //     }
            //     if (gel->Type() != EOned) {
            //         continue;
            //     }
            //     TPZManVector<REAL,3> xco1(3), xco2(3);
            //     gel->Node(0).GetCoordinates(xco1);
            //     gel->Node(1).GetCoordinates(xco2);
            //     if(fabs(xco1[1]-0) < 1.e-6 && fabs(xco2[1]-0) < 1.e-6)
            //     {
            //         TPZGeoElBC(gel,2,Ebc2);
            //     }
            // }
            // TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
            // SBFem->SetDefaultOrder(POrder);
            
#ifdef _AUTODIFF
            // TPZCompMesh *SBFem = SetupCrackedOneElement(irefskeleton, POrder, hasexact, elastic);

            // if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Emat1);
                mat->SetForcingFunction(autodummy);
                // mat->SetForcingFunction(ElastExactLower.ForcingFunction());
            }
            // if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Emat2);
                mat->SetForcingFunction(autodummy);
                // mat->SetForcingFunction(ElastExact.ForcingFunction());
            }
            // if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Emat3);
                mat->SetForcingFunction(autodummy);
                // mat->SetForcingFunction(ElastExactUpper.ForcingFunction());
            }
            // if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(1);
                // bc->SetForcingFunction(0,autodummy);
                // mat->SetForcingFunction(BodyLoads);

            }
            // if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(0);
                // bc->SetForcingFunction(0,autodummy);
                // mat->SetForcingFunction(ElastExact.TensorFunction());
            }
            // if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(1);
                // bc->SetForcingFunction(0,autodummy);
                // mat->SetForcingFunction(ElastExactUpper.TensorFunction());

                // TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Elasticity_exact_upper, 6);
                // TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
                // bc->SetForcingFunction(0,autodummy);
            }
#endif
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                std::ofstream gout("gmesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(SBFem->Reference(), gout,true);
                std::stringstream sout;
                SBFem->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
                std::ofstream gout("gmesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(SBFem->Reference(), gout,true);
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(Analysis, SBFem);
            
            
            
            
            std::cout << "Post processing\n";
            //        ElasticAnalysis->Solution().Print("Solution");
            //        mphysics->Solution().Print("expandec");
#ifdef _AUTODIFF
            Analysis->SetExact(Elasticity_exact);
#endif
            //                ElasticAnalysis->SetExact(Singular_exact);
            
            TPZManVector<REAL> errors(3,0.);
            
            int64_t neq = SBFem->Solution().Rows();
            
            if(!scalarproblem)
            {
                std::stringstream filename;
                filename << "SquareRootOneElement_NR_" << irefskeleton << "_P_" << POrder << ".vtk";
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                vecnames.Push("Displacement");
                scalnames.Push("SigmaX");
                scalnames.Push("SigmaY");
                scalnames.Push("TauXY");
                scalnames.Push("EpsX");
                scalnames.Push("EpsY");
                scalnames.Push("EpsXY");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, filename.str());
                Analysis->PostProcess(3);
            }

            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            if(hasexact)
            {
            
                std::cout << "Compute errors\n";
                
                Analysis->PostProcessError(errors);
                
    //                VerifyShapeFunctionIntegrity(Analysis->Mesh());
                
    //                IntegrateDirect(Analysis->Mesh());
                
                std::stringstream sout;
                sout << "../CrackRestrainedShape";
                if (scalarproblem) {
                    sout << "Scalar.txt";
                }
                else
                    sout << "Elastic2D.txt";
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(*  numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "Errmat[[" << irefskeleton+1 << "]][[" << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            }
            if(0)
            {
                std::cout << "Plotting shape functions\n";
                int numshape = 25;
                if (numshape > SBFem->NEquations()) {
                    numshape = SBFem->NEquations();
                }
                TPZVec<int64_t> eqindex(numshape);
                for (int i=0; i<numshape; i++) {
                    eqindex[i] = i;
                }
                Analysis->ShowShape("OneElementCracked.vtk", eqindex);
            }
            
            delete Analysis;
            delete SBFem;
            //                exit(-1);
        }
        //            exit(-1);
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}






void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

void IntegrateDirect(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr) {
            TPZStack<TPZCompEl *,5> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            TPZElementMatrix ekvol, efvol, ekgrp, efgrp;
            elgr->CalcStiff(ekgrp, efgrp);
            for (int iv=0; iv<nvol; iv++) {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                TPZElementMatrix ek,ef;
                elvol->CalcStiff(ek, ef);
                if (iv==0) {
                    ekvol = ek;
                    efvol = ef;
                }
                else
                {
                    ekvol.fMat += ek.fMat;
                    efvol.fMat += ef.fMat;
                }
            }
            ekgrp.fMat.Print("EKGRP = ",std::cout,EMathematicaInput);
            ekvol.fMat.Print("EKVOL = ",std::cout,EMathematicaInput);
            break;
        }
    }

    
}

using namespace std;
/// Read a UNSWSBFem file
TPZGeoMesh *ReadUNSWSBGeoFile(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices)
{
    
    int maxvol = -1;
    
    std::ifstream file(filename);
    
    map<set<int64_t> , int64_t> midnode;
    string buf;
    getline(file,buf);
    if(!file) DebugStop();
    int64_t nnodes, nvolumes;
    file >> nnodes >> nvolumes;
    elpartition.Resize(nvolumes*6, -1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        for (int i=0; i<3; i++) {
            file >> xco[i];
        }
        gmesh->NodeVec()[in].Initialize(xco, *gmesh);
    }
#ifdef PZDEBUG
    std::set<int64_t> badvolumes;
#endif
    int64_t nothing;
    file >> nothing; // now it represents the number of faces
    for (int64_t iv=0; iv<nvolumes; iv++)
    {
#ifdef PZDEBUG
        map<set<int64_t>,int64_t> nodepairs;
#endif
        int nfaces;
        file >> nfaces;
        for (int face = 0; face < nfaces; face++) {
            int elnnodes;
            file >> elnnodes;

            TPZManVector<int64_t,10> nodes(elnnodes);
            for (int i=0; i<elnnodes; i++) {
                file >> nodes[i];
                nodes[i]--;
#ifdef PZDEBUG
                if (i>0) {
                    set<int64_t> edge;
                    edge.insert(nodes[i-1]);
                    edge.insert(nodes[i]);
                    nodepairs[edge]++;
                }
                if (i==elnnodes-1) {
                    set<int64_t> edge;
                    edge.insert(nodes[0]);
                    edge.insert(nodes[i]);
                    nodepairs[edge]++;
                }
#endif
            }
            
            // tototototo
            if (maxvol != -1 && iv >= maxvol) {
                continue;
            }
            if (elnnodes == 1)
            {
                int64_t index;
                MElementType eltype = EPoint;
                gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                elpartition[index] = iv;
                
            }
            else if (elnnodes == 2)
            {
                int64_t index;
                MElementType eltype = EOned;
                gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                elpartition[index] = iv;
                
            }
            else if (elnnodes == 3 || elnnodes == 4)
            {
                int64_t index;
                MElementType eltype = ETriangle;
                if (elnnodes == 4) {
                    eltype = EQuadrilateral;
                }
                gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                elpartition[index] = iv;
            }
            else if(elnnodes == 8)
            {
                int64_t index;
                // new TPZGeoElRefPattern<pzgeom::TPZQuadraticQuad> (nodes, ESkeleton, *gmesh,  index);
                elpartition[index] = iv;

            }
            else if(elnnodes > 4)
            {
                set<int64_t>  elnodes;
                TPZManVector<REAL,3> midxco(3,0.);
                for (int i=0; i<elnnodes; i++) {
                    elnodes.insert(nodes[i]);
                    TPZManVector<REAL,3> x(3);
                    gmesh->NodeVec()[nodes[i]].GetCoordinates(x);
                    //                    std::cout << "x " << x << endl;
                    for(int j=0; j<3; j++) midxco[j] += x[j]/elnnodes;
                }
                int64_t midindex = -1;
                if (midnode.find(elnodes) == midnode.end()) {
                    midindex = gmesh->NodeVec().AllocateNewElement();
                    gmesh->NodeVec()[midindex].Initialize(midxco, *gmesh);
                    midnode[elnodes] = midindex;
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
                    gmesh->CreateGeoElement(ETriangle, nodeindices, ESkeleton, index);
                    elpartition[index] = iv;
                }
            }
            else
            {
                DebugStop();
            }
        }
#ifdef PZDEBUG
        bool suspicious = false;
        for (auto it = nodepairs.begin(); it != nodepairs.end(); it++) {
            if(it->second != 2) suspicious = true;
        }
        if (suspicious == true) {
            std::cout << "volume " << iv << " has no closure\n";
            badvolumes.insert(iv);
        }
#endif
        if (elpartition.size() < gmesh->NElements()+100) {
            elpartition.Resize(elpartition.size()*2, -1);
        }
    }
    // totototototo
    if (maxvol != -1) {
        nvolumes = maxvol;
    }
    int64_t nmidnodes = midnode.size();
    gmesh->NodeVec().Resize(nvolumes+nmidnodes+nnodes);
    scalingcenterindices.Resize(nvolumes, -1);
    for (int64_t in=0; in<nvolumes; in++) {
        TPZManVector<REAL,3> xco(3);
        for (int i=0; i<3; i++) {
            file >> xco[i];
        }
        gmesh->NodeVec()[nnodes+nmidnodes+in].Initialize(xco, *gmesh);
        scalingcenterindices[in] = nnodes+nmidnodes+in;
    }
    {
        ofstream mirror("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, mirror);
    }
#ifdef PZDEBUG
    if (badvolumes.size()) {
        int64_t nel = gmesh->NElements();
        TPZManVector<REAL> elval(nel,0);
        for (int64_t el=0; el<nel; el++) {
            if (badvolumes.find(elpartition[el]) != badvolumes.end()) {
                elval[el] = 10.;
            }
        }
        {
            ofstream badel("gmesh_bad.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, badel, elval);
        }
    }
#endif
    elpartition.Resize(gmesh->NElements(), -1);
    std::cout << "Building element connectivity\n";
    gmesh->BuildConnectivity();
    return gmesh;
}

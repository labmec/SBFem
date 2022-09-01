#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"

#include "pzgeoelbc.h"

#include "Elasticity/TPZElasticity2D.h"
#include "TPZBndCond.h"

#include "TPZVTKGeoMesh.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"

#include "tpzintpoints.h"
#include "pzgeoelrefless.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"

void AddBoundaryElements(TPZAutoPointer<TPZGeoMesh> & gmesh);

void AdjustElementOrientation(TPZGeoMesh &gmesh, TPZVec<int64_t> &elpartitions, TPZVec<int64_t> &scalingcenterindices);

void SolveSistPolygons(TPZLinearAnalysis *an, TPZCompMesh *Cmesh, int numthreads);

int main(int argc, char *argv[])
{
#if PZ_LOG
    TPZLogger::InitializePZLOG();
#endif // PZ_LOG
    // Initial data
    int minnelxcount = 1, maxnelxcount = 2;
    int minporder = 1, maxporder = 7;
    int numthreads = 0;
    bool scalarproblem = true; // false for elasticity 2D problems
    bool usesbfem = true; // false for FEM simulations
    bool useexact = true;
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
    int counter = 1;
    
    std::string filename;
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int nelx = 1; nelx <=4 ; nelx++){
            if(nelx == 1) filename = "polygon1.txt";
            else if (nelx == 2) filename = "polygon2.txt";
            else if (nelx == 3) filename = "polygon3.txt";
            else if (nelx == 4) filename = "polygon4.txt";
            else if (nelx == 5) filename = "polygon5.txt";
            else DebugStop();
#ifdef MACOSX
            filename = "../"+filename;
#endif
            std::string vtkfilename;
            std::string vtkfilegeom;
            std::string vtkfilecmesh;
            std::string rootname;
            std::string boundaryname;
            {
                int pos = filename.find(".txt");
                std::string truncate = filename;
                truncate.erase(pos);
                rootname = truncate;
                std::stringstream sout;
                sout << truncate << "_p" << POrder << ".vtk";
                vtkfilename = sout.str();
                std::stringstream boundstr;
                boundstr << truncate << "_boundary";
                boundaryname = boundstr.str();
                std::stringstream vtkgeom;
                vtkgeom << truncate << "_geom.vtk";
                vtkfilegeom = vtkgeom.str();
                std::stringstream vtkcmesh;
                vtkcmesh << truncate << "_cmesh.vtk";
                vtkfilecmesh = vtkcmesh.str();
            }
            
            TPZVec<int64_t> elpartition;
            TPZVec<int64_t> scalingcenterindices;
            TPZAutoPointer<TPZGeoMesh> gmesh = ReadUNSWQuadtreeMesh(filename, elpartition, scalingcenterindices);
            
            AdjustElementOrientation(*gmesh, elpartition, scalingcenterindices);
            AddBoundaryElements(gmesh);
            elpartition.Resize(gmesh->NElements(), -1);
            scalingcenterindices.Resize(gmesh->NElements(), -1);

            std::cout << "Building computational mesh\n";
            std::map<int,int> matmap;
            matmap[ESkeleton] = Emat1;
            int EPoly = 100;
            matmap[EPoly] = Emat2;

            TPZBuildSBFem build(gmesh,ESkeleton,matmap);
            build.SetPartitions(elpartition, scalingcenterindices);
            build.DivideSkeleton(0);

            TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
            SBFem->SetDefaultOrder(POrder);
            InsertMaterialObjects(SBFem, scalarproblem, useexact);
            
            // elemento 194 esta com problema...
            build.BuildComputationalMeshFromSkeleton(*SBFem);
            if(1) {
                std::cout << "Plotting the geometric mesh\n";
                std::ofstream outtxt("gmesh.txt");
                gmesh->Print(outtxt);
                std::ofstream out(vtkfilegeom);
                TPZVTKGeoMesh vtk;
                vtk.PrintGMeshVTK(gmesh, out,true);
            }
            if(1) {
                std::cout << "Plotting the computational mesh\n";
                std::ofstream out(vtkfilecmesh);
                TPZVTKGeoMesh vtk;
                vtk.PrintCMeshVTK(SBFem, out,true);
                std::ofstream outc("cmesh.txt");
                SBFem->Print(outc);
            }

            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);

            SolveSistPolygons(Analysis, SBFem, numthreads);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;

            std::cout << "Post processing\n";
            std::string sout;
            sout.append("../PolygonsSolution_SBFem");
            
            PostProcessing(*Analysis, sout, scalarproblem, numthreads, POrder, nelx, 0);
            
            delete SBFem;
            delete Analysis;
        }

    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void AddBoundaryElements(TPZAutoPointer<TPZGeoMesh> & gmesh)
{
    std::set<int64_t> setbottom,setright,settop,setleft;
    int64_t nnodes = gmesh->NNodes();
    int dim = gmesh->Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[in].GetCoordinates(xco);
        if (fabs(xco[1]+1.) < 1.e-3)
        {
            setbottom.insert(in);
        }
        if (fabs(xco[0]-1.) < 1.e-3)
        {
            setright.insert(in);
        }
        if (fabs(xco[1]-1.) < 1.e-3)
        {
            settop.insert(in);
        }
        if (fabs(xco[0]+1.) < 1.e-3)
        {
            setleft.insert(in);
        }
    }
    int64_t nelem = gmesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            int nsidenodes = gel->NSideNodes(is);
            int nfoundbottom = 0;
            int nfoundright = 0;
            int nfoundtop = 0;
            int nfoundleft = 0;
            for (int in=0; in<nsidenodes; in++) {
                int64_t nodeindex = gel->SideNodeIndex(is, in);
                if (setbottom.find(nodeindex) != setbottom.end()) {
                    nfoundbottom++;
                }
                if (setright.find(nodeindex) != setright.end()) {
                    nfoundright++;
                }
                if (settop.find(nodeindex) != settop.end()) {
                    nfoundtop++;
                }
                if (setleft.find(nodeindex) != setleft.end()) {
                    nfoundleft++;
                }
            }
            if (nfoundbottom == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc1);
            }
            if (nfoundright == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc2);
            }
            if (nfoundtop == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc3);
            }
            if (nfoundleft == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc4);
            }
            else
            {
                TPZGeoElSide gelside(gel,is);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if (neighbour == gelside) {
                    int EPoly = 100;
                    TPZGeoElBC(gelside,EPoly);
                }
            }
        }
    }
}

// reorient the elements such that the boundary elements have outward pointing normals
void AdjustElementOrientation(TPZGeoMesh &gmesh, TPZVec<int64_t> &elpartitions, TPZVec<int64_t> &scalingcenterindices)
{
    std::cout << "Adjusting element orientations\n";
    int64_t nel = gmesh.NElements();
    TPZManVector<int64_t,8> nodeindices(8,0);
    int64_t index;
    int matid = 1;
    int64_t numelgood = 0;
    int64_t numelbad = 0;
    TPZGeoEl *hexa = new TPZGeoElRefLess<pzgeom::TPZGeoQuad>(nodeindices, matid,gmesh, index);
    TPZGeoEl *prism = new TPZGeoElRefLess<pzgeom::TPZGeoTriangle>(nodeindices, matid,gmesh, index);
    for (int64_t el=0; el<nel; el++) {
        if (elpartitions[el] == -1) {
            continue;
        }
        int64_t elpartition = elpartitions[el];
        TPZGeoEl *gel = gmesh.Element(el);
        if (gel->Type() == EQuadrilateral) {
            TPZVec<TPZGeoNode> nodes(3);
            nodes[0] = gmesh.NodeVec()[gel->NodeIndex(0)];
            nodes[1] = gmesh.NodeVec()[gel->NodeIndex(1)];
            nodes[2] = gmesh.NodeVec()[gel->NodeIndex(2)];

            REAL axby = nodes[0].Coord(0)*nodes[1].Coord(1);
            REAL aycx = nodes[0].Coord(1)*nodes[2].Coord(0);
            REAL bxcy = nodes[1].Coord(0)*nodes[2].Coord(1);
            REAL bycx = nodes[1].Coord(1)*nodes[2].Coord(0);
            REAL aybx = nodes[0].Coord(1)*nodes[1].Coord(0);
            REAL axcy = nodes[0].Coord(0)*nodes[2].Coord(1);
            
            REAL area = axby + aycx + bxcy - bycx - aybx - axcy;
            
            // REAL area = bax*cay - cax*bay;
            if (area<0.) {
                numelbad++;
                gel->RemoveConnectivities();
                for(int i=0; i<4; i++)
                {
                    gel->SetNodeIndex(i, hexa->NodeIndex(3-i));
                }
            }
            else
            {
                numelgood++;
            }
                
        }
        else if(gel->Type() == ETriangle)
        {

            TPZVec<TPZGeoNode> nodes(3);
            nodes[0] = gmesh.NodeVec()[gel->NodeIndex(0)];
            nodes[1] = gmesh.NodeVec()[gel->NodeIndex(1)];
            nodes[2] = gmesh.NodeVec()[gel->NodeIndex(2)];

            REAL axby = nodes[0].Coord(0)*nodes[1].Coord(1);
            REAL aycx = nodes[0].Coord(1)*nodes[2].Coord(0);
            REAL bxcy = nodes[1].Coord(0)*nodes[2].Coord(1);
            REAL bycx = nodes[1].Coord(1)*nodes[2].Coord(0);
            REAL aybx = nodes[0].Coord(1)*nodes[1].Coord(0);
            REAL axcy = nodes[0].Coord(0)*nodes[2].Coord(1);
            
            REAL area = axby + aycx + bxcy - bycx - aybx - axcy;
            if (area<0.) {
                numelbad++;
                gel->RemoveConnectivities();
                for(int i=0; i<3; i++)
                {
                    gel->SetNodeIndex(i, prism->NodeIndex(2-i));
                }
            }
            else{
                numelgood++;
            }

        }
    }
    delete hexa;
    delete prism;
    std::cout << "Number of elements with original orientation " << numelgood << " number of elements inverted " <<
        numelbad << std::endl;
    gmesh.BuildConnectivity();
}

void SolveSistPolygons(TPZLinearAnalysis *an, TPZCompMesh *Cmesh, int numthreads)
{
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(Cmesh);
    strmat.SetNumThreads(numthreads);
    an->SetStructuralMatrix(strmat);
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    an->Assemble();
//    try {
//        an->Assemble();
//    } catch (...) {
//        exit(-1);
//    }
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
}

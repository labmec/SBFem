#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common3D.h"
#include "TPZBuildSBFem.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzaxestools.h"

#include "pzgeoelbc.h"
#include "TPZBndCond.h"
#include "Elasticity/TPZElasticity3D.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZVTKGeoMesh.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"

#include "tpzintpoints.h"
#include "pzgeoelrefless.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"

#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"

#include <algorithm>

TPZGeoMesh * GenerateSBElementsFromGmsh(TPZGeoMesh *gmsh,TPZManVector<int64_t,1000> &elpartition,TPZManVector<int64_t,1000> &scalingcentreindexes);

void InsertMaterialObjects3DPolygons(TPZCompMesh * SBFem);

void AddBoundaryElements(TPZGeoMesh &gmesh);

int main(int argc, char *argv[])
{
    ExactLaplace.fExact = TLaplaceExample1::EHarmonic2;

    int minrefskeleton = 1, maxrefskeleton = 4;
    int minporder = 1, maxporder = 4;
    int counter = 1;
    int numthreads = 32;
    for ( int POrder = minporder; POrder < maxporder; POrder++)
    {
        for ( int nref = minrefskeleton; nref < maxrefskeleton; nref++)
        {
            std::string filename;
            switch (nref)
            {
            case 1:
                filename = "n8-id1-2.msh";
                break;
            case 2:
                filename = "n64-id1-3.msh";
                break;
            case 3:
                filename = "n512-id1-3.msh";
                break;
            default:
                break;
            }
            
            std::string vtkfilename;
            std::string rootname;
            std::string boundaryname;
            std::string vtkfilegeom;
            {
                int pos = filename.find(".msh");
                std::string truncate = filename;
                truncate.erase(pos);
                rootname = truncate;
                std::stringstream sout;
                sout << truncate << "_t" << numthreads << "_p" << POrder << "_href" << 1 << ".vtk";
                vtkfilename = sout.str();
                std::stringstream boundstr;
                boundstr << truncate << "_boundary";
                boundaryname = boundstr.str();
                std::stringstream vtkgeom;
                vtkgeom << truncate << "_geom.vtk";
                vtkfilegeom = vtkgeom.str();
            }

            std::cout << "Reading " << filename << std::endl;
            TPZManVector<int64_t,1000> elpartitions;
            TPZManVector<int64_t,1000> scalingcenterindices;

            TPZGmshReader gmshrdr;
            TPZGeoMesh *gmsh = new TPZGeoMesh;
            gmsh = gmshrdr.GeometricGmshMesh(filename,gmsh);

            std::ofstream file("polygons.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmsh, file);

            TPZAutoPointer<TPZGeoMesh> gmsh2 = GenerateSBElementsFromGmsh(gmsh,elpartitions,scalingcenterindices);
            AddBoundaryElements(*gmsh2);
            elpartitions.Resize(gmsh2->NElements(), -1);

            std::cout << "Building computational mesh \n";
            std::map<int,int> matidtranslation;
            matidtranslation[ESkeleton] = Emat1;
            TPZBuildSBFem build(gmsh2, ESkeleton, matidtranslation);
            build.SetPartitions(elpartitions, scalingcenterindices);
            build.DivideSkeleton(0);

            TPZCompMesh *SBFem = new TPZCompMesh(gmsh2);
            SBFem->SetDefaultOrder(POrder);
            InsertMaterialObjects3DPolygons(SBFem);
            build.BuildComputationalMeshFromSkeleton(*SBFem);

            if(0){
                TPZVTKGeoMesh vtk;
                std::ofstream outvtk("GMeshVTK.vtk");
                TPZGeoMesh * ref = SBFem->Reference();
                vtk.PrintGMeshVTK(ref,outvtk);
                std::ofstream outtxt("GMeshVTK.txt");
                SBFem->Reference()->Print(outtxt);
            }
            
            // Visualization of computational meshes
            std::cout << "Analysis\n";
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(Analysis, SBFem, numthreads);

            std::cout << "Plotting\n";
            
            Analysis->SetExact(Laplace_exact);

            int64_t neq = SBFem->Solution().Rows();
            if(0)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                {
                    scalnames.Push("State");
                }
                Analysis->DefineGraphMesh(3, scalnames, vecnames, vtkfilename);
                Analysis->PostProcess(3);
            }

            std::cout << "Post processing\n";

            TPZManVector<REAL> errors(3,0.);
            Analysis->SetThreadsForError(numthreads);
            Analysis->PostProcessError(errors);
                
            std::stringstream sout;
            {
                sout << "../Scalar3DSolutionPolygons.txt";
            }
            
            {
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(* POrder " << POrder << " neq " << neq << " *)" << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "ErrmatPoly[[" << nref << "," << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            }

            delete Analysis;
            delete SBFem;
        }
    }

    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

TPZGeoMesh * GenerateSBElementsFromGmsh(TPZGeoMesh *gmsh,TPZManVector<int64_t,1000> &elpartition,TPZManVector<int64_t,1000> &scalingcentreindexes)
{
    TPZGeoMesh * sbgmsh = new TPZGeoMesh;
    sbgmsh->SetDimension(3);
    auto nnodes = gmsh->NNodes();
    sbgmsh->NodeVec() = gmsh->NodeVec();
    
    int64_t nelem = gmsh->NElements();
    elpartition.Resize(nelem,-1);
    TPZFNMatrix<1000,REAL> elscalcentre(nelem,3,0);
    int64_t count = 0;

    for (int64_t iel = 0; iel < nelem; iel++)
    {
        TPZManVector<REAL,3> xsc(3,0);
        TPZGeoEl * gel = gmsh->Element(iel);
        if (!gel)
        {
            continue;
        }
        int64_t matid = gel->MaterialId();
        int nsides = gel->NSides();
        for (int isides = 0; isides < nsides; isides++)
        {
            TPZGeoElSide gelside = gel->Neighbour(isides);
            if (gelside.Dimension() != 2)
            {
                continue;
            }
            int matid2 = gelside.Element()->MaterialId();
            if (matid == matid2 && gelside != gelside.Neighbour())
            {
                continue;
            }
            else
            {
                int64_t index = count;
                TPZManVector<REAL,3> xco(3,0);
                auto nsidenodes = gelside.NSideNodes();
                TPZManVector<int64_t,3> cornerindexes(nsidenodes,0);
                for (int i = 0; i < nsidenodes; i++)
                {
                    cornerindexes[i] = gel->SideNodeIndex(isides, i);
                    elscalcentre(index,0) += sbgmsh->NodeVec()[cornerindexes[i]].Coord(0);
                    elscalcentre(index,1) += sbgmsh->NodeVec()[cornerindexes[i]].Coord(1);
                    elscalcentre(index,2) += sbgmsh->NodeVec()[cornerindexes[i]].Coord(2);
                }
                elscalcentre(index,0) = elscalcentre(index,0)/3;
                elscalcentre(index,1) = elscalcentre(index,1)/3;
                elscalcentre(index,2) = elscalcentre(index,2)/3;
                if (nsidenodes == 3)
                {
                    sbgmsh->CreateGeoElement(ETriangle, cornerindexes, ESkeleton, index);
                } else if(nsidenodes == 4)
                {
                    sbgmsh->CreateGeoElement(EQuadrilateral, cornerindexes, ESkeleton, index);
                }
                elpartition[index] = matid;
                count++;
            }
        }
    }

    nelem = sbgmsh->NElements();
    auto numpart = std::max_element(elpartition.begin(),elpartition.end());
    int64_t n = *numpart+1;
    TPZFNMatrix<100,REAL> scalingcentreco(n,3,0);
    TPZManVector<REAL,100> countpart(n,0);
    for (int64_t iel = 0; iel < nelem; iel++)
    {
        int64_t idpart = elpartition[iel];
        scalingcentreco(idpart,0) += elscalcentre(iel,0);
        scalingcentreco(idpart,1) += elscalcentre(iel,1);
        scalingcentreco(idpart,2) += elscalcentre(iel,2);
        countpart[idpart] += 1;
    }
    sbgmsh->NodeVec().Resize(nnodes+n);
    scalingcentreindexes.Resize(elpartition.size(),-1);
    for (auto i = 0; i < n; i++)
    {
        TPZManVector<REAL,3> xco(3,0);
        xco[0] = scalingcentreco(i,0)/countpart[i];
        xco[1] = scalingcentreco(i,1)/countpart[i];
        xco[2] = scalingcentreco(i,2)/countpart[i];
        // auto index = sbgmsh->NodeVec().AllocateNewElement();
        sbgmsh->NodeVec()[i+nnodes].Initialize(i+nnodes,xco, *sbgmsh);
        countpart[i] = i+nnodes;
        scalingcentreindexes[i] = i+nnodes;
    }
    sbgmsh->BuildConnectivity();
    return sbgmsh;
}

void InsertMaterialObjects3DPolygons(TPZCompMesh * cmesh){

    // Getting mesh dimension
    TPZMatPoisson<STATE> *matloc = new TPZMatPoisson<STATE>(Emat1, 3);
    int nstate = 1;
    cmesh->InsertMaterialObject(matloc);

    TPZFMatrix<STATE> val1(nstate,nstate,0.);
    TPZManVector<STATE> val2(nstate,0.);

    auto BCond1 = matloc->CreateBC(matloc,Ebc1,0, val1, val2);
    BCond1->SetForcingFunctionBC(ExactLaplace.ExactSolution());
    
    auto BCond2 = matloc->CreateBC(matloc,Ebc2,0, val1, val2);
    BCond2->SetForcingFunctionBC(ExactLaplace.ExactSolution());
    
    auto BCond3 = matloc->CreateBC(matloc,Ebc3, 0, val1, val2);
    BCond3->SetForcingFunctionBC(ExactLaplace.ExactSolution());
    
    auto BCond4 = matloc->CreateBC(matloc, Ebc4, 0, val1, val2);
    BCond4->SetForcingFunctionBC(ExactLaplace.ExactSolution());
    
    auto BSkeleton = matloc->CreateBC(matloc, ESkeleton,1, val1, val2);

    cmesh->InsertMaterialObject(BSkeleton);

    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
}

void AddBoundaryElements(TPZGeoMesh &gmesh)
{
    std::set<int64_t> setbottom, settop, setright, setleft, setfront, setback;
    int64_t nnodes = gmesh.NNodes();
    int dim = gmesh.Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        if (xco[1] < -1.e-3) {
            setbottom.insert(in);
        }
    }
    int64_t nelem = gmesh.NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh.Element(el);
        if (gel->Dimension() != dim-1) {
            DebugStop();
        }
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            int nsidenodes = gel->NSideNodes(is);
            int nfoundbottom = 0;
            for (int in=0; in<nsidenodes; in++) {
                int64_t nodeindex = gel->SideNodeIndex(is, in);
                if (setbottom.find(nodeindex) != setbottom.end()) {
                    nfoundbottom++;
                }
            }
            if (nfoundbottom == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc1);
            }
            else
            {
                TPZGeoElSide gelside(gel,is);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if (neighbour == gelside) {
                    TPZGeoElBC(gelside,Ebc2);
                }
            }
        }
    }
}

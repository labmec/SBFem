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
#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void AddBoundaryElementsCook(TPZGeoMesh &gmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 0;
    int maxrefskeleton = 1;
    int minporder = 1;
    int maxporder = 6;
    int counter = 1;
    int numthreads = 20;
    ExactElast.fE = 1000;
    ExactElast.fPoisson = 0.33;
    
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            
//            std::string filename("../CooksMembrane_poly_16_1_1.txt");
//            std::string filename("CooksMembrane_sbfemesh_64_2_1.txt");
//            std::string filename("../CooksMembrane_sbfemesh_128_4_1.txt");
            std::string filename("CooksMembrane_sbfemesh_32_2_1.txt");
            std::string vtkfilename;

            {
                int pos = filename.find(".txt");
                std::string truncate = filename;
                truncate.erase(pos);
                std::stringstream sout;
                sout << truncate << "_t" << numthreads << "_p" << POrder << "_href" << irefskeleton << ".vtk";
                vtkfilename = sout.str();
                
            }

            TPZManVector<int64_t,1000> elpartitions;
            TPZVec<int64_t> scalingcenterindices;
            TPZAutoPointer<TPZGeoMesh> gmesh = ReadUNSWSBGeoFile(filename, elpartitions, scalingcenterindices);
            AddBoundaryElementsCook(gmesh);
            elpartitions.Resize(gmesh->NElements(), -1);

            std::map<int,int> matidtranslation;
            matidtranslation[ESkeleton] = Emat1;
            TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
            build.SetPartitions(elpartitions, scalingcenterindices);
            TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
            SBFem->SetDefaultOrder(POrder);
            bool scalarproblem = false;
            InsertMaterialObjects3D(SBFem, scalarproblem);
            {
                TPZElasticity3D *mat = dynamic_cast<TPZElasticity3D *>(SBFem->FindMaterial(Emat1));
                mat->SetMaterialDataHook(1000., 0.49999);
            }
            {
                TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc2));
                bnd->SetType(0);
                TPZManVector<STATE> val2(3,0.);
                bnd->SetVal2(val2);
            }
            {
                TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *>(SBFem->FindMaterial(Ebc1));
                bnd->SetType(1);
                TPZManVector<STATE> val2(3,0.);
                val2[1] = 10.;
                bnd->SetVal2(val2);
            }

            build.BuildComputationalMeshFromSkeleton(*SBFem);

            TPZVTKGeoMesh vtk;
            std::ofstream out("CMeshVTKCook.txt");
            SBFem->Print(out);
            std::ofstream outvtk("CMeshVTKCook.vtk");
            vtk.PrintCMeshVTK(SBFem,outvtk);
            
            int64_t nelx = SBFem->NElements();
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                std::stringstream sout;
                SBFem->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            if(1)
            {
                std::ofstream outg("GMesh3D.txt");
                gmesh->Print(outg);
                std::ofstream out("Geometry3D.vtk");
                TPZVTKGeoMesh vtk;
                vtk.PrintGMeshVTK(gmesh, out,true);
            }

            std::cout << "nelx = " << nelx << std::endl;
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            
		    SolveSist(Analysis, SBFem, numthreads);
            
            std::cout << "Post processing\n";
            
            int64_t neq = SBFem->Solution().Rows();
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                vecnames.Push("Displacement");
                scalnames.Push("StressX");
                scalnames.Push("StressY");
                scalnames.Push("StressZ");
                Analysis->DefineGraphMesh(3, scalnames, vecnames, vtkfilename);
                Analysis->PostProcess(1);
            }
            
            delete Analysis;
            delete SBFem;
        }
    }
    
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void AddBoundaryElementsCook(TPZGeoMesh &gmesh)
{
    std::set<int64_t> leftset, rightset;
    int64_t nnodes = gmesh.NNodes();
    int dim = gmesh.Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        if (abs(xco[0]) < 1.e-3) {
            leftset.insert(in);
        }
        if (abs(xco[0]-48.) < 1.e-3) {
            rightset.insert(in);
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
            int nfoundleft = 0;
            int nfoundright = 0;
            for (int in=0; in<nsidenodes; in++) {
                int64_t nodeindex = gel->SideNodeIndex(is, in);
                if (leftset.find(nodeindex) != leftset.end()) {
                    nfoundleft++;
                }
                if (rightset.find(nodeindex) != rightset.end()) {
                    nfoundright++;
                }
            }
            if (nfoundright == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc1);
            }
            else if (nfoundleft == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc2);
            }
            else
            {
                TPZGeoElSide gelside(gel,is);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if (neighbour == gelside) {
                    TPZGeoElBC(gelside,Ebc3);
                }
            }
        }
    }
}
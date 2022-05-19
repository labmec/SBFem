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

void AddBoundaryElementsSphere(TPZGeoMesh &gmesh);

void SubstituteBoundaryConditionsSphere(TPZCompMesh &cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 0;
    int maxrefskeleton = 2;
    int minporder = 1;
    int maxporder = 3;
    int counter = 1;
    int numthreads = 32;
    ExactElast.fE = 200;
    ExactElast.fPoisson = 0.3;
    ExactElast.fProblemType = TElasticity3DAnalytic::ESphere;
    
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            for (int nref = 0; nref < 3; nref++)
            {
                std::string filename;
                switch (nref)
                {
                case 0:
                    filename = "spheres_10_50_sbfemesh_32_8_1.txt";
                    break;
                case 1:
                    filename = "spheres_10_50_sbfemesh_64_8_1.txt";
                    break;
                case 2:
                    filename = "spheres_10_50_sbfemesh_128_8_1.txt";
                    break;
                default:
                    std::cout << "Could not read the file\n";
                    DebugStop();
                    break;
                }
                std::string vtkfilename;
                std::string vtkgeofilename;
                std::string rootname;

                {
                    int pos = filename.find(".txt");
                    std::string truncate = filename;
                    truncate.erase(pos);
                    rootname = truncate;
                    {
                        std::stringstream sout;
                        sout << truncate << "_t" << numthreads << "_p" << POrder << "_href" << irefskeleton << ".vtk";
                        vtkfilename = sout.str();
                    }
                    {
                        std::stringstream sout;
                        sout << truncate <<  "_geo.vtk";
                        vtkgeofilename = sout.str();
                    }
                }

                TPZManVector<int64_t,1000> elpartitions;
                TPZVec<int64_t> scalingcenterindices;
                TPZAutoPointer<TPZGeoMesh> gmesh =ReadUNSWSBGeoFile(filename, elpartitions, scalingcenterindices);
                AddBoundaryElementsSphere(gmesh);
                elpartitions.Resize(gmesh->NElements(), -1);
                std::cout << "Done reading the file\n";
                
                std::map<int,int> matidtranslation;
                matidtranslation[ESkeleton] = Emat1;
                TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
                build.SetPartitions(elpartitions, scalingcenterindices);
                build.DivideSkeleton(irefskeleton);
                
                std::cout << "Creating the computational mesh\n";
                TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
                SBFem->SetDefaultOrder(POrder);
                bool scalarproblem = false;
                InsertMaterialObjects3D(SBFem, scalarproblem);
                SubstituteBoundaryConditionsSphere(*SBFem);
                build.BuildComputationalMeshFromSkeleton(*SBFem);
                
                int64_t nelx = SBFem->NElements();
                if(1)
                {
                    std::ofstream out(vtkgeofilename);
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

                int64_t neq = SBFem->Solution().Rows();
                
                if(0)
                {
                    std::cout << "Plotting\n";
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    vecnames.Push("State");
                    scalnames.Push("StressX");
                    scalnames.Push("StressY");
                    scalnames.Push("StressZ");
                    Analysis->DefineGraphMesh(3, scalnames, vecnames, vtkfilename);
                    Analysis->PostProcess(1);
                }
                
                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                std::cout << "Post processing\n";

                TPZManVector<REAL> errors(3,0.);
                Analysis->SetExact(Elasticity_exact);
                Analysis->SetThreadsForError(numthreads);
                Analysis->PostProcessError(errors);
                
                std::stringstream sout("spheres_Error.txt");
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(* nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "Errmat[[" << nref+1 << "," << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                
                delete Analysis;
                delete SBFem;
            }
        }
    }

    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void AddBoundaryElementsSphere(TPZGeoMesh &gmesh)
{
    std::set<int64_t> xset, yset, zset, inner, outer;
    REAL inner_radius = 10.;
    REAL outer_radius = 50.;
    int64_t nnodes = gmesh.NNodes();
    int dim = gmesh.Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        for (int i=0; i<3; i++) {
            xco[i] -= 100.;
        }
        gmesh.NodeVec()[in].SetCoord(xco);
    }
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        REAL radius = sqrt(xco[0]*xco[0]+xco[1]*xco[1]+xco[2]*xco[2]);
        if (abs(xco[0]) < 2.5e-2) {
            xset.insert(in);
            xco[0] = 0.;
            gmesh.NodeVec()[in].SetCoord(xco);
        }
        if (abs(xco[1]) < 2.5e-2) {
            yset.insert(in);
            xco[1] = 0.;
            gmesh.NodeVec()[in].SetCoord(xco);
        }
        if (abs(xco[2]) < 2.5e-2) {
            xco[2] = 0.;
            gmesh.NodeVec()[in].SetCoord(xco);
            zset.insert(in);
        }
        
        if (abs(radius-inner_radius) < 1.e-1) {
            REAL radius = sqrt(xco[0]*xco[0]+xco[1]*xco[1]+xco[2]*xco[2]);
            for (int i=0; i<3; i++) {
                xco[i] *= inner_radius/radius;
            }
            gmesh.NodeVec()[in].SetCoord(xco);
            inner.insert(in);
        }
        if (abs(radius-outer_radius) < 1.e-1) {
            REAL radius = sqrt(xco[0]*xco[0]+xco[1]*xco[1]+xco[2]*xco[2]);
            for (int i=0; i<3; i++) {
                xco[i] *=  outer_radius/radius;
            }
            gmesh.NodeVec()[in].SetCoord(xco);
            outer.insert(in);
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
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int nsidenodes = gel->NSideNodes(is);
            int nfoundx = 0;
            int nfoundy = 0;
            int nfoundz = 0;
            int nfoundinner = 0;
            int nfoundouter = 0;
            for (int in=0; in<nsidenodes; in++) {
                int64_t nodeindex = gel->SideNodeIndex(is, in);
                if (xset.find(nodeindex) != xset.end()) {
                    nfoundx++;
                }
                if (yset.find(nodeindex) != yset.end()) {
                    nfoundy++;
                }
                if (zset.find(nodeindex) != zset.end()) {
                    nfoundz++;
                }
                if (inner.find(nodeindex) != inner.end()) {
                    nfoundinner++;
                }
                if (outer.find(nodeindex) != outer.end()) {
                    nfoundouter++;
                }
            }
            if (nfoundx == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc1);
            }
            else if (nfoundy == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc2);
            }
            else if (nfoundz == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc3);
            }
            else if (nfoundinner == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc4);
            }
            else if (nfoundouter == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc5);
            }
            else if(gelside == neighbour)
            {
                for (int in = 0; in < nsidenodes; in++) {
                    int64_t index = gel->SideNodeIndex(is, in);
                    TPZManVector<REAL,3> xco(3);
                    gmesh.NodeVec()[index].GetCoordinates(xco);
                    REAL radius = sqrt(xco[0]*xco[0]+xco[1]*xco[1]+xco[2]*xco[2]);

                    std::cout << "in " << in << " xco " << xco << " radius " << radius <<  std::endl;
                }
                DebugStop();
            }
        }
    }
}

void SubstituteBoundaryConditionsSphere(TPZCompMesh &cmesh)
{
//    Ebc1 -> x
//    Ebc2 -> y
//    Ebc3 -> z
//    Ebc4 -> inner
//    Ebc5 -> outer
    {
        auto bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc1));
        bc->SetType(2);
        TPZFNMatrix<9,STATE> val1(3,3,0.);
        TPZManVector<REAL> val2(3,0.);
        val1(0,0) = 1.e12;
        bc->SetVal1(val1);
        bc->SetVal2(val2);
    }

    {
        auto bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc2));
        bc->SetType(2);
        TPZFNMatrix<9,STATE> val1(3,3,0.);
        TPZManVector<STATE> val2(3,0.);
        val1(1,1) = 1.e12;
        bc->SetVal1(val1);
        bc->SetVal2(val2);
    }
    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc3));
        bc->SetType(2);
        TPZFNMatrix<9,STATE> val1(3,3,0.);
        TPZManVector<STATE> val2(3,0.);
        val1(2,2) = 1.e12;
        bc->SetVal1(val1);
        bc->SetVal2(val2);
    }
    auto exactsol = [](const TPZVec<REAL>&x, TPZVec<STATE>&u,
                               TPZFMatrix<STATE>&du){
        ExactElast.Exact()->Execute(x, u, du);
    };
    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc4));
        bc->SetType(0);
        bc->SetForcingFunctionBC(exactsol,2);
    }
    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc5));
        bc->SetType(0);
        bc->SetForcingFunctionBC(exactsol,2);
    }

}

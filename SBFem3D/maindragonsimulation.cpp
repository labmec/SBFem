#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common3D.h"
#include "TPZBuildSBFem.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzaxestools.h"

#include "pzgeoelbc.h"
#include "TPZBndCondT.h"
#include "Elasticity/TPZElasticity3D.h"
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

// reorient the elements such that the boundary elements have outward pointing normals
void AdjustElementOrientation(TPZGeoMesh &gmesh, TPZVec<int64_t> &elpartitions, TPZVec<int64_t> &scalingcenterindices);

void AddBoundaryElements(TPZGeoMesh &gmesh, int boundarymatid);


void AddBoundaryElementsDragon(TPZGeoMesh &gmesh, int boundarymatid, int bottommatid);

void SubstituteBoundaryConditionsDragon(TPZCompMesh &cmesh);


void ComputeLoadVector(TPZCompMesh &cmesh, TPZFMatrix<STATE> &rhs);

void SolveSistDragon(TPZLinearAnalysis *an, TPZCompMesh *Cmesh, int numthreads);

// boundary group group index of each boundary element
void BuildBoundaryGroups(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup);

// generate a plot file for each boundary group
void PlotBoundaryGroups(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup, const std::string &rootname);

// output the volume within each boundary
void IntegrateVolumes(TPZGeoMesh &gmesh, TPZVec<int> &boundarygroup);

// print boundary group neighbours
void PrintBoundaryGroupNeighbourPartitions(TPZGeoMesh &gmesh, TPZVec<int> &boundarygroup, TPZVec<int64_t> &elpartitions
                                           , TPZVec<int64_t> &scalingcenterindices);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 0;
    int maxrefskeleton = 1;
    int minporder = 1;
    int maxporder = 3;
    int counter = 1;
    int numthreads = 32;
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            
//            std::string filename("../spheres_10_50_sbfemesh_128_8_1.txt");
//            std::string filename("../dragon_sbfemesh_256.txt");
//            std::string filename("../sphinx_sbfemesh_512.txt");
//            std::string filename("../bell_sbfemesh_512.txt");
            // std::string filename("dragon_sbfemesh_256.txt");
            std::string filename("dragon_remesh_sbfemesh_128.txt");
//            std::string filename("../dolphin_sbfemesh_128.txt");
//            std::string filename("../spheres_10_50_sbfemesh_64_8_1.txt");
            std::string vtkfilename;
            std::string vtkfilegeom;
            std::string rootname;
            std::string boundaryname;

            {
                int pos = filename.find(".txt");
                std::string truncate = filename;
                truncate.erase(pos);
                rootname = truncate;
                std::stringstream sout;
                sout << truncate << "_t" << numthreads << "_p" << POrder << "_href" << irefskeleton << ".vtk";
                vtkfilename = sout.str();
                std::stringstream boundstr;
                boundstr << truncate << "_boundary";
                boundaryname = boundstr.str();
                std::stringstream vtkgeom;
                vtkgeom << truncate << "_geom.vtk";
                vtkfilegeom = vtkgeom.str();
            }

            TPZManVector<int64_t,1000> elpartitions;
            TPZVec<int64_t> scalingcenterindices;
            std::cout << "Reading " << filename << std::endl;
            TPZAutoPointer<TPZGeoMesh> gmesh =ReadUNSWSBGeoFile(filename, elpartitions, scalingcenterindices);
            
            AdjustElementOrientation(gmesh, elpartitions, scalingcenterindices);
            
            std::cout << "Adding boundary conditions\n";
            AddBoundaryElements(gmesh,Ebc3);
            
            AddBoundaryElementsDragon(gmesh, Ebc1, Ebc3);
            elpartitions.Resize(gmesh->NElements(), -1);

            // change if you want to check the mesh
            if(0)
            {
                std::cout << "Checking the mesh\n";
                TPZVec<int> boundarygroups;
                BuildBoundaryGroups(gmesh, Ebc3, boundarygroups);
                // print boundary group neighbours
                PrintBoundaryGroupNeighbourPartitions(gmesh, boundarygroups, elpartitions
                                                      , scalingcenterindices);

                IntegrateVolumes(gmesh, boundarygroups);
                PlotBoundaryGroups(gmesh, Ebc3, boundarygroups, boundaryname);
                elpartitions.Resize(gmesh->NElements(), -1);
            }
            if(0)
            {
                std::cout << "Plotting the geometric mesh\n";
                std::ofstream out(vtkfilegeom);
                TPZVTKGeoMesh vtk;
                vtk.PrintGMeshVTK(gmesh, out,true);
            }
            std::cout << "Building the computational mesh\n";
            std::map<int,int> matidtranslation;
            matidtranslation[ESkeleton] = Emat1;
            TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
            build.SetPartitions(elpartitions, scalingcenterindices);
            build.DivideSkeleton(irefskeleton);
            TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
            SBFem->SetDefaultOrder(POrder);
            bool scalarproblem = false;
            InsertMaterialObjects3D(SBFem, scalarproblem);
            SubstituteBoundaryConditionsDragon(*SBFem);
            build.BuildComputationalMeshFromSkeleton(*SBFem);
            
            int64_t nelx = SBFem->NElements();
            if(1)
            {
                std::cout << "Plotting the geometric mesh\n";
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
            SolveSistDragon(Analysis, SBFem, numthreads);
            
            std::cout << "Post processing\n";
            
            if(1)
            {
                std::cout << "Plotting\n";
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

int64_t SBFemGroup(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *grp = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(grp) return el;
    }
    return -1;
}

void AddBoundaryElementsDragon(TPZGeoMesh &gmesh, int boundarymatid, int bottommatid)
{
    std::set<int64_t> bottomset;
    int64_t nnodes = gmesh.NNodes();
    int dim = gmesh.Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
//        if (xco[1] < 57) {
//        if (xco[1] < 65) {
        if (xco[1] < 61) {
            bottomset.insert(in);
        }
    }
    int64_t nelem = gmesh.NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh.Element(el);
        if (gel->Dimension() != dim-1) {
            DebugStop();
        }
        if (gel->MaterialId() != boundarymatid) {
            continue;
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
                if (bottomset.find(nodeindex) != bottomset.end()) {
                    nfoundbottom++;
                }
            }
            if (nfoundbottom == nsidenodes) {
                gel->SetMaterialId(bottommatid);
            }
        }
    }
}

void AddBoundaryElements(TPZGeoMesh &gmesh, int boundarymatid)
{
    std::set<int64_t> bottomset;
    int64_t nnodes = gmesh.NNodes();
    int dim = gmesh.Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        //        if (xco[1] < 57) {
        //        if (xco[1] < 65) {
        if (xco[1] < 61) {
            bottomset.insert(in);
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
            if (neighbour == gelside) {
                TPZGeoElBC(gelside,boundarymatid);
            }
        }
    }
}

void SubstituteBoundaryConditionsDragon(TPZCompMesh &cmesh)
{
//    Ebc1 -> x
//    Ebc2 -> y
//    Ebc3 -> z
//    Ebc4 -> inner
//    Ebc5 -> outer
    
    {
        TPZElasticity3D *elast = dynamic_cast<TPZElasticity3D *>(cmesh.FindMaterial(Emat1));
        elast->SetMaterialDataHook(30., 0.2);
    }
    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc1));
        bc->SetType(1);
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
        bc->SetVal1(val1);
    }

    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc2));
        bc->SetType(1);
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
        bc->SetVal1(val1);
    }
    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc3));
        bc->SetType(0);
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
        bc->SetVal1(val1);
    }
    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc4));
        bc->SetType(1);
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
        bc->SetVal1(val1);
    }
    {
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(Ebc5));
        bc->SetType(1);
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
        bc->SetVal1(val1);
    }

}

void CornerEquations(TPZSBFemElementGroup *elgr, TPZVec<int64_t> &indices)
{
    TPZVec<TPZCompEl *> elvol;
    TPZCompMesh *cmesh = elgr->Mesh();
    elvol = elgr->GetElGroup();
    int nvol = elvol.size();
    std::set<int64_t> nodeindices;
    for (int iv=0; iv<nvol; iv++) {
        TPZCompEl *cel = elvol[iv];
        TPZSBFemVolume *sbvol = dynamic_cast<TPZSBFemVolume *>(cel);
        int64_t skeleton = sbvol->SkeletonIndex();
        TPZCompEl *cskel = cmesh->Element(skeleton);
        TPZGeoEl *gskel = cskel->Reference();
        int ncorner = gskel->NCornerNodes();
        for (int ic=0; ic<ncorner; ic++) {
            nodeindices.insert(cskel->ConnectIndex(ic));
        }
    }
    int nc = elgr->NConnects();
    int neq = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = elgr->Connect(ic);
        neq += c.NShape()*c.NState();
    }
    indices.Resize(neq, 0);
    indices.Fill(0);
    int count = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = elgr->Connect(ic);
        int neqcon = c.NShape()*c.NState();
        int64_t cindex = elgr->ConnectIndex(ic);
        if (nodeindices.find(cindex) == nodeindices.end()) {
            for (int eq = 0; eq<neqcon; eq++) {
                indices[count++] = 0;
            }
        }
        else
        {
            for (int eq = 0; eq<neqcon; eq++) {
                indices[count++] = 1;
            }
        }
    }

}
void ComputeLoadVector(TPZCompMesh &cmesh, TPZFMatrix<STATE> &rhs)
{
    int64_t nel = cmesh.NElements();
    int64_t neq = cmesh.NEquations();
    rhs.Redim(neq, 1);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!elgr) {
            continue;
        }
        TPZFMatrix<STATE> &mass = elgr->MassMatrix();
        std::cout << "Norm of mass matrix el " << el << " = " << Norm(mass) << std::endl;
        int64_t nrow = mass.Rows();
        TPZManVector<int64_t> indices;
        CornerEquations(elgr, indices);
        if (indices.size() != nrow) {
            DebugStop();
        }
        TPZManVector<STATE> elrhs(nrow,0.);

        if(nrow%3 != 0) DebugStop();
        for (int ir=2; ir<nrow; ir+=3) {
            for (int ic = 0; ic<nrow; ic++) {
                if (indices[ic]  == 0) {
                    continue;
                }
                elrhs[ir] += mass(ir,ic)*2400.*9.81;
            }
        }
        int nc = elgr->NConnects();
        int count = 0;
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = elgr->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int blsize = c.NShape()*c.NState();
            int64_t pos = cmesh.Block().Position(seqnum);
            for (int eq=0; eq<blsize; eq++) {
                rhs(pos+eq,0) += elrhs[count++];
            }
        }
    }
}
#ifdef USING_BOOST
#include "boost/crc.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

static void printvec(const std::string &name, TPZVec<boost::crc_32_type::value_type> &vec)
{
    std::ofstream out(name);
    int64_t nel = vec.size();
    for (int64_t el=0; el<nel; el++) {
        if(vec[el] != 0)
        {
            out << el << " " << vec[el] << std::endl;
        }
    }
}

#endif

void SolveSistDragon(TPZLinearAnalysis *an, TPZCompMesh *Cmesh, int numthreads)
{
    int gnumthreads = numthreads;
    
    int64_t nel = Cmesh->NElements();
#ifdef USING_BOOST
#include "boost/crc.hpp"

extern TPZVec<boost::crc_32_type::value_type> matglobcrc, eigveccrc, stiffcrc, matEcrc, matEInvcrc;

    matglobcrc.Resize(nel, 0);
    eigveccrc.Resize(nel, 0);
    stiffcrc.Resize(nel, 0);
    matEcrc.Resize(nel, 0);
    matEInvcrc.Resize(nel, 0);
    std::stringstream matglob,eigvec,stiff,sol,matE,matEInv;
    matglob << "matglob_" << gnumthreads << "_" << nel << ".txt";
    eigvec << "eigvec_" << gnumthreads << "_" << nel << ".txt";
    stiff << "stiff_" << gnumthreads << "_" << nel << ".txt";
    sol << "sol_" << gnumthreads << "_" << nel << ".txt";
    matE << "matE_" << gnumthreads << "_" << nel << ".txt";
    matEInv << "matEInv_" << gnumthreads << "_" << nel << ".txt";
#endif
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
#ifdef USING_MKL
    //    TPZSkylineStructMatrix strmat(Cmesh);
    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(Cmesh);
#endif
    strmat.SetNumThreads(gnumthreads);
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
    step.SetDirect(ECholesky);
    an->SetSolver(step);
    
    
    try {
        an->Assemble();
        TPZFMatrix<STATE> rhs;
        ComputeLoadVector(*Cmesh, rhs);
        an->Rhs() = rhs;
    } catch (...) {
#ifdef USING_BOOST
        printvec(matglob.str(), matglobcrc);
        printvec(eigvec.str(), eigveccrc);
        printvec(stiff.str(), stiffcrc);
        printvec(matE.str(), matEcrc);
        printvec(matEInv.str(), matEInvcrc);
#endif
        exit(-1);
    }
    
#ifdef USING_BOOST
    printvec(matglob.str(), matglobcrc);
    printvec(eigvec.str(), eigveccrc);
    printvec(stiff.str(), stiffcrc);
    printvec(matE.str(), matEcrc);
    printvec(matEInv.str(), matEInvcrc);
    
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    std::cout << "rhs norm " << Norm(an->Rhs()) << std::endl;
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
    
    
    {
        std::ofstream out(sol.str());
        an->Solution().Print("sol",out);
    }
#endif
    
}
int NumNeigh(TPZGeoElSide &gelside, int matid)
{
    if(gelside.Element()->MaterialId() != matid) DebugStop();
    int count = 1;
    TPZGeoElSide neighbour = gelside.Neighbour();
    while(neighbour != gelside)
    {
        if(neighbour.Element()->MaterialId() == matid) count++;
        neighbour = neighbour.Neighbour();
    }
    return count;
}


void AddNeighbours(TPZGeoMesh &gmesh, int64_t el, int matid, TPZStack<int64_t> &elstack)
{
    TPZGeoEl *gel = gmesh.Element(el);
    int nsides = gel->NSides();
    int dim = gel->Dimension();
    if (dim != gmesh.Dimension()-1) {
        DebugStop();
    }
    for (int is=0; is<nsides; is++) {
        if (gel->SideDimension(is) == dim-1) {
            TPZGeoElSide gelside(gel,is);
            int count = NumNeigh(gelside, matid);
            if (count == 2) {
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->MaterialId() == matid) {
                        elstack.Push(neighbour.Element()->Index());
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            else if(count%2 != 0)
            {
                DebugStop();
            }
        }
    }
}

int64_t ElSeed(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup)
{
    int64_t nel = gmesh.NElements();
    int dim = gmesh.Dimension();
    int64_t elseed = -1;
    for (int64_t el=0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh.Element(el);
        int geldim = gel->Dimension();
        int gelmatid = gel->MaterialId();
        if(boundarygroup[el] == -1 && gelmatid == matid && geldim == dim-1)
        {
            elseed = el;
            break;
        }
    }
    return elseed;
}
// boundary group group index of each boundary element
void BuildBoundaryGroups(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup)
{
    std::cout << "Building boundary groups matid = " << matid << "\n";
    int64_t nel = gmesh.NElements();
    boundarygroup.Resize(nel, -1);
    boundarygroup.Fill(-1);
    int curgroup = 0;
    int64_t elseed = ElSeed(gmesh, matid, boundarygroup);
    while(elseed != -1)
    {
        std::cout << "elseed = " << elseed << " matid " << gmesh.Element(elseed)->MaterialId() << " dim " <<
        gmesh.Element(elseed)->Dimension() << std::endl;
        TPZStack<int64_t> elstack;
        elstack.Push(elseed);
        while (elstack.size()) {
            int64_t el = elstack.Pop();
            if (boundarygroup[el] == -1)
            {
                boundarygroup[el] = curgroup;
                AddNeighbours(gmesh, el, matid, elstack);
            }
        }
        elseed = ElSeed(gmesh, matid, boundarygroup);
        curgroup++;
    }
    std::cout << "Num groups formed " << curgroup << std::endl;
}

// generate a plot file for each boundary group
void PlotBoundaryGroups(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup, const std::string &rootname)
{
    std::cout << "Plotting boundary groups\n";
    int64_t nel = gmesh.NElements();
    std::set<int> groups;
    int maxmat = 0;
    for (int64_t el=0; el<nel; el++) {
        int bgroup = boundarygroup[el];
        if (bgroup != -1) {
            groups.insert(bgroup);
        }
        if (gmesh.Element(el)->MaterialId() > maxmat) {
            maxmat = gmesh.Element(el)->MaterialId();
        }
    }
    int count = 0;
    for (auto it = groups.begin(); it != groups.end(); it++) {
        int64_t numel = 0;
        for (int64_t el=0; el<nel; el++) {
            if(boundarygroup[el] == *it)
            {
                gmesh.Element(el)->SetMaterialId(maxmat+10+count);
                numel++;
            }
        }
        std::cout << "boundary group " << *it << " number of elements " << numel << std::endl;
        std::stringstream sout;
        sout << rootname << "." << count << ".vtk";
        std::ofstream file(sout.str());
        std::set<int> matids;
        matids.insert(maxmat+10+count);
        TPZVTKGeoMesh::PrintGMeshVTKmy_material(&gmesh, file, matids);
        count++;
    }
    for (int64_t el=0; el<nel; el++) {
        if(boundarygroup[el] != -1) gmesh.Element(el)->SetMaterialId(matid);
    }
    std::cout << "done\n";
}

// output the volume within each boundary
void IntegrateVolumes(TPZGeoMesh &gmesh, TPZVec<int> &boundarygroup)
{
    int64_t nel = gmesh.NElements();
    std::map<int, REAL> integral[3];
    for (int64_t el=0; el<nel; el++) {
        int bgroup = boundarygroup[el];
        if (bgroup != -1) {
            TPZGeoEl *gel = gmesh.Element(el);
            int nsides = gel->NSides();
            TPZIntPoints *intrule = gel->CreateSideIntegrationRule(nsides-1, 3);
            int np = intrule->NPoints();
            for (int ip=0; ip<np; ip++) {
                TPZManVector<REAL,3> xi(2), xco(3), normal(3,0.);
                REAL weight;
                intrule->Point(ip,xi,weight);
                TPZFNMatrix<9,REAL> axes(2,3),jac(2,2),jacinv(2,2);
                REAL detjac;
                gel->Jacobian(xi, jac, axes, detjac, jacinv);
                gel->X(xi,xco);
                normal[0] = axes(0,1)*axes(1,2)-axes(1,1)*axes(0,2);
                normal[1] = -(axes(0,0)*axes(1,2)-axes(1,0)*axes(0,2));
                normal[2] = axes(0,0)*axes(1,1)-axes(1,0)*axes(0,1);
                double normalnorm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
                if (abs(normalnorm-1.) > 1.e-6) {
                    std::cout << "normal not computed right normalnorm " << normalnorm << " normal " << normal << std::endl;
                }
                for (int i=0; i<3; i++)
                {
                    integral[i][bgroup] += abs(detjac)*xco[i]*normal[i]*weight;
                }
            }
            delete intrule;
        }
    }
    for (int i=0; i<3; i++)
    {
        for (auto it = integral[i].begin(); it != integral[i].end(); it++) {
            std::cout << "Integral of bgroup " << it->first << "is equal " << it->second << std::endl;
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
    TPZGeoEl *hexa = new TPZGeoElRefLess<pzgeom::TPZGeoCube>(nodeindices, matid,gmesh, index);
    TPZGeoEl *prism = new TPZGeoElRefLess<pzgeom::TPZGeoPrism>(nodeindices, matid,gmesh, index);
    for (int64_t el=0; el<nel; el++) {
        if (elpartitions[el] == -1) {
            continue;
        }
        int64_t elpartition = elpartitions[el];
        TPZGeoEl *gel = gmesh.Element(el);
        if (gel->Type() == EQuadrilateral) {
            for (int i=0; i<4; i++) {
                hexa->SetNodeIndex(i, gel->NodeIndex(i));
            }
            for (int i=4; i<8; i++) {
                hexa->SetNodeIndex(i, scalingcenterindices[elpartition]);
            }
            REAL area = hexa->SideArea(26);
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
            for (int i=0; i<3; i++) {
                prism->SetNodeIndex(i, gel->NodeIndex(i));
            }
            for (int i=3; i<6; i++) {
                prism->SetNodeIndex(i, scalingcenterindices[elpartition]);
            }
            REAL area = prism->SideArea(prism->NSides()-1);
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

// print boundary group neighbours
void PrintBoundaryGroupNeighbourPartitions(TPZGeoMesh &gmesh, TPZVec<int> &boundarygroup, TPZVec<int64_t> &elpartitions
                                           , TPZVec<int64_t> &scalingcenterindices)
{
    int64_t nel = gmesh.NElements();
    for (int64_t el=0; el<nel; el++) {
        if (boundarygroup[el] > 0) {
            TPZGeoEl *gel = gmesh.Element(el);
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int64_t index = neighbour.Element()->Index();
            int partition = elpartitions[index];
            TPZManVector<REAL,3> xco(3);
            gmesh.NodeVec()[scalingcenterindices[partition]].GetCoordinates(xco);
            std::cout << "Boundary group " << boundarygroup[el] << " Element partition " << partition << " with scaling center " << xco << std::endl;
        }
    }
}

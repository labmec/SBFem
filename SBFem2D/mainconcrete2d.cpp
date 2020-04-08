#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"

#include "pzgeoelbc.h"

#include "TPZMatElasticity2D.h"
#include "pzbndcond.h"

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

int Ebr = 0;
int Eagr = 98;
int Ezero = 238;

void BodyLoadsFracture(const TPZVec<REAL> &x, TPZVec<REAL> &val);

void InsertMaterialObjectsConcrete(TPZCompMesh *cmesh);

void AddBoundaryElements(TPZGeoMesh *gmesh);

void BuildBoundaryGroups(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup);

void PrintBoundaryGroupNeighbourPartitions(TPZGeoMesh &gmesh, TPZVec<int> &boundarygroup, TPZVec<int64_t> &elpartitions
                                           , TPZVec<int64_t> &scalingcenterindices);

void IntegrateVolumes(TPZGeoMesh &gmesh, TPZVec<int> &boundarygroup);

void PlotBoundaryGroups(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup, const std::string &rootname);

void AddNeighbours(TPZGeoMesh &gmesh, int64_t el, int matid, TPZStack<int64_t> &elstack);

int NumNeigh(TPZGeoElSide &gelside, int matid);

int64_t ElSeed(TPZGeoMesh &gmesh, int matid, TPZVec<int> &boundarygroup);

void AdjustElementOrientation(TPZGeoMesh &gmesh, TPZVec<int64_t> &elpartitions, TPZVec<int64_t> &scalingcenterindices);

void SolveSistDFN(TPZAnalysis *an, TPZCompMesh *Cmesh, int numthreads);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = false;

    int maxporder = 2;
    int counter = 1;
    
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;

        std::string filename("../../../src/SBFem/SBFem2D/concretemesh.txt");
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
        TPZGeoMesh *gmesh = ReadUNSWSBGeoFile_v3(filename, elpartition, scalingcenterindices);
        AdjustElementOrientation(*gmesh, elpartition, scalingcenterindices);
        AddBoundaryElements(gmesh);
        elpartition.Resize(gmesh->NElements(), -1);

        if(1)
        {
            std::cout << "Checking the mesh\n";
            TPZVec<int> boundarygroups;
            // BuildBoundaryGroups(*gmesh, Ebc1, boundarygroups);
            // BuildBoundaryGroups(*gmesh, Ebc2, boundarygroups);
            // BuildBoundaryGroups(*gmesh, Ebc3, boundarygroups);
            // BuildBoundaryGroups(*gmesh, Ebc4, boundarygroups);
            // BuildBoundaryGroups(*gmesh, ESkeleton, boundarygroups);
            // print boundary group neighbours
            // PrintBoundaryGroupNeighbourPartitions(*gmesh, boundarygroups, elpartition
                    // , scalingcenterindices);

            // elpartition.Resize(gmesh->NElements(), -1);
        }

        elpartition.Resize(gmesh->NElements(), -1);
        scalingcenterindices.Resize(gmesh->NElements(), -1);

        if(1) {
            std::cout << "Plotting the geometric mesh\n";
            std::ofstream outtxt("gmesh.txt");
            gmesh->Print(outtxt);
            std::ofstream out(vtkfilegeom);
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(gmesh, out,true);
        }

        std::cout << "Building computational mesh\n";
        std::map<int,int> matmap;
        matmap[ESkeleton] = Emat1;
        // matmap[ESkeleton] = Ebr;
        // matmap[ESkeleton] = Eagr;
        // matmap[ESkeleton] = Ezero;
        
        TPZBuildSBFem build(gmesh,ESkeleton,matmap);
        build.SetPartitions(elpartition, scalingcenterindices);
        // build.StandardConfiguration();
        // build.Configure(scalingcenterindices);
        build.DivideSkeleton(0);

        TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
        SBFem->SetDefaultOrder(POrder);
        InsertMaterialObjectsConcrete(SBFem);
        
        gmesh = SBFem->Reference();
        elpartition.Resize(gmesh->NElements(), -1);
        scalingcenterindices.Resize(gmesh->NElements(), -1);
        
        // std::set<int> volmatids,boundmatids;
        // volmatids.insert(Eagr);
        // volmatids.insert(Ebr);
        // volmatids.insert(238);
        // boundmatids.insert(Ebc1);
        // boundmatids.insert(Ebc2);
        // boundmatids.insert(Ebc3);
        // boundmatids.insert(Ebc4);
        // boundmatids.insert(ESkeleton);
        // build.BuildComputationMesh(*SBFem,volmatids,boundmatids);

        // std::cout << scalingcenterindices[elpartition[16430]] << std::endl;
        // std::cout << elpartition[1399] << std::endl;
        build.BuildComputationalMeshFromSkeleton(*SBFem);

        if(1) {
            std::cout << "Plotting the computational mesh\n";
            std::ofstream out(vtkfilecmesh);
            TPZVTKGeoMesh vtk;
            vtk.PrintCMeshVTK(SBFem, out,true);
        }

        // Visualization of computational meshes
        bool mustOptimizeBandwidth = true;
        TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
        Analysis->SetStep(counter++);
        std::cout << "neq = " << SBFem->NEquations() << std::endl;

        int numthreads = 0;
        SolveSistDFN(Analysis, SBFem, numthreads);

        TPZStack<std::string> vecnames,scalnames;
        // scalar
        vecnames.Push("Displacement");
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        scalnames.Push("EpsX");
        scalnames.Push("EpsY");
        scalnames.Push("EpsXY");
        std::stringstream sout;
        sout << "../ConcreteElasticity2DSolution";
        sout << "_SBFem.vtk";
        Analysis->DefineGraphMesh(2, scalnames, vecnames,sout.str());
        Analysis->PostProcess(3);

    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void AddBoundaryElements(TPZGeoMesh *gmesh)
{
    std::set<int64_t> setbottom,setright,settop,setleft;
    int64_t nnodes = gmesh->NNodes();
    int dim = gmesh->Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[in].GetCoordinates(xco);
        if (fabs(xco[1]-0.005) < 1.e-3) {
            setbottom.insert(in);
        }
        if (fabs(xco[0]-2.565) < 1.e-3) {
            setright.insert(in);
        }
        if (fabs(xco[1]-2.565) < 1.e-3) {
            settop.insert(in);
        }
        if (fabs(xco[0]-0.005) < 1.e-3) {
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
                    TPZGeoElBC(gelside,Emat1);
                }
            }
        }
    }
}

void InsertMaterialObjectsConcrete(TPZCompMesh *cmesh)
{
    
    // Getting mesh dimension
    int matId1 = Emat1;
    
    TPZMaterial *material;
    int nstate = 2;
    {
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(BodyLoadsFracture, 6);
        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;

        TPZMatElasticity2D *matloc0 = new TPZMatElasticity2D(Emat1);
        matloc0->SetElasticity(20000, 0.2);

        // TPZMatElasticity2D *matloc = new TPZMatElasticity2D(Ebr);
        // matloc->SetElasticity(20000, 0.2);

        // TPZMatElasticity2D *matloc2 = new TPZMatElasticity2D(Eagr);
        // matloc2->SetElasticity(200000, 0.3);

        // TPZMatElasticity2D *matloc3 = new TPZMatElasticity2D(Ezero);
        // matloc2->SetElasticity(200000, 0.3);

        matloc0->SetForcingFunction(autodummy);
        // matloc->SetForcingFunction(autodummy);
        // matloc2->SetForcingFunction(autodummy);
        // matloc3->SetForcingFunction(autodummy);

        cmesh->InsertMaterialObject(matloc0);
        // cmesh->InsertMaterialObject(matloc);
        // cmesh->InsertMaterialObject(matloc2);
        // cmesh->InsertMaterialObject(matloc3);

        material = matloc0;
    }
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
    {
        TPZMaterial * BCond1 = material->CreateBC(material,Ebc1,0, val1, val2);
        cmesh->InsertMaterialObject(BCond1);
        val2.Zero();
    }
    
    {
        val1.Zero();
        val2.Zero();
        TPZMaterial *BCond2 = material->CreateBC(material, Ebc2, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    
    {
        val1.Zero();
        val2.Zero();
        val2(1,0) = -10.0;
        TPZMaterial *BCond2 = material->CreateBC(material, Ebc3, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    
    {
        val1.Zero();
        val2.Zero();
        TPZMaterial *BCond2 = material->CreateBC(material, Ebc4, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    cmesh->InsertMaterialObject(BSkeleton);

}

void BodyLoadsFracture(const TPZVec<REAL> &x, TPZVec<REAL> &val)
{
    val[0] = 0.;
    val[1] = -1.0;
}

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
//            else if(count%2 != 0)
//            {
//                DebugStop();
//            }
        }
    }
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

            // REAL bax = nodes[1].Coord(0) - nodes[0].Coord(0);
            // REAL cay = nodes[2].Coord(1) - nodes[0].Coord(1);
            // REAL cax = nodes[2].Coord(0) - nodes[0].Coord(0);
            // REAL bay = nodes[1].Coord(1) - nodes[0].Coord(1);

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

            // REAL bax = nodes[1].Coord(0) - nodes[0].Coord(0);
            // REAL cay = nodes[2].Coord(1) - nodes[0].Coord(1);
            // REAL cax = nodes[2].Coord(0) - nodes[0].Coord(0);
            // REAL bay = nodes[1].Coord(1) - nodes[0].Coord(1);

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

void SolveSistDFN(TPZAnalysis *an, TPZCompMesh *Cmesh, int numthreads)
{
    int gnumthreads = numthreads;
    
    int64_t nel = Cmesh->NElements();

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix strmat(Cmesh);
#endif
    // strmat.SetNumThreads(gnumthreads);
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
    
    
    try {
        an->Assemble();
    } catch (...) {
        exit(-1);
    }
    
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
#include <chrono>

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common3D.h"
#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "pzgeoelbc.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
void AnalyseSolution(TPZCompMesh *cmesh);
#endif

void AddBoundaryElements(TPZGeoMesh * gmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    auto minporder = 1, maxporder = 4;
	auto numthreads = 8;
    auto scalarproblem = true;
    
#ifdef _AUTODIFF
    ExactLaplace.fExact = TLaplaceExample1::EHarmonic2;
#endif
    
    for (auto POrder = minporder; POrder <= maxporder; POrder += 1)
    {
        TPZGeoMesh * gmesh = new TPZGeoMesh();
        std::string filename = "n512-id1-4.msh";
        std::cout << "Reading " << filename << "and creating geometric mesh...\n";
        TPZGmshReader gmshrdr;
        gmesh = gmshrdr.GeometricGmshMesh4(filename, gmesh);
        AddBoundaryElements(gmesh);
        std::ofstream file("GmeshTetrahedrons.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);

        std::cout << "Creating computational mesh...\n";
        TPZCompMesh cmesh(gmesh);
        cmesh.SetDefaultOrder(POrder);

        int nstate = 1;
        TPZMatLaplacian *matloc = new TPZMatLaplacian(0);
        TPZMaterial * material;
        matloc->SetDimension(gmesh->Dimension());
        matloc->SetSymmetric();
        material = matloc;
        cmesh.InsertMaterialObject(matloc);

        TPZFMatrix<STATE> val1(nstate, nstate, 0.), val2(nstate, 1, 0.);
        TPZMaterial *BCond1 = material->CreateBC(material, Ebc1, 0, val1, val2);
#ifdef _AUTODIFF
        BCond1->SetForcingFunction(ExactLaplace.TensorFunction());
#endif
        cmesh.InsertMaterialObject(BCond1);

        cmesh.AutoBuild();

        std::cout << "Creating analysis...\n";
        TPZAnalysis an(&cmesh);
#ifdef USING_MKL
        TPZSymetricSpStructMatrix strmat(&cmesh);
#else
        TPZSkylineStructMatrix strmat(&cmesh);
#endif
        strmat.SetNumThreads(numthreads);
        an.SetStructuralMatrix(strmat);

        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);

        an.Run();

#ifdef _AUTODIFF
        an.SetExact(Laplace_exact);
#endif

        // Computing errors
        auto start = chrono::steady_clock::now();
        std::cout << "Compute errors\n";

        int64_t neq = cmesh.NEquations();
        TPZManVector<REAL,3> errors(3,0.);
        an.SetThreadsForError(numthreads);
        std::stringstream sout("RegularSolutionFEM3D.txt");
        an.PostProcessError(errors,false);
        
        std::ofstream results(sout.str(),std::ios::app);
        results.precision(15);
        results << "(* POrder " << POrder << " neq " << neq << "*)" << std::endl;
        TPZFMatrix<double> errmat(1,3);
        for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
        std::stringstream varname;
        varname << "ErrmatPolyFEM[[" << POrder << "]] = (1/1000000)*";
        errmat.Print(varname.str().c_str(),results,EMathematicaInput);
        auto end = chrono::steady_clock::now();
        cout << "Elapsed time to compute error (miliseconds): " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << "\n";
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void AddBoundaryElements(TPZGeoMesh * gmesh)
{
    std::set<int64_t> setboundaries;
    for (auto node : gmesh->NodeVec())
    {
        // if (!node) continue;
        TPZManVector<REAL> xco(3,0);
        node.GetCoordinates(xco);
        if (fabs(xco[0]) < 1.e-3) {
            setboundaries.insert(node.Id());
        }
        if (fabs(xco[0]-1.) < 1.e-3) {
            setboundaries.insert(node.Id());
        }
        if (fabs(xco[1]) < 1.e-3) {
            setboundaries.insert(node.Id());
        }
        if (fabs(xco[1]-1.) < 1.e-3) {
            setboundaries.insert(node.Id());
        }
        if (fabs(xco[2]) < 1.e-3) {
            setboundaries.insert(node.Id());
        }
        if (fabs(xco[2]-1.) < 1.e-3) {
            setboundaries.insert(node.Id());
        }
    }
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        int nsides = gel->NSides();
        for (auto isides = 0; isides < nsides; isides++)
        {
            if (gel->SideDimension(isides) != gmesh->Dimension()-1) continue;
            auto nsidenodes = gel->NSideNodes(isides);
            auto nfoundboundaries = 0;
            for (auto in = 0; in < nsidenodes; in++)
            {
                auto nodeindex = gel->SideNodeIndex(isides, in);
                if (setboundaries.find(nodeindex) != setboundaries.end())
                {
                    nfoundboundaries++;
                }
            }
            if (nfoundboundaries == nsidenodes) {
                TPZGeoElBC gelbc(gel,isides,Ebc1);
            }
        }
    }
}
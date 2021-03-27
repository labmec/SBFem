/**
 * @file
 * @brief Contains the methods of the TPZHybridPoissonCollapsed class (multiphysics environment)
 * @author Karolinne Coelho
 * @date 2020/11/17
 */

#include "hybridpoissoncollapsed.h"
#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
#include "pzinterpolationspace.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
static LoggerPtr logerror(Logger::getLogger("pz.mixedpoisson.error"));
#endif

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(): TPZRegisterClassId(&TPZHybridPoissonCollapsed::ClassId), TPZMixedPoisson() {
    fvisc = 1.;
	ff = 0.;
    fPermeabilityFunction = NULL;
}

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(int matid, int dim): TPZRegisterClassId(&TPZHybridPoissonCollapsed::ClassId), TPZMixedPoisson(matid,dim) {
    if (dim < 1) {
        DebugStop();
    }
    fvisc = 1.;
	ff = 0.;
    fPermeabilityFunction = NULL;
}

TPZHybridPoissonCollapsed::~TPZHybridPoissonCollapsed() {
}

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(const TPZHybridPoissonCollapsed &cp) :TPZRegisterClassId(&TPZHybridPoissonCollapsed::ClassId), TPZMixedPoisson(cp) {
    fvisc = cp.fvisc;
    ff = cp.ff;
    fPermeabilityFunction = cp.fPermeabilityFunction;
}

TPZHybridPoissonCollapsed & TPZHybridPoissonCollapsed::operator=(const TPZHybridPoissonCollapsed &copy){
    TPZMatPoisson3d::operator=(copy);
    fvisc = copy.fvisc;
    ff = copy.ff;
    fPermeabilityFunction = copy.fPermeabilityFunction;
    return *this;
} 

int TPZHybridPoissonCollapsed::NStateVariables() const {
	return 1;
}

void TPZHybridPoissonCollapsed::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Base Class properties :";
	TPZMatPoisson3d::Print(out);
	out << "\n";
}

void TPZHybridPoissonCollapsed::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    STATE force = ff;
    if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(datavec[1].x,res);
		force = res[0];
	}

    weight *= 1/datavec[1].detjac;
    
    TPZFNMatrix<9,STATE> PermTensor;
    TPZFNMatrix<9,STATE> InvPermTensor;
    
    GetPermeabilities(datavec[1].x, PermTensor, InvPermTensor);
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
	TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFNMatrix<9,REAL> dphiPXY(3,dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);
    
    
    REAL &faceSize = datavec[0].HSize;
    if(fUseHdois==true){
        fh2 = faceSize*faceSize;
    }else fh2 = 1.;
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
	
    int nactive = 0;
    for (int i=0; i<datavec.size(); i++) {
        if (datavec[i].fActiveApproxSpace) {
            nactive++;
        }
    }
    
	//Calculate the matrix contribution for flux. Matrix A
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
        }

        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fDeformedDirections(id,jvecind);
            }
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<3; id++){
                for(int jd=0; jd<3; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            //jvecZ.Print("mat1 = ");
            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
            
        }
    }
    
	
	// Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        for (int jp=0; jp<phrp; jp++) {
            
            REAL fact = (-1.)*weight*phip(jp,0)*datavec[0].divphi(iq,0);
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
    
}

void TPZHybridPoissonCollapsed::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{

#ifdef PZDEBUG
    int nref =  datavec.size();
	if (nref != 2 ) {
        std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	if (bc.Type() > 2 ) {
        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
		DebugStop();
	}
#endif
	
	TPZFMatrix<REAL> phiQ = datavec[0].phi;
	TPZFMatrix<REAL> phiP = datavec[1].phi;
	int phrq = phiQ.Rows();
	int phrp = phiP.Rows();

	REAL v2;
    if (bc.HasForcingFunction())
    {
		TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(Dimension(),1);
		bc.ForcingFunction()->Execute(datavec[1].x,res,gradu);
		v2 = res[0];
	} else
    {
        v2 = bc.Val2()(0,0);
    }

	switch (bc.Type()) {
		case 0 :		// Dirichlet condition
			//primeira equacao
			for(int iq=0; iq<phrp; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (1.)*v2*phiP(iq,0)*weight*gBigNumber;

                for (int jq = 0; jq < phrp; jq++){
                    ek(iq,jq) += weight * gBigNumber * phiP(iq,0) * phiP(jq,0);
                }
            }
            break;
			
		case 1 :			// Neumann condition
        
			for(int in = 0 ; in < phrp; in++)
            {
				ef(in,0) += v2 *phiP(in,0) * weight;
			}
			break;
        
        case 2 :			// mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
				ef(iq,0) += v2*phiQ(iq,0)*weight;
				for (int jq = 0; jq < phrq; jq++) {
					ek(iq,jq) += weight*bc.Val1()(0,0)*phiQ(iq,0)*phiQ(jq,0);
				}
			}
            
            break;
	}
    
}

/** Returns the variable index associated with the name */
int TPZHybridPoissonCollapsed::VariableIndex(const std::string &name){
	if(!strcmp("Flux",name.c_str()))        return  31;
	if(!strcmp("Pressure",name.c_str()))    return  32;
    if(!strcmp("GradFluxX",name.c_str()))   return  33;
    if(!strcmp("GradFluxY",name.c_str()))   return  34;
    if(!strcmp("DivFlux",name.c_str()))   return  35;
    
    if(!strcmp("ExactPressure",name.c_str()))  return 36;
    if(!strcmp("ExactFlux",name.c_str()))  return 37;
    
    if(!strcmp("POrder",name.c_str()))        return  38;
    if(!strcmp("GradPressure",name.c_str()))        return  39;
    if(!strcmp("Divergence",name.c_str()))      return  40;
    if(!strcmp("ExactDiv",name.c_str()))        return  41;
    if (!strcmp("Derivative",name.c_str())) {
        return 42;
    }
    if (!strcmp("Permeability",name.c_str())) {
        return 43;
    }
    if(!strcmp("g_average",name.c_str()))        return  44;
    if(!strcmp("u_average",name.c_str()))        return  45;

    if(!strcmp("ExactFluxShiftedOrigin",name.c_str()))  return  46;
    return TPZMatPoisson3d::VariableIndex(name);
    
}

int TPZHybridPoissonCollapsed::NSolutionVariables(int var){
	if(var == 31) return 3;
	if(var == 32 || var==8) return 1;
    if(var == 33) return 3;
    if(var == 34) return 3;
    if(var == 35) return 1;
    if(var == 36) return 1;
    if(var == 37) return fDim;
    if(var == 38) return 1;
    if(var == 39) return fDim;
    if(var == 40 || var == 41) return 1;
    if(var == 42) return 3;
    if(var == 43) return 1;
    if(var == 44) return 1;
    if(var == 45) return 1;
    if(var == 46) return fDim;
    return TPZMatPoisson3d::NSolutionVariables(var);
}

void TPZHybridPoissonCollapsed::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    std::cout << "TPZHybridPoissonCollapsed::Errors : Not implemented yet \n";
    DebugStop();
}

void TPZHybridPoissonCollapsed::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    std::cout << "TPZHybridPoissonCollapsed::Errors : Not implemented yet \n";
    DebugStop();
}

void TPZHybridPoissonCollapsed::ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc)
{
    std::cout << "TPZHybridPoissonCollapsed::Errors : Not implemented yet \n";
    DebugStop();
}



int TPZHybridPoissonCollapsed::ClassId() const{
    return Hash("TPZHybridPoissonCollapsed") ^ TPZMatPoisson3d::ClassId() << 1;
}
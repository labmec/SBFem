/**
 * @file
 * @brief Contains the methods of the TPZHybridPoissonCollapsed class (multiphysics environment)
 * @author Karolinne Coelho
 * @date 2020/11/17
 */

#include "hybridpoissoncollapsed.h"
#include "pzlog.h"
#include "TPZBndCondT.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
#include "pzinterpolationspace.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
static LoggerPtr logerror(Logger::getLogger("pz.mixedpoisson.error"));
#endif

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(): TPZMixedDarcyFlow() {
}

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(int matid, int dim): TPZMixedDarcyFlow(matid,dim) {
}
TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(const TPZHybridPoissonCollapsed &cp) : TPZMixedDarcyFlow(cp) {
}

void TPZHybridPoissonCollapsed::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef)
{
    STATE force = 0;
    if (fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction(datavec[1].x, res);
        force = res[0];
    }

    weight *= 1/datavec[1].detjac;
    
    TPZFNMatrix<9,STATE> PermTensor;
    TPZFNMatrix<9,STATE> InvPermTensor;
    
    TPZFNMatrix<1, REAL> K(1, 1), invK(1, 1);
    TPZMixedDarcyFlow::fPermeabilityFunction(datavec[1].x, PermTensor, InvPermTensor);
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
	TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFNMatrix<9,REAL> dphiPXY(3,dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);
    
    
    REAL &faceSize = datavec[0].HSize;
    
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
            ek(iq,jq) += weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
            
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

void TPZHybridPoissonCollapsed::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
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
    if (bc.HasForcingFunctionBC())
    {
		TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(Dimension(),1);
		bc.ForcingFunctionBC()(datavec[1].x,res,gradu);
		v2 = res[0];
	} else
    {
        v2 = bc.Val2()[0];
    }
    auto BigNumber = TPZMixedDarcyFlow::fBigNumber;

	switch (bc.Type()) {
		case 0 :		// Dirichlet condition
			//primeira equacao
			for(int iq=0; iq<phrp; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (1.)*v2*phiP(iq,0)*weight*BigNumber;

                for (int jq = 0; jq < phrp; jq++){
                    ek(iq,jq) += weight * BigNumber * phiP(iq,0) * phiP(jq,0);
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

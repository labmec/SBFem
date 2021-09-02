//
//  hybridpoissoncollapsed.h
//  PZ
//

#pragma once

#include "TPZMatCombinedSpaces.h"
//#include "pzdiscgal.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZMaterial.h"
#include "pzfunction.h"

/**
 * @ingroup material
 * @author Karolinne Coelho
 * @since 11/17/2020
 * @brief Material to solve a mixed poisson problem 2d by multiphysics simulation with HDivCollapsed spaces
 * @brief Pressure(p): uses L2 space.  Velocity (Q): uses Hdiv space
 */
/**
 * \f$ Q = -(k/visc)*grad(p)  ==> Int{Q.q}dx - (k/visc)*Int{p*div(q)}dx + (k/visc)*Int{pD(q.n)}ds = 0  (Eq. 1)  \f$ 
 *
 * \f$ div(Q) = f  ==> Int{div(Q)*v}dx = Int{f*v}dx (Eq. 2) \f$ 
 *
 * \f$ p = pD in Dirichlet boundary and Q.n = qN in Neumann boundary\f$
 */


class TPZHybridPoissonCollapsed : public TPZMixedDarcyFlow {
    
// protected:
// 	/** @brief Forcing function value */
// 	REAL ff;
    
//     /** @brief Fluid viscosity*/
// 	REAL fvisc;
    
//     /** @brief Pointer to forcing function, it is the Permeability and its inverse */
//     TPZAutoPointer<TPZFunction<STATE> > fPermeabilityFunction;

//     /** @brief Lagrange multiplier */
//     REAL fMultiplier;
    
public:
    TPZHybridPoissonCollapsed();

        /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    [[maybe_unused]] TPZHybridPoissonCollapsed(int id, int dim);
	
	[[nodiscard]] virtual std::string Name(){ return "TPZHybridPoissonCollapsed"; }
    
    TPZHybridPoissonCollapsed(const TPZHybridPoissonCollapsed &cp);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) override;
	
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

};

//
//  hybridpoissoncollapsed.h
//  PZ
//

#pragma once

#include "TPZMaterial.h"
#include "pzdiscgal.h"
#include "mixedpoisson.h"
#include "pzpoisson3d.h"
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


class TPZHybridPoissonCollapsed : public TPZMixedPoisson {
    
protected:
	/** @brief Forcing function value */
	REAL ff;
    
    /** @brief Fluid viscosity*/
	REAL fvisc;
    
    /** @brief Pointer to forcing function, it is the Permeability and its inverse */
    TPZAutoPointer<TPZFunction<STATE> > fPermeabilityFunction;

    /** @brief Lagrange multiplier */
    REAL fMultiplier;
    
public:
    TPZHybridPoissonCollapsed();
    
    TPZHybridPoissonCollapsed(int matid, int dim);
    
    virtual ~TPZHybridPoissonCollapsed();
    
    TPZHybridPoissonCollapsed(const TPZHybridPoissonCollapsed &cp);
    
    TPZHybridPoissonCollapsed &operator=(const TPZHybridPoissonCollapsed &copy);
    
    virtual TPZMaterial * NewMaterial(){
        return new TPZHybridPoissonCollapsed(*this);
    }
    
    virtual void Print(std::ostream & out);
	
	virtual std::string Name(){ return "TPZHybridPoissonCollapsed"; }
    
    virtual int NStateVariables() const override;
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */	
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
	
    virtual int NEvalErrors() override {return 5;}

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;
    
    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc) override;

    virtual void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }

    public:
    virtual int ClassId() const  override;

};
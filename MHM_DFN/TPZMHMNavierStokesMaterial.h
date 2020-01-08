/*
 *  TPZMHMNavierStokesMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMatWithMem.h"
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "TPZNavierStokesMaterial.h"
#include "pztrnsform.h"

#ifndef TPZMHMNavierStokesMATERIAL
#define TPZMHMNavierStokesMATERIAL

#ifdef _AUTODIFF
#include "fad.h"
#include "tfad.h"
#endif

class TPZMHMNavierStokesMaterial : public TPZNavierStokesMaterial  {
    
protected:
    
    STATE fMultiplier;
    
    
public:

    /**
     * Empty Constructor
     */
    TPZMHMNavierStokesMaterial() : TPZNavierStokesMaterial(), fMultiplier(1.)
    {
    }
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZMHMNavierStokesMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma) : TPZNavierStokesMaterial(matid,dimension,space,viscosity,theta,Sigma), fMultiplier(1.)
    {

    }
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZMHMNavierStokesMaterial(const TPZNavierStokesMaterial &mat) : TPZNavierStokesMaterial(mat)
    {}
    
    /**
     * Destructor
     */
    ~TPZMHMNavierStokesMaterial()
    {
    }
    
    virtual void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }

    virtual void FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right) override;
    
    virtual TPZMaterial *NewMaterial() override
    {
        return new TPZMHMNavierStokesMaterial(*this);
    }
    
    TPZFMatrix<STATE> Transpose(TPZFMatrix<STATE> &GradU );
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    
    TPZManVector<REAL,3> ComputeNormal(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright);
    
};

#endif

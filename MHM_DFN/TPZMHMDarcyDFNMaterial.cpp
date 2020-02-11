/*
 *  TPZMHMDarcyDFNMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMHMDarcyDFNMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"

using namespace std;

void TPZMHMDarcyDFNMaterial::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
{
    TPZMaterial::FillDataRequirementsInterface(data, datavec_left, datavec_right);
    int nref_left = datavec_left.size();
    datavec_left[0].fNeedsNormal = true;
    datavec_left[0].fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

void TPZMHMDarcyDFNMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
 
    DebugStop();
    
}



void TPZMHMDarcyDFNMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int nshapeLambda = 0;
    nshapeLambda = data.phi.Rows();
    
    //TPZFNMatrix<9, STATE> u_D(dim,1);
    STATE p_D =0., value =0.;
    TPZFMatrix<STATE> * phi;
    phi = &data.phi;
    
    TPZManVector<REAL> x;
    x=data.x;
    
    TPZManVector<STATE> p_bc(1,0.);
    if(bc.HasBCForcingFunction())
    {
        TPZFMatrix<STATE> gradu;
        bc.BCForcingFunction()->Execute(x,p_bc,gradu);
        p_D = p_bc[0];
        value = p_D;
        
        
        switch (bc.Type()) {
            case 0: //Dirichlet for continuous formulation
            {
                for(int j1 = 0; j1 < phi->Rows(); j1++)
                {
                    ef(j1,0) += p_D*(*phi)(j1,0)*weight;
                }
            }
                break;
                
                
            case 1: //Neumann for continuous formulation
            {

                for(int j1 = 0; j1 < phi->Rows(); j1++)
                {
                    ef(j1,0) += gBigNumber*value*(*phi)(j1,0)*weight;
                    
                    for(int i1 = 0; i1 < phi->Rows(); i1++)
                    {
                        ek(i1,j1) += gBigNumber*weight*(*phi)(j1,0)*(*phi)(i1,0);
                    }
                }
                
            }
                break;
            
            case 5: //Ponto pressao
            {

                //return;
                p_D = bc.Val2()(0,0);
                
                
                for(int i = 0; i < phi->Rows(); i++ )
                {
                    ef(i) += 1.0 * p_D * (*phi)(i,0);
                    
                    for(int j = 0; j < phi->Rows(); j++){
        
                        ek(i,j) += 1.0 * ((*phi)(i,0) * (*phi)(j,0));
        
                    }
                    
                }
                
            }
            break;
            
            default:
            {
                std::cout << "Boundary not implemented " << std::endl;
                DebugStop();
            }
            break;
    
        }
    }

}


void TPZMHMDarcyDFNMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int pindex = PIndex();
    int Lindex = LambdaIndex();
    
    int nshapeP = datavec[pindex].phi.Rows();
    int nshapeL = datavec[Lindex].phi.Rows();
    
    //TPZFNMatrix<9, STATE> u_D(dim,1);
    STATE p_D =0., value =0.;
    TPZFMatrix<STATE> * phi;
    phi = &datavec[pindex].phi;
    
    if(nshapeP!=0&&nshapeL!=0){
        DebugStop();
    }
    TPZManVector<REAL> x;
    if (nshapeP!=0) {
        x=datavec[pindex].x;
    }else{
        x=datavec[Lindex].x;
    }

    
    TPZManVector<STATE> p_bc(1,0.);
    if(bc.HasBCForcingFunction())
    {
        TPZFMatrix<STATE> gradu;
        bc.BCForcingFunction()->Execute(x,p_bc,gradu);
        p_D = p_bc[0];
        
    }else{
        p_D = bc.Val1()(0,0);
    }
    value = p_D;
    
        switch (bc.Type()) {
            case 0: //Dirichlet for continuous formulation
            {
                for(int j1 = 0; j1 < phi->Rows(); j1++)
                {
                    ef(j1,0) += p_D*(*phi)(j1,0)*weight;
                }
            }
                break;
                
                
            case 1: //Neumann for continuous formulation
            {
                
                //DebugStop();
                
                for(int j1 = 0; j1 < phi->Rows(); j1++)
                {
                    ef(j1,0) += gBigNumber*value*(*phi)(j1,0)*weight;
                    
                    for(int i1 = 0; i1 < phi->Rows(); i1++)
                    {
                        ek(i1,j1) += gBigNumber*weight*(*phi)(j1,0)*(*phi)(i1,0);
                    }
                }
            }
                break;
                
            case 5: //Ponto pressao
            {
                
                //return;
                p_D = bc.Val2()(0,0);
                
                
                for(int i = 0; i < phi->Rows(); i++ )
                {
                    ef(i) += 1.0 * p_D * (*phi)(i,0);
                    
                    for(int j = 0; j < phi->Rows(); j++){
                        
                        ek(i,j) += 1.0 * ((*phi)(i,0) * (*phi)(j,0));
                        
                    }
                    
                }
                
            }
                break;
                
            default:
            {
                std::cout << "Boundary not implemented " << std::endl;
                DebugStop();
            }
                break;
                
        }
    
    
    
}


void TPZMHMDarcyDFNMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    DebugStop();
    
    
}

TPZFMatrix<STATE> TPZMHMDarcyDFNMaterial::Transpose(TPZFMatrix<STATE> &MatrixU ){

    int dim = Dimension();
    TPZFMatrix<STATE> MatrixUt(dim,dim);
    
    if((MatrixU.Rows()!=dim)&&(MatrixU.Cols()!=dim)){
        DebugStop();
    }
    
    for (int i = 0; i < dim; i++ ) {
        for (int j = 0; j < dim; j++ ) {
            MatrixUt(i,j) = MatrixU(j,i);
         }
    }
    
    return MatrixUt;
    
}

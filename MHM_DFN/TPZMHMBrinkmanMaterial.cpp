/*
 *  TPZMHMBrinkmanMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMHMBrinkmanMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"

using namespace std;

void TPZMHMBrinkmanMaterial::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
{
    TPZMaterial::FillDataRequirementsInterface(data, datavec_left, datavec_right);
    int nref_left = datavec_left.size();
    datavec_left[0].fNeedsNormal = true;
    datavec_left[0].fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

void TPZMHMBrinkmanMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
 
    // DebugStop();
    // Verificar que
    // os termos mistos devem estar sem viscosidade!
    
    //2 = 1 Vel space + 1 Press space for datavecleft
    int nrefleft =  datavecleft.size();
    if (nrefleft != 4 ) {
    //    std::cout << " Erro. The size of the datavec is different from 2 \n";
    //    DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright
    int nrefright =  datavecright.size();
    if (nrefright != 4 ) {
   //     std::cout << " Erro. The size of the datavec is different from 2 \n";
   //     DebugStop();
    }
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    // Setting the phis
    // V - left
    
    TPZFNMatrix<9,REAL>  &tan = datavecright[pindex].axes;
    
    int nshapeV , nshapeP , nshapeLambda;
    
    nshapeV = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeP = datavecleft[pindex].phi.Rows();
    nshapeLambda = datavecright[pindex].phi.Rows();
    
    
    int normvecRows = datavecleft[vindex].fNormalVec.Rows();
    int normvecCols = datavecleft[vindex].fNormalVec.Cols();
    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
    
    if(datavecleft[vindex].fNeedsNormalVecFad){
#ifdef _AUTODIFF
        for (int e = 0; e < normvecRows; e++) {
            for (int s = 0; s < normvecCols; s++) {
                Normalvec(e,s)=datavecleft[vindex].fNormalVecFad(e,s).val();
            }
        }
#else
        DebugStop();
#endif
    }else{
        Normalvec=datavecleft[vindex].fNormalVec;
    }
    
    for(int i1 = 0; i1 < nshapeV; i1++)
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        TPZFNMatrix<9, STATE> phiVi(3,1);
        for (int e=0; e< 3 ; e++) {
            phiVi(e,0)=Normalvec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
        }
        
        // K12 e K21 - (test V left) * (trial Lambda right)
        for(int j1 = 0; j1 < nshapeLambda; j1++)
        {
            // Var. Sigma, Sn :
            TPZFNMatrix<9, STATE> lambda_j(3,1,0.);
            TPZFNMatrix<9, STATE> phiLamb = datavecright[pindex].phi;
            
            // Tangencial comp. vector (t x t)Sn :
            for (int e=0; e< 3 ; e++) {
                lambda_j(e,0) = phiLamb(j1,0)*tan(0,e);
            }
            
            STATE fact = fMultiplier * weight * InnerVec(phiVi,lambda_j);
            ek(i1,j1+nshapeV) += fact;
            ek(j1+nshapeV,i1) += fact;
        }
        
    }
    
    
    
}



void TPZMHMBrinkmanMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    int nshapeV , nshapeLambda;
    
    nshapeV = datavec[0].phi.Rows();
    nshapeLambda = datavec[1].phi.Rows();
    
    if(nshapeV!=0&&nshapeLambda!=0){
        DebugStop();
    }
    
    //TPZFNMatrix<9, STATE> u_D(dim,1);
    TPZVec<STATE> v_Dirichlet(3,0.);
    TPZVec<STATE> v_Neumann(3,0.);
    STATE p_D =0.;
    
    TPZFMatrix<STATE> * phi;
    
    TPZManVector<REAL> x;
    
    if (nshapeV!=0) {
        x=datavec[0].x;
    }else{
        x=datavec[1].x;
    }
    
    TPZFMatrix<STATE> gradu(3,3,0.);
    TPZManVector<STATE> vbc(4,0.);
    TPZFMatrix<STATE> Du(3,3,0.),Dun(3,1,0.);
    if(bc.HasBCForcingFunction())
    {
        
        bc.BCForcingFunction()->Execute(x,vbc,gradu);
        
        v_Dirichlet[0] = vbc[0];
        v_Dirichlet[1] = vbc[1];
        v_Dirichlet[2] = vbc[2];
        p_D = vbc[3];
        
        // Calculo do vetor (Du)n :
        if (nshapeV!=0) {
            gradu.Resize(3,3);
            STATE visco = GetViscosity();
            for (int i = 0; i<3; i++) {
                for (int j = 0; j<3; j++) {
                    Du(i,j)=  visco * (gradu(i,j) + gradu(j, i));
                }
            }
        
            for (int i = 0; i<3; i++) {
                for (int j = 0; j<3; j++) {
                    Dun(i,0) +=  Du(i,j) * datavec[0].normal[j];
                }
            }

            // Neumann vector
            
            for (int i = 0; i<3; i++) {
                v_Neumann[i] = Dun(i,0) - p_D * datavec[0].normal[i];
            }
            
        }
        
        
    }else{
        //DebugStop();
        //std::cout<<"!!!!!!!!! falta esta condição de contorno aqui !!!!!!!!!!!!!"<<std::endl;
    }

    
    STATE value = 0.;
    if(nshapeV!=0){
        
        if (bc.Type()==0) {

            for (int i = 0; i<3; i++) {
                value += v_Dirichlet[i]*datavec[0].normal[i];
            }

        }else if (bc.Type()==1){
            
            for (int j = 0; j<3; j++) {
                value += v_Neumann[j] * datavec[0].normal[j];
            }
        }
        
        phi = &datavec[0].phi;
        
    }else if (nshapeLambda!=0&&(bc.Type()!=5)){
        
        for (int i = 0; i<3; i++) {
            value += v_Dirichlet[i]*datavec[1].axes(0,i);
        }
        
        phi = &datavec[1].phi;
        
    }else if (bc.Type()==5){
        
        phi = &datavec[1].phi;
        
    }else {
        DebugStop();
    }
    
    switch (bc.Type()) {
        case 0: //Dirichlet for continuous formulation
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
            
            
        case 1: //Neumann for continuous formulation
        {
            for(int j1 = 0; j1 < phi->Rows(); j1++)
            {
                ef(j1,0) += value*(*phi)(j1,0)*weight;
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


void TPZMHMBrinkmanMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    
    STATE rhsnorm = Norm(ef);
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavec[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[vindex]);
    }
    
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    // nshapeV = phiV.Rows()*NStateVariables();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    //Dirichlet
    
    TPZFMatrix<STATE> v_2=bc.Val2();
    TPZFMatrix<STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    
    switch (bc.Type()) {
            
            
        case 0: //Dirichlet for continuous formulation
        {
            
            if(bc.HasBCForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                TPZFMatrix<STATE> gradu;
                bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D=vbc[2];
                
            }
            
            //Componente tangencial -> imposta fracamente:
            
            for(int i = 0; i < nshapeV; i++ )
            {
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                TPZFNMatrix<9,STATE> GradVni(fDimension,1,0.),phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                GradVni.Zero();
                
                TPZFNMatrix<4,STATE> GradVi(fDimension,fDimension,0.),GradVit(fDimension,fDimension,0.),Dui(fDimension,fDimension,0.),Duni(fDimension,1,0.);
                
                for (int e=0; e<fDimension; e++) {
                    
                    for (int f=0; f<fDimension; f++) {
                        GradVi(e,f) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                        //termo transposto:
                        GradVit(f,e) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                        
                    }
                }
                
                //Du = 0.5(GradU+GradU^T)
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        Dui(e,f)= (1./2.) * (GradVi(e,f) + GradVit(e,f));
                    }
                }
                
                //Duni
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        Duni(e,0) += Dui(e,f)*normal[f] ;
                    }
                }
                
                //GradVni
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        GradVni(e,0) += GradVi(e,f)*normal[f] ;
                    }
                }
                
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*datavec[vindex].phi(iphi,0);
                    phiVni(0,0)+=phiVi(e,0)*normal[e];
                    
                }
                
                TPZManVector<REAL> n = data.normal;
                TPZManVector<REAL> t(2);
                t[0]=n[1];
                t[1]=n[0];
                
                
                
                phiVti(0,0)= t[0] * phiVi(0,0) + t[1] * phiVi(1,0);
                TPZFNMatrix<9,STATE> phiVtit(fDimension,1,0.);
                phiVtit(0,0)=phiVti(0,0)*t[0];
                phiVtit(1,0)=phiVti(0,0)*t[1];
                
                TPZFNMatrix<9,STATE> phiVnin(fDimension,1,0.);
                phiVnin(0,0)=phiVni(0,0)*n[0];
                phiVnin(1,0)=phiVni(0,0)*n[1];
                
                
                if(fSpace==1||fSpace==3){
                    
                    REAL vh_t = 0.;
                    if(v_h.size()>1){
                        vh_t = v_h[1];
                    }
                    
                    REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
                    
                    TPZManVector<REAL> v_tt(2);
                    v_tt[0]=v_t*t[0];
                    v_tt[1]=v_t*t[1];
                    
                    TPZManVector<REAL> vh_tt(2);
                    vh_tt[0]=vh_t*t[0];
                    vh_tt[1]=vh_t*t[1];
                    
                    TPZFNMatrix<9,STATE> diffvt(fDimension,1,0.);
                    diffvt(0,0)=v_tt[0];
                    diffvt(1,0)=v_tt[1];
                    
                    
                    STATE factef = weight * fSigma * v_t * phiVti(0,0) * fViscosity;
                    
                    ef(i,0) += factef;
                    
                    
                    STATE fact= 2. * weight * fViscosity * InnerVec(diffvt, Duni) ;
                    
                    ef(i,0) += fTheta*fact;
                    
                    
                    for(int j = 0; j < nshapeV; j++){
                        int jphi = datavec[vindex].fVecShapeIndex[j].second;
                        int jvec = datavec[vindex].fVecShapeIndex[j].first;
                        
                        TPZFNMatrix<9,STATE> GradVnj(fDimension,1),phiVtj(1,1,0.),phiVj(fDimension,1);
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*datavec[vindex].phi(jphi,0);
                        }
                        
                        
                        phiVtj(0,0)= t[0] * phiVj(0,0) + t[1] * phiVj(1,0);
                        
                        
                        
                        TPZFNMatrix<4,STATE> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);
                        
                        for (int e=0; e<fDimension; e++) {
                            
                            for (int f=0; f<fDimension; f++) {
                                GradVj(e,f) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                                //termo transposto:
                                GradVjt(f,e) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                                
                            }
                        }
                        
                        //Du = 0.5(GradU+GradU^T)
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Duj(e,f)= (1./2.) * (GradVj(e,f) + GradVjt(e,f));
                            }
                        }
                        
                        //Du2nj
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Dunj(e,0) += Duj(e,f)*normal[f] ;
                            }
                        }
                        
                        STATE factek = weight * fSigma * phiVtj(0,0)* phiVti(0,0) * fViscosity;
                        ek(i,j) +=  factek;
                        
                        
                        STATE fact =(-1.) * weight * 2. * fViscosity * InnerVec(phiVtit, Dunj) ;
                        ek(i,j) += fact ;
                        ek(j,i) += -fTheta*fact;
                        
                        
                    }
                    
                    
                    //Componente normal -> imposta fortemente:
                    if(fSpace==1||fSpace==3){
                        for(int i = 0; i < nshapeV; i++ )
                        {
                            
                            int iphi = datavec[vindex].fVecShapeIndex[i].second;
                            int ivec = datavec[vindex].fVecShapeIndex[i].first;
                            TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                            
                            
                            for (int e=0; e<fDimension; e++) {
                                phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*datavec[vindex].phi(iphi,0);
                                phiVni(0,0)+=phiVi(e,0)*n[e];
                                phiVti(0,0)+=phiVi(e,0)*t[e];
                            }
                            
                            
                            REAL vh_n = v_h[0];
                            REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                            
                            ef(i,0) += -weight * gBigNumber * (vh_n-v_n) * (phiVni(0,0));
                            
                            
                            for(int j = 0; j < nshapeV; j++){
                                
                                int jphi = datavec[vindex].fVecShapeIndex[j].second;
                                int jvec = datavec[vindex].fVecShapeIndex[j].first;
                                
                                TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
                                
                                for (int e=0; e<fDimension; e++) {
                                    phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*datavec[vindex].phi(jphi,0);
                                    phiVnj(0,0)+=phiVj(e,0)*n[e];
                                    phiVtj(0,0)+=phiVj(e,0)*t[e];
                                    
                                }
                                
                                ek(i,j) += weight * gBigNumber * (phiVni(0,0)) * (phiVnj(0,0)) ;
                                
                            }
                            
                        }
                        
                        
                    }
                    
                }
                
                
            }
            break;
            
        case 1: //Neumann for continuous formulation
            {
                
                
                if(bc.HasBCForcingFunction())
                {
                    TPZManVector<STATE> vbc(3);
                    TPZFMatrix<STATE> gradu;
                    bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                    v_2(0,0) = vbc[0];
                    v_2(1,0) = vbc[1];
                    p_D=vbc[2];
                    
                }
                
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
                    TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*phiV(iphi,0);
                    }
                    
                    TPZManVector<REAL> n = data.normal;
                    
                    TPZFNMatrix<9,STATE> pn(fDimension,1);
                    
                    
                    for (int f=0; f<fDimension; f++) {
                        pn(f,0)=n[f]*v_1(0,0);
                    }
                    
                    //Adaptação para Hdiv
                    
                    STATE factef=0.0;
                    
                    factef += InnerVec(pn, phiVi) ;
                    
                    
                    ef(i,0) += weight * factef;
                    
                }
                
                
            }
            
            
            
            break;
            
        case 2: //Penetração com slip for continuous formulation
            {
                
                if(bc.HasBCForcingFunction())
                {
                    TPZManVector<STATE> vbc(3);
                    TPZFMatrix<STATE> gradu;
                    bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                    v_2(0,0) = vbc[0];
                    v_2(1,0) = vbc[1];
                    p_D=vbc[2];
                    
                }
                
                
                if(fSpace==1||fSpace==2){
                    TPZManVector<REAL> n = data.normal;
                    TPZManVector<REAL> t(2);
                    t[0]=-n[1];
                    t[1]=n[0];
                    
                    //Componente normal -> imposta fortemente:
                    
                    for(int i = 0; i < nshapeV; i++ )
                    {
                        
                        int iphi = datavec[vindex].fVecShapeIndex[i].second;
                        int ivec = datavec[vindex].fVecShapeIndex[i].first;
                        TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                        
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*datavec[vindex].phi(iphi,0);
                            phiVni(0,0)+=phiVi(e,0)*n[e];
                            phiVti(0,0)+=phiVi(e,0)*t[e];
                        }
                        
                        REAL vh_n = v_h[0];
                        REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                        
                        ef(i,0) += -weight * gBigNumber * (vh_n-v_n) * (phiVni(0,0));
                        
                        
                        for(int j = 0; j < nshapeV; j++){
                            
                            int jphi = datavec[vindex].fVecShapeIndex[j].second;
                            int jvec = datavec[vindex].fVecShapeIndex[j].first;
                            
                            TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
                            
                            for (int e=0; e<fDimension; e++) {
                                phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*datavec[vindex].phi(jphi,0);
                                phiVnj(0,0)+=phiVj(e,0)*n[e];
                                phiVtj(0,0)+=phiVj(e,0)*t[e];
                                
                            }
                            
                            ek(i,j) += weight * gBigNumber * (phiVni(0,0)) * (phiVnj(0,0)) ;
                            
                        }
                        
                    }
                }
                
                
                if(fSpace==3){
                    
                    for(int i = 0; i < nshapeV; i++ )
                    {
                        int iphi = datavec[vindex].fVecShapeIndex[i].second;
                        int ivec = datavec[vindex].fVecShapeIndex[i].first;
                        TPZFNMatrix<9,STATE> GradVni(fDimension,1,0.),phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                        GradVni.Zero();
                        
                        TPZFNMatrix<4,STATE> GradVi(fDimension,fDimension,0.),GradVit(fDimension,fDimension,0.),Dui(fDimension,fDimension,0.),Duni(fDimension,1,0.);
                        
                        for (int e=0; e<fDimension; e++) {
                            
                            for (int f=0; f<fDimension; f++) {
                                GradVi(e,f) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                                //termo transposto:
                                GradVit(f,e) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                                
                            }
                        }
                        
                        //Du = 0.5(GradU+GradU^T)
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Dui(e,f)= (1./2.) * (GradVi(e,f) + GradVit(e,f));
                            }
                        }
                        
                        //Duni
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Duni(e,0) += Dui(e,f)*normal[f] ;
                            }
                        }
                        
                        //GradVni
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                GradVni(e,0) += GradVi(e,f)*normal[f] ;
                            }
                        }
                        
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*datavec[vindex].phi(iphi,0);
                            phiVni(0,0)+=phiVi(e,0)*normal[e];
                        }
                        
                        TPZManVector<REAL> n = data.normal;
                        TPZManVector<REAL> t(2);
                        t[0]=n[1];
                        t[1]=n[0];
                        
                        phiVti(0,0)= t[0] * phiVi(0,0) + t[1] * phiVi(1,0);
                        TPZFNMatrix<9,STATE> phiVtit(fDimension,1,0.);
                        phiVtit(0,0)=phiVti(0,0)*t[0];
                        phiVtit(1,0)=phiVti(0,0)*t[1];
                        
                        TPZFNMatrix<9,STATE> phiVnin(fDimension,1,0.);
                        phiVnin(0,0)=phiVni(0,0)*n[0];
                        phiVnin(1,0)=phiVni(0,0)*n[1];
                        
                        
                        //Componente normal -> imposta fortemente:
                        
                        for(int i = 0; i < nshapeV; i++ )
                        {
                            
                            int iphi = datavec[vindex].fVecShapeIndex[i].second;
                            int ivec = datavec[vindex].fVecShapeIndex[i].first;
                            TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                            
                            
                            for (int e=0; e<fDimension; e++) {
                                phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*datavec[vindex].phi(iphi,0);
                                phiVni(0,0)+=phiVi(e,0)*n[e];
                                phiVti(0,0)+=phiVi(e,0)*t[e];
                            }
                            
                            
                            REAL vh_n = v_h[0];
                            REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                            
                            ef(i,0) += -weight * gBigNumber * (vh_n-v_n) * (phiVni(0,0));
                            
                            
                            for(int j = 0; j < nshapeV; j++){
                                
                                int jphi = datavec[vindex].fVecShapeIndex[j].second;
                                int jvec = datavec[vindex].fVecShapeIndex[j].first;
                                
                                TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
                                
                                for (int e=0; e<fDimension; e++) {
                                    phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*datavec[vindex].phi(jphi,0);
                                    phiVnj(0,0)+=phiVj(e,0)*n[e];
                                    phiVtj(0,0)+=phiVj(e,0)*t[e];
                                    
                                }
                                
                                ek(i,j) += weight * gBigNumber * (phiVni(0,0)) * (phiVnj(0,0)) ;
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            break;
            
            
            
        case 10: //Penetração com slip for continuous formulation
            {
                
                if(bc.HasBCForcingFunction())
                {
                    TPZManVector<STATE> vbc(3);
                    TPZFMatrix<STATE> gradu;
                    bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                    v_2(0,0) = vbc[0];
                    v_2(1,0) = vbc[1];
                    p_D=vbc[2];
                    
                }
                
                //Componente tangencial -> imposta fracamente:
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
                    TPZFNMatrix<9,STATE> GradVni(fDimension,1,0.),phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                    GradVni.Zero();
                    
                    TPZFNMatrix<4,STATE> GradVi(fDimension,fDimension,0.),GradVit(fDimension,fDimension,0.),Dui(fDimension,fDimension,0.),Duni(fDimension,1,0.);
                    
                    for (int e=0; e<fDimension; e++) {
                        
                        for (int f=0; f<fDimension; f++) {
                            GradVi(e,f) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                            //termo transposto:
                            GradVit(f,e) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                            
                        }
                    }
                    
                    //Du = 0.5(GradU+GradU^T)
                    for (int e=0; e<fDimension; e++) {
                        for (int f=0; f<fDimension; f++) {
                            Dui(e,f)= (1./2.) * (GradVi(e,f) + GradVit(e,f));
                        }
                    }
                    
                    //Duni
                    for (int e=0; e<fDimension; e++) {
                        for (int f=0; f<fDimension; f++) {
                            Duni(e,0) += Dui(e,f)*normal[f] ;
                        }
                    }
                    
                    //GradVni
                    for (int e=0; e<fDimension; e++) {
                        for (int f=0; f<fDimension; f++) {
                            GradVni(e,0) += GradVi(e,f)*normal[f] ;
                        }
                    }
                    
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*datavec[vindex].phi(iphi,0);
                        phiVni(0,0)+=phiVi(e,0)*normal[e];
                        
                    }
                    
                    TPZManVector<REAL> n = data.normal;
                    TPZManVector<REAL> t(2);
                    t[0]=-n[1];
                    t[1]=n[0];
                    
                    
                    
                    phiVti(0,0)= t[0] * phiVi(0,0) + t[1] * phiVi(1,0);
                    TPZFNMatrix<9,STATE> phiVtit(fDimension,1,0.);
                    phiVtit(0,0)=phiVti(0,0)*t[0];
                    phiVtit(1,0)=phiVti(0,0)*t[1];
                    
                    TPZFNMatrix<9,STATE> phiVnin(fDimension,1,0.);
                    phiVnin(0,0)=phiVni(0,0)*n[0];
                    phiVnin(1,0)=phiVni(0,0)*n[1];
                    
                    
                    REAL vh_t = v_h[1];
                    REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
                    TPZManVector<REAL> v_tt(2);
                    v_tt[0]=v_t*t[0];
                    v_tt[1]=v_t*t[1];
                    TPZManVector<REAL> vh_tt(2);
                    vh_tt[0]=vh_t*t[0];
                    vh_tt[1]=vh_t*t[1];
                    TPZFNMatrix<9,STATE> diffvt(fDimension,1,0.);
                    diffvt(0,0)=v_tt[0];
                    diffvt(1,0)=v_tt[1];
                    
                    REAL vh_n = v_h[0];
                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                    TPZManVector<REAL> v_nn(2);
                    v_nn[0]=v_n*n[0];
                    v_nn[1]=v_n*n[1];
                    TPZManVector<REAL> vh_nn(2);
                    vh_nn[0]=vh_n*n[0];
                    vh_nn[1]=vh_n*n[1];
                    TPZFNMatrix<9,STATE> diffvn(fDimension,1,0.);
                    diffvn(0,0)=v_nn[0];
                    diffvn(1,0)=v_nn[1];
                    
                    STATE factefn = weight * fSigma * v_n * phiVni(0,0);
                    ef(i,0) += factefn;
                    STATE factn= 2. * weight * fViscosity * InnerVec(diffvn, Duni);
                    ef(i,0) += fTheta*factn;
                    
                    STATE facteft = weight * fSigma * v_t * phiVti(0,0);
                    ef(i,0) += facteft;
                    STATE factt= 2. * weight * fViscosity * InnerVec(diffvt, Duni);
                    ef(i,0) += fTheta*factt;
                    
                    for(int j = 0; j < nshapeV; j++){
                        int jphi = datavec[vindex].fVecShapeIndex[j].second;
                        int jvec = datavec[vindex].fVecShapeIndex[j].first;
                        TPZFNMatrix<9,STATE> GradVnj(fDimension,1),phiVtj(1,1,0.),phiVnj(1,1,0.),phiVj(fDimension,1);
                        for (int e=0; e<fDimension; e++) {
                            phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*datavec[vindex].phi(jphi,0);
                        }
                        
                        phiVnj(0,0)= n[0] * phiVj(0,0) + n[1] * phiVj(1,0);
                        phiVtj(0,0)= t[0] * phiVj(0,0) + t[1] * phiVj(1,0);
                        
                        TPZFNMatrix<4,STATE> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);
                        for (int e=0; e<fDimension; e++) {
                            
                            for (int f=0; f<fDimension; f++) {
                                GradVj(e,f) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                                //termo transposto:
                                GradVjt(f,e) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                            }
                        }
                        
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Duj(e,f)= (1./2.) * (GradVj(e,f) + GradVjt(e,f));
                            }
                        }
                        
                        //Du2nj
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Dunj(e,0) += Duj(e,f)*normal[f] ;
                            }
                        }
                        
                        STATE factek = weight * fSigma * phiVtj(0,0)* phiVti(0,0);
                        ek(i,j) +=  factek;
                        
                        STATE fact =(-1.) * weight * 2. * fViscosity * InnerVec(phiVtit, Dunj) ;
                        ek(i,j) += fact ;
                        ek(j,i) += -fTheta*fact;
                        
                        
                        STATE factekn = weight * fSigma * phiVnj(0,0)* phiVni(0,0);
                        ek(i,j) +=  factekn;
                        
                        STATE factn =(-1.) * weight * 2. * fViscosity * InnerVec(phiVnin, Dunj) ;
                        ek(i,j) += factn ;
                        ek(j,i) += -fTheta*factn;
                        
                    }
                    
                }
                
                
            }
            break;
            
            
            
        }
            
            
            
    }
    
    
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    
    {
        std::ofstream fileEK("FileEKContributeBCInterf.txt");
        std::ofstream fileEF("FileEFContributeBCInterf.txt");
        ek.Print("MatrizBCint = ",fileEK,EMathematicaInput);
        ef.Print("ForceBCint = ",fileEF,EMathematicaInput);
    }
    
    
    
}


TPZManVector<REAL,3> TPZMHMBrinkmanMaterial::ComputeNormal(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright){

    int vindex = VIndex();
    int pindex = PIndex();
    
    TPZManVector<REAL,3> xcenterL = datavecleft[vindex].XCenter;
    TPZManVector<REAL,3> xcenterR = datavecright[pindex].XCenter;
    TPZManVector<REAL,3> normalV = data.normal;
    
    if(xcenterR[0]>xcenterL[0]){
        normalV[0]=1.;

    }else if(xcenterR[0]>xcenterL[0]){
        normalV[0]=-1.;
    }
    return normalV;
    
    if(xcenterR[1]>xcenterL[1]){
        
    }
    
    return normalV;
    
}

TPZFMatrix<STATE> TPZMHMBrinkmanMaterial::Transpose(TPZFMatrix<STATE> &MatrixU ){

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

/*
 *  TPZNavierStokesMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZNavierStokesMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"
#include "pzlog.h"

using namespace std;


TPZNavierStokesMaterial::TPZNavierStokesMaterial() : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(){
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1;
    
}

////////////////////////////////////////////////////////////////////

TPZNavierStokesMaterial::TPZNavierStokesMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma) : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(matid),fDimension(dimension),fSpace(space),fViscosity(viscosity),fTheta(theta),fSigma(Sigma)
{
    // symmetric version
    //fTheta = -1;
    
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1.;
    

    
}

////////////////////////////////////////////////////////////////////

TPZNavierStokesMaterial::TPZNavierStokesMaterial(const TPZNavierStokesMaterial &mat) : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(mat),fDimension(mat.fDimension),fSpace(mat.fSpace), fViscosity(mat.fViscosity), fTheta(mat.fTheta), fSigma(mat.fSigma)
{
    fk= mat.fk;
    
}

////////////////////////////////////////////////////////////////////

TPZNavierStokesMaterial::~TPZNavierStokesMaterial(){
    
    
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
}


////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
    datavec[0].fNeedsNormalVecFad = NeedsNormalVecFad;
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
        datavec[idata].fNeedsNeighborSol = true;
    }
    datavec[0].fNeedsNormalVecFad = NeedsNormalVecFad;
}

void TPZNavierStokesMaterial::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data){

    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    data.fNeedsNeighborSol = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::FillDataRequirementsInterface(TPZMaterialData &data)
{
    data.fNeedsNormal = true;
    data.fNeedsNeighborCenter = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
{
    int nref_left = datavec_left.size();
    for(int iref = 0; iref<nref_left; iref++){
        datavec_left[iref].SetAllRequirements(false);
        datavec_left[iref].fNeedsNormal = true;
    }
    datavec_left[0].fNeedsNormalVecFad = NeedsNormalVecFad;
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZNavierStokesMaterial::VariableIndex(const std::string &name) {
    
    if (!strcmp("P", name.c_str()))  return 0;
    if (!strcmp("Pressure", name.c_str()))  return 0;
    if (!strcmp("V", name.c_str()))  return 1;
    if (!strcmp("State", name.c_str()))  return 0;
    if (!strcmp("f", name.c_str()))         return 2;
    if (!strcmp("V_exact", name.c_str()))   return 3;
    if (!strcmp("P_exact", name.c_str()))   return 4;
    if (!strcmp("Div", name.c_str()))   return 5;
    if (!strcmp("SymTensorNorm", name.c_str()))   return 6;
    //    if (!strcmp("V_exactBC", name.c_str()))   return 5;
    
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

////////////////////////////////////////////////////////////////////

int TPZNavierStokesMaterial::NSolutionVariables(int var) {
    
    switch(var) {
            
        case 0:
            return 1; // Pressure, Scalar
        case 1:
            return 3; // Velocity, Vector
        case 2:
            return 3; // f, Vector
        case 3:
            return 3; // V_exact, Vector
        case 4:
            return 1; // P_exact, Scalar
        case 5:
            return 1; // Divergente
        case 6:
            return 1; // Symetric tensor norm
            
            //        case 5:
            //            return this->Dimension(); // V_exactBC, Vector
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) {
    
    
    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZManVector<STATE,3> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE,3> p_h = datavec[pindex].sol[0];
    
    TPZFNMatrix<9,STATE> gradu(3,1);
    
    // TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    // TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    TPZFMatrix<STATE> &dsol = datavec[vindex].dsol[0];
   // dsol.Resize(3,3);
    TPZFNMatrix<9,STATE> dsolxy(3,3),dsolxyp(3,1);
    dsolxy = dsol;
    if (fSpace!=1) {
        TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, datavec[vindex].axes);
    }

    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
            
        case 0: //Pressure
        {
            Solout[0] = p_h[0];
        }
            break;
            
        case 1: //Velocity
        {
            Solout[0] = v_h[0]; // Vx
            Solout[1] = v_h[1]; // Vy
            Solout[2] = v_h[2]; // Vz
        }
            break;
        case 2: //f
        {
            TPZVec<STATE> f(3,0.0);
            if(this->HasForcingFunction()){
                TPZVec<STATE> x(3,0.);
                x=datavec[vindex].x;
                this->ForcingFunction()->Execute(x, f, gradu);
                
            }
            
            
            Solout[0] = f[0]; // fx
            Solout[1] = f[1]; // fy
            Solout[2] = f[2]; // fz
        }
            break;
            
        case 3: //v_exact
        {
            TPZVec<STATE> sol(4,0.0);
            if(this->HasForcingFunctionExact()){
                TPZVec<STATE> x(3,0.);
                x=datavec[vindex].x;
                this->fForcingFunctionExact->Execute(x, sol, gradu); // @omar::check it!

            }
            Solout[0] = sol[0]; // vx
            Solout[1] = sol[1]; // vy
            Solout[2] = sol[2]; // vz
         
        }
            break;
            
        case 4: //p_exact
        {
            TPZVec<STATE> sol(4,0.0);
            if(this->HasForcingFunctionExact()){
                TPZVec<STATE> x(3,0.),xrot(3,0.);
                x=datavec[pindex].x;

                this->fForcingFunctionExact->Execute(x, sol, gradu); // @omar::check it!
            }
            Solout[0] = sol[3]; // px
            
        }
            break;
         
        case 5: //div
        {
            STATE Div=0.;
            for(int i=0; i<3; i++) {
                Div+=dsolxy(i,i);
            }
            Solout[0] = Div;
            
        }
            break;

        case 6: //norm of tensor
        {
            TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
            TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
            
            TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
            TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
            
            int nshapeV;
            nshapeV = datavec[vindex].fVecShapeIndex.NElements();
            
            int normvecRows = datavec[vindex].fNormalVec.Rows();
            int normvecCols = datavec[vindex].fNormalVec.Cols();
            TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
            TPZManVector<TPZFNMatrix<4,REAL>,18> GradNormalvec(18);
            
            STATE asd1 = 0., asd2 = 0.,asd3 = 0., asd4 = 0.;
            if (datavec[vindex].fNeedsNormalVecFad) {
#ifdef _AUTODIFF
                for (int e = 0; e < normvecRows; e++) {
                    for (int s = 0; s < normvecCols; s++) {
                        Normalvec(e,s)=datavec[vindex].fNormalVecFad(e,s).val();
                    }
                }
                
                for (int s = 0; s < normvecCols; s++) {
                    TPZFNMatrix<4,REAL> Grad0(3,3,0.); // 2x2
                    Grad0(0,0)=datavec[vindex].fNormalVecFad(0,s).fastAccessDx(0);
                    Grad0(0,1)=datavec[vindex].fNormalVecFad(0,s).fastAccessDx(1);
                    Grad0(1,0)=datavec[vindex].fNormalVecFad(1,s).fastAccessDx(0);
                    Grad0(1,1)=datavec[vindex].fNormalVecFad(1,s).fastAccessDx(1);
                    GradNormalvec[s] = Grad0;
                    //Grad0.Print(std::cout);
                }
            
#else
                DebugStop();
#endif
            }else{
                Normalvec=datavec[vindex].fNormalVec;
            }
            
            TPZFMatrix<STATE> phiVi(3,1,0.0),phiVj(3,1,0.0);
            TPZFNMatrix<4,STATE> GradSol(3,3,0.),GradSolt(3,3,0.),DuSol(3,3,0),GradVi(3,3,0.),GradVit(3,3,0.),Dui(3,3,0.);
            STATE normDu = 0.;
            STATE normDuSol = 0.;
            
            
            GradSol = datavec[vindex].dsol[vindex];
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    GradSolt(e,f) = GradSol(f,e);
                }
            }
            
            DuSol = GradSolt + GradSol;
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    normDuSol += DuSol(e,f)*DuSol(e,f);
                }
            }

            
            
            for(int i = 0; i < nshapeV; i++)
            {
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;

                for (int e=0; e<3; e++) {
                    for (int f=0; f<3; f++) {
                        GradVi(e,f) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
                        GradVit(f,e) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
                    }
                }
                
                
                
                for (int e=0; e<3; e++) {
                    for (int f=0; f<3; f++) {
                        Dui(e,f)= 0.5 * (GradVi(e,f) + GradVit(e,f));
                    }
                }
                
               
                for (int e=0; e<3; e++) {
                    for (int f=0; f<3; f++) {
                        normDu += Dui(e,f)*Dui(e,f);
                    }
                }

            }

//            GradVi.Print(std::cout);
//            Dui.Print(std::cout);
            
//            std::cout<<datavec[0].xParametric<<std::endl;
//            std::cout<<datavec[0].x<<std::endl;
//            Normalvec.Print(std::cout);
//            datavec[0].fNormalVecFad.Print(std::cout);
            //std::cout<<GradNormalvec<<std::endl;
            
            
            Solout[0] = normDuSol;
            
        }
            break;
            
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

////////////////////////////////////////////////////////////////////

// Divergence on master element
void TPZNavierStokesMaterial::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)

{
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFNMatrix<100,REAL> phiuH1        = datavec[ublock].phi;   // For H1  test functions Q
    TPZFNMatrix<300,REAL> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFNMatrix<300,REAL> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFNMatrix<9,STATE> gradu = datavec[ublock].dsol[0];
    TPZFNMatrix<9,STATE> graduMaster;
    gradu.Transpose();
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFNMatrix<9,REAL> Qaxes = datavec[ublock].axes;
    TPZFNMatrix<9,REAL> QaxesT;
    TPZFNMatrix<9,REAL> Jacobian = datavec[ublock].jacobian;
    TPZFNMatrix<9,REAL> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFNMatrix<9,REAL> GradOfX;
    TPZFNMatrix<9,REAL> GradOfXInverse;
    TPZFNMatrix<9,REAL> VectorOnMaster;
    TPZFNMatrix<9,REAL> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    TPZFMatrix<STATE> GradOfXInverseSTATE(GradOfXInverse.Rows(), GradOfXInverse.Cols());
    for (unsigned int i = 0; i < GradOfXInverse.Rows(); ++i) {
        for (unsigned int j = 0; j < GradOfXInverse.Cols(); ++j) {
            GradOfXInverseSTATE(i,j) = GradOfXInverse(i,j);
        }
    }
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            for (int k = 0; k < 3; k++) {
                VectorOnXYZ(k,0) = datavec[ublock].fNormalVec(k,ivectorindex);
            }
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            for (int k = 0; k < fDimension; k++) {
                DivergenceofPhi(iq,0) +=  dphiuH1(k,ishapeindex)*VectorOnMaster(k,0);
            }
            
        }
    }
    
    return;
    
}



////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::Write(TPZStream &buf, int withclassid) const{
    
    TPZMaterial::Write(buf, withclassid);
    
    
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::Read(TPZStream &buf, void *context) {
    
    TPZMaterial::Read(buf, context);
    
}

////////////////////////////////////////////////////////////////////

void TPZNavierStokesMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi){
    
    
    TPZFMatrix<REAL> &dphiV = dataV.dphix;
    
    const int dim = this->Dimension();
    
    GradPhi.clear();
    GradPhi.resize(dim);
    
    //for each shape
    for(int shape = 0; shape < dphiV.Rows(); shape++){
        
        TPZFMatrix<REAL> GPhi(dim,dim,0.);
        
        for(int i = 0; i < dim; i++){
            
            for(int j = 0; j < dim; j++){
                
                GPhi(i,j) = dphiV(j,shape);// itapopo H1 ??
                
            }//j
        }//i
        
        GradPhi[shape] = GPhi;
        
    }//shape
    
}

// Contricucao dos elementos internos
void TPZNavierStokesMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){

    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);

    return;
    
}


// Contricucao dos elementos internos
void TPZNavierStokesMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
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
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);

    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();

    int normvecRows = datavec[vindex].fNormalVec.Rows();
    int normvecCols = datavec[vindex].fNormalVec.Cols();
    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
    TPZManVector<TPZFNMatrix<9,REAL>,18> GradNormalvec(normvecCols);
    for (int i=0; i<normvecRows; i++) {
        GradNormalvec[i].Redim(2,2);
    }
    
    if (datavec[vindex].fNeedsNormalVecFad) {
#ifdef _AUTODIFF
        for (int e = 0; e < normvecRows; e++) {
            for (int s = 0; s < normvecCols; s++) {
                Normalvec(e,s)=datavec[vindex].fNormalVecFad(e,s).val();
            }
        }
        for (int s = 0; s < normvecCols; s++) {
            TPZFNMatrix<9,REAL> Grad0(3,3,0.); // 2x2
            Grad0(0,0)=datavec[vindex].fNormalVecFad(0,s).fastAccessDx(0);
            Grad0(0,1)=datavec[vindex].fNormalVecFad(0,s).fastAccessDx(1);
            Grad0(1,0)=datavec[vindex].fNormalVecFad(1,s).fastAccessDx(0);
            Grad0(1,1)=datavec[vindex].fNormalVecFad(1,s).fastAccessDx(1);
            GradNormalvec[s] = Grad0;
        }
#else
        DebugStop();
#endif
    }else{
        Normalvec=datavec[vindex].fNormalVec;
    }
    
    
    TPZVec<STATE> f(3,0.), f_rot(3,0.);
    for (int e=0; e<3; e++) {
        f[e] = 0.;
    }
    
    TPZFMatrix<STATE> phiVi(3,1,0.0),phiVj(3,1,0.0);

    TPZFNMatrix<100,STATE> divphi;
    TPZFNMatrix<40,STATE> div_on_master;
    TPZFNMatrix<10,STATE> gradUn = datavec[vindex].dsol[0];
    TPZManVector<STATE,3> u_n    = datavec[vindex].sol[0];
    STATE p_n                  = datavec[pindex].sol[0][0];
    
    if (fSpace==1) {
        datavec[0].ComputeFunctionDivergence();
    }
    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<9,STATE> GradVi(3,3,0.),GradVit(3,3,0.),Dui(3,3,0.);
        for (int e=0; e<3; e++) {
            phiVi(e,0) = phiV(iphi,0)*Normalvec(e,ivec);
            for (int f=0; f<3; f++) {
                GradVi(e,f) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
                GradVit(f,e) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
            }
        }
        for (int e=0; e<3; e++) {
            for (int f=0; f<3; f++) {
                Dui(e,f)= 0.5 * (GradVi(e,f) + GradVit(e,f));
            }
        }

        STATE divui = 0.;
        divui = Tr( GradVi );
        
        if(this->HasForcingFunction()){
            TPZFMatrix<STATE> gradu;
            TPZVec<STATE> x(3,0.),xrot(3,0.);
            x=datavec[vindex].x;
            this->ForcingFunction()->Execute(x, f, gradu);
        }
        
        STATE phi_dot_f = 0.0, un_dot_phiV = 0.0; // f - Source term
        for (int e=0; e<3; e++) {
            phi_dot_f += phiVi(e)*f[e];
            un_dot_phiV += phiVi(e)*u_n[e];
        }
        ef(i) += weight * (phi_dot_f-un_dot_phiV*0.);
    
        STATE A_term_f = 0.; // A - Flux term
        TPZFNMatrix<9,STATE> DUn_j(3,3,0.);
        for (int e=0; e<3; e++) {
            for (int f=0; f<3; f++) {
                 DUn_j(e,f)= 0.5 * (gradUn(e,f) + gradUn(f,e));
            }
        }
        A_term_f = Inner(Dui, DUn_j);
        ef(i) += 2. * weight * (-A_term_f);
        
        STATE B_term_f = 0.; // B - Mixed term
        B_term_f = - p_n * datavec[0].divphi(i);
        ef(i) += weight * (-B_term_f);
        
        STATE C_term_f = 0.; // C - Trilinear terms
        TPZFNMatrix<9,STATE> GradUn_phiU(3,1,0.);
        for (int e=0; e<3; e++) {
            for (int f=0; f<3; f++) {
                GradUn_phiU(e,0) += gradUn(e,f)*u_n[f];
            }
        }
        C_term_f += InnerVec(GradUn_phiU, phiVi);
        ef(i) += -weight * C_term_f*0.;
        
        // A, C e D - velocity X velocity
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            
            for (int e=0; e<3; e++) {
                phiVj(e,0) = phiV(jphi,0)*Normalvec(e,jvec);
            }
            
            TPZFNMatrix<9,STATE> GradVj(3,3,0.),GradVjt(3,3,0.),Duj(3,3,0.);
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    GradVj(e,f) = Normalvec(e,jvec)*dphiVx(f,jphi)+GradNormalvec[jvec](e,f)*phiV(jphi,0);
                    GradVjt(f,e) = Normalvec(e,jvec)*dphiVx(f,jphi)+GradNormalvec[jvec](e,f)*phiV(jphi,0);
                }
            }        
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    Duj(e,f)= 0.5 * (GradVj(e,f) + GradVjt(e,f));
                }
            }
            STATE A_term = Inner(Dui, Duj);
            
            STATE C_term = 0.;
            
            TPZFNMatrix<9,STATE> GradU_Un(3,1,0.);
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    GradU_Un(e,0) += GradVj(e,f)*u_n[f];
                }
            }
    
            C_term = InnerVec(GradU_Un, phiVi);

            STATE D_term = 0.;
            
            TPZFNMatrix<9,STATE> GradUn_phiU(3,1,0.);
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    GradUn_phiU(e,0) += gradUn(e,f)*phiVj(f,0);
                }
            }
            
            D_term = InnerVec(GradUn_phiU, phiVi);
            
            ek(i,j) += 2. * weight * fViscosity * A_term;  // A - Bilinear gradV * gradU
        
            ek(i,j) += weight * C_term*0.;  // C - Trilinear terms
            
            ek(i,j) += weight * D_term*0.;  // D - Trilinear terms
            
        }
        

        
        
        // B - pressure and velocity
        for (int j = 0; j < nshapeP; j++) {
            
            TPZManVector<REAL,3> GradPj(3,0.);
            for (int e=0; e<3; e++) {
                GradPj[e] = dphiPx(e,j);
            }

            ///p*div(U)
            STATE fact = 0.;
            if (fSpace==1) {
                fact = (-1.) * weight * phiP(j,0) * datavec[0].divphi(i);
            }else{
                fact = (-1.) * weight * phiP(j,0) * divui;
            }

            // Matrix B
            ek(i, nshapeV+j) += fact;
            // Matrix B^T
            ek(nshapeV+j,i) += fact;
            
        }
        
    }
    
    for (int i = 0; i < nshapeP; i++) {
 
        STATE B_term_f = 0.; // B - Mixed term
        B_term_f = - phiP(i,0)*Tr(gradUn);
        ef(i) += weight * (-B_term_f);
        
    }
    
    
    // Preparação para formulação MHM :
    if (datavec.size()>2) {
        
        TPZFMatrix<REAL> &phigM = datavec[2].phi;
        TPZFMatrix<REAL> &phipM = datavec[3].phi;
        
        // matrix C - pressure and average-pressure
        for (int j = 0; j < nshapeP; j++) {
            
            STATE fact = (1.) * weight * phiP(j,0) * phigM(0,0);
            // Matrix C
            ek(nshapeV+nshapeP, nshapeV+j) += fact;
            // Matrix C^T
            ek(nshapeV+j,nshapeV+nshapeP) += fact;
            
        }
        
        // matrix D - injection and average-pressure

        STATE factG = (1.) * weight * phigM(0,0) * phipM(0,0);
        // Matrix C
        ek(nshapeV+nshapeP+1, nshapeV+nshapeP) += factG;
        // Matrix C^T
        ek(nshapeV+nshapeP,nshapeV+nshapeP+1) += factG;
        
    }

    if (datavec.size()>4) {
        
        TPZFMatrix<REAL> &phigM0 = datavec[4].phi;
        TPZFMatrix<REAL> &phipM0 = datavec[5].phi;
        
        // matrix C0 - pressure and average-pressure
        for (int j = 0; j < nshapeP; j++) {
            
            STATE fact0 = (1.) * weight * phiP(j,0) * phigM0(0,0);
            // Matrix C
            ek(nshapeV+nshapeP+2, nshapeV+j) += fact0;
            // Matrix C^T
            ek(nshapeV+j,nshapeV+nshapeP+2) += fact0;
            
        }
        
        // matrix D0 - injection and average-pressure
        
        STATE factG0 = (1.) * weight * phigM0(0,0) * phipM0(0,0);
        // Matrix C
        ek(nshapeV+nshapeP+3, nshapeV+nshapeP+2) += factG0;
        // Matrix C^T
        ek(nshapeV+nshapeP+2,nshapeV+nshapeP+3) += factG0;
        
    }
    

#ifdef PZDEBUG
    {
        std::ofstream fileEK("FileEKContribute.txt");
        ek.Print("stiff = ",fileEK,EMathematicaInput);
        
        std::ofstream fileEF("FileEFContribute.txt");
        ef.Print("rhs = ",fileEF,EMathematicaInput);
    }
#endif
    
}


void TPZNavierStokesMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
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
    
    // Getting the linear combination or finite element approximations
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
 //   nshapeV = phiV.Rows()*NStateVariables();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    int normvecRows = datavec[vindex].fNormalVec.Rows();
    int normvecCols = datavec[vindex].fNormalVec.Cols();
    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
    
//    if (datavec[vindex].fNeedsNormalVecFad==false) {
        Normalvec=datavec[vindex].fNormalVec;
//    }else{
//        for (int e = 0; e < normvecRows; e++) {
//            for (int s = 0; s < normvecCols; s++) {
//                Normalvec(e,s)=datavec[vindex].fNormalVecFad(e,s).val();
//            }
//        }
//    }
    
//    Normalvec.Print(std::cout);
    
    if (fSpace==1) {
        nshapeV = nshapeV/2.;
    }


    int gy=v_h.size();
    
    TPZFNMatrix<9,STATE> phiVi(fDimension,1,0.),phiVni(1,1,0.), phiVj(fDimension,1,0.),phiVnj(1,1,0.), phiPi(fDimension,1),phiPj(fDimension,1);
    
    TPZFNMatrix<3,STATE> v_2=bc.Val2();
    TPZFNMatrix<3,STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    switch (bc.Type()) {
        case 0: //Dirichlet for continuous formulation
        {
            
            TPZFMatrix<STATE> gradu(3,3,0.);
            TPZManVector<STATE> vbc(4,0.);
            TPZFMatrix<STATE> Du(3,3,0.),Dun(3,1,0.);

            if(bc.HasBCForcingFunction())
            {
                TPZManVector<STATE> vbc(4,0.);
                TPZFMatrix<STATE> gradu(3,3,0.);
                bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                v_2(2,0) = vbc[2];
                p_D = vbc[3];
            }
            
            if(fSpace==1){

                
                for(int i = 0; i < nshapeV; i++ )
                {

                    //Adaptação para Hdiv

                    TPZManVector<REAL> n = datavec[0].normal;

                    REAL vh_n = v_h[0];
                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];

                    ef(i,0) += -weight * gBigNumber * (vh_n - v_n) * phiV(i,0);

                    for(int j = 0; j < nshapeV; j++){

                        ek(i,j) += weight * gBigNumber * phiV(j,0) * phiV(i,0);
                        
                    }
                    


                }
                
                
                
            }else{
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=Normalvec(e,ivec)*phiV(iphi,0);
                    }
                    
                    
                    //Adaptação para Hdiv
                    
                    STATE factef=0.0;
                    for(int is=0; is<gy ; is++){
                        factef += -1.0*(v_h[is] - v_2(is,0)) * phiVi(is,0);
                    }
                    
                    ef(i,0) += weight * gBigNumber * factef;
                    
                    for(int j = 0; j < nshapeV; j++){
                        int jphi = datavec[vindex].fVecShapeIndex[j].second;
                        int jvec = datavec[vindex].fVecShapeIndex[j].first;
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVj(e,0)=Normalvec(e,jvec)*phiV(jphi,0);
                        }
                        
                        //Adaptação para Hdiv
                        
                        STATE factek = 0.0;
                        for(int is=0; is<gy ; is++){
                            factek += phiVj(is,0) * phiVi(is,0);
                        }
                        
                        ek(i,j) += weight * gBigNumber * factek;
                        
                    }
                    
                }
            }
            
        }
            break;
            
        case 1: //Neumann for continuous formulation
        {
            
            
            if(bc.HasBCForcingFunction())
            {
                TPZManVector<STATE> vbc(4,0.);
                TPZFMatrix<STATE> gradu(3,3,0.);
                bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                v_2(2,0) = vbc[2];
                p_D = vbc[3];
            }
            
            
            for(int i = 0; i < nshapeV; i++ )
            {
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=Normalvec(e,ivec)*phiV(iphi,0);
                }
                
                TPZManVector<REAL> n = datavec[vindex].normal;
                
                TPZFNMatrix<9,STATE> pn(fDimension,1);
                
                
                for (int f=0; f<fDimension; f++) {
                    pn(f,0)=n[f]*v_1(0,0);
                }
                
                //Adaptação para Hdiv
                
                STATE factef=0.0;
                for(int is=0; is<gy ; is++){
                    factef += (pn(is,0))* phiVi(is,0);
                }
                
                ef(i,0) += weight * factef;
                
            }
            
            
        }
            
            
            
            break;
            
            
        case 2: //Condição Penetração
        {
            
            if(bc.HasBCForcingFunction())
            {
                TPZManVector<STATE> vbc(4,0.);
                TPZFMatrix<STATE> gradu(3,3,0.);
                bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                v_2(2,0) = vbc[2];
                p_D = vbc[3];
            }

            
            if(fSpace==1){

                
                for(int i = 0; i < nshapeV; i++ )
                {
                    
                    //Adaptação para Hdiv
                    
                    TPZManVector<REAL> n = datavec[0].normal;
                    
                    REAL vh_n = v_h[0];
                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                    
                    ef(i,0) += -weight * gBigNumber * (vh_n - v_n) * phiV(i,0);
                    
                    for(int j = 0; j < nshapeV; j++){
                        
                        ek(i,j) += weight * gBigNumber * phiV(j,0) * phiV(i,0);
                        
                    }
                    
                }
                
                
            }else{
                
                
                
                
                TPZManVector<REAL> n = datavec[0].normal;
                TPZManVector<REAL> t(2);
                t[0]=-n[1];  //oioioio
                t[1]=n[0];
                
                
                
                //Componente normal -> imposta fortemente:
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    
                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
                    TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                    
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=Normalvec(e,ivec)*datavec[vindex].phi(iphi,0);
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
                            phiVj(e,0)=Normalvec(e,jvec)*datavec[vindex].phi(jphi,0);
                            phiVnj(0,0)+=phiVj(e,0)*n[e];
                            phiVtj(0,0)+=phiVj(e,0)*t[e];
                            
                        }
                        
                        ek(i,j) += weight * gBigNumber * phiVni(0,0) * phiVnj(0,0);
                        
                    }
                    
                }
                
                //Componente tangencial -> imposta fortemente:

//                for(int i = 0; i < nshapeV; i++ )
//                {
//
//                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
//                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
//                    TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
//
//
//                    for (int e=0; e<fDimension; e++) {
//                        phiVi(e,0)=Normalvec(e,ivec)*datavec[vindex].phi(iphi,0);
//                        phiVni(0,0)+=phiVi(e,0)*n[e];
//                        phiVti(0,0)+=phiVi(e,0)*t[e];
//                    }
//
//                    REAL vh_t = v_h[1];
//                    REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
//
//                    ef(i,0) += -weight * gBigNumber * (vh_t-v_t) * (phiVti(0,0));
//
//
//                    for(int j = 0; j < nshapeV; j++){
//
//                        int jphi = datavec[vindex].fVecShapeIndex[j].second;
//                        int jvec = datavec[vindex].fVecShapeIndex[j].first;
//
//                        TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
//
//                        for (int e=0; e<fDimension; e++) {
//                            phiVj(e,0)=Normalvec(e,jvec)*datavec[vindex].phi(jphi,0);
//                            phiVnj(0,0)+=phiVj(e,0)*n[e];
//                            phiVtj(0,0)+=phiVj(e,0)*t[e];
//
//                        }
//
//                        ek(i,j) += weight * gBigNumber * phiVti(0,0) * phiVtj(0,0);
//
//                    }
//
//                }

                
                
            }
            
        }
            break;
            
            
        case 3: //Contribuicao ponto no x
        {
            
            REAL p_D = v_2(0,0);
            
            if(bc.HasBCForcingFunction())
            {
                TPZManVector<STATE> pbc(1);
                bc.BCForcingFunction()->Execute(datavec[vindex].x,pbc);
                p_D = pbc[0];
                
            }
            
            TPZManVector<REAL> n = datavec[0].normal;
            
            REAL phiVi_n;
            
            
            for(int i = 0; i < 3; i++ )
            {
                phiVi_n = 0.0;
                
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=Normalvec(e,ivec)*phiV(iphi,0);
                    phiVi_n += phiVi(e,0)*n[e];
                }
                
                ef(i*2,0) += -1.0*weight * p_D * phiVi_n;
                
                
            }
            
            
        }
            break;
            
        case 4: //Contribuicao ponto no y
        {
            
            REAL p_D = v_2(0,0);
            
            if(bc.HasBCForcingFunction())
            {
                TPZManVector<STATE> pbc(1);
                bc.BCForcingFunction()->Execute(datavec[vindex].x,pbc);
                p_D = pbc[0];
                
            }
            
            TPZManVector<REAL> n = datavec[0].normal;
            
            REAL phiVi_n;
            
            
            for(int i = 0; i < 3; i++ )
            {
                phiVi_n = 0.0;
                
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=Normalvec(e,ivec)*phiV(iphi,0);
                    phiVi_n += phiVi(e,0)*n[e];
                }
                
                ef(i*2+1,0) += -1.0*weight * p_D * phiVi_n;
                
                
            }
            
            
        }
            break;
            
        case 5: //Ponto pressao
        {
           
            //return;
            p_D = bc.Val2()(0,0);
            
            
            for(int i = 0; i < nshapeP; i++ )
            {
                
                
                ef(i) += 1.0 * p_D * phiP(i,0);
                
                for(int j = 0; j < nshapeP; j++){
                    
                    ek(i,j) += 1.0 * (phiP(i,0) * phiP(j,0));
                    
                }
                
            }
            
        }
            break;
            
        case 6: //Pressao Dirichlet
        {
            
            if(bc.HasBCForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                bc.BCForcingFunction()->Execute(datavec[pindex].x,vbc);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D  = vbc[2]*0.;
                
                
            }
            
            //pressao
            
            for(int i = 0; i < nshapeP; i++ )
            {
                
                
                ef(i) += -weight * gBigNumber * (p_h[0] - p_D) * phiP(i,0);
                
                for(int j = 0; j < nshapeP; j++){
                    
                    ek(i,j) += weight * gBigNumber * (phiP(i,0) * phiP(j,0));
                    
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
    
    
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }

    
#ifdef PZDEBUG
    {
        std::ofstream fileEK("FileEKContributeBC.txt");
        std::ofstream fileEF("FileEFContributeBC.txt");
        ek.Print("stiff = ",fileEK,EMathematicaInput);
        ef.Print("force = ",fileEF,EMathematicaInput);
    }
#endif

}

////////////////////////////////////////////////////////////////////
template <typename TVar>
TVar TPZNavierStokesMaterial::Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
    
    //inner product of two tensors
    
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    TVar Val = 0;
    
    for(int i = 0; i < S.Cols(); i++){
        for(int j = 0; j < S.Cols(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}


////////////////////////////////////////////////////////////////////
STATE TPZNavierStokesMaterial::InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
    //inner product of two vectors
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0;
    
    for(int j = 0; j < S.Cols(); j++){
        for(int i = 0; i < S.Rows(); i++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}



////////////////////////////////////////////////////////////////////

STATE TPZNavierStokesMaterial::Tr( TPZFMatrix<REAL> &GradU ){
    
#ifdef DEBUG
    if( GradU.Rows() != GradU.Cols() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0.;
    
    for(int i = 0; i < GradU.Rows(); i++){
        Val += GradU(i,i);
    }
    
    return Val;
}


/// transform a H1 data structure to a vector data structure
void TPZNavierStokesMaterial::FillVecShapeIndex(TPZMaterialData &data)
{
    data.fNormalVec.Resize(fDimension,fDimension);
    data.fNormalVec.Identity();
    data.fVecShapeIndex.Resize(fDimension*data.phi.Rows());
    for (int d=0; d<fDimension; d++) {
        for (int i=0; i<data.phi.Rows(); i++) {
            data.fVecShapeIndex[i*fDimension+d].first = d;
            data.fVecShapeIndex[i*fDimension+d].second = i;
        }
    }
}



void TPZNavierStokesMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &sol_exact, TPZFMatrix<STATE> &dsol_exact, TPZVec<REAL> &errors)
{
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    TPZManVector<STATE> Velocity(3,0.), Pressure(3,0.);
    
    this->Solution(data,VariableIndex("V"), Velocity);
    this->Solution(data,VariableIndex("P"), Pressure);
    
    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZFMatrix<REAL> dudx(3,3);
    TPZFMatrix<STATE> &dsol = data[vindex].dsol[0];
    TPZFMatrix<STATE> &dsolp = data[pindex].dsol[0];
    //std::cout<<dsol<<std::endl;
    
    //Adaptação feita para Hdiv
    dsol.Resize(3,3);
    
    TPZFNMatrix<2,STATE> dsolxy(3,3,0.), dsolxyp(3,1,0.);
    dsolxy = dsol;
    if (fSpace!=1) {
        TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, data[vindex].axes);
    }
   // TPZAxesTools<STATE>::Axes2XYZ(dsolp, dsolxyp, data[pindex].axes);
    
    dsolxyp = dsolp;
    
//    std::cout<<Velocity<<std::endl;
//    std::cout<<sol_exact<<std::endl;
    
    int shift = 3;
    // velocity
    
    //values[2] : erro norma L2
    STATE diff, diffp;
    errors[1] = 0.;
    for(int i=0; i<3; i++) {
        diff = Velocity[i] - sol_exact[i];
        errors[1]  += diff*diff;
    }
    
    ////////////////////////////////////////////////// H1 / GD
    
    if(fSpace==2){
        
//        //values[2] : erro em semi norma H1
//        errors[2] = 0.;
//        TPZFMatrix<STATE> S(Dimension(),Dimension(),0.0);
//        for(int i=0; i<Dimension(); i++) {
//            for(int j=0; j<Dimension(); j++) {
//                S(i,j) = dsolxy(i,j) - du_exact(i,j);
//            }
//        }
//
//        diff = Inner(S, S);
//        errors[2]  += diff;
//
//        //values[0] : erro em norma H1 <=> norma Energia
//        errors[0]  = errors[1]+errors[2];

        /// erro norma HDiv
        
        STATE Div_exact=0., Div=0.;
        for(int i=0; i<3; i++) {
            Div_exact+=dsol_exact(i,i);
            Div+=dsolxy(i,i);
        }
        
        diff = Div-Div_exact;
        
        errors[2]  = diff*diff;
        
        // errors[0]  = errors[1]+errors[2];
        
    }
    
    if(fSpace==3){
        
//        //values[2] : erro em semi norma H1
//        errors[2] = 0.;
//        TPZFMatrix<STATE> S(Dimension(),Dimension(),0.0);
//        for(int i=0; i<Dimension(); i++) {
//            for(int j=0; j<Dimension(); j++) {
//                S(i,j) = dsolxy(i,j) - du_exact(i,j);
//            }
//        }
//
//        diff = Inner(S, S);
//        errors[2]  += diff;
//
//        //values[0] : erro em norma H1 <=> norma Energia
//        errors[0]  = errors[1]+errors[2];

        /// erro norma HDiv
        
        STATE Div_exact=0., Div=0.;
        for(int i=0; i<3; i++) {
            Div_exact+=dsol_exact(i,i);
            Div+=dsolxy(i,i);
        }
        
        diff = Div-Div_exact;
        
        errors[2]  = diff*diff;
        
   //     errors[0]  = errors[1]+errors[2];
        
    }
    
    
    ////////////////////////////////////////////////// H1 / GD
    
    // pressure
    
    /// values[1] : eror em norma L2
    diffp = Pressure[0] - sol_exact[3];
    errors[shift+1]  = diffp*diffp;
    

    
    ////////////////////////////////////////////////// HDIV
    
    if(fSpace==1){
        /// erro norma HDiv
        
        STATE Div_exact=0., Div=0.;
        for(int i=0; i<3; i++) {
            Div_exact+=dsol_exact(i,i);
            Div+=dsolxy(i,i);
        }
        
        diff = Div-Div_exact;
        
        errors[2]  = diff*diff;
        
    //    errors[0]  = errors[1]+errors[2];
        
    }
    
    ////////////////////////////////////////////////// HDIV
    
}


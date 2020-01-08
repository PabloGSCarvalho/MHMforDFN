/*
 *  TPZBrinkmanMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZBrinkmanMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"
#include "pzlog.h"

using namespace std;


TPZBrinkmanMaterial::TPZBrinkmanMaterial() : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(){
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1;
    
}

////////////////////////////////////////////////////////////////////

TPZBrinkmanMaterial::TPZBrinkmanMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma) : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(matid),fDimension(dimension),fSpace(space),fViscosity(viscosity),fTheta(theta),fSigma(Sigma)
{
    // symmetric version
    //fTheta = -1;
    
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1.;
    

    
}

////////////////////////////////////////////////////////////////////

TPZBrinkmanMaterial::TPZBrinkmanMaterial(const TPZBrinkmanMaterial &mat) : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(mat),fDimension(mat.fDimension),fSpace(mat.fSpace), fViscosity(mat.fViscosity), fTheta(mat.fTheta), fSigma(mat.fSigma)
{
    fk= mat.fk;
    
}

////////////////////////////////////////////////////////////////////

TPZBrinkmanMaterial::~TPZBrinkmanMaterial(){
    
    
}

////////////////////////////////////////////////////////////////////

void TPZBrinkmanMaterial::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
}


////////////////////////////////////////////////////////////////////

void TPZBrinkmanMaterial::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
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

void TPZBrinkmanMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
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

void TPZBrinkmanMaterial::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data){

    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    data.fNeedsNeighborSol = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

////////////////////////////////////////////////////////////////////

void TPZBrinkmanMaterial::FillDataRequirementsInterface(TPZMaterialData &data)
{
    data.fNeedsNormal = true;
    data.fNeedsNeighborCenter = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

////////////////////////////////////////////////////////////////////

void TPZBrinkmanMaterial::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
{
    int nref_left = datavec_left.size();
    for(int iref = 0; iref<nref_left; iref++){
        datavec_left[iref].SetAllRequirements(false);
        datavec_left[iref].fNeedsNormal = true;
    }
    datavec_left[0].fNeedsNormalVecFad = NeedsNormalVecFad;
}

////////////////////////////////////////////////////////////////////

void TPZBrinkmanMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZBrinkmanMaterial::VariableIndex(const std::string &name) {
    
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

int TPZBrinkmanMaterial::NSolutionVariables(int var) {
    
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

void TPZBrinkmanMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) {
    
    
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
            TPZFNMatrix<10,STATE> gradV_axes = datavec[vindex].dsol[0];
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
void TPZBrinkmanMaterial::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)

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

void TPZBrinkmanMaterial::Write(TPZStream &buf, int withclassid) const{
    
    TPZMaterial::Write(buf, withclassid);
    
    
}

////////////////////////////////////////////////////////////////////

void TPZBrinkmanMaterial::Read(TPZStream &buf, void *context) {
    
    TPZMaterial::Read(buf, context);
    
}

////////////////////////////////////////////////////////////////////

void TPZBrinkmanMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi){
    
    
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

void TPZBrinkmanMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
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
    
//    Normalvec.Print(std::cout);
//    std::cout<<GradNormalvec<<std::endl;
    
    TPZVec<STATE> f(3,0.), f_rot(3,0.);
    for (int e=0; e<3; e++) {
        f[e] = 0.;
    }
    
    TPZFMatrix<STATE> phiVi(3,1,0.0),phiVj(3,1,0.0);

    TPZFNMatrix<100,STATE> divphi;
    TPZFNMatrix<40,STATE> div_on_master;
    TPZFNMatrix<10,STATE> gradV_axes = datavec[vindex].dsol[0];
    
    STATE jac_det;
    if (fSpace==1) {
        datavec[0].ComputeFunctionDivergence();
//        this->ComputeDivergenceOnMaster(datavec, div_on_master);
    }
 
  //  phiV.Print(std::cout);
    
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

//        Normalvec.Print(std::cout);
//        dphiVx.Print(std::cout);
//       GradVi.Print(std::cout);
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<3; e++) {
            for (int f=0; f<3; f++) {
                Dui(e,f)= 0.5 * (GradVi(e,f) + GradVit(e,f));
            }
        }
        
        //Dui.Print(std::cout);
        //Divergente (Calculo em H1)
        STATE divui = 0.;
        divui = Tr( GradVi );
        
        
        if(this->HasForcingFunction()){
            TPZFMatrix<STATE> gradu;
            TPZVec<STATE> x(3,0.),xrot(3,0.);
            x=datavec[vindex].x;
            this->ForcingFunction()->Execute(x, f, gradu);
        }
        
        STATE phi_dot_f = 0.0;
        for (int e=0; e<3; e++) {
            phi_dot_f += phiVi(e)*f[e];
        }
        
        ef(i) += weight * phi_dot_f;
        
        // matrix A - gradV
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
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    Duj(e,f)= 0.5 * (GradVj(e,f) + GradVjt(e,f));
                }
            }
        
            //Duj.Print(std::cout);
            
            STATE val = Inner(Dui, Duj);
            STATE val1 = InnerVec(phiVi, phiVj);

            ek(i,j) += 2. * weight * fViscosity * val + 0.*weight * val1; ///Visc*(GradU+GradU^T):GradPhi
        }
        
        // matrix B - pressure and velocity
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
        STATE factf= -weight * phiP(i,0)*f[2]*0.;
        ef(nshapeV+i,0) += factf;
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
    }
#endif
    
}


void TPZBrinkmanMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
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

void TPZBrinkmanMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    
    // Verificar que
    // os termos mistos devem estar sem viscosidade!
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    // Setting forcing function
    /*STATE force = 0.;
     if(this->fForcingFunction) {
     TPZManVector<STATE> res(1);
     fForcingFunction->Execute(datavec[pindex].x,res);
     force = res[0];
     }*/
    
    //Gravity
    STATE rhoi = 900.; //itapopo
    STATE g = 9.81; //itapopo
    STATE force = rhoi*g;
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    // V - right
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    // P - right
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    data.fNeedsNormal = true;
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    //Detjac
    REAL Detjac=fabs(data.detjac);
    
    REAL sigmaConst = fSigma;
    
    //Triangle Verification:
    if (fabs(normal[0])<1.&&fabs(normal[1])<1.) {
        sigmaConst = fSigma/(sqrt(2.));
    }
    

    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);

    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);

    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);

    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);

    
    TPZManVector<REAL, 3> tangent(fDimension,0.);
    TPZFNMatrix<3,STATE> tangentV(fDimension,1,0.);
    for(int i=0; i<fDimension; i++) tangent[i] = data.axes(0,i);
    for(int i=0; i<fDimension; i++) tangentV(i,0) = data.axes(0,i);
    
    
    //TPZManVector<REAL,3> normalx(fDimension,phiP2.Cols());
    //TPZAxesTools<REAL>::Axes2XYZ(normal, normalx, data.axes);
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    //NormalVec Left
    int normvecRowsL = datavecleft[vindex].fNormalVec.Rows();
    int normvecColsL = datavecleft[vindex].fNormalVec.Cols();
    TPZFNMatrix<3,REAL> NormalvecLeft(normvecRowsL,normvecColsL,0.);
    TPZManVector<TPZFNMatrix<4,REAL>,18> GradNormalvecLeft(18);
    
    if (datavecleft[vindex].fNeedsNormalVecFad){
#ifdef _AUTODIFF
        for (int e = 0; e < normvecRowsL; e++) {
            for (int s = 0; s < normvecColsL; s++) {
                NormalvecLeft(e,s)=datavecleft[vindex].fNormalVecFad(e,s).val();
            }
        }
        
        for (int s = 0; s < normvecColsL; s++) {
            TPZFNMatrix<4,REAL> Grad0(3,3,0.); // 2x2
            Grad0(0,0)=datavecleft[vindex].fNormalVecFad(0,s).fastAccessDx(0);
            Grad0(0,1)=datavecleft[vindex].fNormalVecFad(0,s).fastAccessDx(1);
            Grad0(1,0)=datavecleft[vindex].fNormalVecFad(1,s).fastAccessDx(0);
            Grad0(1,1)=datavecleft[vindex].fNormalVecFad(1,s).fastAccessDx(1);
            GradNormalvecLeft[s] = Grad0;
        }
        
#else
        DebugStop();
#endif
    }else{
        NormalvecLeft=datavecleft[vindex].fNormalVec;
    }
    
    //NormalVec Right
    int normvecRowsR = datavecright[vindex].fNormalVec.Rows();
    int normvecColsR = datavecright[vindex].fNormalVec.Cols();
    TPZFNMatrix<3,REAL> NormalvecRight(normvecRowsR,normvecColsR,0.);
    TPZManVector<TPZFNMatrix<4,REAL>,18> GradNormalvecRight(18);
    
    if (datavecright[vindex].fNeedsNormalVecFad) {
#ifdef _AUTODIFF
        for (int e = 0; e < normvecRowsR; e++) {
            for (int s = 0; s < normvecColsR; s++) {
                NormalvecRight(e,s)=datavecright[vindex].fNormalVecFad(e,s).val();

            }
        }
        for (int s = 0; s < normvecColsR; s++) {
            TPZFNMatrix<4,REAL> Grad0(3,3,0.); // 2x2
            Grad0(0,0)=datavecright[vindex].fNormalVecFad(0,s).fastAccessDx(0);
            Grad0(0,1)=datavecright[vindex].fNormalVecFad(0,s).fastAccessDx(1);
            Grad0(1,0)=datavecright[vindex].fNormalVecFad(1,s).fastAccessDx(0);
            Grad0(1,1)=datavecright[vindex].fNormalVecFad(1,s).fastAccessDx(1);
            GradNormalvecRight[s] = Grad0;
        }
        
#else
        DebugStop();
#endif
    }else{
        NormalvecRight=datavecright[vindex].fNormalVec;
    }
    
    for(int i1 = 0; i1 < nshapeV1; i1++ )
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        TPZFNMatrix<9,STATE> GradV1ni(fDimension,1,0.),phiV1i(fDimension,1),phiV1ni(1,1,0.), phiV1ti(fDimension,1,0.);
        TPZFNMatrix<4,STATE> GradV1i(fDimension,fDimension,0.),GradV1it(fDimension,fDimension,0.),Du1i(fDimension,fDimension,0.),Du1ni(fDimension,1,0.),  Du1ti(fDimension,1,0.);
        STATE phiit = 0.;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=NormalvecLeft(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV1i(e,f) = NormalvecLeft(e,ivec1)*dphiVx1(f,iphi1);
                //termo transposto:
                GradV1it(f,e) = NormalvecLeft(e,ivec1)*dphiVx1(f,iphi1);
            }
        }
        
//                dphiV1.Print("dphiV = ",cout);
//                dphiVx1.Print("dphiVx = ",cout);
//                datavecleft[vindex].fNormalVec.Print("normalvec = ",cout);
//                GradV1i.Print("GradV1i = ",cout);
//                phiV1i.Print("phiV1i = ",cout);
        
        phiit = InnerVec(phiV1i,tangentV);
        for(int e = 0; e<fDimension; e++)
        {
            phiV1ti(e,0) += phiit*tangent[e];
        }
                
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du1i(e,f)= (1./2.) * (GradV1i(e,f) + GradV1it(e,f));
            }
        }
        
        //Du1ni e Du1ti
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du1ni(e,0) += Du1i(e,f)*normal[f];
                Du1ti(e,0) += Du1i(e,f)*tangent[f];
            }
        }
        
        TPZFNMatrix<9,STATE> GradV1nj(fDimension,1,0.),phiV1j(fDimension,1),phiV1nj(1,1,0.);
        
        // K11 - (trial V left) * (test V left)
        for(int j1 = 0; j1 < nshapeV1; j1++){
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            
            TPZFNMatrix<3,STATE> phiV1j(fDimension,1,0.),phiV1tj(fDimension,1,0.), phiV1nj(1,1,0.);
            TPZFNMatrix<4,STATE> GradV1j(fDimension,fDimension,0.),GradV1jt(fDimension,fDimension,0.),Du1j(fDimension,fDimension,0.),Du1nj(fDimension,1,0.),Du1tj(fDimension,1,0.);
            STATE phijt = 0.;
            
            for (int e=0; e<fDimension; e++) {
                
                phiV1j(e,0)=NormalvecLeft(e,jvec1)*datavecleft[vindex].phi(jphi1,0);
                phiV1nj(0,0)+=phiV1j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV1j(e,f) = NormalvecLeft(e,jvec1)*dphiVx1(f,jphi1);
                    //termo transposto:
                    GradV1jt(f,e) = NormalvecLeft(e,jvec1)*dphiVx1(f,jphi1);
                    
                }
            }
            
            phijt = InnerVec(phiV1j,tangentV);
            for(int e = 0; e<fDimension; e++)
            {
                phiV1tj(e,0) += phijt*tangent[e];
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1j(e,f)= (1./2.) * (GradV1j(e,f) + GradV1jt(e,f));
                }
            }
            
            //Du1nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1nj(e,0) += Du1j(e,f)*normal[f];
                    Du1tj(e,0) += Du1j(e,f)*tangent[f];
                }
            }
            
            // 1-1 : Inter : phiVi, Dunj
            
            STATE fact = (-1./2.) * weight * 2.* fViscosity * InnerVec(phiV1i, Du1nj);
            
            ek(i1,j1) +=fact;
            ek(j1,i1) +=-fact*fTheta;
            
            
            // 1-1 : Penalidade : phiVi, phiVj
            
            STATE penalty = sigmaConst * weight * fViscosity * InnerVec(phiV1i, phiV1j);
            ek(i1,j1) +=penalty;
            
        }
        
        // K12 e K21 - (trial V left) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            
            TPZFNMatrix<9,STATE> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            // 1-1 : Pressure : phiVni, phiPj
            
            STATE fact = (1./2.) * weight * Inner(phiV1ni,phiP1j);
            
            ek(i1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i1) += fact;
            
        }
        
        
        // K13 - (trial V left) * (test V right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9,STATE> GradV2nj(fDimension,1),phiV2j(fDimension,1),phiV2nj(1,1,0.);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            TPZFNMatrix<4,STATE> GradV2j(fDimension,fDimension,0.),GradV2jt(fDimension,fDimension,0.),Du2j(fDimension,fDimension,0.),Du2nj(fDimension,1,0.);
            
            for (int e=0; e<fDimension; e++) {
                
                phiV2j(e,0)=NormalvecRight(e,jvec2)*datavecright[vindex].phi(jphi2,0);
                phiV2nj(0,0)+=phiV2j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV2j(e,f) = NormalvecRight(e,jvec2)*dphiVx2(f,jphi2);
                    //termo transposto:
                    GradV2jt(f,e) = NormalvecRight(e,jvec2)*dphiVx2(f,jphi2);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2j(e,f)= (1./2.) * (GradV2j(e,f) + GradV2jt(e,f));
                }
            }
            
            //Du2nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2nj(e,0) += Du2j(e,f)*normal[f] ;
                }
            }
            
            // 1-2 : Inter : phiVi, Dunj
            
            STATE fact = (-1./2.) * weight * 2. * fViscosity * InnerVec(phiV1i,Du2nj);
            
            ek(i1,j2+nshapeV1+nshapeP1) += fact;
            ek(j2+nshapeV1+nshapeP1,i1) += -fact*fTheta;
            
            // 1-2 : Penalidade : phiVi, phiVj
            
            STATE penalty = sigmaConst * weight * fViscosity * InnerVec(phiV1i, phiV2j);
            ek(i1,j2+nshapeV1+nshapeP1) += -penalty;
            
        }
        
        // K14 e K41 - (trial V left) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9,STATE> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            // 1-2 : Pressure : phiVni, phiPj
            
            STATE fact = (1./2.) * weight * InnerVec(phiV1ni,phiP2j);
            
            ek(i1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i1) += fact;
            
        }
        
    }
    
    
    for(int i2 = 0; i2 < nshapeV2; i2++ ){
        
        TPZFNMatrix<9,STATE> GradV2ni(fDimension,1),phiV2i(fDimension,1),phiV2ni(1,1,0.);
        TPZFNMatrix<4,STATE> GradV2i(fDimension,fDimension,0.),GradV2it(fDimension,fDimension,0.),Du2i(fDimension,fDimension,0.),Du2ni(fDimension,1,0.);
        
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV2i(e,0)=NormalvecRight(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV2i(e,f) = NormalvecRight(e,ivec2)*dphiVx2(f,iphi2);
                //termo transposto:
                GradV2it(f,e) = NormalvecRight(e,ivec2)*dphiVx2(f,iphi2);
            }
        }
        
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du2i(e,f)= (1./2.) * (GradV2i(e,f) + GradV2it(e,f));
            }
        }
        
        //Du2ni
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du2ni(e,0) += Du2i(e,f)*normal[f] ;
            }
        }
        
        
        
        // K31 - (trial V right) * (test V left)
        for(int j1 = 0; j1 < nshapeV1; j1++){
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            
            TPZFNMatrix<4,STATE> GradV1j(fDimension,fDimension,0.),GradV1jt(fDimension,fDimension,0.),Du1j(fDimension,fDimension,0.),Du1nj(fDimension,1,0.);
            
            TPZFNMatrix<9,STATE> phiV1j(fDimension,1),phiV1nj(1,1,0.);
            
            for (int e=0; e<fDimension; e++) {
                
                phiV1j(e,0)=NormalvecLeft(e,jvec1)*datavecleft[vindex].phi(jphi1,0);
                phiV1nj(0,0)+=phiV1j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV1j(e,f) = NormalvecLeft(e,jvec1)*dphiVx1(f,jphi1);
                    //termo transposto:
                    GradV1jt(f,e) = NormalvecLeft(e,jvec1)*dphiVx1(f,jphi1);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1j(e,f)= (1./2.) * (GradV1j(e,f) + GradV1jt(e,f));
                }
            }
            
            //Du1nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1nj(e,0) += Du1j(e,f)*normal[f] ;
                }
            }
            
            // 2-1 : Inter : phiVi, Dunj
            
            STATE fact = (1./2.) * weight * 2. * fViscosity * InnerVec(phiV2i, Du1nj);
            
            ek(i2+nshapeV1+nshapeP1,j1) += fact;
            ek(j1,i2+nshapeV1+nshapeP1) += -fact*fTheta;
            
            // 2-1 : Penalidade : phiVi, phiVj
            
            STATE penalty = sigmaConst * weight * fViscosity * InnerVec(phiV2i, phiV1j);
            ek(i2+nshapeV1+nshapeP1,j1) += -penalty;
            
            
        }
        
        // K32 e K23 - (trial V right) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            TPZFNMatrix<9,STATE> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            // 2-1 : Pressure : phiVni, phiPj
            
            STATE fact = (-1./2.) * weight * InnerVec(phiV2ni,phiP1j);
            
            ek(i2+nshapeV1+nshapeP1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i2+nshapeV1+nshapeP1) += fact;
            
        }
        
        
        // K33 - (trial V right) * (test V right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9,STATE> GradV2nj(fDimension,1);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            TPZFNMatrix<4,STATE> GradV2j(fDimension,fDimension,0.),GradV2jt(fDimension,fDimension,0.),Du2j(fDimension,fDimension,0.),Du2nj(fDimension,1,0.);
            
            TPZFNMatrix<9,STATE> phiV2j(fDimension,1),phiV2nj(1,1,0.);
            
            
            for (int e=0; e<fDimension; e++) {
                
                phiV2j(e,0)=NormalvecRight(e,jvec2)*datavecright[vindex].phi(jphi2,0);
                phiV2nj(0,0)+=phiV2j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV2j(e,f) = NormalvecRight(e,jvec2)*dphiVx2(f,jphi2);
                    //termo transposto:
                    GradV2jt(f,e) = NormalvecRight(e,jvec2)*dphiVx2(f,jphi2);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2j(e,f)= (1./2.) * (GradV2j(e,f) + GradV2jt(e,f));
                }
            }
            
            //Du2nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2nj(e,0) += Du2j(e,f)*normal[f] ;
                }
            }
            
            // 2-2 : Inter : phiVi, Dunj
            
            STATE fact = (1./2.) * weight *  2. * fViscosity * InnerVec(phiV2i,Du2nj);
            
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact;
            ek(j2+nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += -fact*fTheta;
            
            // 2-2 : Penalidade : phiVi, phiVj
            
            STATE penalty = sigmaConst * weight * fViscosity * InnerVec(phiV2i, phiV2j);
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) +=penalty;
            
            
        }
        
        // K34 e K43- (trial V right) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9,STATE> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            // 2-2 : Pressure : phiVni, phiPj
            
            STATE fact = (-1./2.) * weight * InnerVec(phiV2ni,phiP2j);
            
            ek(i2+nshapeV1+nshapeP1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact;
        }
        
    }
    
#ifdef PZDEBUG
    std::ofstream plotfileM("ekInterface.txt");
    ek.Print("Kint = ",plotfileM,EMathematicaInput);
#endif
    
}


void TPZBrinkmanMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
   

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
    
    int normvecRows = datavec[vindex].fNormalVec.Rows();
    int normvecCols = datavec[vindex].fNormalVec.Cols();
    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
    TPZManVector<TPZFNMatrix<4,REAL>,18> GradNormalvec(18);
    
    if (datavec[vindex].fNeedsNormalVecFad) {
#ifdef _AUTODIFF
        for (int e = 0; e < normvecRows; e++) {
            for (int s = 0; s < normvecCols; s++) {
                //datavec[vindex].fNormalVecFad.Print(std::cout);
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
        }
        
#else
        DebugStop();
#endif
    }else{
        Normalvec=datavec[vindex].fNormalVec;
    }
    
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
                TPZManVector<STATE> vbc(4,0.);
                TPZFMatrix<STATE> gradu(3,3,0.);
                bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                v_2(2,0) = vbc[2];
                p_D = vbc[3];
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
                        GradVi(e,f) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
                        GradVit(f,e) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
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
                    phiVi(e,0)=Normalvec(e,ivec)*datavec[vindex].phi(iphi,0);
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
                            phiVj(e,0)=Normalvec(e,jvec)*datavec[vindex].phi(jphi,0);
                        }
                        
                        
                        phiVtj(0,0)= t[0] * phiVj(0,0) + t[1] * phiVj(1,0);
                        
                        
                        
                        TPZFNMatrix<4,STATE> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);
                        
                        for (int e=0; e<fDimension; e++) {
                            
                            for (int f=0; f<fDimension; f++) {
                                GradVj(e,f) = Normalvec(e,jvec)*dphiVx(f,jphi)+GradNormalvec[jvec](e,f)*phiV(jphi,0);
                                GradVjt(f,e) = Normalvec(e,jvec)*dphiVx(f,jphi)+GradNormalvec[jvec](e,f)*phiV(jphi,0);
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
                    if(0){
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
                    TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=Normalvec(e,ivec)*phiV(iphi,0);
                    }
                    
                    TPZManVector<REAL> n = data.normal;
                    TPZManVector<REAL> n2 = datavec[0].normal;
                    
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
                    TPZManVector<STATE> vbc(4,0.);
                    TPZFMatrix<STATE> gradu(3,3,0.);
                    bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                    v_2(0,0) = vbc[0];
                    v_2(1,0) = vbc[1];
                    v_2(2,0) = vbc[2];
                    p_D = vbc[3];
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
                                GradVi(e,f) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
                                GradVit(f,e) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
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
                            phiVi(e,0)=Normalvec(e,ivec)*datavec[vindex].phi(iphi,0);
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
                    TPZManVector<STATE> vbc(4,0.);
                    TPZFMatrix<STATE> gradu(3,3,0.);
                    bc.BCForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                    v_2(0,0) = vbc[0];
                    v_2(1,0) = vbc[1];
                    v_2(2,0) = vbc[2];
                    p_D = vbc[3];
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
                            GradVi(e,f) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
                            GradVit(f,e) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
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
                        phiVi(e,0)=Normalvec(e,ivec)*datavec[vindex].phi(iphi,0);
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
                                phiVj(e,0)=Normalvec(e,jvec)*datavec[vindex].phi(jphi,0);
                            }
                            
                            phiVnj(0,0)= n[0] * phiVj(0,0) + n[1] * phiVj(1,0);
                            phiVtj(0,0)= t[0] * phiVj(0,0) + t[1] * phiVj(1,0);
                            
                            TPZFNMatrix<4,STATE> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);
                            for (int e=0; e<fDimension; e++) {
                                
                                for (int f=0; f<fDimension; f++) {
                                    GradVj(e,f) = Normalvec(e,jvec)*dphiVx(f,jphi)+GradNormalvec[jvec](e,f)*phiV(jphi,0);
                                    GradVjt(f,e) = Normalvec(e,jvec)*dphiVx(f,jphi)+GradNormalvec[jvec](e,f)*phiV(jphi,0);
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



////////////////////////////////////////////////////////////////////
template <typename TVar>
TVar TPZBrinkmanMaterial::Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
    
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
STATE TPZBrinkmanMaterial::InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
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

STATE TPZBrinkmanMaterial::Tr( TPZFMatrix<REAL> &GradU ){
    
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
void TPZBrinkmanMaterial::FillVecShapeIndex(TPZMaterialData &data)
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



void TPZBrinkmanMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &sol_exact, TPZFMatrix<STATE> &dsol_exact, TPZVec<REAL> &errors)
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


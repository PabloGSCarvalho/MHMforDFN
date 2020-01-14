/*
 *  TPZDarcyDFNMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZDarcyDFNMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"
#include "pzlog.h"

using namespace std;


TPZDarcyDFNMaterial::TPZDarcyDFNMaterial() : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(){
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1;
    
}

////////////////////////////////////////////////////////////////////

TPZDarcyDFNMaterial::TPZDarcyDFNMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma) : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(matid),fDimension(dimension),fSpace(space),fViscosity(viscosity),fTheta(theta),fSigma(Sigma)
{
    // symmetric version
    //fTheta = -1;
    
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1.;
    

}

////////////////////////////////////////////////////////////////////

TPZDarcyDFNMaterial::TPZDarcyDFNMaterial(const TPZDarcyDFNMaterial &mat) : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(mat),fDimension(mat.fDimension),fSpace(mat.fSpace), fViscosity(mat.fViscosity), fTheta(mat.fTheta), fSigma(mat.fSigma)
{
    fk= mat.fk;
    
}

////////////////////////////////////////////////////////////////////

TPZDarcyDFNMaterial::~TPZDarcyDFNMaterial(){
    
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
}


////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
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

void TPZDarcyDFNMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
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

void TPZDarcyDFNMaterial::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data){

    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    data.fNeedsNeighborSol = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::FillDataRequirementsInterface(TPZMaterialData &data)
{
    data.fNeedsNormal = true;
    data.fNeedsNeighborCenter = true;
    data.fNeedsNormalVecFad = NeedsNormalVecFad;
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
{
    int nref_left = datavec_left.size();
    for(int iref = 0; iref<nref_left; iref++){
        datavec_left[iref].SetAllRequirements(false);
        datavec_left[iref].fNeedsNormal = true;
    }
    datavec_left[0].fNeedsNormalVecFad = NeedsNormalVecFad;
}

////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZDarcyDFNMaterial::VariableIndex(const std::string &name) {
    
    if (!strcmp("P", name.c_str()))  return 0;
    if (!strcmp("Pressure", name.c_str()))  return 0;
    if (!strcmp("f", name.c_str()))         return 1;
    if (!strcmp("P_exact", name.c_str()))   return 2;
    if (!strcmp("Div", name.c_str()))   return 3;
    
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

////////////////////////////////////////////////////////////////////

int TPZDarcyDFNMaterial::NSolutionVariables(int var) {
    
    switch(var) {
            
        case 0:
            return 1; // Pressure, Scalar
        case 1:
            return 1; // Force, Scalar
        case 2:
            return 1; // P_exact, Scalar
        case 3:
            return 1; // Div, Scalar

        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) {
    
    TPZManVector<STATE,3> p_h = data.sol[0];
    TPZFMatrix<STATE> &dsol = data.dsol[0];
    
    TPZFNMatrix<9,STATE> gradu(3,1);

    TPZFNMatrix<9,STATE> dsolxy(3,3),dsolxyp(3,1);
    dsolxy = dsol;
    if (fSpace!=1) {
        TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, data.axes);
    }

    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
            
        case 0: //Pressure
        {
            Solout[0] = p_h[0];
        }
            break;
            
        case 1: //f
        {
            TPZVec<STATE> f(1,0.0);
            if(this->HasForcingFunction()){
                TPZVec<STATE> x(3,0.);
                x=data.x;
                this->ForcingFunction()->Execute(x, f, gradu);
            }
            
            Solout[0] = f[0]; // fx
         }
            break;
            
        case 2: //p_exact
        {
            TPZVec<STATE> sol(1,0.0);
            if(this->HasForcingFunctionExact()){
                TPZVec<STATE> x(3,0.),xrot(3,0.);
                x=data.x;

                this->fForcingFunctionExact->Execute(x, sol, gradu); // @omar::check it!
            }
            Solout[0] = sol[0]; // px
            
        }
            break;
         
        case 3: //div
        {
            STATE Div=0.;
            for(int i=0; i<3; i++) {
                Div+=dsolxy(i,i);
            }
            Solout[0] = Div;
            
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

void TPZDarcyDFNMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) {
    
    int pindex = PIndex();
    TPZManVector<STATE,3> p_h = datavec[pindex].sol[0];
    TPZFMatrix<STATE> &dsol = datavec[pindex].dsol[0];
    
    TPZFNMatrix<9,STATE> gradu(3,1);
    
    TPZFNMatrix<9,STATE> dsolxy(3,3),dsolxyp(3,1);
    dsolxy = dsol;
    if (fSpace!=1) {
        TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, datavec[pindex].axes);
    }
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
            
        case 0: //Pressure
        {
            Solout[0] = p_h[0];
        }
            break;
            
        case 1: //f
        {
            TPZVec<STATE> f(1,0.0);
            if(this->HasForcingFunction()){
                TPZVec<STATE> x(3,0.);
                x=datavec[pindex].x;
                this->ForcingFunction()->Execute(x, f, gradu);
            }
            
            Solout[0] = f[0]; // fx
        }
            break;
            
        case 2: //p_exact
        {
            TPZVec<STATE> sol(1,0.0);
            if(this->HasForcingFunctionExact()){
                TPZVec<STATE> x(3,0.),xrot(3,0.);
                x=datavec[pindex].x;
                
                this->fForcingFunctionExact->Execute(x, sol, gradu); // @omar::check it!
            }
            Solout[0] = sol[0]; // px
            
        }
            break;
            
        case 3: //div
        {
            STATE Div=0.;
            for(int i=0; i<3; i++) {
                Div+=dsolxy(i,i);
            }
            Solout[0] = Div;
            
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
void TPZDarcyDFNMaterial::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)

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

void TPZDarcyDFNMaterial::Write(TPZStream &buf, int withclassid) const{
    
    TPZMaterial::Write(buf, withclassid);
    
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::Read(TPZStream &buf, void *context) {
    
    TPZMaterial::Read(buf, context);
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyDFNMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi){
    
    
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

void TPZDarcyDFNMaterial::Contribute(TPZMaterialData &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    if (datavec.fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec);
    }
    
    // P
    TPZFMatrix<REAL> &phiP = datavec.phi;
    TPZFMatrix<REAL> &dphiP = datavec.dphix;
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec.axes);
    
    int nshapeP = 0;
    nshapeP = phiP.Rows();

    // Grad pressure * Grad pressure
    for (int i = 0; i < nshapeP; i++) {
            
        TPZFNMatrix<100,STATE> GradPi(3,1,0.);
        for (int e=0; e<3; e++) {
            GradPi(e,0) = dphiPx(e,i);
        }

        STATE termA = 0.;
        for (int j = 0; j < nshapeP; j++) {
            TPZFNMatrix<100,STATE> GradPj(3,1,0.);
            for (int e=0; e<3; e++) {
                GradPj(e,0) = dphiPx(e,j);
            }
        
            termA =  weight * InnerVec(GradPi, GradPj);
            ek(i, j) += termA;
        
        }

        STATE termF = 0.0;
        TPZFMatrix<STATE> gradu;
        TPZVec<STATE> f(1,0.);
        
        if(this->HasForcingFunction()){
            this->ForcingFunction()->Execute(datavec.x, f, gradu);
        }
        termF = weight * f[0]*phiP(i,0);
        ef(i,0) += termF;
    }
    
    
    

#ifdef PZDEBUG
    {
       std::ofstream fileEK("FileEKContribute.txt");
        ek.Print("stiff = ",fileEK,EMathematicaInput);
    }
#endif
    
}


////////////////////////////////////////////////////////////////////


void TPZDarcyDFNMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int pindex = PIndex();
    int Lindex = LambdaIndex();
    int Apindex = AveragePIndex();
    
    if (datavec[pindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[pindex]);
    }

    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeP = 0, nshapeL = 0, nshapeAp = 0;
    nshapeP = phiP.Rows();
    nshapeL = datavec[Lindex].phi.Rows();
    nshapeAp = datavec[Apindex].phi.Rows();
    
    // Grad pressure * Grad pressure
    for (int i = 0; i < nshapeP; i++) {
        
        TPZFNMatrix<100,STATE> GradPi(3,1,0.);
        for (int e=0; e<3; e++) {
            GradPi(e,0) = dphiPx(e,i);
        }
        
        STATE termA = 0.;
        for (int j = 0; j < nshapeP; j++) {
            TPZFNMatrix<100,STATE> GradPj(3,1,0.);
            for (int e=0; e<3; e++) {
                GradPj(e,0) = dphiPx(e,j);
            }
            
            termA =  weight * InnerVec(GradPi, GradPj);
            ek(i, j) += termA;
            
        }
        
        STATE termF = 0.0;
        TPZFMatrix<STATE> gradu;
        TPZVec<STATE> f(1,0.);
        
        if(this->HasForcingFunction()){
            this->ForcingFunction()->Execute(datavec[pindex].x, f, gradu);
        }
        termF = weight * f[0]*phiP(i,0);
        ef(i,0) += termF;
    }
    
    
    
    // Preparação para formulação MHM :
    if (datavec.size()>2) {

        TPZFMatrix<REAL> &phiLM = datavec[1].phi;

        // matrix - pressure and average-pressure
        for (int j = 0; j < nshapeP; j++) {

            STATE fact = (1.) * weight * phiP(j,0) * phiLM(0,0);
            // Matrix C
            ek(nshapeP, j) += fact;
            // Matrix C^T
            ek(j,nshapeP) += fact;

        }
        
        TPZFMatrix<REAL> &phipM = datavec[2].phi;
    
        // matrix - pressure and average-pressure
        for (int j = 0; j < nshapeL; j++) {
    
            STATE fact = (-1.) * weight * phiLM(0,0) * phipM(0,0);
                // Matrix C
            ek(nshapeP+nshapeL,nshapeP+j) += fact;
                // Matrix C^T
            ek(nshapeP+j,nshapeP+nshapeL) += fact;
    
        }

        
    
    }
    
    
    
#ifdef PZDEBUG
    {
        std::ofstream fileEK("FileEKContribute.txt");
        ek.Print("stiff = ",fileEK,EMathematicaInput);
    }
#endif
    
}

////////////////////////////////////////////////////////////////////
template <typename TVar>
TVar TPZDarcyDFNMaterial::Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
    
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
STATE TPZDarcyDFNMaterial::InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
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

STATE TPZDarcyDFNMaterial::Tr( TPZFMatrix<REAL> &GradU ){
    
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
void TPZDarcyDFNMaterial::FillVecShapeIndex(TPZMaterialData &data)
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



void TPZDarcyDFNMaterial::Errors(TPZMaterialData &data, TPZVec<STATE> &sol_exact, TPZFMatrix<STATE> &dsol_exact, TPZVec<REAL> &errors)
{
    DebugStop();
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    TPZManVector<STATE> Velocity(1,0.), Pressure(1,0.);
    
    this->Solution(data,VariableIndex("P"), Pressure);
    
    TPZFMatrix<REAL> dudx(3,3);
    TPZFMatrix<STATE> &dsol = data.dsol[0];
    
    TPZFNMatrix<2,STATE> dsolxy(3,3,0.), dsolxyp(3,1,0.);
    dsolxy = dsol;
    if (fSpace!=1) {
        TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, data.axes);
    }

    STATE diff_sol=0.,diff_dsol=0.;
    
    // Pressure : eror L2
    diff_sol = Pressure[0] - sol_exact[3];
    errors[0]  = diff_sol*diff_sol;

    
    // Grad Pressure : : eror L2
    for(int i=0; i<3; i++) {
        diff_dsol = dsol(i,0) - dsol_exact(i,0);
        errors[1]  += diff_dsol*diff_dsol;
    }
    
}

void TPZDarcyDFNMaterial::Errors(TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
       TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
       TPZVec<STATE> &uexact, TPZFMatrix<STATE> &duexact,
                                 TPZVec<REAL> &errors){
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    STATE diff_sol=0.,diff_dsol=0.;
    
    // Pressure : eror L2
    diff_sol = sol[0] - uexact[0];
    errors[0]  = diff_sol*diff_sol;

    // Grad Pressure : : eror L2
    for(int i=0; i<2; i++) {
        diff_dsol = dsol(i,0) - duexact(i,0);
        errors[1]  += diff_dsol*diff_dsol;
    }
    
}

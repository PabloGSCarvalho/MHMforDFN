//
//  TPZMHMixedMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMDarcyDFNMeshControl.h"
#include "TPZMHMeshControl.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
#include "TPZMatLaplacian.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZCompElLagrange.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZFractureInsertion.h"
#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"

#include "pzmat1dlin.h"
#include "pzinterpolationspace.h"
#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"
#include "pzl2projection.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mhmdarcyDFNmeshcontrol"));
#endif

TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices), fFractureMatId(),  fFractureFlowDim1MatId(), fBCLagrangeFluxMatIds()
{
    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCMeshLagrangeLocal = new TPZCompMesh(fGMesh);
    fCMeshConstantPressureLocal = new TPZCompMesh(fGMesh);
    
}

TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(int dimension) : TPZMHMeshControl(dimension), fFractureMatId(10.), fFractureFlowDim1MatId(), fBCLagrangeFluxMatIds() {
    
    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCMeshLagrangeLocal = new TPZCompMesh(fGMesh);
    fCMeshConstantPressureLocal = new TPZCompMesh(fGMesh);
    fmatFracPointL = -7;
    fmatFracPointR = -8;
}

TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMeshControl(gmesh), fFractureMatId(10.), fFractureFlowDim1MatId(), fBCLagrangeFluxMatIds() {
    
    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCMeshLagrangeLocal = new TPZCompMesh(fGMesh);
    fCMeshConstantPressureLocal = new TPZCompMesh(fGMesh);
    fmatFracPointL = -7;
    fmatFracPointR = -8;
}


TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(const TPZMHMDarcyDFNMeshControl &copy) : TPZMHMeshControl(copy){
    
    this->operator=(copy);
}

TPZMHMDarcyDFNMeshControl &TPZMHMDarcyDFNMeshControl::operator=(const TPZMHMDarcyDFNMeshControl &cp){
    
    fDistrFluxMesh = cp.fDistrFluxMesh;
    fAveragePressMesh = cp.fAveragePressMesh;
    fCMeshLagrangeLocal = cp.fCMeshLagrangeLocal;
    fCMeshConstantPressureLocal = cp.fCMeshConstantPressureLocal;

    return *this;
}


TPZMHMDarcyDFNMeshControl::~TPZMHMDarcyDFNMeshControl()
{
    
}

void TPZMHMDarcyDFNMeshControl::BuildComputationalMesh(bool usersubstructure)
{

    int nstate = fNState;
    int dim = fGMesh->Dimension();
    fCMesh->SetDimModel(dim);

    fFractureFlowDim1MatId.insert(fFractureMatId);
    
    InsertH1Fracture(fCMesh);

    
    InsertPeriferalMaterialObjects();
    CreateInternalElements();
    
    CreateFractureBC();
    CreateH1Wrappers();    
    CreateHDivWrappers();


    //DivideSkeletonElements(0);
    //InsertBCMultipliers();
    //InsertInternalMultipliers();
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo_test.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk1("MalhaGeo_test.vtk"); //Impressão da malha geométrica (formato vtk)
    fGMesh->Print(fileg1);
    TPZVTKGeoMesh::PrintGMeshVTK(fGMesh, filegvtk1,true);
#endif
    
//    InsertPeriferalLagrangeFluxObjects();
//    CreateLagrangeFluxElements();
    CreateSkeleton();
//    CreateSkeletonElements2();
    
    CreateInterfaceElements();
    CreateFractureInterfaces();
    
    
    //    AddBoundaryInterfaceElements();
    fCMesh->ExpandSolution();
    fCMesh->CleanUpUnconnectedNodes();
    
    
    // Distributed Flux and Average pressure for internal elements
    
    //SetFracSubDomain(fFractureMatId);
    
    fCMeshLagrangeLocal = new TPZCompMesh(fGMesh);
    this->CreateLagrangeMultiplierMeshLocal(fGMesh->Dimension());
//    SetFracSubDomain(fFractureMatId);
    this->CreateLagrangeMultiplierMeshLocal(fGMesh->Dimension()-1);
    
    
//    InsertDistributedFluxMaterialObjects();
//    CreateDistributedFluxMHMMesh();
//    InsertPeriferalAveragePressMaterialObjects();
//    CreateAveragePressMHMMesh();
    

    if (fLagrangeAveragePressure) {
        fCMeshLagrange = new TPZCompMesh(fGMesh);
        this->CreateLagrangeMultiplierMesh(fGMesh->Dimension());
        
        this->CreateLagrangeMultiplierMesh(fGMesh->Dimension()-1);
//
//        std::ofstream fileCData("ConnectsData.txt");
//        this->Print(fileCData);

        std::ofstream filec2("MalhaC_Cmesh.txt");
        fCMesh->Print(filec2);
        
        this->TransferToMultiphysics();
        

        std::ofstream fileCData("ConnectsData.txt");
        this->Print(fileCData);
        
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "*********** BEFORE SUBSTRUCTURING *************\n";
        fCMesh->Print(sout);
        int64_t nel = fCMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMesh->Element(el);
            if(!cel) continue;
            sout << "el index " << el << " subdomain " << WhichSubdomain(cel) << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();
    if(usersubstructure==true){
        this->SubStructure();
    }
    fNumeq = fCMesh->NEquations();
}



void TPZMHMDarcyDFNMeshControl::InsertDistributedFluxMaterialObjects(){
    
    TPZCompMesh *cmeshDistributedFlux = fDistrFluxMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshDistributedFlux->MaterialVec().find(matid) == cmeshDistributedFlux->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshDistributedFlux->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
}


void TPZMHMDarcyDFNMeshControl::CreateDistributedFluxMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshDistributedFlux = fDistrFluxMesh.operator->();
    gmesh->ResetReference();
    cmeshDistributedFlux->SetName("DistributedFluxMesh");
    cmeshDistributedFlux->SetDimModel(gmesh->Dimension());
    cmeshDistributedFlux->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshDistributedFlux->SetDefaultOrder(porder);
    
    int meshdim = cmeshDistributedFlux->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshDistributedFlux->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshDistributedFlux->AutoBuild(matids);
    fDistrFluxMesh->ExpandSolution();
    
    
    int64_t nel = fDistrFluxMesh->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmeshDistributedFlux->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) DebugStop();
//        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
//        {
//            DebugStop();
//        }
        celdisc->SetTotalOrderShape();
        celdisc->SetFalseUseQsiEta();
    }
    
    int64_t nc = cmeshDistributedFlux->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshDistributedFlux->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshDistributedFlux->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            continue;
        }
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
    }
    
    
    if(1)
    {
        std::ofstream out("DistributedFluxMesh.txt");
        fDistrFluxMesh->Print(out);
    }
    
}

void TPZMHMDarcyDFNMeshControl::InsertPeriferalAveragePressMaterialObjects(){
    
    TPZCompMesh *cmeshAverPressure = fAveragePressMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshAverPressure->MaterialVec().find(matid) == cmeshAverPressure->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshAverPressure->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
    
}

void TPZMHMDarcyDFNMeshControl::CreateAveragePressMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshAverPressute = fAveragePressMesh.operator->();
    gmesh->ResetReference();
    cmeshAverPressute->SetName("AveragePressureMesh");
    cmeshAverPressute->SetDimModel(gmesh->Dimension());
    cmeshAverPressute->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshAverPressute->SetDefaultOrder(porder);
    
    int meshdim = cmeshAverPressute->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshAverPressute->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshAverPressute->AutoBuild(matids);
    fAveragePressMesh->ExpandSolution();
    
    
    int64_t nel = fAveragePressMesh->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmeshAverPressute->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) DebugStop();
        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
        {
            DebugStop();
        }
        celdisc->SetTotalOrderShape();
        celdisc->SetFalseUseQsiEta();
    }
    
    int64_t nc = cmeshAverPressute->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshAverPressute->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshAverPressute->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            continue;
        }
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
    }
    
    if(1)
    {
        std::ofstream out("AveragePressureMesh.txt");
        fAveragePressMesh->Print(out);
    }
    
    
}







void TPZMHMDarcyDFNMeshControl::Print(std::ostream & out) const {
    
    //ComputeNodElCon();
//    fCMesh->LoadReferences();
//    fGMesh->Reference()->LoadReferences();

    out << "number of connects            = " << fCMesh->NConnects() << std::endl;
    out << "number of elements            = " << fCMesh->NElements() << std::endl;

    out << "\n\t Connect Information:\n\n";
    int64_t i, nelem = 0.;
//
//    int64_t nel = fGMesh->NElements();
//    for (int64_t el=0; el<nel; el++) {
//        TPZGeoEl *gel = fGMesh->Element(el);
//        if (gel->MaterialId()!=fFractureMatId) {
//            continue;
//        }
//
//        out << "\n\t Computable Element Information:\n\n";
//
//        TPZCompEl *cel = gel->Reference();
//        if(!cel){
//            out << " Comp El is NULL" << std::endl;
//            out << " Geo Index = " << gel->Index();
//            out << " - MatID  " << gel->MaterialId() << std::endl;;
//
//            continue;
//        }else{
//
//            out << " Geo Index = " << cel->Reference()->Index();
//            out << " - Comp Index =  = " << cel->Index();
//            out << " - MatID = " << cel->Reference()->MaterialId();
//            out << " - Center coordinate: ";
//            for (int i = 0; i < cel->Reference()->NCornerNodes(); i++) {
//                TPZVec< REAL > center(3, 0.);
//                for (int j = 0; j < 3; j++) center[j] = cel->Reference()->NodePtr(i)->Coord(j);
//                out << "[" << i << "]" << center << " ";
//            }
//            out << std::endl;
//            int nconects = cel->NConnects();
//            out << " Number of connects = " << cel->NConnects() << " Connect indexes : ";
//            int nod;
//            for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
//                TPZConnect &c = cel->Connect(nod);
//                out << cel->ConnectIndex(nod) << ' ';
//            }
//            out << std::endl;
//            out << " Connect / number of connected elem : ";
//            for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
//                TPZConnect &c = cel->Connect(nod);
//                out << cel->ConnectIndex(nod) << "/" << c.NElConnected()  << "  " ;
//            }
//
//            out << std::endl;
//
//
//        }
//
//        out << " ------- Neighbours by sides ------- " << std::endl;
//        int nsides = gel->NSides();
//        out << " Nsides : " << nsides << std::endl;
//
//        out << " ----------------------------------- " << std::endl;
//
//        for (int side = 0 ; side < nsides; side++) {
//
//            if (gel->SideDimension(side) != fGMesh->Dimension()-1) {
//                       continue;
//            }
//
//            TPZGeoElSide gelside(gel,side);
//            TPZGeoElSide neighbour = gelside.Neighbour();
//            while(neighbour != gelside)
//            {
//                TPZCompEl *neigh = neighbour.Element()->Reference();
//
//                out << " Neighbour by side : " << side << std::endl;
//
//                if(!neigh){
//                    out << " Neighbour is NULL" << std::endl;
//                    out << " Geo Index = " << neighbour.Element()->Index();
//                    out << " - MatID  " << neighbour.Element()->MaterialId() << std::endl;
//
//                }else{
//
//                    out << " Geo Index = " << neigh->Reference()->Index();
//                    out << " - Comp Index = " << neigh->Index();
//                    out << " - MatID = " << neigh->Reference()->MaterialId();
//                    out << " - Center coordinate: ";
//                    for (int i = 0; i < neigh->Reference()->NCornerNodes(); i++) {
//                        TPZVec< REAL > center(3, 0.);
//                        for (int j = 0; j < 3; j++) center[j] = neigh->Reference()->NodePtr(i)->Coord(j);
//                        out << "[" << i << "]" << center << " ";
//                    }
//                    out << std::endl;
//                    int nconects = neigh->NConnects();
//                    out << " Number of connects = " << neigh->NConnects() << " Connect indexes : ";
//                    int nod;
//                    for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
//                        TPZConnect &c = neigh->Connect(nod);
//                        out << neigh->ConnectIndex(nod) << ' ';
//                    }
//                    out << std::endl;
//                    out << " Connect / number of connected elem : ";
//                    for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
//                        TPZConnect &c = neigh->Connect(nod);
//                        out << neigh->ConnectIndex(nod) << "/" << c.NElConnected()  << "  " ;
//                    }
//                    out << std::endl;
//                }
//
//                out << " ----------------------------------- " << std::endl;
//                neighbour = neighbour.Neighbour();
//            }
//        }
//            out << endl;
//            out << endl;
//
//    }

    out << "\n\t Computable Element Information:\n\n";
    nelem = fCMesh->NElements();
    for(i=0; i<nelem; i++) {
        if(!fCMesh->Element(i)) continue;
        TPZCompEl *cel = fCMesh->Element(i);
        out << "CompEl Index = " << cel->Index();
        out << " - GeoEl Index = " << cel->Reference()->Index();
        out << " - MatID = " << cel->Reference()->MaterialId();
        out << " - Center coordinate: ";
        for (int i = 0; i < cel->Reference()->NCornerNodes(); i++) {
            TPZVec< REAL > center(3, 0.);
            for (int j = 0; j < 3; j++) center[j] = cel->Reference()->NodePtr(i)->Coord(j);
            out << "[" << i << "]" << center << " ";
        }
        out << std::endl;
        int nconects = cel->NConnects();
        out << "Number of connects = " << cel->NConnects() << " Connect indexes : ";
        int nod;
        for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
            TPZConnect &c = cel->Connect(nod);
            out << cel->ConnectIndex(nod) << ' ';
        }
        out << std::endl;
        out << "Connect / number of connected elem : ";
        for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
            TPZConnect &c = cel->Connect(nod);
            out << cel->ConnectIndex(nod) << "/" << c.NElConnected()  << "  " ;
        }

        out << endl;
        out << endl;
    }

}

//Insert fracture and break H1 connectivity
void TPZMHMDarcyDFNMeshControl::InsertH1Fracture(TPZCompMesh &cmesh)
{
    
    for(auto it = fFractureFlowDim1MatId.begin(); it != fFractureFlowDim1MatId.end(); it++)
    {
        int fracture_id = *it;
        TPZFractureInsertion fracture(cmesh.Reference(),fracture_id,fMaterialBCIds);
        fracture.ClassifyNeighboursofPivots();
        fracture.OpenFractureOnH1(&cmesh); // (ok)
        fracture.SetDiscontinuosFrac(&cmesh); // (ok)
        //SetFracSubDomain(fracture_id);
        //fracture.AddWraps(&cmesh, fLagrangeFluxMatIdL);
        //fracture.SetInterfaces(&cmesh, fmatInterfaceLeft, fmatInterfaceRight);
        //fractureInsert=fracture;
    }
    cmesh.ComputeNodElCon();
    
}

void TPZMHMDarcyDFNMeshControl::SetFracSubDomain(int fracMatId){
    
    
    int nelem = fGMesh->NElements();
    
    for (int64_t el=0; el<nelem ; el++) {
        
        if (fGMesh->Element(el)->MaterialId()!=fracMatId) {
            continue;
        }

        if(fGMesh->Element(el)!=fGMesh->Element(el)->LowestFather()){
            continue;
        }
            
        
        TPZGeoEl *gel = fGMesh->Element(el);

        int ngeoToMHM = fGeoToMHMDomain.size();
        
        int elindex = gel->Index();
        
        ngeoToMHM++;
        fGeoToMHMDomain[elindex] = fGeoToMHMDomain[elindex];
        fMHMtoSubCMesh[gel->Index()] = -1;
        
    }

//
//
//
//#ifdef PZDEBUG
//    {
//        int64_t nel = fGMesh->NElements();
//        for (int64_t el = 0; el<nel; el++) {
//            TPZGeoEl *gel = fGMesh->Element(el);
//            if (gel && gel->MaterialId() == 1) {
//                if (gel->Dimension() != fGMesh->Dimension()) {
//                    DebugStop();
//                }
//            }
//        }
//    }
//#endif
//
//    int64_t nconnects = fCMesh->NConnects();
//
//    for (int64_t ic=nconnects; ic<nconnects; ic++) {
//        fCMesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
//    }
//    fGMesh->ResetReference();
//
//    // associate the connects with the proper subdomain
//    int64_t nel = fCMesh->NElements();
//
//    for (int64_t el=0; el<nel; el++)
//    {
//        TPZCompEl *cel = fCMesh->Element(el);
//#ifdef PZDEBUG
//        if (! cel) {
//            DebugStop();
//        }
//#endif
//        TPZGeoEl *gel = cel->Reference();
//
//        int dim = fGMesh->Dimension();
//
//        if (gel->Dimension()==dim) {
//            continue;
//        }
//        int gelindex = gel->Index();
//
//        int nsides = gel->NSides();
//
//        TPZGeoElSide gelside(gel,nsides-1);
//        TPZGeoElSide neighbour = gel->Neighbour(nsides-1);
//
//        int notintern = 0;
//
//        while(neighbour != gelside)
//        {
//            if(neighbour.Element()->MaterialId()==fInternalWrapMatId) {
//                notintern = 1;
//                break;
//            }
//            neighbour = neighbour.Neighbour();
//        }
//
//        if (notintern==1) {
//
//            neighbour = gel->Neighbour(nsides-1);
//
//            while(neighbour != gelside)
//            {
//                if (neighbour.Element()->Dimension() == dim) {
//                    break;
//                }
//                neighbour = neighbour.Neighbour();
//            }
//
//            int domain = fGeoToMHMDomain[neighbour.Element()->Index()];
//#ifdef PZDEBUG
//            if (domain == -1) {
//                DebugStop();
//            }
//#endif
//
//            SetSubdomain(cel, -1);
//        }else{
//            SetSubdomain(cel, -1);
//        }
//
//
//    }
//
//    std::ofstream out("PressureAndTractionFineMesh.txt");
//    fPressureFineMesh->Print(out);
    
    
}



// Insert internal multiplier by fracture elements
void TPZMHMDarcyDFNMeshControl::InsertInternalMultipliers(){

   // fGMesh->ResetReference();
    fCMesh->LoadReferences();
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if(!gel){
            continue;
        }
        if(gel->HasSubElement())
        {
            continue;
        }
        if ((gel->Dimension() != fGMesh->Dimension()-1)||(gel->MaterialId()!=fFractureMatId)) {
            continue;
        }
        int nsides = gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZGeoElSide> allneigh;
        gelside.AllNeighbours(allneigh);

        TPZGeoElSide neighbour = gelside.Neighbour();
        if(neighbour.Element()->HasSubElement())
        {
            continue;
        }
        //TPZStack<TPZGeoElSide> allneigh;
            
        gelside.AllNeighbours(allneigh);
            
        int count =0;
            while (neighbour != gelside) {
                
                if(neighbour.Element()->Dimension()==fGMesh->Dimension()){
                    
                    TPZCompElSide neigh(neighbour.Reference());
                    
                    if(neighbour.Side()==5) {
                        
                        CreateWrapElement(neigh,fLagrangeFluxMatIdL);
                        //TPZGeoElBC(gelside, fLagrangeFluxMatIdL);
                        count++;
                    }else if (neighbour.Side()==7) {
//                        TPZCompElSide left1(neighbour.Reference());
//                        TPZCompElSide right1(gelside.Reference());

                        CreateWrapElement(neigh,fLagrangeFluxMatIdR);
                        //TPZGeoElBC(gelside, fLagrangeFluxMatIdR);
                        count++;
                    }
                    //TPZGeoElBC(gelside, fLagrangeFluxMatIdL);

                }
                if (count == 2) {
                    break;
                }

                neighbour = neighbour.Neighbour();
            }
            
            if (neighbour == gelside) {
                continue;
            }
        
    }
    
}

void TPZMHMDarcyDFNMeshControl::CreateWrapElement(const TPZCompElSide &left, int multId){

//    TPZCompMesh *multiMesh = mfcel->Mesh();
//    TPZInterpolationSpace *hdivel = dynamic_cast<TPZInterpolationSpace *> (left.Element());
//    TPZCompEl *cel = dynamic_cast<TPZCompEl *>(left.Element());
//    MElementType celType = cel->Type();
//
//    TPZGeoEl *gel = left.Element()->Reference();
//
//    int dimMesh = left.Element()->Mesh()->Dimension();
//    if (!hdivel || !cel || gel->Dimension() != dimMesh) {
//        DebugStop();
//    }
//
//    //wrap element
//    TPZStack<TPZCompEl *, 7> wrapEl;
//    wrapEl.push_back(left.Element());
//
////    for (int side = 0; side < gel->NSides(); side++)
////    {
//    int side = left.Side();
//
////    if (gel->SideDimension(side) != gel->Dimension()-1) {
////            continue;
////        }
//
//    TPZGeoEl *gelbound = gel->CreateBCGeoEl(side, multId);
//    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(hdivel);
//    int loccon = intel->SideConnectLocId(0,side);
//
//    //Impose that the skeleton element has the same polynomial order  to the element of side.
//    TPZConnect &conside = intel->Connect(loccon);
//    int sideorder = conside.Order();
//    intel->Mesh()->SetDefaultOrder(sideorder);
//
//    int64_t index;
//    TPZInterpolationSpace *bound;
//    MElementType elType = gel->Type(side);
//    switch(elType)
//    {
//        case(EOned)://line
//        {
//            bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(* intel->Mesh(),gelbound,index);
//            int sideorient = intel->GetSideOrient(side);
//            TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
//            hdivbound->SetSideOrient(pzshape::TPZShapeLinear::NSides-1,sideorient);
//            break;
//        }
//        case(ETriangle)://triangle
//        {
//            bound = new TPZCompElHDivBound2<pzshape::TPZShapeTriang>(* intel->Mesh(),gelbound,index);
//            int sideorient = intel->GetSideOrient(side);
//            TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(bound);
//            hdivbound->SetSideOrient(pzshape::TPZShapeTriang::NSides-1,sideorient);
//            break;
//        }
//        case(EQuadrilateral)://quadrilateral
//        {
//            bound = new TPZCompElHDivBound2<pzshape::TPZShapeQuad>(* intel->Mesh(),gelbound,index);
//            int sideorient = intel->GetSideOrient(side);
//            TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(bound);
//            hdivbound->SetSideOrient(pzshape::TPZShapeQuad::NSides-1,sideorient);
//                break;
//
//        }
//
//        default:
//        {
//            bound=0;
//            std::cout << "ElementType not found!";
//            DebugStop();
//            break;
//        }
//    }
//
//    int64_t sideconnectindex = intel->ConnectIndex(loccon);
//
//    TPZConnect &co = bound->Connect(0);
//    if(co.HasDependency()){
//        if(bound->NConnects()!=1) DebugStop();
//        //int64_t cindex_bound = bound->ConnectIndex(0);
//        co.RemoveDepend();
//    }
//    bound->SetConnectIndex(0, sideconnectindex);
//    //bound->Print(std::cout);
//    gelbound->ResetReference();
//
//    TPZCompEl *newMFBound = fCMesh->CreateCompEl(gelbound, index);
    
//    TPZMultiphysicsElement *locMF = dynamic_cast<TPZMultiphysicsElement *>(newMFBound);
//
//    locMF->AddElement(bound, 0);
//    //locMF->Print(std::cout);
//
//    if(celType==EDiscontinuous){
//        TPZCompElDisc *discel = dynamic_cast<TPZCompElDisc *>(left.Element());
//        locMF->AddElement(TPZCompElSide(discel,side), 1);
//    }else{
//        TPZInterpolationSpace *contcel = dynamic_cast<TPZInterpolationSpace *>(left.Element());
//        locMF->AddElement(TPZCompElSide(contcel,side), 1);
//    }
//
//    wrapEl.push_back(locMF);
//    newMFBound->Reference()->ResetReference();
//    newMFBound->LoadElementReference();
    
   // ListGroupEl.push_back(wrapEl);
    
    
    
    TPZGeoElSide gleft(left.Reference());

    fCMesh->ApproxSpace().SetAllCreateFunctionsHDiv(fCMesh->Dimension());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());

    //intelleft->SetSideOrient(left.Side(), 1);
    //intelright->SetSideOrient(right.Side(), 1);

    TPZConnect &cleft = intelleft->SideConnect(0, left.Side());


    int64_t index = fCMesh->AllocateNewConnect(cleft);
    TPZConnect &newcon = fCMesh->ConnectVec()[index];
    cleft.DecrementElConnected();
    newcon.ResetElConnected();
    newcon.IncrementElConnected();
    newcon.SetSequenceNumber(fCMesh->NConnects() - 1);


    int sideorder = cleft.Order();
    fCMesh->SetDefaultOrder(sideorder);

    // create HDivBound on the sides of the elements
    TPZCompEl *wrap1;
    {
        intelleft->LoadElementReference();
        intelleft->SetPreferredOrder(sideorder);
        TPZGeoElBC gbc(gleft, multId);
        int64_t index;
        wrap1 = fCMesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fCMesh, index);
        if(cleft.Order() != sideorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap1);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        
        
        intel->SetSideOrient(wrapside, 1);
        intelleft->Reference()->ResetReference();
        wrap1->Reference()->ResetReference();
    }

    wrap1->LoadElementReference();

    intelleft->LoadElementReference();

    
}

void TPZMHMDarcyDFNMeshControl::InsertBCMultipliers(){
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        
        if(gel->HasSubElement()){
            continue;
        }
        
        for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
        {
            if (gel->MaterialId()==*it) {
                TPZGeoElBC(gelside, fBCLagrangeFluxMatIds[*it]);
            }
        }
    }
    
}


bool TPZMHMDarcyDFNMeshControl::HasFracNeighbour(TPZGeoElSide &gelside){
    

    TPZStack<TPZGeoElSide> allneigh;
    gelside.AllNeighbours(allneigh);
    for (auto ineigh : allneigh) {
        if (ineigh.Element()->HasSubElement()) {
            continue;
        }
        
        if (ineigh.Element()->MaterialId() == fFractureMatId) {
            return true;
        }
    }
    
    return false;

}

void TPZMHMDarcyDFNMeshControl::CreateFractureBC()
{

//    fGMesh->ResetReference();
//    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
//    {
//        int64_t nel = fCMesh->NElements();
//        for (int64_t el=0; el<nel; el++) {
//            TPZCompEl *cel = fCMesh->Element(el);
//            if (cel && cel->Reference()) {
//                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
//            }
//        }
//    }
//
//    //Create fracture boundaries
//    int nel = fGMesh->NElements();
//    for (int64_t el = 0; el<nel; el++) {
//        TPZGeoEl *gel = fGMesh->Element(el);
//        if(!gel){
//            continue;
//        }
//
//        int64_t index;
//
//        if (gel->MaterialId() == fmatFracPointL) {
//
////            TPZGeoElSide gelsideL(gel,0);
////            TPZInterpolatedElement *intelL = fluxelements[gelsideL.Element()->Index()];
////            if (intelL) {
////                intelL->LoadElementReference();
//                TPZCompEl *sidecelL = fCMesh->CreateCompEl(gel, index);
////                gel->ResetReference();
////            }
//        }
//
//        if (gel->MaterialId() == fmatFracPointR) {
////            TPZGeoElSide gelsideR(gel,1);
////            TPZInterpolatedElement *intelR = fluxelements[gelsideR.Element()->Index()];
////            if (intelR) {
////                intelR->LoadElementReference();
//                TPZCompEl *sidecelR = fCMesh->CreateCompEl(gel, index);
////                gel->ResetReference();
////            }
//        }
//
//    }
//    fCMesh->ExpandSolution();

    
///
    fGMesh->ResetReference();
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        int64_t nel = fCMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMesh->Element(el);
            if (cel && cel->Reference()) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    
    //Create fracture boundaries
    int nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if(!gel){
            continue;
        }
        if (gel->MaterialId() != fFractureMatId) {
            continue;
        }
        if (gel->HasSubElement()) {
            continue;
        }
        
        TPZGeoElSide gelsideL(gel,0);

        bool FractureQL = HasFracNeighbour(gelsideL);
        
        int64_t index;
        if(!FractureQL){
            int leftindex = gelsideL.Element()->Index();
            
            TPZGeoElBC gelbcL(gelsideL, fmatFracPointL);
            int dim = gelbcL.CreatedElement()->Dimension();

            
            TPZInterpolatedElement *intelL = fluxelements[gelsideL.Element()->Index()];
            if (intelL) {
                intelL->LoadElementReference();
                TPZCompEl *sidecelL = fCMesh->CreateCompEl(gelbcL.CreatedElement(), index);
                gelbcL.CreatedElement()->ResetReference();
            }
        }

        TPZGeoElSide gelsideR(gel,1);
        
        bool FractureQR = HasFracNeighbour(gelsideR);
        
        if(!FractureQR){
            TPZGeoElBC gelbcR(gelsideR, fmatFracPointR);
            
            TPZInterpolatedElement *intelR = fluxelements[gelsideR.Element()->Index()];
            if (intelR) {
                intelR->LoadElementReference();
                TPZCompEl *sidecelR = fCMesh->CreateCompEl(gelbcR.CreatedElement(), index);
                gelbcR.CreatedElement()->ResetReference();
            }
        }
        
    }
    fCMesh->ExpandSolution();
}
    
void TPZMHMDarcyDFNMeshControl::CreateH1Wrappers()
{
    fGMesh->ResetReference();
    int meshdim = fGMesh->Dimension();
    // build a vector which for each geoel index points to the flux element
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        int64_t nel = fCMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMesh->Element(el);
            if (cel && cel->Reference()) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    //fCMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->MaterialId() != fFractureMatId) {
            continue;
        }
        // now we have a geometric element of material id fInternalWrapMatId
        int nsides = gel->NSides();
        int side = nsides-1;
        if(gel->SideDimension(side) != meshdim-1) continue;
        // look for a neighbouring flux element, determine its side order
        TPZGeoElSide gelside(gel,side);
        if (!NeedsHDivWrapper(gelside)) {
            continue;
        }
        // create the HDiv wrapper element for all connected flux elements
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
            if (intel) {
                intel->LoadElementReference();
                
#ifdef PZDEBUG
                int midsideconnectlocid = intel->MidSideConnectLocId(neighbour.Side());
                int cindex = intel->ConnectIndex(midsideconnectlocid);
                // check whether a wrapper with this connect already exists
                // if affirmative - DebugStop
#endif
                TPZConnect &midsideconnect = intel->MidSideConnect(neighbour.Side());
                int nshapeC = midsideconnect.NShape();
                
                int order = midsideconnect.Order();
                fCMesh->SetDefaultOrder(order);
                TPZGeoElBC gelbc(gelside, fWrapH1MatIdL);
                int64_t index;
                TPZCompEl *sidecel = fCMesh->CreateCompEl(gelbc.CreatedElement(), index);
                
                // all the elements have normal pointing outwards
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Creating cap on element " << intel->Index() << " side " << neighbour.Side() << " element created " << sidecel->Index();
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                // intel->SetSideOrient(neighbour.Side(), 1);
                TPZInterpolatedElement *bound = dynamic_cast<TPZInterpolatedElement *>(sidecel);
                //bound->SetSideOrient(side, 1);
                gelbc.CreatedElement()->ResetReference();
                //neighbour.Element()->ResetReference();
                
            }
            ResetAllNeigReferences(gelside);
            neighbour = neighbour.Neighbour();
        }
    }
    
    fCMesh->ExpandSolution();
    
}


//void TPZMHMDarcyDFNMeshControl::CreateH1Wrappers()
//{
//    fGMesh->ResetReference();
//    int meshdim = fGMesh->Dimension();
//    // build a vector which for each geoel index points to the flux element
//    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
//    {
//        int64_t nel = fCMesh->NElements();
//        for (int64_t el=0; el<nel; el++) {
//            TPZCompEl *cel = fCMesh->Element(el);
//            if (cel && cel->Reference()) {
//                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
//            }
//        }
//    }
//    fCMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
//
//    int64_t nel = fGMesh->NElements();
//    for (int64_t el=0; el<nel; el++) {
//        TPZGeoEl *gel = fGMesh->Element(el);
//        if (!gel || gel->MaterialId() != fFractureMatId) {
//            continue;
//        }
//        // now we have a geometric element of material id fInternalWrapMatId
//        int nsides = gel->NSides();
//        int side = nsides-1;
//        if(gel->SideDimension(side) != meshdim-1) continue;
//        // look for a neighbouring flux element, determine its side order
//        TPZGeoElSide gelside(gel,side);
//        if (!NeedsHDivWrapper(gelside)) {
//            continue;
//        }
//        // create the HDiv wrapper element for all connected flux elements
//        TPZGeoElSide neighbour = gelside.Neighbour();
//        while(neighbour != gelside)
//        {
//            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
//            if (intel) {
//                intel->LoadElementReference();
//
//#ifdef PZDEBUG
//                int midsideconnectlocid = intel->MidSideConnectLocId(neighbour.Side());
//                int cindex = intel->ConnectIndex(midsideconnectlocid);
//                // check whether a wrapper with this connect already exists
//                // if affirmative - DebugStop
//#endif
//                TPZConnect &midsideconnect = intel->MidSideConnect(neighbour.Side());
//                int order = midsideconnect.Order();
//                fCMesh->SetDefaultOrder(order);
//                TPZGeoElBC gelbc(gelside, fWrapH1MatIdL);
//                int64_t index;
//                TPZCompEl *sidecel = fCMesh->CreateCompEl(gelbc.CreatedElement(), index);
//                // all the elements have normal pointing outwards
//#ifdef LOG4CXX
//                if (logger->isDebugEnabled())
//                {
//                    std::stringstream sout;
//                    sout << "Creating cap on element " << intel->Index() << " side " << neighbour.Side() << " element created " << sidecel->Index();
//                    LOGPZ_DEBUG(logger, sout.str())
//                }
//#endif
//                // intel->SetSideOrient(neighbour.Side(), 1);
//                TPZInterpolatedElement *bound = dynamic_cast<TPZInterpolatedElement *>(sidecel);
//                //bound->SetSideOrient(side, 1);
//                gelbc.CreatedElement()->ResetReference();
//                neighbour.Element()->ResetReference();
//            }
//            neighbour = neighbour.Neighbour();
//        }
//    }
//}


void TPZMHMDarcyDFNMeshControl::CreateHDivWrappers()
{
    fGMesh->ResetReference();
    int indw = 0;
    int meshdim = fGMesh->Dimension();
    // build a vector which for each geoel index points to the flux element
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        int64_t nel = fCMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMesh->Element(el);
            if (cel && cel->Reference()) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    fCMesh->ApproxSpace().SetAllCreateFunctionsHDiv(fCMesh->Dimension());
    
    int64_t nel = fGMesh->NElements();
    std::map<int64_t, std::pair<int64_t,int64_t> > mapdivided;
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->MaterialId() != fFractureMatId) {
            continue;
        }
        int gelindex = gel->Index(); // Fracture index
        // now we have a geometric element of material id fInternalWrapMatId
        int nsides = gel->NSides();
        int side = nsides-1;
        if(gel->SideDimension(side) != meshdim-1) continue;
        // look for a neighbouring flux element, determine its side order
        TPZGeoElSide gelside(gel,side);
        if (!NeedsHDivWrapper(gelside)) {
            continue;
        }
        // create the HDiv wrapper element for all connected flux elements
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            
            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
            if (intel){
                intel->LoadElementReference();
            }
            if (intel && neighbour.Element()->MaterialId()==fWrapH1MatIdL) { //Force to find wrapID
                
#ifdef PZDEBUG
                int midsideconnectlocid = intel->MidSideConnectLocId(neighbour.Side());
                int cindex = intel->ConnectIndex(midsideconnectlocid);
                // check whether a wrapper with this connect already exists
                // if affirmative - DebugStop
#endif
                TPZConnect &midsideconnect = intel->MidSideConnect(neighbour.Side());
                int nconnectedEl = midsideconnect.NElConnected();
                int order = midsideconnect.Order();

                int64_t index;
                
                int nshape = midsideconnect.NShape();
                
                //Duplica um connect
                index = fCMesh->AllocateNewConnect(midsideconnect.NShape(), midsideconnect.NState(), midsideconnect.Order());
                
                //intel->SetConnectIndex(2, index);
                
                intel->Reference()->ResetReference();
                ResetAllNeigReferences(gelside);
                // Create Geo and Comp Element
                fCMesh->SetDefaultOrder(order);
                TPZGeoElBC gelbc(gelside, fLagrangeFluxMatIdL);
                
                TPZCompEl *sidecel = fCMesh->CreateCompEl(gelbc.CreatedElement(), index);
                
                int64_t neigh_index = neighbour.Element()->Index(); //BC for H1
                int64_t gelbc_index = gelbc.CreatedElement()->Index(); // Hdiv 1D
                

                mapdivided[gelindex] = std::make_pair(neigh_index, gelbc_index);
                fFracInterfaces.insert(mapdivided);
                
                mapdivided.clear();
                
                //sidecel->SetConnectIndex(0, index);
                
                //midsideconnect.DecrementElConnected();
                //fCMesh->ConnectVec()[index].IncrementElConnected();
                
                // all the elements have normal pointing outwards
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Creating cap on element " << intel->Index() << " side " << neighbour.Side() << " element created " << sidecel->Index();
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                //intel->SetSideOrient(neighbour.Side(), 1);
                TPZInterpolatedElement *bound = dynamic_cast<TPZInterpolatedElement *>(sidecel);
                bound->SetSideOrient(side, 1);
                //gelbc.CreatedElement()->ResetReference();
                //neighbour.Element()->ResetReference();
                sidecel->Reference()->ResetReference();
                intel->Reference()->ResetReference();
            }
            
            ResetAllNeigReferences(gelside);
            neighbour = neighbour.Neighbour();
        }
        
        
        
    }
    fCMesh->ExpandSolution();
}

void TPZMHMDarcyDFNMeshControl::ResetAllNeigReferences(TPZGeoElSide gelside)
{

    TPZGeoElSide neighbour = gelside.Neighbour();
    while(neighbour != gelside)
    {
        
        neighbour.Element()->ResetReference();
        
        neighbour = neighbour.Neighbour();
    }
}


/// verifies if a HDiv wrapper needs to be created for a given element/side
// the method will look if there is a neighbouring geometric element with material id fPressureDim1MatId
bool TPZMHMDarcyDFNMeshControl::NeedsHDivWrapper(TPZGeoElSide gelside)
{
    TPZGeoElSide wrap = gelside.Neighbour();
    while(wrap != gelside && wrap.Element()->MaterialId() != fInternalWrapMatId)
    {
        wrap = wrap.Neighbour();
    }
    if (wrap.Element()->MaterialId() != fInternalWrapMatId) {
        DebugStop();
    }
    // we go all the way up
    while (wrap.Father2()) {
        wrap = wrap.Father2();
    }
    // now we look for a neighbour one level at a time
    TPZStack<TPZGeoElSide> elstack;
    elstack.Push(wrap);
    while (elstack.size()) {
        TPZGeoElSide search = elstack.Pop();
        TPZGeoElSide neighbour = search.Neighbour();
        while (neighbour != search) {
            if (fFractureFlowDim1MatId.find(neighbour.Element()->MaterialId()) != fFractureFlowDim1MatId.end()) {
                return true;
            }
            neighbour = neighbour.Neighbour();
        }
        if (search.HasSubElement()) {
            int nsub = search.Element()->NSubElements();
            for (int isub = 0; isub < nsub; isub++)
            {
                TPZGeoEl *subel = search.Element()->SubElement(isub);
                elstack.Push(TPZGeoElSide(subel,subel->NSides()-1));
            }
        }
    }
    return false;
}

/// will create the interface elements between the internal elements and the skeleton
void TPZMHMDarcyDFNMeshControl::CreateInterfaceElements()
{
    fCMesh->LoadReferences();
    int dim = fGMesh->Dimension();
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
    /// loop over the skeleton elements
    for (it=fInterfaces.begin(); it != fInterfaces.end(); it++) {
        // index of the skeleton element
        int64_t elindex = it->first;
        // left and right indexes in the coarse mesh
        int64_t leftelindex = it->second.first;
        int64_t rightelindex = it->second.second;
        // skip boundary elements
        int matid = 0, matidleft = 0, matidright = 0;
        
        // second condition indicates a boundary element
        if (leftelindex < rightelindex || rightelindex == elindex)
        {
            matidleft = fLagrangeMatIdLeft;
            matidright = fLagrangeMatIdRight;
        }
        else
        {
            matidleft = fLagrangeMatIdRight;
            matidright = fLagrangeMatIdLeft;
        }
        int numlado = 2;
        if (rightelindex == elindex) {
            numlado = 1;
        }
        
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZCompElSide celskeleton = gelside.Reference();
        // the skeleton element must exist
        if(!celskeleton.Element()) DebugStop();
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            int neighmatid = neighbour.Element()->MaterialId();
            if(neighmatid == fSkeletonWrapMatId || neighmatid == fBoundaryWrapMatId || neighmatid == fInternalWrapMatId)
            {
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        if (neighbour == gelside) {
            DebugStop();
        }
        TPZStack<TPZGeoElSide> gelstack;
        gelstack.Push(neighbour);
        while (gelstack.size())
        {
            TPZGeoElSide smallGeoElSide = gelstack.Pop();
            // the smaller elements returned by GetSubElements include element/sides of lower dimension
            if (smallGeoElSide.Dimension() != gel->Dimension()) {
                continue;
            }
            // look for the neighbours of smallGeoElSide for elements which are of dimension of the fGMesh
            // and which beint64_t to the subdomains
            TPZGeoElSide neighbour = smallGeoElSide.Neighbour();
            while (neighbour != smallGeoElSide)
            {
                if (neighbour.Element()->Dimension() != fGMesh->Dimension() ) {
                    neighbour = neighbour.Neighbour();
                    continue;
                }
                TPZCompElSide csmall = neighbour.Reference();
                if (csmall)
                {
                    int matid = -1;
                    if(WhichSubdomain(csmall.Element()) == leftelindex) {
                        matid = matidleft;
                    }
                    else if(WhichSubdomain(csmall.Element()) == rightelindex)
                    {
                        matid = matidright;
                    }
                    if (matid == -1) {
                      //  DebugStop();
                    }
                    // create an interface between the finer element and the MHM flux
                    int64_t index;
                    TPZGeoEl *gelnew = smallGeoElSide.Element()->CreateBCGeoEl(smallGeoElSide.Side(), matid);
                    new TPZInterfaceElement(fCMesh, gelnew , index, csmall, celskeleton);
#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "New interface left " << smallGeoElSide.Element()->Index() << " right " << gel->Index() << " matid " << matid;
                        sout << " interface index " << gelnew->Reference()->Index() << " beint64_ts to subdomain " << WhichSubdomain(gelnew->Reference());
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                neighbour = neighbour.Neighbour();
            }
            
            smallGeoElSide.GetSubElements2(gelstack);
        }
    }
}


void TPZMHMDarcyDFNMeshControl::InsertPeriferalMaterialObjects(){
    
    int matid = *fMaterialIds.begin();
    TPZMaterial *mat = fCMesh->FindMaterial(matid);
    if (!mat) {
        DebugStop();
    }
    
    TPZFNMatrix<1,STATE> xk(fNState,fNState,0.),xb(fNState,fNState,0.),xc(fNState,fNState,0.),xf(fNState,1,0.);
    TPZFNMatrix<4,STATE> val1(fNState,fNState,0.), val2Flux(fNState,1,0.);
    TPZMat1dLin *matPerif = NULL;
    
    if (fCMesh->FindMaterial(fSkeletonMatId)) {
        DebugStop();
    }
    matPerif = new TPZMat1dLin(fSkeletonMatId);
    matPerif->SetMaterial(xk, xc, xb, xf);
    fCMesh->InsertMaterialObject(matPerif);
    
    if (1) {
        if (fCMesh->FindMaterial(fPressureSkeletonMatId)) {
            DebugStop();
        }
        matPerif = new TPZMat1dLin(fPressureSkeletonMatId);
        matPerif->SetMaterial(xk, xc, xb, xf);
        fCMesh->InsertMaterialObject(matPerif);
        
        if (fCMesh->FindMaterial(fSecondSkeletonMatId)) {
            DebugStop();
        }
        matPerif = new TPZMat1dLin(fSecondSkeletonMatId);
        matPerif->SetMaterial(xk, xc, xb, xf);
        
        fCMesh->InsertMaterialObject(matPerif);
        
        
        int LagrangeMatIdLeft = 50;
        int LagrangeMatIdRight = 51;
        int nstate = fNState;
        int dim = fGMesh->Dimension();
        
        if (fCMesh->FindMaterial(fLagrangeMatIdLeft)) {
            DebugStop();
        }
        if (fCMesh->FindMaterial(fLagrangeMatIdRight)) {
            DebugStop();
        }
        if (fCMesh->FindMaterial(fLagrangeMatIdFrac)) {
            DebugStop();
        }
        TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(fLagrangeMatIdLeft,dim,nstate);
        TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(fLagrangeMatIdRight,dim,nstate);
        TPZLagrangeMultiplier *matFracInterfaces = new TPZLagrangeMultiplier(fLagrangeMatIdFrac,dim,nstate);
        if (fSwitchLagrangeSign) {
            matleft->SetMultiplier(-1.);
            matright->SetMultiplier(1.);
            matFracInterfaces->SetMultiplier(1.);
        }
        else
        {
            matleft->SetMultiplier(1.);
            matright->SetMultiplier(-1.);
            matFracInterfaces->SetMultiplier(1.);
        }
        fCMesh->InsertMaterialObject(matleft);
        fCMesh->InsertMaterialObject(matright);
        fCMesh->InsertMaterialObject(matFracInterfaces);
    }

}


/// will create the internal elements, one coarse element at a time
void TPZMHMDarcyDFNMeshControl::CreateInternalElements()
{
    TPZCompEl::SetgOrder(fpOrderInternal);
    
    fCMesh->ApproxSpace().SetCreateLagrange(false);
    fCMesh->SetAllCreateFunctionsContinuous();
    //fCMesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Criar elementos computacionais malha MHM
    
    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    TPZCompMesh *cmesh = fCMesh.operator->();
    fConnectToSubDomainIdentifier[cmesh].Expand(10000);
    
    int64_t nel = fGMesh->NElements();
    for (auto itMHM = fMHMtoSubCMesh.begin(); itMHM != fMHMtoSubCMesh.end(); itMHM++)
    {
        fGMesh->ResetReference();
        bool LagrangeCreated = false;
        for (int64_t el=0; el< nel; el++)
        {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel->MaterialId()==fFractureMatId) {
                continue;
            }
            
            if (!gel || gel->HasSubElement() || fGeoToMHMDomain[el] != itMHM->first) {
                continue;
            }
            int64_t index;
            // create the flux element
            fCMesh->CreateCompEl(gel, index);
            TPZCompEl *cel = fCMesh->Element(index);
            /// associate the connects with the subdomain
            SetSubdomain(cel, itMHM->first);
            // we need to create a lagrange multiplier element in order to delay decomposition of an equation
            if (!LagrangeCreated)
            {
                LagrangeCreated = true;
                TPZCompEl *cel = fCMesh->ElementVec()[index];
                int64_t cindex = cel->ConnectIndex(0);
                int nshape(1), nvar(1), order(1);
                int lagrangelevel = 0;
                if (this->fLagrangeAveragePressure) {
                    lagrangelevel = 1;
                }
                else
                {
                    lagrangelevel = 3;
                }
                fCMesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                if (fProblemType == EElasticity2D) {
                    cindex = cel->ConnectIndex(2);
                    fCMesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                }
                if (fProblemType == EElasticity3D) {
                    cindex = cel->ConnectIndex(6);
                    fCMesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                }
            }
        }
    }
    fGMesh->ResetReference();
}


void TPZMHMDarcyDFNMeshControl::InsertPeriferalLagrangeFluxObjects()
{
    
    TPZCompMesh * cmeshLagrangeFlux = fCMesh.operator->();
    
    // Material for interior traction:
    
    for (auto it:fMaterialIds){
        if (cmeshLagrangeFlux->MaterialVec().find(fLagrangeFluxMatIdL) == cmeshLagrangeFlux->MaterialVec().end())
        {
            TPZVecL2 *matLagrangeFluxL = new TPZVecL2(fLagrangeFluxMatIdL);
            matLagrangeFluxL->SetDimension(fGMesh->Dimension()-1);
            cmeshLagrangeFlux->InsertMaterialObject(matLagrangeFluxL);
        }
        if (cmeshLagrangeFlux->MaterialVec().find(fLagrangeFluxMatIdR) == cmeshLagrangeFlux->MaterialVec().end())
        {
            TPZVecL2 *matLagrangeFluxR = new TPZVecL2(fLagrangeFluxMatIdR);
            matLagrangeFluxR->SetDimension(fGMesh->Dimension()-1);
            cmeshLagrangeFlux->InsertMaterialObject(matLagrangeFluxR);
        }
    }
    

    for (auto it:fMaterialBCIds)
    {
        if (fBCLagrangeFluxMatIds.size()!=fMaterialBCIds.size()) {
            DebugStop();
        }
        
        int matid= fBCLagrangeFluxMatIds[it];
        if (cmeshLagrangeFlux->MaterialVec().find(matid) == cmeshLagrangeFlux->MaterialVec().end())
        {
            TPZVecL2 *matBCLagrangeFlux = new TPZVecL2(matid);
            matBCLagrangeFlux->SetDimension(fGMesh->Dimension()-1);
            cmeshLagrangeFlux->InsertMaterialObject(matBCLagrangeFlux);
        }
        
    }

}



/// will create the internal elements, one coarse element at a time
void TPZMHMDarcyDFNMeshControl::CreateLagrangeFluxElements()
{
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    
    int64_t nskeletonconnects = fCMesh->NConnects();
    int porder = fpOrderInternal;
    TPZCompMesh * cmeshLFlux = fCMesh.operator->();
    gmesh->ResetReference();
    cmeshLFlux->SetName("LagrangeFluxMesh");
    cmeshLFlux->SetDimModel(gmesh->Dimension()-1);
    cmeshLFlux->SetAllCreateFunctionsHDiv();
    cmeshLFlux->SetDefaultOrder(porder-1);
    int meshdim = cmeshLFlux->Dimension();
    
    std::set<int> matids;
    TPZMaterial *matL = cmeshLFlux->FindMaterial(fLagrangeFluxMatIdL);
    TPZMaterial *matR = cmeshLFlux->FindMaterial(fWrapH1MatIdL);
    if (matL && matR && matR->Dimension() == meshdim) {
        matids.insert(fLagrangeFluxMatIdL);
//        matids.insert(fWrapH1MatIdL);
    }
    
    for (auto it:fMaterialBCIds) {
        int dsmatid = fBCLagrangeFluxMatIds[it];
        TPZMaterial *mat = cmeshLFlux->FindMaterial(fBCLagrangeFluxMatIds[it]);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(fBCLagrangeFluxMatIds[it]);
        }
    }
    
    cmeshLFlux->AutoBuild(matids);
    fCMesh->ExpandSolution();
    
    
    
#ifdef PZDEBUG
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    int64_t nc = cmeshLFlux->NConnects();
    if(nskeletonconnects != 0){
        //      DebugStop();
    }
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshLFlux->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
    
    // associate the connects with the proper subdomain
    int64_t nel = cmeshLFlux->NElements();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshLFlux->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        
        int dim = gmesh->Dimension();
        
        if (gel->Dimension()==dim) {
            continue;
        }
        int gelindex = gel->Index();
        
        int nsides = gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gel->Neighbour(nsides-1);
        
        int notintern = 0;
        
        while(neighbour != gelside)
        {
            if(neighbour.Element()->MaterialId()==fInternalWrapMatId) {
                notintern = 1;
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        
        if (notintern==1) {
            
            neighbour = gel->Neighbour(nsides-1);
            
            while(neighbour != gelside)
            {
                if (neighbour.Element()->Dimension() == dim) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            
            int domain = fGeoToMHMDomain[neighbour.Element()->Index()];
#ifdef PZDEBUG
            if (domain == -1) {
                DebugStop();
            }
#endif
            
            SetSubdomain(cel, domain);
        }else{
            SetSubdomain(cel, -1);
        }
    }
    
    std::ofstream out("PressureAndLagrangeFluxMesh.txt");
    fCMesh->Print(out);
    
}

/// will create the elements on the skeleton
void TPZMHMDarcyDFNMeshControl::CreateSkeleton()
{
    fGMesh->ResetReference();
    // comment this line or not to switch the type of skeleton elements
    int meshdim = fCMesh->Dimension();
    fCMesh->SetDimModel(meshdim);
    //fCMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    fCMesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    int order = fpOrderSkeleton;
    if (order < 0) {
        order = 0;
    }
    fCMesh->SetDefaultOrder(order);
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        int64_t elindex = it->first;
        // skip the boundary elements
        //        if (elindex == it->second.second) {
        //            it++;
        //            continue;
        //        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        int64_t index;
        // create a discontinuous element to model the flux
        fCMesh->CreateCompEl(gel, index);
        TPZCompEl *cel = fCMesh->ElementVec()[index];
        int Side = gel->NSides()-1;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int lagrangelevel = 0;
            if (this->fLagrangeAveragePressure) {
                lagrangelevel = 2;
            }
            else
            {
                lagrangelevel = 2;
            }
            
            if(gel->MaterialId()==fFractureMatId){
                int cindex = cel->ConnectIndex(ic);
                intel->SetConnectIndex(0, cindex);
            }
            cel->Connect(ic).SetLagrangeMultiplier(lagrangelevel);
        }
        SetSubdomain(cel, -1);
        
        if (elindex == it->second.second) {
            // set the side orientation of the boundary elements
            intel->SetSideOrient(Side, 1);
            SetSubdomain(cel, it->second.first);
        }
        else
        {
            if (it->second.first < it->second.second) {
                // set the flux orientation depending on the relative value of the element ids
                intel->SetSideOrient(Side, 1);
            }
            else
            {
                intel->SetSideOrient(Side, -1);
            }
            SetSubdomain(cel, -1);
        }
        gel->ResetReference();
        it++;
    }
    fCMesh->SetDimModel(meshdim);
}

/// will create 1D elements on the interfaces between the coarse element indices
void TPZMHMDarcyDFNMeshControl::CreateSkeletonElements()
{
    if (fInterfaces.size()) {
        DebugStop();
    }
    
    if(fSkeletonMatId < 0) DebugStop();
    
    TPZCompMesh *cmesh = CriaMalhaTemporaria();
    
    int64_t nel = fGMesh->NElements();
    int dimension = fGMesh->Dimension();
    int ninterf;
    
    for(int64_t iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = fGMesh->ElementVec()[iel];
        if(!gel) continue;
        
        ninterf = gel->NumInterfaces();
        if(ninterf > 1) DebugStop();
        if (ninterf==1 && gel->MaterialId()!= fFractureMatId) // Não deve inserir interfaces na fratura
        {
            TPZCompEl *cel = gel->Reference();
            TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
            if (!intface) {
                DebugStop();
            }
            int interfacematid = fSkeletonMatId;
            TPZCompEl *left = intface->LeftElement();
            TPZCompEl *right = intface->RightElement();
            int64_t leftind = left->Reference()->Index();
            int64_t rightind = right->Reference()->Index();
            if (left->Reference()->Dimension() == dimension-1) {
                // the boundary element is the skeleton element
                fInterfaces[leftind]=std::make_pair(rightind, leftind);
                continue;
            }
            if (right->Reference()->Dimension() == dimension-1) {
                // the boundary element is the skeleton element
                fInterfaces[rightind]=std::make_pair(leftind, rightind);
                continue;
            }
            fInterfaces[iel] = std::make_pair(leftind, rightind);
            gel->SetMaterialId(interfacematid);
            // in order to prevent the element from being deleted
            gel->DecrementNumInterfaces();
        }
    }
    
    fGMesh->ResetReference();
    delete cmesh;
    
    //    BuildWrapMesh(fGMesh->Dimension());
    //    BuildWrapMesh(fGMesh->Dimension()-1);
    
    fGeoToMHMDomain.Resize(fGMesh->NElements(), -1);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it = fInterfaces.begin();
        while (it != fInterfaces.end()) {
            int leftdim = fGMesh->ElementVec()[it->second.first]->Dimension();
            int rightdim = fGMesh->ElementVec()[it->second.second]->Dimension();
            sout << "Interface index " << it->first << " Left Element " << it->second.first << "/" << leftdim
            << " Right Element " << it->second.second << "/" << rightdim << std::endl;
            it++;
        }
        sout << "Geometric mesh after creating the skeleton\n";
        fGMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// will create 1D elements on the interfaces between the coarse element indices
void TPZMHMDarcyDFNMeshControl::CreateFractureInterfaces()
{
    
    fCMesh->LoadReferences();
    int dim = fGMesh->Dimension();
    
    
    
//    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
    /// loop over the skeleton elements
    
    
    
//    for (it=fInterfaces.begin(); it != fInterfaces.end(); it++) {
//        // index of the skeleton element
//        int64_t elindex = it->first;
//        // left and right indexes in the coarse mesh
//        int64_t leftelindex = it->second.first;
//        int64_t rightelindex = it->second.second;
//        // skip boundary elements
//        int matid = 0, matidleft = 0, matidright = 0;
//
//        // second condition indicates a boundary element
//        if (leftelindex < rightelindex || rightelindex == elindex)
//        {
//            matidleft = fLagrangeMatIdLeft;
//            matidright = fLagrangeMatIdRight;
//        }
//        else
//        {
//            matidleft = fLagrangeMatIdRight;
//            matidright = fLagrangeMatIdLeft;
//        }
//        int numlado = 2;
//        if (rightelindex == elindex) {
//            numlado = 1;
//        }

    int nelem = fGMesh->NElements();
//
//    for (int iel = 0; iel < nelem; iel++) {
//        TPZGeoEl *gel = fGMesh->Element(iel);
//        if (gel->MaterialId()!=fGMesh->Dimension()) {
//            continue;
//        }
//        if (gel->HasSubElement()) {
//            continue;
//        }
//    }
//
//
//
//    int64_t leftelindex = it->second.first;
//    int64_t rightelindex = it->second.second;
    
    for (auto mapfrac : fFracInterfaces) {

        std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it = mapfrac.begin();
        while (it != mapfrac.end()) {
            int64_t elindex = it->first;
            //    for (int iel = 0; iel< nelem; iel++) {
            
            TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
            if (gel->MaterialId()!=fFractureMatId) {
                DebugStop();
            }
            
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZCompElSide celfrac = gelside.Reference();
            TPZCompElSide celskeleton;
            // the skeleton element must exist
            if(!celfrac.Element()) DebugStop();
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside)
            {
                int neighmatid = neighbour.Element()->MaterialId();
                if(neighbour.Element()->Index()==it->second.first)
                {
                    celskeleton = neighbour.Reference();
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                DebugStop();
            }
            TPZStack<TPZGeoElSide> gelstack;
            gelstack.Push(neighbour);
            while (gelstack.size())
            {
                TPZGeoElSide smallGeoElSide = gelstack.Pop();
                // the smaller elements returned by GetSubElements include element/sides of lower dimension
                if (smallGeoElSide.Dimension() != gel->Dimension()) {
                    continue;
                }
                // look for the neighbours of smallGeoElSide for elements which are of dimension of the fGMesh
                // and which beint64_t to the subdomains
                TPZGeoElSide neighbour = smallGeoElSide.Neighbour();
                while (neighbour != smallGeoElSide)
                {
                    //                if (neighbour.Element()->Dimension() != fGMesh->Dimension() ) {
                    //                    neighbour = neighbour.Neighbour();
                    //                    continue;
                    //                }
                    if (neighbour.Element()->Index()!=it->second.second) {
                        neighbour = neighbour.Neighbour();
                        continue;
                    }
                    
                    TPZCompElSide csmall = neighbour.Reference();
                    if (csmall)
                    {
                        int matid = fLagrangeMatIdFrac;
                        //                    if(WhichSubdomain(csmall.Element()) == leftelindex) {
                        //                        matid = matidleft;
                        //                    }
                        //                    else if(WhichSubdomain(csmall.Element()) == rightelindex)
                        //                    {
                        //                        matid = matidright;
                        //                    }
                        //                    if (matid == -1) {
                        //                        //  DebugStop();
                        //                    }
                        // create an interface between the finer element and the MHM flux
                        int64_t index;
                        TPZGeoEl *gelnew = smallGeoElSide.Element()->CreateBCGeoEl(smallGeoElSide.Side(), matid);
                        new TPZInterfaceElement(fCMesh, gelnew , index, csmall, celskeleton);

                        TPZGeoEl *gelnew2 = smallGeoElSide.Element()->CreateBCGeoEl(smallGeoElSide.Side(), matid);
                        new TPZInterfaceElement(fCMesh, gelnew2 , index, csmall, celfrac);
                        
                        
                        //                    std::cout << "New interface left " << smallGeoElSide.Element()->Index() << " right " << gel->Index() << " matid " << matid << std::endl;
                        //                    std::cout  << " interface index " << gelnew->Reference()->Index() << " beint64_ts to subdomain " << WhichSubdomain(gelnew->Reference())<< std::endl;
                        
#ifdef LOG4CXX
                        if (logger->isDebugEnabled()) {
                            std::stringstream sout;
                            sout << "New interface left " << smallGeoElSide.Element()->Index() << " right " << gel->Index() << " matid " << matid;
                            sout << " interface index " << gelnew->Reference()->Index() << " beint64_ts to subdomain " << WhichSubdomain(gelnew->Reference());
                            LOGPZ_DEBUG(logger, sout.str())
                        }
#endif
                    }
                    neighbour = neighbour.Neighbour();
                }
                
                smallGeoElSide.GetSubElements2(gelstack);
            }
            it++;
            
        }

    }
    
}

/// create the lagrange multiplier mesh, one element for each subdomain
void TPZMHMDarcyDFNMeshControl::CreateLagrangeMultiplierMesh(int dim)
{
    
    //int dim = fGMesh->Dimension();
    fCMeshLagrange->SetDimModel(dim);
    fCMeshLagrange->SetAllCreateFunctionsDiscontinuous();
    if (fProblemType == EScalar) {
        fCMeshLagrange->SetDefaultOrder(0);
    }
    else if(fProblemType == EElasticity2D)
    {
        fCMeshLagrange->SetDefaultOrder(1);
    }
    fGMesh->ResetReference();
    int64_t connectcounter = fCMesh->NConnects();
    /// criar materiais
    std::set<int> matids;
    TPZGeoMesh &gmesh = fGMesh;
    int64_t nel = gmesh.NElements();
    // this code needs to be modified to create lagrange computational elements which share a connect
    // between each other
    //DebugStop();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(el) == fMHMtoSubCMesh.end()) {
            continue;
        }
        int materialid = gel->MaterialId();
        matids.insert(materialid);
    }
    
//    if(dim==1){
//        matids.insert(fFractureMatId);   // Insert matID seting fracture for Darcy
//    }

    
    TPZManVector<STATE,1> sol(1,0.);
    int nstate = 1;
    std::set<int>::iterator it = matids.begin();
    TPZMaterial *meshmat = 0;
    while (it != matids.end()) {
        TPZNullMaterial *material = new TPZNullMaterial(*it);
        fCMeshLagrange->InsertMaterialObject(material);
        if (!meshmat) {
            meshmat = material;
        }
        it++;
        
    }
    if (!meshmat) {
        DebugStop();
    }
    TPZCompElDisc::SetTotalOrderShape(fCMeshLagrange.operator->());
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(el) == fMHMtoSubCMesh.end()) {
            continue;
        }
        int64_t index;
        TPZCompElDisc *disc = new TPZCompElDisc(fCMeshLagrange,gel,index);
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        int64_t cindex = disc->ConnectIndex(0);
#ifdef PZDEBUG
        static int count = 0;
        if (count == 0)
        {
            TPZConnect &c = disc->Connect(0);
            std::cout << "Number of shape functions of discontinuous element " << c.NShape() << std::endl;
            count++;
        }
#endif
        SetSubdomain(disc, el);
        
        //        fCMeshConstantStates->CreateCompEl(gel, index);
    }
    fCMeshLagrange->ExpandSolution();
    fGMesh->ResetReference();
    
    fCMeshConstantPressure = new TPZCompMesh(fCMeshLagrange);
    {
        int64_t nel = fCMeshConstantPressure->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMeshConstantPressure->Element(el);
            TPZCompEl *cel2 = fCMeshLagrange->Element(el);
            int subdomain = WhichSubdomain(cel2);
            SetSubdomain(cel, subdomain);
        }
    }

    if (dim==1) {
        std::ofstream filec1("MalhaC_Lagrange.txt");
        fCMeshLagrange->Print(filec1);
        
        std::ofstream filec2("MalhaC_ConstPress.txt");
        fCMeshConstantPressure->Print(filec2);
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        fCMeshLagrange->Print(sout);
        fCMeshConstantPressure->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fGMesh->ResetReference();
}

/// create the lagrange multiplier mesh, one element for each subdomain
void TPZMHMDarcyDFNMeshControl::CreateLagrangeMultiplierMeshLocal(int dim)
{
    
    fCMeshLagrangeLocal->SetDimModel(dim);
    fCMeshLagrangeLocal->SetAllCreateFunctionsDiscontinuous();
    if (fProblemType == EScalar) {
        fCMeshLagrangeLocal->SetDefaultOrder(0);
    }
    else if(fProblemType == EElasticity2D)
    {
        fCMeshLagrangeLocal->SetDefaultOrder(1);
    }
    fGMesh->ResetReference();
    int64_t connectcounter = fCMesh->NConnects();
    /// criar materiais
    std::set<int> matids;
    TPZGeoMesh &gmesh = fGMesh;
    int64_t nel = gmesh.NElements();
    // this code needs to be modified to create lagrange computational elements which share a connect
    // between each other
    //DebugStop();
//    for (int64_t el=0; el<nel; el++) {
//        TPZGeoEl *gel = gmesh.ElementVec()[el];
//        if (!gel) {
//            continue;
//        }
//        if (gel->HasSubElement()) {
//            continue;
//        }
//        if (gel->Dimension() != dim) {
//            continue;
//        }
//        int materialid = gel->MaterialId();
//        matids.insert(materialid);
//    }
    
    int materialid =0;
    
    if (dim == fGMesh->Dimension()) {
        materialid = 1;
    }else if(dim == fGMesh->Dimension()-1){
        materialid = fFractureMatId;
    }
    matids.insert(materialid);
    
    
    TPZManVector<STATE,1> sol(1,0.);
    int nstate = 1;
    std::set<int>::iterator it = matids.begin();
    TPZMaterial *meshmat = 0;
    while (it != matids.end()) {
        TPZNullMaterial *material = new TPZNullMaterial(*it);
        fCMeshLagrangeLocal->InsertMaterialObject(material);
        if (!meshmat) {
            meshmat = material;
        }
        it++;
        
    }
    if (!meshmat) {
        DebugStop();
    }
    TPZCompElDisc::SetTotalOrderShape(fCMeshLagrangeLocal.operator->());
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) continue;
        if (gel->MaterialId() != materialid) continue;
        if (gel->Dimension() != dim) continue;
        if (gel->HasSubElement()) continue;

        int64_t index;
        TPZCompElDisc *disc = new TPZCompElDisc(fCMeshLagrangeLocal,gel,index);
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        int64_t cindex = disc->ConnectIndex(0);
#ifdef PZDEBUG
        static int count = 0;
        if (count == 0)
        {
            TPZConnect &c = disc->Connect(0);
            std::cout << "Number of shape functions of discontinuous element " << c.NShape() << std::endl;
            count++;
        }
#endif
//        SetSubdomain(disc, el);
        
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(disc, fGeoToMHMDomain[gel->Index()]);

    }
    fCMeshLagrangeLocal->ExpandSolution();
    fGMesh->ResetReference();
    
    fCMeshConstantPressureLocal = new TPZCompMesh(fCMeshLagrangeLocal);
    {
        int64_t nel = fCMeshConstantPressureLocal->NElements();
        for (int64_t el=0; el<nel; el++) {
//
//            TPZGeoEl *gel = gmesh.ElementVec()[el];
//            if (!gel) continue;
//            if (gel->MaterialId() != materialid) continue;
//            if (gel->Dimension() != dim) continue;
//            if (gel->HasSubElement()) continue;

            TPZCompEl *cel = fCMeshConstantPressureLocal->Element(el);
            TPZCompEl *cel2 = fCMeshLagrangeLocal->Element(el);
            int subdomain = WhichSubdomain(cel2);
            SetSubdomain(cel, subdomain);
        }
    }
    
//    if (dim==1) {
        std::ofstream filec1("MalhaC_Lagrange_Local.txt");
        fCMeshLagrangeLocal->Print(filec1);
        
        std::ofstream filec2("MalhaC_ConstPress_Local.txt");
        fCMeshConstantPressureLocal->Print(filec2);
//    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        fCMeshLagrangeLocal->Print(sout);
        fCMeshConstantPressureLocal->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fGMesh->ResetReference();
}

/// transform the computational mesh into a multiphysics mesh
void TPZMHMDarcyDFNMeshControl::TransferToMultiphysics()
{
    fGMesh->ResetReference();
    this->fCMesh = new TPZCompMesh(fGMesh);
    this->fCMesh->SetDimModel(fGMesh->Dimension());
    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    
    // copy the material objects
    std::map<int,TPZMaterial *>::iterator it = fPressureFineMesh->MaterialVec().begin();
    while (it != fPressureFineMesh->MaterialVec().end()) {
        it->second->Clone(fCMesh->MaterialVec());
        it++;
    }
    
    int64_t nel = fPressureFineMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
        if (intface) {
            continue;
        }
        TPZCompElLagrange *lagr = dynamic_cast<TPZCompElLagrange *>(cel);
        if (lagr) {
            int64_t index;
            new TPZCompElLagrange(fCMesh,*cel,index);
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        int64_t index;
        fCMesh->CreateCompEl(gel, index);
        TPZCompEl *celnew = fCMesh->ElementVec()[index];
        TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celnew);
        if (!mult) {
            DebugStop();
        }
        mult->AddElement(cel, 0);
    }
    fCMesh->LoadReferences();
    

    
    nel = fCMeshLagrangeLocal->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMeshLagrangeLocal->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
//        gelside.ConnectedCompElementList(celstack, 0, 0);
        if (gel->Reference()) {
            celstack.Push(gelside.Reference());
        }
        int nstack = celstack.size();
        if(nstack != 1) DebugStop();
        for (int ist=0; ist<nstack; ist++) {
            TPZCompEl *celmult = celstack[ist].Element();
            TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celmult);
            mult->AddElement(cel, 1);
        }
    }
    fGMesh->ResetReference();
    fCMesh->LoadReferences();
   
    
//    nel = fCMeshLagrange->NElements();
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = fCMeshLagrange->ElementVec()[el];
//        TPZGeoEl *gel = cel->Reference();
//        int nsides = gel->NSides();
//        TPZGeoElSide gelside(gel,nsides-1);
//        TPZStack<TPZCompElSide> celstack;
//        //        gelside.ConnectedCompElementList(celstack, 0, 0);
//        if (gel->Reference()) {
//            celstack.Push(gelside.Reference());
//        }
//        int nstack = celstack.size();
//        if(nstack != 1) DebugStop();
//        for (int ist=0; ist<nstack; ist++) {
//            TPZCompEl *celmult = celstack[ist].Element();
//            TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celmult);
//            mult->AddElement(cel, 3);
//        }
//    }
//    fGMesh->ResetReference();
//    fCMesh->LoadReferences();
    
    nel = fCMeshConstantPressureLocal->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMeshConstantPressureLocal->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
      //  gelside.ConnectedCompElementList(celstack, 0, 0);
        if (gel->Reference()) {
            celstack.Push(gelside.Reference());
        }
        int nstack = celstack.size();
        if(nstack != 1) DebugStop();
        for (int ist=0; ist<nstack; ist++) {
            TPZCompEl *celmult = celstack[ist].Element();
            TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celmult);
            mult->AddElement(cel, 2);
        }
    }
    
    nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZMultiphysicsElement *multel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!multel) {
            continue;
        }
        int nmeshes = multel->NMeshes();
        
        for (int im=nmeshes; im<3; im++) {
            multel->AddElement(0, im);
        }
    }
    
    //void TPZBuildMultiphysicsMesh::AddConnects(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh)
    TPZManVector<TPZCompMesh *,3> cmeshvec(3,0);
    cmeshvec[0] = fPressureFineMesh.operator->();
    cmeshvec[1] = fCMeshLagrangeLocal.operator->();
    cmeshvec[2] = fCMeshConstantPressureLocal.operator->();
    TPZCompMesh *cmesh = fCMesh.operator->();
    TPZBuildMultiphysicsMesh::AddConnects(cmeshvec,cmesh);
    
    nel = fPressureFineMesh->NElements();
    
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!intface) {
            continue;
        }
        //        TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right);
        TPZCompElSide pressleft = intface->LeftElementSide();
        TPZCompElSide pressright = intface->RightElementSide();
        TPZGeoElSide gleft = pressleft.Reference();
        TPZGeoElSide gright = pressright.Reference();
        TPZCompElSide multleft = gleft.Reference();
        TPZCompElSide multright = gright.Reference();
        int64_t index;
        new TPZMultiphysicsInterfaceElement(fCMesh,intface->Reference(),index,multleft,multright);
        
    }
    
    
    nel = fCMeshConstantPressureLocal->NElements();
    int64_t npressconnect = fPressureFineMesh->NConnects();
    int64_t nlagrangeconnect = fCMeshLagrangeLocal->NConnects();
    // nel numero de dominios MHM, tem um connect associado a cada um e os mesmos estao no final
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = this->fCMeshConstantPressureLocal->Element(el);
        int64_t pressureconnect = cel->ConnectIndex(0);
        int64_t cindex = npressconnect+nlagrangeconnect+pressureconnect;
        fCMesh->ConnectVec()[cindex].SetLagrangeMultiplier(3);
    }
    fCMesh->ExpandSolution();
    //fCMesh->SaddlePermute();

    {
        std::ofstream filecp("Malha_FinePressure.txt");
        fPressureFineMesh->Print(filecp);
    
        std::ofstream filecL("Malha_LagrangeMult.txt");
        fCMeshLagrangeLocal->Print(filecL);
        
        std::ofstream filecm("MalhaC2_MHM.txt");
        fCMesh->Print(filecm);
    }
    
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        fCMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

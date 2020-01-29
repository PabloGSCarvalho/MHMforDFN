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

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mhmdarcyDFNmeshcontrol"));
#endif

TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices), fFractureMatId(),  fFractureFlowDim1MatId(), fBCLagrangeFluxMatIds()
{
    
}

TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(int dimension) : TPZMHMeshControl(dimension), fFractureMatId(), fFractureFlowDim1MatId(), fBCLagrangeFluxMatIds() {
    
}

TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMeshControl(gmesh), fFractureMatId(), fFractureFlowDim1MatId(), fBCLagrangeFluxMatIds() {
    
}


TPZMHMDarcyDFNMeshControl::TPZMHMDarcyDFNMeshControl(const TPZMHMDarcyDFNMeshControl &copy) : TPZMHMeshControl(copy){
    
    this->operator=(copy);
}

TPZMHMDarcyDFNMeshControl &TPZMHMDarcyDFNMeshControl::operator=(const TPZMHMDarcyDFNMeshControl &cp){
    
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
    
    CreateH1Wrappers();
    CreateHDivWrappers();

    std::ofstream fileCData("ConnectsData.txt");
    this->Print(fileCData);
    
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
    
    std::ofstream filecm("MalhaC2_MHM.txt");
    CMesh()->Print(filecm);

    //CreateSkeleton();
    

    
    CreateInterfaceElements();
    //    AddBoundaryInterfaceElements();
    fCMesh->ExpandSolution();
    fCMesh->CleanUpUnconnectedNodes();
    if (fLagrangeAveragePressure) {
        this->CreateLagrangeMultiplierMesh();
        this->TransferToMultiphysics();
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

void TPZMHMDarcyDFNMeshControl::Print(std::ostream & out) const {
    
    //ComputeNodElCon();
    out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
    
    out << "number of connects            = " << fCMesh->NConnects() << std::endl;
    out << "number of elements            = " << fCMesh->NElements() << std::endl;

    
    out << "\n\t Connect Information:\n\n";
    int64_t i, nelem = 0.;

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
//#ifdef PZDEBUG
//            if (c.NShape() != cel->NConnectShapeF(nod, c.Order())) {
//                DebugStop();
//            }
//#endif
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
//        out << "Side orders = ";
//        int NSides = cel->Reference()->NSides();
//        for (int is = 0; is < NSides; is++) {
//            int nsconnect = NSideConnects(is);
//            out << " side " << is << " orders ";
//            for (int ic = 0; ic < nsconnect; ic++) {
//                int sloc = SideConnectLocId(ic, is);
//                int order = cel->Connect(sloc).Order();
//                out << order << " ";
//            }
//            out << std::endl;
//        }
        
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
    
    int64_t nconnects = fCMesh->NConnects();

    for (int64_t ic=nconnects; ic<nconnects; ic++) {
        fCMesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    fGMesh->ResetReference();
    
    // associate the connects with the proper subdomain
    int64_t nel = fCMesh->NElements();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = fCMesh->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        
        int dim = fGMesh->Dimension();
        
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
            
            SetSubdomain(cel, -1);
        }else{
            SetSubdomain(cel, -1);
        }
        
        
    }
    
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
                ResetAllNeigReferences(gelside);
            }
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
//                sidecel->Reference()->ResetReference();
                intel->Reference()->ResetReference();
                ResetAllNeigReferences(gelside);
                
            }
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
                if (neighbour.Element()->Dimension() != fGMesh->Dimension()) {
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
                        DebugStop();
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
    
    if (0) {
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
        TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(fLagrangeMatIdLeft,dim,nstate);
        TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(fLagrangeMatIdRight,dim,nstate);
        if (fSwitchLagrangeSign) {
            matleft->SetMultiplier(-1.);
            matright->SetMultiplier(1.);
        }
        else
        {
            matleft->SetMultiplier(1.);
            matright->SetMultiplier(-1.);
        }
        fCMesh->InsertMaterialObject(matleft);
        fCMesh->InsertMaterialObject(matright);
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

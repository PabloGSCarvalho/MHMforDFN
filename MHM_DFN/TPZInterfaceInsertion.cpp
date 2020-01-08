
//  TPZInterfaceInsertion.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 02/08/18.
//

#include "TPZInterfaceInsertion.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "pzcompelwithmem.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "TPZMultiphysicsInterfaceEl.h"

/// Default constructor
TPZInterfaceInsertion::TPZInterfaceInsertion(){
    m_id_flux_wrap = 0;
    m_interface_id = 0;
    m_interfaceVector_ids.Resize(0);
    m_multiplier_id = 0;
    m_multiplierBC_id = 0;
    m_cmesh = NULL;
    m_geometry = NULL;
    m_boundaries_ids.clear();
    m_Eltype = EQuadrilateral;
}

/// Default desconstructor
TPZInterfaceInsertion::~TPZInterfaceInsertion(){
    
}

/// Copy constructor
TPZInterfaceInsertion::TPZInterfaceInsertion(TPZInterfaceInsertion & other){
    m_interface_id       = other.m_interface_id;
    m_interfaceVector_ids     = other.m_interfaceVector_ids;
    m_id_flux_wrap       = other.m_id_flux_wrap;
    m_multiplier_id      = other.m_multiplier_id;
    m_multiplierBC_id      = other.m_multiplierBC_id;
    m_geometry           = other.m_geometry;
    m_boundaries_ids     = other.m_boundaries_ids;
    m_Eltype            = other.m_Eltype;
}

/// Constructor based on a computational mesh and fracture material id
TPZInterfaceInsertion::TPZInterfaceInsertion(TPZCompMesh *cmesh, int mat_multiplier, std::set<int> & boundaries_ids, MElementType eltype){
    m_cmesh = cmesh;
    m_geometry = cmesh->Reference();
    m_multiplier_id = mat_multiplier;
    m_boundaries_ids = boundaries_ids;
    m_interface_id = 0;
    m_interfaceVector_ids.Resize(0);
    m_Eltype            = eltype;
    
}

/// Set Interface Identifier (interfaces between multiplier and volumetric materials)
void TPZInterfaceInsertion::SetInterfaceId(int Interface_id){
    m_interface_id = Interface_id;
}

/// Get Interface Identifier
int & TPZInterfaceInsertion::GetInterfaceId(){
    return m_interface_id;
}

/// Set Interface Identifier (interfaces between multiplier and volumetric materials)
void TPZInterfaceInsertion::SetInterfaceVectorId(TPZManVector<int64_t,3> interfaceVector_ids){
    m_interfaceVector_ids = interfaceVector_ids;
}

/// Get Interface Identifier
TPZManVector<int64_t,3> & TPZInterfaceInsertion::GetInterfaceVectorId(){
    return m_interfaceVector_ids;
}

/// Set wrap Identifier
void TPZInterfaceInsertion::SetWrapFluxIdentifier(int wrapFlux){
    m_id_flux_wrap = wrapFlux;
}

/// Get wrap Identifier
int & TPZInterfaceInsertion::GetWrapFluxId(){
    return m_id_flux_wrap;
}

/// Set multiplier material id
void TPZInterfaceInsertion::SetMultiplierMatId(int multiId){
    m_multiplier_id  = multiId;
}

/// Get multiplier material id
int & TPZInterfaceInsertion::GetMultiplierMatId(){
    return m_multiplier_id;
}

/// Set multiplier material id (BC)
void TPZInterfaceInsertion::SetMultiplierBCMatId(int multiBCId){
    m_multiplierBC_id  = multiBCId;
}

/// Get multiplier material id (BC)
int & TPZInterfaceInsertion::GetMultiplierBCMatId(){
    return m_multiplierBC_id;
}

/// Open the connects of a Interface, create dim-1 Interface elements (Hdiv version)
void TPZInterfaceInsertion::InsertHdivBound(int mat_id_flux_wrap){
    
    TPZCompMesh *cmesh = m_cmesh;
#ifdef PZDEBUG
    if (!m_geometry) {
        DebugStop();
    }
#endif
    SetWrapFluxIdentifier(mat_id_flux_wrap);
    TPZGeoMesh *gmesh = m_geometry;
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    
    int64_t nel = m_geometry->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = m_geometry->Element(el);
        int matid = gel->MaterialId();
        
        if (matid != m_multiplier_id) {
            continue;
        }
        
        TPZStack<TPZCompElSide> neigh;
        int nsides = gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        
        gelside.EqualLevelCompElementList(neigh, 0, 0);
        
        if(neigh.size()!=2){
            DebugStop();
        }
        gel->ResetReference();
        neigh[0].Element()->Reference()->ResetReference();
        neigh[1].Element()->Reference()->ResetReference();
        
        //working on element 0
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neigh[0].Element());
            
            if(!intel){
                DebugStop();
            }
            
            intel->LoadElementReference();
            
            int locindex = intel->MidSideConnectLocId(neigh[0].Side());
            TPZConnect &midsideconnect = intel->MidSideConnect(neigh[0].Side());
            if(midsideconnect.NElConnected() != 2)
            {
                DebugStop();
            }
            
            //Duplica um connect
            int64_t index = cmesh->AllocateNewConnect(midsideconnect.NShape(), midsideconnect.NState(), midsideconnect.Order());
            
            intel->SetConnectIndex(locindex, index);
            midsideconnect.DecrementElConnected();
            cmesh->ConnectVec()[index].IncrementElConnected();
            intel->SetSideOrient(neigh[0].Side(), 1);
            
            
            TPZGeoElBC bc(intel->Reference(),neigh[0].Side(),mat_id_flux_wrap);
            cmesh->CreateCompEl(bc.CreatedElement(), index);
            
            TPZCompEl *var = cmesh->Element(index);
            var->Reference()->ResetReference();
            intel->Reference()->ResetReference();
            
            
        }
        
        // working on element 1
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neigh[1].Element());
            
            if(!intel){
                DebugStop();
            }
            
            intel->LoadElementReference();
            
            intel->SetSideOrient(neigh[1].Side(), 1);
            
            int64_t index;
            
            TPZGeoElBC bc(intel->Reference(),neigh[1].Side(),mat_id_flux_wrap);
            cmesh->CreateCompEl(bc.CreatedElement(), index);
            TPZCompEl *var = cmesh->Element(index);
            var->Reference()->ResetReference();
            
            intel->Reference()->ResetReference();
            
        }
    }
        cmesh->ExpandSolution();
    
    
}
    
void TPZInterfaceInsertion::AddMultiphysicsInterfacesLeftNRight(int matfrom)
{
    // create interface elements between the tangent velocity element and the volumetric elements
    
    TPZCompMesh *cmesh = m_cmesh;
    int dim = cmesh->Dimension();
    
    if (! (m_geometry) ) {
        DebugStop();
    }
    
    if (m_interfaceVector_ids.size()==0) {
        DebugStop();
    }
    
    int64_t nel = m_geometry->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = m_geometry->Element(el);
        int matid = gel->MaterialId();
        
        if (matid != matfrom) {
            continue;
        }
        
        if (gel->HasSubElement() == 1) {
            continue;
        }
        
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZCompElSide celside = gelside.Reference();
        
        TPZStack<TPZGeoElSide> neighbourset;
        gelside.AllNeighbours(neighbourset);
        
        gelside.LowerLevelCompElementList2(1);
    
        int nneighs = neighbourset.size();
        
        TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
        LeftElIndices[0]=0;
        RightElIndices[0]=1;
        
        for(int stack_i=0; stack_i <nneighs; stack_i++){
            TPZGeoElSide neigh = neighbourset[stack_i];
            TPZCompElSide celneigh = neigh.Reference();
            if (!celside || !celneigh) {
                DebugStop();
            }
            int64_t neigh_index = neigh.Element()->Index();
            if (neigh.Element()->Dimension()!=2){
                continue;
            }
            
            TPZGeoElBC gbc(gelside,m_interfaceVector_ids[stack_i]);
            int64_t index;
            
            TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
            elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
            
            std::cout << "Created an interface element between volumetric element " << neigh.Element()->Index() <<
            " side " << neigh.Side() <<
            " and interior 1D element " << gelside.Element()->Index() << std::endl;
            
        }
        
        if (nneighs==1) {
            
            LeftElIndices[0]=0;
            RightElIndices[0]=1;
            
            TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
            if(!clarge) DebugStop();
            TPZGeoElSide glarge = clarge.Reference();
            
            TPZGeoElBC gbc(gelside, m_interfaceVector_ids[1]);
            
            int64_t index;
            TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, clarge, celside);
            intface->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
            
            std::cout << "Created an interface element between volumetric element " << glarge.Element()->Index() <<
            " side " << glarge.Side() <<
            " and interior 1D element " << gelside.Element()->Index() << std::endl;
            nneighs++;
        }
        
        if(nneighs!=2){
            DebugStop();
        }

    }
    
}

void TPZInterfaceInsertion::AddMultiphysicsInterfacesLeftNRight2(int matfrom)
{
 
    TPZCompMesh *cmesh = m_cmesh;
    int dim = cmesh->Dimension();
    
    if (! (m_geometry) ) {
        DebugStop();
    }
    
    if (m_interfaceVector_ids.size()==0) {
        DebugStop();
    }
    
    int64_t nel = m_geometry->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = m_geometry->Element(el);
        int meshdim = m_geometry->Dimension();
        int matid = gel->MaterialId();
        
        if (matid != matfrom) {
            continue;
        }
        
        if (gel->HasSubElement() == 1) {
            continue;
        }
        
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZCompElSide celside = gelside.Reference();
        
        TPZStack<TPZGeoElSide> neighbourset;
        gelside.AllNeighbours(neighbourset);
        
        gelside.LowerLevelCompElementList2(1);
        
        int nneighs = neighbourset.size();
        
        TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
        LeftElIndices[0]=0;
        RightElIndices[0]=1;
        
        for(int stack_i=0; stack_i <nneighs; stack_i++){
            TPZGeoElSide neigh = neighbourset[stack_i];
            TPZCompElSide celneigh = neigh.Reference();
            if (!celside ) {
                DebugStop();
            }
            
            
            
            if (neigh.Element()->HasSubElement()) {
              //  DebugStop();
                TPZStack<TPZGeoElSide> subelements;

                TPZStack<TPZGeoElSide> subel;
                neigh.GetAllSiblings(subel);
                
                for (int i_sub =0; i_sub<subel.size(); i_sub++) {
                    
                    TPZCompElSide cel_sub_neigh = subel[i_sub].Reference();
                    
                    TPZGeoElBC gbc_sub(subel[i_sub],m_interfaceVector_ids[stack_i]);
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc_sub.CreatedElement(),index,cel_sub_neigh,celside);
                    elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                    
                    std::cout << "****Created an interface element between volumetric element " << subel[i_sub].Element()->Index() <<
                    " side " << subel[i_sub].Side() <<
                    " and interior 1D element " << gelside.Element()->Index() << std::endl;

                    
                }
                
                
            }else{
               
                int64_t neigh_index = neigh.Element()->Index();
                if (neigh.Element()->Dimension()!=meshdim){
                    continue;
                }
                
                TPZGeoElBC gbc(gelside,m_interfaceVector_ids[stack_i]);
                int64_t index;
                
                TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
                elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                
                std::cout << "Created an interface element between volumetric element " << neigh.Element()->Index() <<
                " side " << neigh.Side() <<
                " and interior 1D element " << gelside.Element()->Index() << std::endl;
                
            }
            
            
            
            
            
        }
        
        
        
//        if (nneighs==1) {
//
//            LeftElIndices[0]=0;
//            RightElIndices[0]=1;
//
//            TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
//            if(!clarge) DebugStop();
//            TPZGeoElSide glarge = clarge.Reference();
//
//            TPZGeoElBC gbc(gelside, m_interfaceVector_ids[1]);
//
//            int64_t index;
//            TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, clarge, celside);
//            intface->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
//
//            std::cout << "Created an interface element between volumetric element " << glarge.Element()->Index() <<
//            " side " << glarge.Side() <<
//            " and interior 1D element " << gelside.Element()->Index() << std::endl;
//            nneighs++;
//        }
//
//        if(nneighs!=2){
//            DebugStop();
//        }
        
    }
    
}



void TPZInterfaceInsertion::AddMultiphysicsInterfaces()
{
    // create interface elements between the tangent velocity element and the volumetric elements
    {
        
        TPZCompMesh *cmesh = m_geometry->Reference();
        TPZGeoMesh *gmesh = m_geometry;
        std::set<int> velmatid;
        
        velmatid = m_boundaries_ids;
        velmatid.insert(m_multiplier_id);

        int64_t nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            int matid = gel->MaterialId();
            if(velmatid.find(matid) != velmatid.end())
            {
                int nsides = gel->NSides();
                TPZGeoElSide gelside(gel,nsides-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    
                    TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
                    LeftElIndices[0]=0;
                    RightElIndices[0]=0;
                    
                    if (neighbour.Element()->Dimension() == 2 && gelside.Element()->Dimension() == 1) { //oioioi IDFlux -> ID
                        // create an interface element
                        TPZCompElSide celside = gelside.Reference();
                        TPZCompElSide celneigh = neighbour.Reference();
                        if (!celside || !celneigh) {
                            DebugStop();
                        }
                        std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
                        " side " << neighbour.Side() <<
                        " and interface element " << gelside.Element()->Index() << std::endl;
                        TPZGeoElBC gelbc(gelside,m_interface_id);
                        int64_t index;
                        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(*cmesh,gelbc.CreatedElement(),index,celneigh,celside);
                        intf->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                        
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
        }
    }
    
}


void TPZInterfaceInsertion::AddMultiphysicsInterfaces(int matfrom, int mattarget)
{
    if (!(m_geometry&&m_cmesh)) {
        DebugStop();
    }
    
    TPZGeoMesh *gmesh = m_geometry;
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != matfrom) {
            continue;
        }
        
        int nsides= gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 2) {
            DebugStop();
        }
        gel->SetMaterialId(mattarget);
        int64_t index;
        new TPZMultiphysicsInterfaceElement(*m_cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}

void TPZInterfaceInsertion::AddMultiphysicsBCInterface(int matfrom, int matBCinterface)
{

    TPZCompMesh *cmesh = m_cmesh;
    
    if (! (m_geometry) ) {
        DebugStop();
    }
    
    if (m_interfaceVector_ids.size()==0) {
        DebugStop();
    }
    
    int64_t nel = m_geometry->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = m_geometry->Element(el);
        int matid = gel->MaterialId();
        
        if (matid != matfrom) {
            continue;
        }
        
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZCompElSide celside = gelside.Reference();
        
        TPZStack<TPZGeoElSide> neighbourset;
        gelside.AllNeighbours(neighbourset);
        
        int nneighs = neighbourset.size();
        if(nneighs!=2&&m_Eltype==EQuadrilateral){
            DebugStop();
        }
        if(nneighs!=3&&m_Eltype==ETriangle){
      //      DebugStop();
        }
        
        TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
        LeftElIndices[0]=0;
        RightElIndices[0]=1;
        
        for(int stack_i=0; stack_i <nneighs; stack_i++){
            TPZGeoElSide neigh = neighbourset[stack_i];
            
            TPZCompElSide celneigh = neigh.Reference();
            if (!celside || !celneigh) {
            //    DebugStop();
            }
            int64_t neigh_index = neigh.Element()->Index();
            if (neigh.Element()->Dimension()!=2){
                continue;
            }
            
            TPZGeoElBC gbc(gelside,matBCinterface);
            int64_t index;
            
            TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
            elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
            
            
            std::cout << "Created an BC interface element between volumetric element " << neigh.Element()->Index() <<
            " side " << neigh.Side() <<
            " and boundary 1D element " << gelside.Element()->Index() << std::endl;
            
        }
    }
    
    
}

void TPZInterfaceInsertion::AddMultiphysicsBCInterface2(int matfrom, int matBCinterface)
{
    
    TPZCompMesh *cmesh = m_cmesh;
    
    if (! (m_geometry) ) {
        DebugStop();
    }
    
    if (m_interfaceVector_ids.size()==0) {
        DebugStop();
    }
    
    int64_t nel = m_geometry->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = m_geometry->Element(el);
        int meshdim = m_geometry->Dimension();
        int matid = gel->MaterialId();
        
        if (matid != matfrom) {
            continue;
        }
        
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZCompElSide celside = gelside.Reference();
        
        TPZStack<TPZGeoElSide> neighbourset;
        gelside.AllNeighbours(neighbourset);
        
        int nneighs = neighbourset.size();
        if(nneighs!=2&&m_Eltype==EQuadrilateral){
        //    DebugStop();
        }
        if(nneighs!=3&&m_Eltype==ETriangle){
            //      DebugStop();
        }
        
        TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
        LeftElIndices[0]=0;
        RightElIndices[0]=1;
        
        for(int stack_i=0; stack_i <nneighs; stack_i++){
            TPZGeoElSide neigh = neighbourset[stack_i];
            
            TPZCompElSide celneigh = neigh.Reference();
            if (!celside || !celneigh) {
                //    DebugStop();
            }
            int64_t neigh_index = neigh.Element()->Index();
            if (neigh.Element()->Dimension()!=meshdim){
                continue;
            }
            
            if (neigh.Element()->HasSubElement()) {
                
                TPZStack<TPZGeoElSide > subelside;
                neigh.GetAllSiblings(subelside);
                
                for (int i_sub = 0; i_sub<subelside.size(); i_sub++) {
 
                    TPZCompElSide cel_sub_neigh = subelside[i_sub].Reference();
                    
                    TPZGeoElBC gbc_sub(subelside[i_sub],matBCinterface);
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc_sub.CreatedElement(),index,cel_sub_neigh,celside);
                    elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                    
                    
                    std::cout << "****Created an BC interface element between volumetric element " << subelside[i_sub].Element()->Index() <<
                    " side " << subelside[i_sub].Side() <<
                    " and boundary 1D element " << gelside.Element()->Index() << std::endl;
                }
                
            }else{
                
                TPZGeoElBC gbc(gelside,matBCinterface);
                int64_t index;
                
                TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
                elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                
                
                std::cout << "Created an BC interface element between volumetric element " << neigh.Element()->Index() <<
                " side " << neigh.Side() <<
                " and boundary 1D element " << gelside.Element()->Index() << std::endl;
                
            }
            

            
        }
    }
    
    
}



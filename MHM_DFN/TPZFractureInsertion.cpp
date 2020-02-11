
//  TPZFractureInsertion.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 02/08/18.
//

#include "TPZFractureInsertion.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "pzcompelwithmem.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"


/// Default constructor
TPZFractureInsertion::TPZFractureInsertion(){
    
}

/// Default desconstructor
TPZFractureInsertion::~TPZFractureInsertion(){
    
}

/// Copy constructor
TPZFractureInsertion::TPZFractureInsertion(TPZFractureInsertion & other){
    m_fracture_id               = other.m_fracture_id;
    m_geometry                  = other.m_geometry;
    m_boundaries_material_ids   = other.m_boundaries_material_ids;
    m_pivot_indexes             = other.m_pivot_indexes;
    m_non_pivot_indexes         = other.m_non_pivot_indexes;
    m_gel_left_indexes          = other.m_gel_left_indexes;
    m_gel_right_indexes         = other.m_gel_right_indexes;
    mat_id_flux_wrap            = other.mat_id_flux_wrap;
}

/// Set fracture Identifier
void TPZFractureInsertion::SetFractureIdentifier(int fracture_id){
    m_fracture_id = fracture_id;
}

/// Get fracture Identifier
int & TPZFractureInsertion::GetFractureMaterialId(){
    return m_fracture_id;
}

/// Set wrap Identifier
void TPZFractureInsertion::SetWrapFluxIdentifier(int wrapFlux){
    mat_id_flux_wrap = wrapFlux;
}

/// Get wrap Identifier
int & TPZFractureInsertion::GetWrapFluxId(){
    return mat_id_flux_wrap;
}

/// Get node pivots
std::vector<TPZGeoElSide> & TPZFractureInsertion::GetPivotIndexes(){
    return m_pivot_indexes;
}

/// Get node non pivots
std::vector<TPZGeoElSide> & TPZFractureInsertion::GetNonPivotIndexes(){
    return m_non_pivot_indexes;
}

/// Get geometric fracture indexes
std::vector<int64_t> & TPZFractureInsertion::GetFractureIndexes(){
    return m_fracture_indexes;
}

/// Get geometric indexes for left
std::set<int64_t> & TPZFractureInsertion::GetLeftIndexes(){
    return m_gel_left_indexes;
}

/// Get geometric indexes for right
std::set<int64_t> & TPZFractureInsertion::GetRightIndexes(){
    return m_gel_right_indexes;
}

/// Constructor based on a computational mesh and fracture material id
TPZFractureInsertion::TPZFractureInsertion(TPZGeoMesh * geometry, int fracture_id, std::set<int> & boundaries_ids){
    m_geometry = geometry;
    m_fracture_id = fracture_id;
    m_boundaries_material_ids = boundaries_ids;
    BuildFractureElements();
}

void TPZFractureInsertion::BuildFractureElements(){
    
    int64_t n_gel = m_geometry->NElements();
    for (int64_t iel=0; iel < n_gel; iel++) {

        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(iel);
        if(gel->MaterialId() != m_fracture_id || gel->Dimension()!=m_geometry->Dimension()-1 || gel->HasSubElement() == 1){
            continue;
        }
        m_fracture_indexes.push_back(gel->Index());
    }
    
    if(m_fracture_indexes.size() == 0){
        std::cout << "No fracture elements with id = " << m_fracture_id << "  were detected." << std::endl;
        DebugStop();
    }
    
}

void TPZFractureInsertion::BuildPivotDataStructure(){
    
    for (auto ifrac : m_fracture_indexes) {
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(ifrac);
        unsigned int n_corner_sides = gel->NCornerNodes();
        for (unsigned int i_side = 0; i_side < n_corner_sides; i_side++) {
            
            TPZStack<TPZGeoElSide> all_neighbors;
            TPZGeoElSide gel_corner_side(gel,i_side);
            gel_corner_side.AllNeighbours(all_neighbors);
            
            unsigned int d_minus_bc_neighbors_counter = 0;
            unsigned int d_minus_same_fracture_neighbors_counter = 0;
            
            bool has_point_neighbor_Q = false;
            unsigned int n_neighbors = all_neighbors.size();
            for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
                TPZGeoEl * gel_neighbor = all_neighbors[i_neighbor].Element();
#ifdef PZDEBUG
                if (!gel_neighbor) {
                    DebugStop();
                }
#endif
                
                if (gel_neighbor->HasSubElement() == 1) {
                    continue;
                }
                
                if (gel_neighbor->Dimension() == m_geometry->Dimension()-1) {
                    int material_id = gel_neighbor->MaterialId();
                    bool is_boundary_memeber_Q = m_boundaries_material_ids.find(material_id) != m_boundaries_material_ids.end();
                    if (is_boundary_memeber_Q) {
                        d_minus_bc_neighbors_counter++;
                    }else{
                        if (material_id == m_fracture_id) {
                            d_minus_same_fracture_neighbors_counter++;
                        }
                    }
                }
                
                if (gel_neighbor->Dimension() == m_geometry->Dimension()-2) {
                    

                    has_point_neighbor_Q = true;
                }
                
            }
            
            bool is_boundary_neighbour_Q = (d_minus_bc_neighbors_counter != 0);
            bool is_the_end_of_fractrue_Q = d_minus_same_fracture_neighbors_counter == 0 ;
            bool is_internal_side_Q = !is_the_end_of_fractrue_Q || d_minus_same_fracture_neighbors_counter != 0;
            
            bool is_pivot_Q = (is_boundary_neighbour_Q && is_the_end_of_fractrue_Q) || is_internal_side_Q;
            
            if (is_pivot_Q) {
                bool is_already_inserted_Q = HasInsertedPivot(gel_corner_side);
                if (!is_already_inserted_Q) {
                    m_pivot_indexes.push_back(gel_corner_side);
                }
            }else{
                bool is_already_inserted_Q = HasInsertedNonPivot(gel_corner_side);
                if (!is_already_inserted_Q) {
                    m_non_pivot_indexes.push_back(gel_corner_side);
                }
            }
        }
    }
    
    if(m_pivot_indexes.size() == 0){
        std::cout << "No pivots were detected." << std::endl;
        DebugStop();
    }
    
}

/// Verify if the non pivot has been inserted on non pivot structure
bool TPZFractureInsertion::HasInsertedNonPivot(TPZGeoElSide pivot_side){
    
    TPZGeoNode node_candidate = pivot_side.Element()->Node(pivot_side.Side());
    
    bool has_inserted_non_pivot_Q = false;
    for (auto non_pivot: m_non_pivot_indexes) {
        TPZGeoNode node_target = non_pivot.Element()->Node(non_pivot.Side());
        if (node_candidate.Id() == node_target.Id()) {
            has_inserted_non_pivot_Q = true;
        }
    }

    return has_inserted_non_pivot_Q;
    
}

/// Verify if the pivot has been inserted on pivot structure
bool TPZFractureInsertion::HasInsertedPivot(TPZGeoElSide pivot_side){
    
    TPZGeoNode node_candidate = pivot_side.Element()->Node(pivot_side.Side());
    
    bool has_inserted_pivot_Q = false;
    for (auto pivot: m_pivot_indexes) {
        TPZGeoNode node_target = pivot.Element()->Node(pivot.Side());
        if (node_candidate.Id() == node_target.Id()) {
            has_inserted_pivot_Q = true;
        }
    }
    
    return has_inserted_pivot_Q;
    
}

std::set<int64_t> TPZFractureInsertion::PivotNeighbours(TPZGeoElSide pivotside){
    
    TPZStack<TPZGeoElSide> all_neighbors;
    pivotside.AllNeighbours(all_neighbors);
    std::set<int64_t> neigh_indexes;
    for (int i = 0; i < all_neighbors.size(); i++) {
        
        if (all_neighbors[i].HasSubElement() == 1) {
            continue;
        }
//        if (all_neighbors[i].Element()!=all_neighbors[i].Element()->LowestFather()) {
//            continue;
//        }
        
        if(all_neighbors[i].Element()->Dimension() == m_geometry->Dimension()){
            neigh_indexes.insert(all_neighbors[i].Element()->Index());
        }
        int matid = all_neighbors[i].Element()->MaterialId();
        if(m_boundaries_material_ids.find(matid) != m_boundaries_material_ids.end())
        {
            neigh_indexes.insert(all_neighbors[i].Element()->Index());
        }
        
    }
    return neigh_indexes;
}

/// Find the 2 elements that will initiate the element sets
// m_gel_left_indexes and m_gel_right_indexes have to be empty
// take a neighbouring fracture element - insert its two neighbours along the face
void TPZFractureInsertion::InsertFractureNeighbours(std::set<int64_t> pivot_neighbours){
    
    
    // this is the index of the fracture element
    int64_t fractture_index = m_fracture_indexes[0];
    TPZGeoEl * fracture_element = m_geometry->Element(fractture_index);
    int count = 0;
    /// identify left and right seeds
    TPZStack<TPZGeoElSide> all_neigh_t,all_neigh;
    TPZGeoElSide fracture_side(fracture_element,fracture_element->NSides()-1);
    fracture_side.AllNeighbours(all_neigh_t);
    if (all_neigh_t.size() != 2) {
        for (int i = 0; i < all_neigh_t.size(); i++) {
            if (all_neigh_t[i].Element()->MaterialId()!=501) {
//                break;
//            }else{
                all_neigh.resize(count+1);
                all_neigh[count] = all_neigh_t[i];
                count++;
            }
        }
        if (all_neigh.size() != 2) {
            DebugStop();
        }
    }
    
    // verify if the center of the neighbour is to the left or right of the element
    
    // Fist, define fracture reference vector:
    int frac_dim = fracture_side.Element()->Dimension();
    
    TPZVec<REAL> coX_frac_side0(frac_dim,0.), coord_frac_side0(frac_dim,0);
    TPZVec<REAL> coX_frac_side1(frac_dim,0.), coord_frac_side1(frac_dim,0);

    fracture_element->CenterPoint(0, coX_frac_side0);
    fracture_element->CenterPoint(1, coX_frac_side1);
    fracture_element->X(coX_frac_side0, coord_frac_side0);
    fracture_element->X(coX_frac_side1, coord_frac_side1);
    
    
    // Second, find center coordinates of fracture neighbour
    TPZGeoEl * gel_neigh0 = all_neigh[0].Element();
    TPZGeoEl * gel_neigh1 = all_neigh[1].Element();
    
    int neigh0_dim = gel_neigh0->Dimension();
    int neigh1_dim = gel_neigh1->Dimension();
    
    TPZVec<REAL> coX_neigh0(neigh0_dim,0.), coord_neigh0(3,0);
    int gel_neigh0_side = gel_neigh0->NSides()-1;
    gel_neigh0->CenterPoint(gel_neigh0_side, coX_neigh0);
    gel_neigh0->X(coX_neigh0, coord_neigh0);
    
    // classify left right using fracture_neighbour - gel_left + gel_right
    if( ( (coord_frac_side1[0]-coord_neigh0[0])>0 && (coord_frac_side0[1]-coord_neigh0[1]) < 0) || ( (coord_frac_side1[0]-coord_neigh0[0]) <= 0 && (coord_frac_side0[1]-coord_neigh0[1]) >= 0 ) ) {
        neigh_is_left = true;
    }else{
        neigh_is_left = false;
    }
    
    if(neigh_is_left){
        m_gel_left_indexes.insert(gel_neigh0->Index());
        m_gel_right_indexes.insert(gel_neigh1->Index());
    }else{
        m_gel_left_indexes.insert(gel_neigh1->Index());
        m_gel_right_indexes.insert(gel_neigh0->Index());
    }
    
}


void TPZFractureInsertion::ClassifyElements(std::set<int64_t> pivotneighbours)
{
    while(pivotneighbours.size())
    {
        std::set<int64_t> found;
        for (auto elindex : pivotneighbours) {
            bool has_left_index_neighbor_Q = HasLeftIndexNeighbour(elindex);
            bool has_right_index_neighbor_Q = HasRightIndexNeighbour(elindex);
            if(has_left_index_neighbor_Q) m_gel_left_indexes.insert(elindex);
            if(has_right_index_neighbor_Q) m_gel_right_indexes.insert(elindex);
            if (has_left_index_neighbor_Q || has_right_index_neighbor_Q) {
                found.insert(elindex);
            }
        }
        for (auto elindex: found) {
            pivotneighbours.erase(elindex);
        }
    }
}


/// look if there is a fracture element neighbour (ok)
int64_t TPZFractureInsertion::FractureNeighbourIndex(TPZGeoElSide gelside){

#ifdef PZDEBUG
    {
        int meshdim = m_geometry->Dimension();
        if (gelside.Dimension() != meshdim -1) {
            DebugStop();
        }
    }
#endif
    
    if (gelside.Element()->HasSubElement() == 1) {
        return -1;
    }
    
    TPZGeoElSide neighbour = gelside.Neighbour();
    while(neighbour != gelside)
    {
        if(neighbour.Element()->MaterialId() == m_fracture_id) return neighbour.Element()->Index();
        neighbour = neighbour.Neighbour();
    }
    return -1;
    
}

bool TPZFractureInsertion::HasLeftIndexNeighbour(int64_t gel_index){
    
    TPZGeoEl * gel = m_geometry->Element(gel_index);
    int meshdim = m_geometry->Dimension();
    unsigned int n_corner_sides = gel->NCornerNodes();
    unsigned int gel_candidate_side = gel->NSides();
    
    if (gel->HasSubElement() == 1) {
        return true;
    }
    
    if (m_gel_left_indexes.find(gel_index) != m_gel_left_indexes.end()) {
        return true;
    }
    
    for(int side=n_corner_sides; side<gel_candidate_side; side++)
    {
        TPZGeoElSide gelside(gel,side);
        if(gelside.Dimension() != meshdim-1) continue;
        if (FractureNeighbourIndex(gelside) != -1) {
            continue;
        }
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            int64_t neighindex = neighbour.Element()->Index();
            if (m_gel_left_indexes.find(neighindex) != m_gel_left_indexes.end()) {
                return true;
            }
            if (neighbour.Element()->HasSubElement()) {  // procurando entre os subelementos
                TPZStack<TPZGeoEl *> unrefinedSons;
                neighbour.Element()->YoungestChildren(unrefinedSons);
                for (auto sonindex : unrefinedSons) {
                    if (m_gel_left_indexes.find(sonindex->Index()) != m_gel_left_indexes.end()) {
                        return true;
                    }
                }
            }
            neighbour = neighbour.Neighbour();
        }
    }
    return false;
    
}

bool TPZFractureInsertion::HasRightIndexNeighbour(int64_t gel_index){
    
    TPZGeoEl * gel = m_geometry->Element(gel_index);
    unsigned int n_corner_sides = gel->NCornerNodes();
    unsigned int gel_candidate_side = gel->NSides() ;
    int meshdim = m_geometry->Dimension();
    if (gel->HasSubElement() == 1) {
        return true;
    }
    
    if (m_gel_right_indexes.find(gel_index) != m_gel_right_indexes.end()) {
        return true;
    }
    
    for(int side=n_corner_sides; side<gel_candidate_side; side++)
    {
        TPZGeoElSide gelside(gel,side);
        if(gelside.Dimension() != meshdim-1) continue;
        if (FractureNeighbourIndex(gelside) != -1) {
            continue;
        }
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            int64_t neighindex = neighbour.Element()->Index();
            if (m_gel_right_indexes.find(neighindex) != m_gel_right_indexes.end()) {
                return true;
            }
            if (neighbour.Element()->HasSubElement()) {  // procurando entre os subelementos
                TPZStack<TPZGeoEl *> unrefinedSons;
                neighbour.Element()->YoungestChildren(unrefinedSons);
                for (auto sonindex : unrefinedSons) {
                    if (m_gel_right_indexes.find(sonindex->Index()) != m_gel_right_indexes.end()) {
                        return true;
                    }
                }
            }
            
            neighbour = neighbour.Neighbour();
        }
    }
    return false;
    
}

bool TPZFractureInsertion::HasClassifiedNeighbour(TPZGeoElSide & pivot_side){
    
    std::set<int64_t> neighbors = PivotNeighbours(pivot_side);
    
    bool has_classified_Q = false;
    bool has_left_index_Q;
    bool has_right_index_Q;
    for (auto index: neighbors) {
        if(HasLeftIndexNeighbour(index)){
            has_left_index_Q  = true;
        };
        if(HasRightIndexNeighbour(index)){
            has_right_index_Q  = true;
        };
    }
    if (has_left_index_Q && has_right_index_Q) {
        has_classified_Q = true;
    }
    
    return has_classified_Q;
    
}


/// Classify the neighbouring elements of the pivots
void TPZFractureInsertion::ClassifyNeighboursofPivots(){

    std::ofstream filegvtk_b("Geometry_labels.vtk"); //Impressão da malha geométrica (formato vtk)
    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry, filegvtk_b,true);
    
    BuildPivotDataStructure();
    
    // Computes the seeds
    TPZGeoElSide pivot = m_pivot_indexes[0];
    std::set<int64_t> neighbours =  PivotNeighbours(pivot);
    InsertFractureNeighbours(neighbours);
    
    // Auxiliary list for account pivots that are processed
    std::set<TPZGeoElSide> pivots;
    for (auto gel_side: m_pivot_indexes) {
        pivots.insert(gel_side);
    }
    
    while (pivots.size()) {
        std::set<TPZGeoElSide> to_be_deleted;
        for (auto pivot: pivots) {
            // return true if their is at least one element neighbour that is classified
            bool has_classified_neighbor_Q = HasClassifiedNeighbour(pivot);
            if (has_classified_neighbor_Q) {
                std::set<int64_t> neighbours =  PivotNeighbours(pivot);
                ClassifyElements(neighbours);
                to_be_deleted.insert(pivot);
            }
        }
        for (auto pivot: to_be_deleted) {
            pivots.erase(pivot);
        }
    }
    
    //Definição dos Ids dos materiais vizinhos(Print)
    
//    int left_id = 100;
//    for (auto it=m_gel_left_indexes.begin(); it != m_gel_left_indexes.end(); ++it){
//        int64_t iel = *it;
//        m_geometry->Element(iel)->SetMaterialId(left_id);
//    }
//
//    int right_id = 200;
//    for (auto it=m_gel_right_indexes.begin(); it != m_gel_right_indexes.end(); ++it){
//        int64_t iel = *it;
//        m_geometry->Element(iel)->SetMaterialId(right_id);
//    }
//
//    std::ofstream filegvtk("Geometry_labels.vtk"); //Impressão da malha geométrica (formato vtk)
//    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry, filegvtk,true);
    
}

/// Open the connects of a fracture, create dim-1 fracture elements
void  TPZFractureInsertion::OpenFractureOnH1(TPZCompMesh *cmesh){
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    int meshdim = gmesh->Dimension();
    cmesh->LoadReferences();
    std::map<int64_t,int64_t> connectmap_to_be_duplicated;
    std::set<int64_t> geonodes;
    for (auto gelside : m_pivot_indexes) {
        geonodes.insert(gelside.SideNodeIndex(0));
    }
    for (auto ifrac : m_fracture_indexes) {
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(ifrac);
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        unsigned int n_sides = gel->NSides();
        for (unsigned int i_side = 0; i_side < n_sides; i_side++) {
            
            int64_t nodeindex = gel->NodeIndex(i_side);
            
            bool is_non_pivot_Q = geonodes.find(nodeindex) == geonodes.end();
            if (i_side < gel->NCornerNodes() && is_non_pivot_Q) { // Filtering just non pivots sides
                continue;
            }
            
            // For all pivots and fracture itself side
            TPZStack<TPZGeoElSide> all_neighbors;
            TPZGeoElSide gel_corner_side(gel,i_side);
            gel_corner_side.AllNeighbours(all_neighbors);
            unsigned int n_neighbors = all_neighbors.size();
            for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
                TPZGeoEl * gel_neighbor = all_neighbors[i_neighbor].Element();
                if(gel_neighbor->Dimension() != meshdim) continue;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(gel_neighbor->Reference());
                if(!intel) continue;
                
                bool is_not_neigh_bound = m_boundaries_material_ids.find(gel_neighbor->MaterialId()) == m_boundaries_material_ids.end();
                
                if (is_not_neigh_bound && gel_neighbor->Dimension() != meshdim) {
                    continue;
                }
                
//                int neigh_side = all_neighbors[i_neighbor].Side();
//                int nconnec = intel->NConnects();
//                if (nconnec==1) { //aquipapap
//                    continue;
//                }
                int64_t conindex = intel->ConnectIndex(all_neighbors[i_neighbor].Side());
                connectmap_to_be_duplicated[conindex] = -1;
            }
        
        }
    }
    for (auto indexpair :connectmap_to_be_duplicated) {
        int64_t conindex = indexpair.first;
        TPZConnect &cnew (cmesh->ConnectVec()[conindex]);
        int64_t newindex = cmesh->AllocateNewConnect(cnew);
        connectmap_to_be_duplicated[conindex] = newindex;
    }
    
    for (auto ifrac : m_fracture_indexes) {
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(ifrac);
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        unsigned int n_corner_sides = gel->NSides();
        for (unsigned int i_side = 0; i_side < n_corner_sides; i_side++) {
            
            int64_t nodeindex = gel->NodeIndex(i_side);
            
            bool is_non_pivot_Q = geonodes.find(nodeindex) == geonodes.end();
            if (i_side < gel->NCornerNodes() && is_non_pivot_Q) { // Filtering just non pivots sides
                continue;
            }
            
            TPZStack<TPZGeoElSide> all_neighbors;
            TPZGeoElSide gel_corner_side(gel,i_side);
            gel_corner_side.AllNeighbours(all_neighbors);
            unsigned int n_neighbors = all_neighbors.size();
            for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
                TPZGeoEl * gel_neighbor = all_neighbors[i_neighbor].Element();
                int64_t gelindex = gel_neighbor->Index();
                bool is_not_left_element_Q = m_gel_left_indexes.find(gelindex) == m_gel_left_indexes.end();
                if(is_not_left_element_Q)
                {
                    continue;
                }
                
                // Set duplicated connect for left side elements
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(gel_neighbor->Reference());
                if(!intel) continue;
                int64_t conindex = intel->ConnectIndex(all_neighbors[i_neighbor].Side());
                if(connectmap_to_be_duplicated.find(conindex) == connectmap_to_be_duplicated.end())
                {
                    continue;
                }
                int64_t newindex = connectmap_to_be_duplicated[conindex];
                intel->SetConnectIndex(all_neighbors[i_neighbor].Side(), newindex);
            }            
        }
    }
    
    cmesh->InitializeBlock();
    cmesh->ExpandSolution();

}

/// Open the connects of a fracture, create dim-1 fracture elements (Hdiv version)
void TPZFractureInsertion::AddWraps(TPZCompMesh *cmesh, int mat_id_flux_wrap){
    
#ifdef PZDEBUG
    if (!cmesh) {
        DebugStop();
    }
#endif
    SetWrapFluxIdentifier(mat_id_flux_wrap);
    TPZGeoMesh *gmesh = cmesh->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    
    for (auto ifrac : m_fracture_indexes) {
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(ifrac);
            
        if (gel->HasSubElement()) {
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
            //    DebugStop();
            }
            
            //Duplica um connect
            int64_t index = cmesh->AllocateNewConnect(midsideconnect.NShape(), midsideconnect.NState(), midsideconnect.Order());
            
            intel->SetConnectIndex(locindex, index);
            midsideconnect.DecrementElConnected();
            cmesh->ConnectVec()[index].IncrementElConnected();
            
            
            
           // intel->SetSideOrient(neigh[0].Side(), 1);
            
        
            TPZGeoElBC bc(intel->Reference(),neigh[0].Side(),mat_id_flux_wrap);
            cmesh->CreateCompEl(bc.CreatedElement(), index);
            
            TPZInterpolationSpace *bound;
            bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(*intel->Mesh(),bc.CreatedElement(),index);
            int sideorient = neigh[0].Side();
            TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
            hdivbound->SetSideOrient(pzshape::TPZShapeLinear::NSides-1,sideorient);

            
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
            
//            intel->SetSideOrient(neigh[1].Side(), 1);
            
            int64_t index;
            
            TPZGeoElBC bc(intel->Reference(),neigh[1].Side(),mat_id_flux_wrap+1);
            cmesh->CreateCompEl(bc.CreatedElement(), index);
            
            
            TPZInterpolationSpace *bound;
            bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(*intel->Mesh(),bc.CreatedElement(),index);
            int sideorient = neigh[1].Side();
            TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
            hdivbound->SetSideOrient(pzshape::TPZShapeLinear::NSides-1,sideorient);
            
            
            
            TPZCompEl *var = cmesh->Element(index);
            var->Reference()->ResetReference();
            
            intel->Reference()->ResetReference();
            
        }
        
    }
    
    cmesh->ExpandSolution();
    
}
    

/// Set Discontinuous elements on fractures
void TPZFractureInsertion::SetDiscontinuosFrac(TPZCompMesh *cmesh){
    
    int meshdim = cmesh->Dimension();
    cmesh->SetDimModel(meshdim);
    //cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    //cmesh->ApproxSpace().CreateDisconnectedElements(true);
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();

    int cmesh_order = cmesh->GetDefaultOrder();
    int fracture_cel_order = cmesh->GetDefaultOrder()-1;
    cmesh->SetDefaultOrder(fracture_cel_order); // Because we need p - 1 order on fractures
    
    // Created fracture elements based on fracture material with id  GetFractureMaterialId()
    for (auto ifrac : m_fracture_indexes) {
            
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(ifrac);
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        int64_t index;
        cmesh->CreateCompEl(gel, index);
        TPZCompEl *cel = cmesh->Element(index);
        cel->SetIntegrationRule(4);
        cel->PrepareIntPtIndices();
        int64_t cindex = cel->ConnectIndex(0);
        cel->SetgOrder(fracture_cel_order);
        int lagrangelevel = 1;
        cmesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
        cmesh->ConnectVec()[cindex].SetOrder(fracture_cel_order, cindex);
    }
    
    cmesh->Reference()->ResetReference();
    cmesh->ExpandSolution();
    cmesh->SetDefaultOrder(cmesh_order); // Because we desire to preserve the original p order for the cmesh
    
}

/// Set interfaces elements between fracture and volumetric elements
void TPZFractureInsertion::SetInterfaces(TPZCompMesh *cmesh, int matInterfaceLeft, int matInterfaceRight){
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    std::vector<int64_t> fracture_index = GetFractureIndexes();
    
    int fracture_id = GetFractureMaterialId();
    TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh->Reference()->ElementVec();
    int meshdim = cmesh->Dimension();
    int64_t nelem = elvec.NElements();
    
    int64_t index;
    
    std::set<int64_t> left_el_indexes = GetLeftIndexes();
    std::set<int64_t> right_el_indexes = GetRightIndexes();
    
    for (int iel = 0; iel < fracture_index.size(); iel++) {
        TPZGeoEl *gel = elvec[fracture_index[iel]];
        
        const int gelMatId = gel->MaterialId();
#ifdef PZDEBUG
        if (gelMatId != fracture_id)
            DebugStop();
        if (gel->HasSubElement())
            DebugStop();
#endif
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        if (!celside) {
            DebugStop();
        }
        //        TPZGeoElSide neigh = gelside.Neighbour();
        //        int64_t neigh_index = neigh.Element()->Index();
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 2) {
            DebugStop();
        }
        for(int stack_i=0; stack_i <celstack.size(); stack_i++){
            TPZGeoElSide neigh = celstack[stack_i].Reference();
            int64_t neigh_index = neigh.Element()->Index();
            
            if(left_el_indexes.find(neigh_index)!=left_el_indexes.end()){
                
                TPZGeoElBC gbcleft(gelside,matInterfaceLeft);
                int64_t index;
                auto elem_left  = new TPZCompElWithMem<TPZInterfaceElement>(*cmesh, gbcleft.CreatedElement(),index);
                elem_left->SetLeftRightElements(celside, celstack[stack_i]);
                //new TPZInterfaceElement(*cmesh,gbcleft.CreatedElement(),index,celside,celstack[stack_i]);
                
            }else if((right_el_indexes.find(neigh_index) != right_el_indexes.end())){
                
                TPZGeoElBC gbcright(gelside,matInterfaceRight);
                int64_t index;
                auto elem_right  = new TPZCompElWithMem<TPZInterfaceElement>(*cmesh, gbcright.CreatedElement(),index);
                elem_right->SetLeftRightElements(celside, celstack[stack_i]);
                //new TPZInterfaceElement(*cmesh,gbcright.CreatedElement(),index,celside,celstack[stack_i]);
                
            }else{
                DebugStop();
            }
        }
        
    }

    
    
}



/// Set interfaces elements between fracture and volumetric elements
void TPZFractureInsertion::SetMultiphysicsInterfaces(TPZCompMesh *cmesh, int matInterfaceLeft, int matInterfaceRight, int mat_id_flux_wrap){
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    std::vector<int64_t> fracture_index = GetFractureIndexes();
    
    int fracture_id = GetFractureMaterialId();
    TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh->Reference()->ElementVec();
    int meshdim = cmesh->Dimension();
    int64_t nelem = elvec.NElements();
    
    int64_t index;
    
    std::set<int64_t> left_el_indexes = GetLeftIndexes();
    std::set<int64_t> right_el_indexes = GetRightIndexes();
    
    for (int iel = 0; iel < fracture_index.size(); iel++) {
        TPZGeoEl *gel = elvec[fracture_index[iel]];
        
        const int gelMatId = gel->MaterialId();
#ifdef PZDEBUG
        if (gelMatId != fracture_id)
            DebugStop();
        if (gel->HasSubElement())
            DebugStop();
#endif
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        TPZGeoElSide neighbour = gelside.Neighbour();
        if (!celside) {
            DebugStop();
        }
        //        TPZGeoElSide neigh = gelside.Neighbour();
        //        int64_t neigh_index = neigh.Element()->Index();
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 4) {
            DebugStop();
        }
        
        TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
        LeftElIndices[0]=0;
        RightElIndices[0]=1;
        
        
        for(int stack_i=0; stack_i <celstack.size(); stack_i++){
            TPZGeoElSide neigh = celstack[stack_i].Reference();
            int64_t neigh_index = neigh.Element()->Index();
            if (neigh.Element()->Dimension()!=2){
                continue;
            }
            
            TPZGeoElSide wrap_neigh = neigh.Neighbour();
            if (wrap_neigh.Element()->Dimension() != 1 && wrap_neigh.Element()->MaterialId() != mat_id_flux_wrap){
                continue;
            }

            TPZCompElSide Comp_wrapneigh = wrap_neigh.Reference();
            
            if(left_el_indexes.find(neigh_index)!=left_el_indexes.end()){
                
                TPZGeoElBC gbcleft(gelside,matInterfaceLeft);
                int64_t index;

                TPZMultiphysicsInterfaceElement *elem_left = new TPZMultiphysicsInterfaceElement(*cmesh,gbcleft.CreatedElement(),index,Comp_wrapneigh,celside);
                elem_left->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                
                
            }else if((right_el_indexes.find(neigh_index) != right_el_indexes.end())){
                
                TPZGeoElBC gbcright(gelside,matInterfaceRight);
                int64_t index;

                TPZMultiphysicsInterfaceElement *elem_right = new TPZMultiphysicsInterfaceElement(*cmesh,gbcright.CreatedElement(),index,Comp_wrapneigh,celside);
                elem_right->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                
            }else{
                DebugStop();
            }
        }
        
    }
    
    
    
}



void TPZFractureInsertion::AdjustSideOrient(TPZCompMesh *cmesh){
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    std::vector<int64_t> fracture_index = GetFractureIndexes();
    
    int fracture_id = GetFractureMaterialId();
    TPZAdmChunkVector<TPZGeoEl *> &elvec = gmesh->ElementVec();
    
    int meshdim = cmesh->Dimension();
    int64_t nelem = elvec.NElements();
    
    int64_t index;
    
    for (int iel = 0; iel < fracture_index.size(); iel++) {
        TPZGeoEl *gel = elvec[fracture_index[iel]];
        
        const int gelMatId = gel->MaterialId();
#ifdef PZDEBUG
        if (gelMatId != fracture_id)
            DebugStop();
        if (gel->HasSubElement())
            DebugStop();
#endif
        TPZCompEl *cel = gel->Reference();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        
        
        TPZStack<TPZCompElSide> neigh0, neigh1;
        int nsides = gel->NSides();
        
        TPZGeoElSide gelside0(gel,0);
        TPZGeoElSide gelside1(gel,1);
        TPZGeoElSide neighbour0 = gelside0.Neighbour();
        TPZGeoElSide neighbour1 = gelside1.Neighbour();
        
        gelside0.EqualLevelCompElementList(neigh0, 0, 0);
        gelside1.EqualLevelCompElementList(neigh1, 0, 0);
        TPZCompElSide celside0 = gelside0.Reference();
        TPZCompElSide celside1 = gelside1.Reference();
        
        TPZCompElSide neigh_elside;
        
        for (int i =0; i<neigh1.size(); i++) {
            if (neigh1[i].Element()->Material()->Id()==m_fracture_id) {
                neigh_elside=neigh1[i];
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(celside1.Element());
                intel->SetSideOrient(1,- 1);
            }
        
            bool isboundary = m_boundaries_material_ids.find(neigh1[i].Element()->Material()->Id()) != m_boundaries_material_ids.end();
            bool ispoint = neigh1[i].Element()->Dimension() == 0;
            if (isboundary && ispoint) {
                neigh_elside=neigh1[i];
            }
        }
        
        
    }
    
    VerifySideOrientation(cmesh);
    
}

void TPZFractureInsertion::VerifySideOrientation(TPZCompMesh *cmesh){
    
    TPZGeoMesh *gmesh = cmesh->Reference();

    int fracture_id = GetFractureMaterialId();
    TPZAdmChunkVector<TPZGeoEl *> &elvec = gmesh->ElementVec();
    
    int Wrap_Flux_ID = GetWrapFluxId();
    
    for (int iel = 0 ; iel< gmesh->NElements(); iel++) {
        TPZGeoEl *gel = elvec[iel];
        int dim_gel = gel->Dimension();
        if (dim_gel!=gmesh->Dimension()) {
            continue;
        }
        
        for (int i_side = 0 ; i_side < gel->NSides() ; i_side++) {
            
            TPZGeoElSide gelside(gel, i_side);
            
            if (gelside.Dimension()==dim_gel-1) {
                TPZStack<TPZCompElSide> neigh;
                gelside.EqualLevelCompElementList(neigh, 0, 0);
                for (int i_n =0; i_n<neigh.size(); i_n++) {

                    if (neigh[i_n].Element()->Material()->Id()!=Wrap_Flux_ID) {
                        continue;
                    }
                    
                    if (neigh[i_n].Element()->NConnects()!=1) {
                        DebugStop();
                    }
                    
                    TPZCompElSide celside =  gelside.Reference();
                    
                    TPZInterpolatedElement *intel_vol = dynamic_cast<TPZInterpolatedElement*>(celside.Element());
                    int SideOrient_vol = intel_vol->GetSideOrient(celside.Side()) ;
                    
                    TPZInterpolatedElement *intel_wrap = dynamic_cast<TPZInterpolatedElement*>(neigh[i_n].Element());
                    int SideOrient_wrap = intel_wrap->GetSideOrient(neigh[i_n].Side());
                    
                    if (SideOrient_vol != SideOrient_wrap) {
                        DebugStop();
                    }
                    
                    
                }

            
            }
        

        }
        
    }
    
    
    
}


REAL TPZFractureInsertion::VerifyLeftRightPoints(int LeftP, int RightP){
    
    int64_t n_gel = m_geometry->NElements();
    TPZGeoEl *gelPleft, *gelPright;
    TPZFMatrix<REAL > cooL(3,1,0.), cooR(3,1,0.);
    
    for (int64_t iel=0; iel < n_gel; iel++) {
        TPZGeoEl *gel = m_geometry->Element(iel);
        if (gel->Dimension()!=0) {
            continue;
        }
        if (gel->MaterialId()==LeftP) {
            gelPleft = gel;
            gelPleft->NodesCoordinates(cooL);
        }
        if(gel->MaterialId()==RightP){
            gelPright = gel;
            gelPright->NodesCoordinates(cooR);
        }
    }
    
    if (cooL(0,0)>cooR(0,0)) {
        DebugStop();
    }
    //Define frature orientation:
    if(cooL(1,0)>cooR(1,0)){
        return -1.;
    }else if (cooL(1,0)<cooR(1,0)){
        return 1.;
    }else{
        return 0.;
    }
    
    
    
}

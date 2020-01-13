//
//  TPZFractureInsertion.h
//  Benchmark0a
//  Class that stores fracture neighbor data in terms of geometric element indexes
//  Created by Pablo Carvalho on 02/08/18.
//

#ifndef TPZFractureInsertion_h
#define TPZFractureInsertion_h

#include <stdio.h>
#include <iostream>
#include <set>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"

class TPZFractureInsertion {
  
private:
    
    int m_fracture_id;
    
    TPZGeoMesh * m_geometry;
    
    std::set<int> m_boundaries_material_ids;
    
    std::vector<TPZGeoElSide> m_pivot_indexes;
    
    std::vector<TPZGeoElSide> m_non_pivot_indexes;
        
    std::vector<int64_t> m_fracture_indexes;
    
    std::set<int64_t> m_gel_left_indexes;
    
    std::set<int64_t> m_gel_right_indexes;
    
    /// look if there is a fracture element neighbour (ok)
    int64_t FractureNeighbourIndex(TPZGeoElSide gel_side);
    
    /// look if the element has a neighbour in left indices (ok)
    bool HasLeftIndexNeighbour(int64_t gel_index);
    
    /// look if the element has a neighbour in right indices (ok)
    bool HasRightIndexNeighbour(int64_t gel_index);
    
    /// Find all elements neighbour of a pivot (excluding lower dimension elements) (ok)
    std::set<int64_t> PivotNeighbours(TPZGeoElSide pivot_side);
    
    /// Verify if the non pivot has been inserted on non pivot structure
    bool HasInsertedNonPivot(TPZGeoElSide pivot_side);
    
    /// Verify if the pivot has been inserted on pivot structure
    bool HasInsertedPivot(TPZGeoElSide pivot_side);
    
    /// verify if the pivot has at least one element left and one element right as neighbour (ok)
    bool HasClassifiedNeighbour(TPZGeoElSide & pivot_side);
    
    /// Find the 2 elements that will initiate the element sets
    // m_gel_left_indexes and m_gel_right_indexes have to be empty
    // take a neighbouring fracture element - insert its two neighbours along the face (ok)
    void InsertFractureNeighbours(std::set<int64_t> pivot_neighbours);
    
    /// Classify elements into left and right elements (ok)
    void ClassifyElements(std::set<int64_t> pivot_neighbours);
    
    /// look for all fracture elements (ok)
    void BuildFractureElements();
    
    /// build the pivot and non pivot data structure (ok)
    void BuildPivotDataStructure();
    
    /// Return true if neighboour is at left position
    bool neigh_is_left;
    
    /// Wrap flux id
    int mat_id_flux_wrap;
    
public:

    /// @TODO:: OD, Rename TPZFractureInsertion -> TPZFractureDescription
    /// @TODO:: OD, Refactor and rename the methods dependent on the approximation space
    
    /// Default constructor
    TPZFractureInsertion();
    
    /// Default desconstructor
    ~TPZFractureInsertion();
    
    /// Copy constructor
    TPZFractureInsertion(TPZFractureInsertion & other);
    
    /// Constructor based on a computational mesh and fracture material id
    TPZFractureInsertion(TPZGeoMesh * geometry, int fracture_id, std::set<int> & boundaries_ids);

    /// Set fracture Identifier
    void SetFractureIdentifier(int fracture_id);
    
    /// Get fracture material Identifier
    int & GetFractureMaterialId();
    
    /// Set wrap Identifier
    void SetWrapFluxIdentifier(int wrapFlux);
    
    /// Get wrap Identifier
    int & GetWrapFluxId();

    /// Get node pivots
    std::vector<TPZGeoElSide> & GetPivotIndexes();

    /// Get node non pivots
    std::vector<TPZGeoElSide> & GetNonPivotIndexes();

    /// Get geometric fracture indexes
    std::vector<int64_t> & GetFractureIndexes();
    
    /// Get geometric element indexes for left
    std::set<int64_t> & GetLeftIndexes();
    
    /// Get geometric element indexes for right
    std::set<int64_t> & GetRightIndexes();
    
    /// Classify the neighbouring elements of the pivots
    void ClassifyNeighboursofPivots();
    
    /// Open the connects of a fracture, create dim-1 fracture elements (H1 version)
    void OpenFractureOnH1(TPZCompMesh *cmesh);
    
    /// Open the connects of a fracture, create dim-1 fracture elements (Hdiv version)
    void OpenFractureOnHdiv(TPZCompMesh *cmesh, int mat_id_flux_wrap);
    
    /// Set Discontinuous elements on fractures
    void SetDiscontinuosFrac(TPZCompMesh *cmesh);
    
    /// Set interfaces elements between fracture and volumetric elements
    void SetInterfaces(TPZCompMesh *cmesh, int matInterfaceLeft, int matInterfaceRight);
    
    /// Set multiphysics interfaces elements between fracture and volumetric elements
    void SetMultiphysicsInterfaces(TPZCompMesh *cmesh, int matInterfaceLeft, int matInterfaceRight, int mat_id_flux_wrap);
    
    void AdjustSideOrient(TPZCompMesh *cmesh);
    
    void VerifySideOrientation(TPZCompMesh *cmesh);
    
    REAL VerifyLeftRightPoints(int LeftP, int RightP);
    
};


#endif /* TPZFractureInsertion_h */

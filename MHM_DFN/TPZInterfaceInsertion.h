//
//  TPZInterfaceInsertion.h
//  Benchmark0a
//  Class that stores Interface neighbor data in terms of geometric element indexes
//  Created by Pablo Carvalho on 02/08/18.
//

#ifndef TPZInterfaceInsertion_h
#define TPZInterfaceInsertion_h

#include <stdio.h>
#include <iostream>
#include <set>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMultiphysicsCompMesh.h"

class TPZInterfaceInsertion {
    
private:

    /// Mesh geometry
    TPZCompMesh *m_cmesh;
    
    TPZGeoMesh * m_geometry;
        
    /// Wrap flux id
    int m_id_flux_wrap;
    
    /// Interfaces ids

    int m_interface_id;
    
    TPZManVector<int64_t,3> m_interfaceVector_ids;
    
    /// Multiplier material id (interior)
    int m_multiplier_id;
    
    /// Multiplier material id (BC)
    int m_multiplierBC_id;

    /// Set of boundary material ids
    std::set<int> m_boundaries_ids;
    
    /// Set mesh type
    MElementType m_Eltype;

public:
    
    /// @TODO:: OD, Rename TPZInterfaceInsertion -> TPZInterfaceDescription
    /// @TODO:: OD, Refactor and rename the methods dependent on the approximation space
    
    /// Default constructor
    TPZInterfaceInsertion();
    
    /// Default desconstructor
    ~TPZInterfaceInsertion();
    
    /// Copy constructor
    TPZInterfaceInsertion(TPZInterfaceInsertion & other);
    
    /// Constructor based on a computational mesh and Interface material id
    TPZInterfaceInsertion(TPZCompMesh *m_cmesh, int Interface_id, std::set<int> & boundaries_ids, MElementType eltype);
    
    /// Set Interface Identifier
    void SetInterfaceId(int Interface_id);
    
    /// Get Interface material Identifier
    int & GetInterfaceId();
    
    /// Set Interface Vector Ids
    void SetInterfaceVectorId(TPZManVector<int64_t,3> interfaceVector_ids);

    /// Get Interface Vector Ids
    TPZManVector<int64_t,3> & GetInterfaceVectorId();
    
    /// Set wrap Identifier
    void SetWrapFluxIdentifier(int wrapFlux);
    
    /// Get wrap Identifier
    int & GetWrapFluxId();

    /// Set multiplier material id
    void SetMultiplierMatId(int multiplier);
    
    /// Get multiplier material id
    int & GetMultiplierMatId();
    
    /// Set multiplier material id - BC
    void SetMultiplierBCMatId(int multiplier_BC);
    
    /// Get multiplier material id - BC
    int & GetMultiplierBCMatId();
    
    /// Add multiphysics interfaces for all boundaries and internal elements
    void AddMultiphysicsInterfaces();
    
    /// Add interface from a reference material, only one side
    void AddMultiphysicsInterfaces(int matfrom, int mattarget);
    
    /// Add BC interface from a reference material, only one side
    void AddMultiphysicsBCInterface(int matfrom, int mattarget);
    
    /// Add BC interface from a reference material, only one side
    void AddMultiphysicsBCInterface2(int matfrom, int mattarget);
    
    /// Add interfaces in left and right sides of a reference material, it should be given an interface vector (left and right materials)
    void AddMultiphysicsInterfacesLeftNRight(int matfrom);

    /// Add interfaces in left and right sides of a reference material, it should be given an interface vector (left and right materials)
    /// This method consider the lowest fathers mesh skeleton
    void AddMultiphysicsInterfacesLeftNRight2(int matfrom);
    
    /// Open the connects of a Interface, create dim-1 Interface elements (Hdiv version)
    void InsertHdivBound(int mat_id_flux_wrap);
    
};


#endif /* TPZInterfaceInsertion_h */


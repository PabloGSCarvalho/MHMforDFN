//
//  TPZMHMStokesMeshControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#ifndef TPZMHMStokesMeshControl_hpp
#define TPZMHMStokesMeshControl_hpp

#include <stdio.h>

#include "TPZMHMixedMeshControl.h"

/// class for creating TPZMHMM with Stokes problem Meshes
class TPZMHMStokesMeshControl : public TPZMHMixedMeshControl
{
    
protected:
    
    /// Computational mesh to contain the avarege pressure elements
    TPZAutoPointer<TPZCompMesh> fAveragePressMesh;

    /// Computational mesh to contain the distributed flux elements
    TPZAutoPointer<TPZCompMesh> fDistrFluxMesh;

    /// Computational mesh to contain the coarse avarege pressure elements
    TPZAutoPointer<TPZCompMesh> fCoarseAveragePressMesh;
    
    /// Computational mesh to contain the coarse distributed flux elements
    TPZAutoPointer<TPZCompMesh> fCoarseDistrFluxMesh;
    
    /// Set coarse average pressure and distributed flux
    bool fsetCoarseAverageMultipliers;


public:
    
    TPZMHMStokesMeshControl() : TPZMHMixedMeshControl()
    {
        
    }
    
    TPZMHMStokesMeshControl(int dimension);
    
    TPZMHMStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices);
    
    TPZMHMStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh);
    
    TPZMHMStokesMeshControl(const TPZMHMStokesMeshControl &copy);
    
    TPZMHMStokesMeshControl &operator=(const TPZMHMStokesMeshControl &cp);
    
    virtual ~TPZMHMStokesMeshControl();
    
    /// Create all data structures for the computational mesh
    virtual void BuildComputationalMesh(bool usersubstructure);
    
    void InsertInternalSkeleton();
    
    void InsertBCSkeleton();
    
    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        meshvec.Resize(4);
        meshvec[0] = fFluxMesh.operator->();
        meshvec[1] = fPressureFineMesh.operator->();
        meshvec[2] = fAveragePressMesh.operator->();
        meshvec[3] = fDistrFluxMesh.operator->();
    }
    
    TPZVec<TPZAutoPointer<TPZCompMesh> > GetMeshes()
    {
        TPZManVector<TPZAutoPointer<TPZCompMesh>,4> result(4);
        result[0] = fFluxMesh;
        result[1] = fPressureFineMesh;
        result[2] = fAveragePressMesh;
        result[3] = fDistrFluxMesh;
        
        return result;
    }
    
    /// Set the flag for creating Lagrange Dofs for the average pressure and distributed flux on coarse mesh
    void SetCoarseAverageMultipliers(bool flag)
    {
        fsetCoarseAverageMultipliers = flag;
    }
    
protected:

    
    /// Insert the necessary pressure material objects to create the pressure mesh
    virtual void InsertPeriferalPressureMaterialObjects();
    
    /// Create the pressure mesh which contain traction variable
    void CreatePressureAndTractionMHMMesh();

    /// Insert the necessary average pressure material objects to create the average pressure mesh
    void InsertPeriferalAveragePressMaterialObjects();

    /// Create the average pressure mesh
    void CreateAveragePressMHMMesh();

    /// Insert the necessary coarse average pressure material objects to create the average pressure mesh
    void InsertPeriferalCoarseAveragePressMaterialObjects();
    
    /// Create the coarse average pressure mesh
    void CreateCoarseAveragePressMHMMesh();
    
    /// Insert the necessary distributed flux material objects to create the distributed flux material pressure mesh
    void InsertDistributedFluxMaterialObjects();
    
    /// Create the distributed flux mesh
    void CreateDistributedFluxMHMMesh();
    
    /// Insert the necessary distributed flux material objects to create the distributed flux material pressure mesh
    void InsertCoarseDistributedFluxMaterialObjects();
    
    /// Create the distributed flux mesh
    void CreateCoarseDistributedFluxMHMMesh();
    
    /// Create multiphysics MHM mesh
    void CreateMultiPhysicsMHMMesh();
    
//    /// build the multi physics mesh (not at the finest geometric mesh level)
//    virtual void BuildMultiPhysicsMesh();

    /// Create the multiphysics interface elements between elements of traction material id
    void CreateMultiPhysicsInterfaceElements();

    /// Create the multiphysics BC interface elements between elements of BC traction material id
    void CreateMultiPhysicsBCInterfaceElements();

    void GroupandCondenseSubMeshes();
    
    void GroupAndCondense(TPZCompMesh *cmesh_m);

    void BuildSubMeshes();
    
    
public:
    
    /// material id associated with the internal skeleton elements
    int64_t fTractionMatId = 3;

    /// material ids associated with the BC skeleton elements
    std::map<int64_t,int64_t> fBCTractionMatIds;
    
};

#endif /* TPZMHMixedMeshControl_hpp */


//
//  TPZMHMDarcyDFNMeshControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#ifndef TPZMHMDarcyDFNMeshControl_hpp
#define TPZMHMDarcyDFNMeshControl_hpp

#include <stdio.h>

#include "TPZMHMeshControl.h"

/// class for creating TPZMHMM with DarcyDFN problem Meshes
class TPZMHMDarcyDFNMeshControl : public TPZMHMeshControl
{
    
protected:
    
    /// Computational mesh to contain the avarege pressure elements
    TPZAutoPointer<TPZCompMesh> fAveragePressMesh;
    
    /// Computational mesh to contain the distributed flux elements
    TPZAutoPointer<TPZCompMesh> fDistrFluxMesh;
    
    
public:
    
    TPZMHMDarcyDFNMeshControl() : TPZMHMeshControl()
    {
        
    }
    
    TPZMHMDarcyDFNMeshControl(int dimension);
    
    TPZMHMDarcyDFNMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices);
    
    TPZMHMDarcyDFNMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh);
    
    TPZMHMDarcyDFNMeshControl(const TPZMHMDarcyDFNMeshControl &copy);
    
    TPZMHMDarcyDFNMeshControl &operator=(const TPZMHMDarcyDFNMeshControl &cp);
    
    virtual ~TPZMHMDarcyDFNMeshControl();
    
    /// Create all data structures for the computational mesh
    virtual void BuildComputationalMesh(bool usersubstructure);
    
    void InsertH1Fracture(TPZCompMesh &cmesh);

    virtual void CreateInterfaceElements();
    
    virtual void CreateSkeleton();
    
    void InsertInternalMultipliers();
    
    void InsertBCMultipliers();
    
    /// Insert material objects that do not perform any actual computation
    void InsertPeriferalMaterialObjects();
    
    void CreateInternalElements();
    
    void InsertPeriferalLagrangeFluxObjects();
    
    void CreateLagrangeFluxElements();

    /// Insert the necessary average pressure material objects to create the average pressure mesh
    void InsertPeriferalAveragePressMaterialObjects();
    
    /// Create the average pressure mesh
    void CreateAveragePressMHMMesh();
    
    /// Insert the necessary distributed flux material objects to create the distributed flux material pressure mesh
    void InsertDistributedFluxMaterialObjects();
    
    /// Create the distributed flux mesh
    void CreateDistributedFluxMHMMesh();
    
    void CreateWrapElement(const TPZCompElSide &left, int multId);

    void CreateFractureBC();
    
    void CreateH1Wrappers();
    
    void CreateHDivWrappers();
    
    void ResetAllNeigReferences(TPZGeoElSide gelside);
    
    bool NeedsHDivWrapper(TPZGeoElSide gelside);
    
    void Print(std::ostream & out = std::cout) const;
    
    void SetFracSubDomain(int fracMatId);
    
    void CreateSkeletonElements();
    
    void CreateFractureInterfaces();
    
    void CreateLagrangeMultiplierMesh(int dim);
    
    void TransferToMultiphysics();
    
    bool HasFracNeighbour(TPZGeoElSide &gelside);
    
public:
    
    /// material id associated with the internal skeleton elements
    int64_t fLagrangeFluxMatIdL = 11;
    
    /// material id associated with the internal skeleton elements
    int64_t fLagrangeFluxMatIdR = 12;

    /// wrap material id for H1 elements
    int64_t fWrapH1MatIdL = 13;

    /// wrap material id for H1 elements
    int64_t fWrapH1MatIdR = 14;

    /// material id associated with the fractue elements
    int64_t fFractureMatId = 10;

    /// material id associated with the fractue boundaries
    int64_t fmatFracPointL = -7;
    int64_t fmatFracPointR = -8;
    
    /// material id of the flow elements of dimension fDim-1
    std::set<int> fFractureFlowDim1MatId;
    
    /// material ids associated with the BC skeleton elements
    std::map<int64_t,int64_t> fBCLagrangeFluxMatIds;
    
    /// material id associated with the lagrange multiplier elements
    int64_t  fLagrangeMatIdFrac = 552;
    
    // geoindex of fractures and next pair of bounds
    std::set<std::map<int64_t, std::pair<int64_t,int64_t> > > fFracInterfaces;
    
};

#endif /* TPZMHMixedMeshControl_hpp */



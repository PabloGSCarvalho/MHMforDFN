//
//  TPZSegregatedAnalysis.h
//
//  Created by Omar Durán on 9/13/18.
//

#ifndef TPZNSAnalysis_h
#define TPZNSAnalysis_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "pzanalysis.h"
#include "TPZSimulationData.h"
#include "pzpostprocanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMHMNavierStokesMaterial.h"


class TPZNSAnalysis : public TPZAnalysis {
    
private:
    
    /// Pointer to Simulation data object
    TPZSimulationData * m_simulation_data;
    
    /// Solution at n+1 state
    TPZFMatrix<STATE> m_U_Plus;
    
    /// Solution at n (past) state
    TPZFMatrix<STATE> m_U_n;
    
    /// Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> m_mesh_vec;
    
    /// Residue error
    STATE m_res_error;
    
    
    /// number of Newton iterations
    int m_k_iterations;
    
    /// Post-processor object
    TPZPostProcAnalysis * m_post_processor;
    
    /// Variables being postprocessed as scalar
    TPZStack<std::string> m_var_names;
    
    /// Variables being postprocessed as vector on standard post-process
    TPZStack<std::string> m_vec_var_names;
    

    /** @brief Part of residue at n+1 */
    TPZFMatrix<REAL> m_R_Plus;
    
    /** @brief Part of residue at n  */
    TPZFMatrix<REAL> m_R_n;

    /// Correction variation
    STATE m_dU_norm;

    
    
public:
    
    /// Default constructor
    TPZNSAnalysis();
    
    /// Destructor
    ~TPZNSAnalysis();
    
    /// Copy constructor
    TPZNSAnalysis(const TPZNSAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }

    /// Configurate internal analysis objects and linked them through the memory shared pointer
    void ConfigurateAnalysis(DecomposeType decompose_geo, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_M, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZVec<std::string> & var_names);
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep();
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string & res_file);
    
    /// Execute the transient evolution using Fixed Stress Split Iteration
    void ExecuteTimeEvolution();
    
    // Set initial parameters
    void SetInitialParameters();
    
    // Update parameters
    void UpdateParameters();
    
    // Adjust integration orders
    void AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d);
    
    /// Update the memory with the converged time step solution
    void AcceptTimeStepSolution();
    
    /// Execute a single newton iteration
    void ExecuteNewtonInteration();
    
    /// Load the current state for the hdiv and 2 meshes
    void LoadCurrentState();
    
    /// Load the last state for the hdiv and 2 meshes
    void LoadLastState();
    
    /// Update solution state x = x_n
    void UpdateState(){
        m_U_Plus = m_U_n;
    }
    
    /** @brief Set Residue error */
    void Set_error(STATE error)
    {
        m_res_error = error;
    }
    
    /** @brief Get Residue error */
    STATE Get_error()
    {
        return m_res_error;
    }
    
    /** @brief Set Correction variation */
    void Set_dU_norm(STATE dUnorm)
    {
        m_dU_norm = dUnorm;
    }
    
    /** @brief Get Correction variation */
    STATE Get_dU_norm()
    {
        return m_dU_norm;
    }
    
    /** @brief Set number of Newton iterations */
    void Set_k_iterations(int kiterations)
    {
        m_k_iterations = kiterations;
    }
    
    /** @brief Get number of Newton iterations */
    int Get_k_iterations()
    {
        return m_k_iterations;
    }
    
    /** @brief Set variables being postprocessed as vector on standard post-process */
    void Set_vec_var_names(TPZStack<std::string> & vec_var_names)
    {
        m_vec_var_names = vec_var_names;
    }
    
    /** @brief Get variables being postprocessed as vector on standard post-process */
    TPZStack<std::string> & Get_vec_var_names()
    {
        return m_vec_var_names;
    }

};

#endif /* TPZNSAnalysis.h */

//
//  TPZNSAnalysis.cpp
//  PMRS
//
//  Created by Omar Durán on 9/13/18.
//

#include "TPZNSAnalysis.h"
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("NavierStokes.Analysis"));
#endif

TPZNSAnalysis::TPZNSAnalysis() : TPZAnalysis() {
    m_simulation_data = NULL;
    m_U_Plus.Resize(0, 0);
    m_U_n.Resize(0, 0);
    m_mesh_vec.Resize(0);
    m_res_error = 0;
    m_dU_norm = 0;
    m_k_iterations = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    m_vec_var_names.resize(0);
    m_R_Plus.Resize(0,0);
    m_R_n.Resize(0,0);
}

TPZNSAnalysis::~TPZNSAnalysis(){
    
}

TPZNSAnalysis::TPZNSAnalysis(const TPZNSAnalysis & other){
    m_simulation_data   = other.m_simulation_data;
    m_U_Plus            = other.m_U_Plus;
    m_U_n               = other.m_U_n;
    m_mesh_vec          = other.m_mesh_vec;
    m_res_error         = other.m_res_error;
    m_dU_norm           = other.m_dU_norm;
    m_k_iterations      = other.m_k_iterations;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    m_vec_var_names     = other.m_vec_var_names;
    m_R_Plus            = other.m_R_Plus;
    m_R_n               = other.m_R_n;
}


void TPZNSAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_M, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZVec<std::string> & var_names){
    
    SetSimulationData(simulation_data);
    bool mustOptimizeBandwidth = simulation_data->GetOptimizeBandwidthQ();
    this->SetCompMesh(cmesh_M,mustOptimizeBandwidth);
    TPZStepSolver<STATE> step;
    unsigned int n_threads = m_simulation_data->GetNthreads();
    
    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    
    m_mesh_vec = mesh_vec;
    switch (decomposition) {
        case ELU:
        {
            //#ifdef USING_MKL
            //            TPZSpStructMatrix struct_mat(Mesh());
            //            struct_mat.SetNumThreads(n_threads);
            //            this->SetStructuralMatrix(struct_mat);
            //#else
            TPZSkylineNSymStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(n_threads);
            this->SetStructuralMatrix(struct_mat);
            //#endif
        }
            break;
        case ELDLt:
        {
            TPZSymetricSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(n_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    step.SetDirect(decomposition);
    this->SetSolver(step);
    this->Solution().Resize(Mesh()->Solution().Rows(), 1);
    m_U_n.Resize(Mesh()->Solution().Rows(), 1);
    m_U_Plus.Resize(Mesh()->Solution().Rows(), 1);
    
    m_post_processor = new TPZPostProcAnalysis;
    m_post_processor->SetCompMesh(Mesh());
    
    int n_vols = m_simulation_data->Get_volumetric_material_id().size();
    TPZManVector<int,10> post_mat_id(n_vols);
    for (int ivol = 0; ivol < n_vols; ivol++)
    {
        int matid = m_simulation_data->Get_volumetric_material_id()[ivol];
        post_mat_id[ivol] = matid;
    }
    
    for (auto i : var_names) {
        m_var_names.Push(i);
    }
    
    m_post_processor->SetPostProcessVariables(post_mat_id, m_var_names);
    int dim = Mesh()->Dimension();
    int div = 0;
    TPZStack< std::string> vecnames;
    std::string plotfile("NavierStokesPostProcess.vtk");
    
    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,plotfile);
    
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);
    
}


void TPZNSAnalysis::ExecuteOneTimeStep(){

    // m_X means the solution at the previous time step
    if (m_simulation_data->IsInitialStateQ()) {
        m_U_n = Solution();
    }
    // Set current state false means overwriting p of the memory
    //m_simulation_data->SetCurrentStateQ(false);
    // Accect time solution means writing one of the vectors of this object in the memory
    //AcceptTimeStepSolution();

    //    // Initial guess
    //    m_X_n = m_X;
    // Set current state true means overwriting p_n of the memory object
    //m_simulation_data->SetCurrentStateQ(true);

    // Accept time solution here means writing one of the vectors of the object into the memory
    //AcceptTimeStepSolution();

    
    std::ofstream plotNavierEK("NavierStiffness.txt");
    std::ofstream plotNavierEF("NavierRhs.txt");
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fRhs.Print("Rhs =",sout);
        //EFormatted, EInputFormat, EMathematicaInput, EMatlabNonZeros, EMatrixMarket
        //  fSolver->Matrix()->Print("ek = ",plotDarcyEK,EMathematicaInput);
        fRhs.Print("ef = ",plotNavierEF,EMathematicaInput);

        PrintVectorByElement(sout, fRhs);
        PrintVectorByElement(sout, fSolution);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZFMatrix<STATE> dU;
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dU;
    REAL res_norm = m_simulation_data->Get_epsilon_res();
    REAL dU_norm = m_simulation_data->Get_epsilon_cor();
    int n_it = m_simulation_data->Get_n_iterations();
    
    for (int i = 1; i <= n_it; i++) {
        this->ExecuteNewtonInteration();
        
        dU = Solution();
        
//        std::cout<<dU<<std::endl;
        
        norm_dU  = Norm(dU);
        m_U_Plus += dU;
        
//        std::cout<<m_U_Plus<<std::endl;
        
        //LoadCurrentState();
        LoadLastState();
        AssembleResidual();
        m_R_n = this->Rhs();
//        std::cout<<this->Rhs()<<std::endl;
        LoadCurrentState();
        AssembleResidual();
        m_R_Plus = this->Rhs();
//        std::cout<<this->Rhs()<<std::endl;
//        REAL test_RHS_norm_n = Norm(m_R_n);
//        REAL test_RHS_norm = Norm(m_R_Plus);
        
        //m_R_Plus += m_R_n; // total residue
        m_res_error =  Norm(m_R_Plus); // residue error
        
        
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            fRhs.Print("Rhs =",sout);
            {
                std::ofstream plotDarcyEK("DarcyStiffness.txt");
                std::ofstream plotDarcyEF("DarcyRhs.txt");
                fSolver->Matrix()->Print("ek = ",plotDarcyEK,EMathematicaInput);
                fRhs.Print("ef = ",plotDarcyEF,EMathematicaInput);
            }
            PrintVectorByElement(sout, fRhs);
            PrintVectorByElement(sout, fSolution);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        
        norm_res = Norm(Rhs());
        residual_stop_criterion_Q   = norm_res < res_norm;
        correction_stop_criterion_Q = norm_dU  < dU_norm;
        
        m_k_iterations = i;
        m_res_error = norm_res;
        m_dU_norm = norm_dU;
        
        
        if (residual_stop_criterion_Q ||  correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPZNSAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPZNSAnalysis:: Number of iterations = " << i << std::endl;
            std::cout << "TPZNSAnalysis:: Correction norm = " << norm_dU << std::endl;
#endif
            m_simulation_data->SetCurrentStateQ(true);
            this->AcceptTimeStepSolution();
            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        std::cout << "TPZNSAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
    

}

void TPZNSAnalysis::PostProcessTimeStep(std::string & res_file){

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = 2;
    TPZStack<std::string> scalnames, vecnames;
//    if (fSimulationData->IsInitialStateQ()) {
//        plotfile =  "DualSegregatedDarcyOnBox_I.vtk";
//        return;
//    }
//    else{
//        plotfile =  "DualSegregatedDarcyOnBox.vtk";
//    }
    
    //Pós-processamento (paraview):

    scalnames.Push("P");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    scalnames.Push("Div");
    
    this->DefineGraphMesh(dim, scalnames, vecnames, res_file);
    this->PostProcess(div,dim);
    

}

void TPZNSAnalysis::ExecuteTimeEvolution(){

    std::string file_NavierStokes("NavierStokes.vtk");
    //Testes
    std::string file_NavierStokes_test("NavierStokes_test.vtk");

    int n_max_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    int n_enforced_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    int n_time_steps = 1;
    REAL res_norm = m_simulation_data->Get_epsilon_res();
    REAL dU_norm = m_simulation_data->Get_epsilon_cor();
    bool error_stop_criterion_Q = false;
    bool dU_stop_criterion_Q = false;
    this->SetInitialParameters();
    
    {
   //     std::ofstream filecE("MalhaC_E_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
   //     m_elastoplast_analysis->Mesh()->Print(filecE);
    }
    {
        std::ofstream filecM("MalhaC_M_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
        Mesh()->Print(filecM);
    }
    
    
    for (int it = 0; it < n_time_steps; it++) { //??
        for (int k = 1; k <= n_max_fss_iterations; k++) {
            this->ExecuteOneTimeStep();
            
            error_stop_criterion_Q = this->Get_error() < res_norm;
            dU_stop_criterion_Q = this->Get_dU_norm() < dU_norm;
            
            this->PostProcessTimeStep(file_NavierStokes_test);
            
            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) && dU_stop_criterion_Q) {
                this->PostProcessTimeStep(file_NavierStokes);
                std::cout << "TPZNSAnalysis:: Iterative process converged with residue norm  = " << this->Get_error() << std::endl;
                UpdateState();
                break;
            }
        }
    }

}

void TPZNSAnalysis::SetInitialParameters(){
    
    // Updating volumetric parameters :
    
//    int matid = m_simulation_data->Get_elasticity_matid();
//
//    TPZMaterial * material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(matid);
//    TPZMaterial * material_darcy = m_darcy_analysis->Mesh()->FindMaterial(matid);
//
//    if (!material_elastoplast || !material_darcy) {
//        DebugStop();
//    }
//
//    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_darcy);
//
//    long N_ipoints = mat_with_memory_darcy->GetMemory().get()->NElements();
//    for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
//        TPZMemoryDFN & memory = mat_with_memory_darcy->GetMemory().get()->operator[](ip_index);
//        TPZFNMatrix<9,REAL>  k_0 = m_simulation_data->Get_PermeabilityTensor_0();
//        REAL phi0 = m_simulation_data->Get_Porosity_0();
//        memory.Setkappa_0(k_0);
//        memory.Setkappa(k_0);
//        memory.Setphi_0(phi0);
//    }
//
//    // Updating Fracture parameters :
//
//    for (int imat_frac = 0; imat_frac < m_simulation_data->Get_fracture_material_id().size(); imat_frac++) {
//        int frac_matid = m_simulation_data->Get_fracture_material_id()[imat_frac];
//
//        TPZMaterial * frac_material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(frac_matid);
//        TPZMaterial * frac_material_darcy = m_darcy_analysis->Mesh()->FindMaterial(frac_matid);
//
//
//
//        if (!frac_material_elastoplast || !frac_material_darcy) {
//            DebugStop();
//        }
//
//        TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_darcy);
//
//        long N_ipoints = matfrac_with_memory_darcy->GetMemory().get()->NElements();
//
//        for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
//            TPZMemoryFracDFN & memory_frac = matfrac_with_memory_darcy->GetMemory().get()->operator[](ip_index);
//            TPZFNMatrix<9,REAL>  k_0(3,3,0.);
//            k_0(0,0) = m_simulation_data->Get_Permeability_0().find(frac_matid)->second;
//            k_0(1,1) = m_simulation_data->Get_Permeability_0().find(frac_matid)->second;
//            k_0(2,2) = m_simulation_data->Get_Permeability_0().find(frac_matid)->second;
//
//            REAL Vm = m_simulation_data->Get_Vm().find(frac_matid)->second; //Max fracture closure
//            REAL a0 = m_simulation_data->Get_a0().find(frac_matid)->second;
//            REAL Kni = m_simulation_data->Get_Kni();
//
//            REAL Du_0 = Vm-a0; // initial closure
//
//           // Du_0 = (-17.2*Vm)/(-17.2-Kni*Vm);
//
//            memory_frac.SetVm(Vm);
//           // memory_frac.SetDu_0(Du_0);
//            memory_frac.Setkappa_0(k_0);
//            memory_frac.Setkappa(k_0);
//        }
//
//
//    }
    
}



void TPZNSAnalysis::UpdateParameters(){
    
    // Updating volumetric parameters :
    
//    int matid = m_simulation_data->Get_elasticity_matid();
//
//    TPZMaterial * material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(matid);
//    TPZMaterial * material_darcy = m_darcy_analysis->Mesh()->FindMaterial(matid);
//
//    if (!material_elastoplast || !material_darcy) {
//        DebugStop();
//    }
//
//    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_elastoplast = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_elastoplast);
//    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_darcy);
//
//    long N_ipoints = mat_with_memory_darcy->GetMemory().get()->NElements();
//
//    for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
//        TPZMemoryDFN & memory = mat_with_memory_darcy->GetMemory().get()->operator[](ip_index);
//        TPZFNMatrix<9,REAL>  k_0 = memory.kappa_0();
//        REAL phi_0 = memory.phi_0();
//        REAL nu = m_simulation_data->Get_Poisson();
//        REAL E = m_simulation_data->Get_Eyoung();
//        TPZFNMatrix<9,REAL> k_n(3,3,0.);
//        for(int i = 0; i < k_0.Rows(); i++){
//            k_n(i,i) = memory.Permeability(k_0(i,i),phi_0,nu,E,m_simulation_data->Get_Stress_Vol0()[ip_index]);
//        }
//        memory.Setkappa(k_n);
//    }
//
//    // Updating Fracture parameters :
//
//    for (int imat_frac = 0; imat_frac < m_simulation_data->Get_fracture_material_id().size(); imat_frac++) {
//        int frac_matid = m_simulation_data->Get_fracture_material_id()[imat_frac];
//
//        TPZMaterial * frac_material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(frac_matid);
//        TPZMaterial * frac_material_darcy = m_darcy_analysis->Mesh()->FindMaterial(frac_matid);
//
//        if (!frac_material_elastoplast || !frac_material_darcy) {
//            DebugStop();
//        }
//
//        TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_elastoplast = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_elastoplast);
//        TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_darcy);
//
//        long N_ipoints = matfrac_with_memory_darcy->GetMemory().get()->NElements();
//
//        for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
//            TPZMemoryFracDFN & memory_frac = matfrac_with_memory_darcy->GetMemory().get()->operator[](ip_index);
//            TPZFNMatrix<9,REAL>  k_0 = memory_frac.kappa_0();
//
//            TPZFNMatrix<9,REAL> k_n(k_0.Rows(),k_0.Cols(),0.);
//            for(int i = 0; i < k_0.Rows(); i++){
//                k_n(i,i) = memory_frac.Permeability(k_0(i,i));
//            }
//            memory_frac.Setkappa(k_n);
//        }
//
//
//    }
    
}




void TPZNSAnalysis::AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    
    // Assuming the cmesh_o as directive.
    
//    cmesh_d->LoadReferences();
//    int nel_o = cmesh_o->NElements();
//    int nel_d = cmesh_d->NElements();
//
//    int nvol_el_o = 0;
//    int nvol_el_d = 0;
//
//    for (long el = 0; el < nel_o; el++) {
//        TPZCompEl *cel_o = cmesh_o->Element(el);
//        int matid_o = cel_o->Material()->Id();
//        for (int imat = 0; imat < m_simulation_data->Get_volumetric_material_id().size(); imat++) {
//            int matid = m_simulation_data->Get_volumetric_material_id()[imat];
//            if (matid==matid_o) {
//                nvol_el_o ++;
//                break;
//            }
//        }
//    }
//
//    for (long el = 0; el < nel_d; el++) {
//        TPZCompEl *cel_d = cmesh_d->Element(el);
//        int matid_d = cel_d->Material()->Id();
//        for (int imat = 0; imat < m_simulation_data->Get_volumetric_material_id().size(); imat++) {
//            int matid = m_simulation_data->Get_volumetric_material_id()[imat];
//            if (matid==matid_d) {
//                nvol_el_d ++;
//                break;
//            }
//        }
//    }
//
//    if (nvol_el_o != nvol_el_d) {
//        std::cout << "The geometrical partitions are not the same." << std::endl;
//        DebugStop();
//    }
//
//    for (long el = 0; el < nel_o; el++) {
//        TPZCompEl *cel_o = cmesh_o->Element(el);
//        if (!cel_o) {
//            continue;
//        }
//
//        TPZGeoEl * gel = cel_o->Reference();
//        if (!gel) {
//            continue;
//        }
//
//        // Finding the other computational element
//        TPZCompEl * cel_d = gel->Reference();
//        if (!cel_d) {
//            continue;
//        }
//        cel_o->SetFreeIntPtIndices();
//        cel_o->ForcePrepareIntPtIndices();
//        const TPZIntPoints & rule = cel_o->GetIntegrationRule();
//        TPZIntPoints * cloned_rule = rule.Clone();
//
//        TPZManVector<int64_t,20> indices;
//        cel_o->GetMemoryIndices(indices);
//        cel_d->SetFreeIntPtIndices();
//        cel_d->SetMemoryIndices(indices);
//        cel_d->SetIntegrationRule(cloned_rule);
//    }
//
//#ifdef PZDEBUG
//    std::ofstream out_geo("Cmesh_origin_adjusted.txt");
//    cmesh_o->Print(out_geo);
//#endif
//
//
//#ifdef PZDEBUG
//    std::ofstream out_res("Cmesh_destination_adjusted.txt");
//    cmesh_d->Print(out_res);
//#endif
}

void TPZNSAnalysis::ExecuteNewtonInteration(){
    this->Assemble();
    this->Rhs() *= 1.0;
    this->Solve();
}

void TPZNSAnalysis::LoadCurrentState(){
     LoadSolution(m_U_Plus);
     TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());

}

void TPZNSAnalysis::LoadLastState(){
    LoadSolution(m_U_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
}

void TPZNSAnalysis::AcceptTimeStepSolution(){
    
    bool state = m_simulation_data->IsCurrentStateQ();
    if (state) {
        // must accept solution changes a global data structure shared by the material objects
        // which indicates the solution should be overwritten in memory
        m_simulation_data->Set_must_accept_solution_Q(true);
        // load current state copies m_X_n into the solution vector
        LoadCurrentState();
        // puts the solution vector into a variable depending on yet another global variable
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }else{
        // m_simulation_data is pointer shared by the material object
        // this call forces the solution to be loaded into the memory object
        m_simulation_data->Set_must_accept_solution_Q(true);
        // put m_X in the mesh solution
        LoadLastState();
        // copy the state vector into the memory because must_accept_solution_Q in the m_simulation_data is true
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }
}



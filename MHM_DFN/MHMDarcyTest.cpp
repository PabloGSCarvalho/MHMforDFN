/*
 *  MHMDarcyTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "MHMDarcyTest.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "TPZInterfaceInsertion.h"
#include "pzinterpolationspace.h"
#include "pzcompel.h"
#include "TPZVecL2.h"
#include "pzintel.h"
#include "TPZNullMaterial.h"
#include "pzgengrid.h"
#include "TPZLagrangeMultiplier.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZExtendGridDimension.h"
#include "tpzgeoelrefpattern.h"
//#include "TPZMHMDarcyDFNMeshControl.h"
#include "TPZMHMeshControl.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"

using namespace std;

const REAL Pi=M_PI;

const REAL phi_r = 0.;

TPZTransform<REAL> MHMDarcyTest::f_T(3,3);

TPZTransform<REAL> MHMDarcyTest::f_InvT(3,3);

MHMDarcyTest::MHMDarcyTest()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    
    //Materiais das condições de contorno
    fmatBCbott=-1;
    fmatBCtop=-2;
    fmatBCleft=-3;
    fmatBCright=-4;
    fmatBCtop_z=-5; //3D
    fmatBCbott_z=-6; //3D normal negativa
    
    //Material do elemento de interface
    fmatLambda=9; // Multiplier material
//    fmatLambdaBC=3;

    fmatFrac=10;
    
    fmatLambdaBC_bott=-11;
    fmatLambdaBC_top=-12;
    fmatLambdaBC_left=-13;
    fmatLambdaBC_right=-14;
    fmatLambdaBC_top_z=-15;
    fmatLambdaBC_bott_z=-16;
    
    fmatInterfaceLeft=5;
    fmatInterfaceRight=6;
    fmatWrap = 7;
    
    //Materia de um ponto
    fmatPoint=-15;
    
    //Condições de contorno do problema
    fdirichlet=0;
    fneumann=1;
    fpenetration=2;
    fpointtype=5;
    fdirichletvar=4;
    
    
    fquadmat1=1; //Parte inferior do quadrado
    fquadmat2=2; //Parte superior do quadrado
    fquadmat3=3; //Material de interface
    
    fviscosity=1.;
    fpermeability=1.;
    ftheta=-1.;
    
    fphi_r=0;
    
    f_is_hdivFull = false;
    
    f_hdivPlus = false;
    
    feltype = EQuadrilateral;
    
    f_mesh_vector.resize(4);
    
    f_T = TPZTransform<>(3,3);
    f_InvT = TPZTransform<>(3,3);
    
}

MHMDarcyTest::~MHMDarcyTest()
{
    
}

void MHMDarcyTest::Run()
{
    int int_order = fsimData.GetInternalOrder();
    int skeleton_order = fsimData.GetSkeletonOrder();
    TPZVec<int> n_s = fsimData.GetCoarseDivisions();
    TPZVec<REAL> h_s = fsimData.GetDomainSize();
    int nrefs = fsimData.GetNInterRefs();
    if(feltype==ECube||feltype==EPrisma||feltype==ETetraedro){
        Set3Dmesh();
    }
    //Gerando malha geométrica com elementos coarse, os quais serão subdomínios MHM
    TPZGeoMesh *gmesh;
    
    if (f_3Dmesh) {
        //DebugStop();
        gmesh = CreateGMesh3D(n_s, h_s);
    }else{
    //    gmesh = CreateGMeshCurve();
        gmesh = CreateGMesh(n_s, h_s);
    }
    
    //Vetor com os indices dos elementos coarse
    TPZVec<int64_t> coarseindex;
    GetElIndexCoarseMesh(gmesh, coarseindex);
    
    //Fracture
    //InsertFractureMaterial(gmesh);
    
    //Refinamento de subelemntos
    SubdomainRefine(nrefs,gmesh,coarseindex);
    
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo_0.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo_0.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Criando objeto para gerenciar a malha MHM
    TPZAutoPointer<TPZGeoMesh> gmeshpointer(gmesh);
    TPZMHMeshControl *DarcyControl;
    DarcyControl = new TPZMHMeshControl(gmeshpointer);
    DarcyControl->DefinePartitionbyCoarseIndices(coarseindex); //Define the MHM partition by the coarse element indices
    
    std::set<int> matids;
    matids.insert(fmatID);
    DarcyControl->fMaterialIds = matids;
    matids.clear();
    matids.insert(fmatBCtop);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCleft);
    matids.insert(fmatBCright);
    if (f_3Dmesh) {
        matids.insert(fmatBCtop_z);
        matids.insert(fmatBCbott_z);
    }
    
    DarcyControl->fMaterialBCIds = matids;
    
    InsertMaterialObjects(DarcyControl);
    
    DarcyControl->SetInternalPOrder(int_order);
    DarcyControl->SetSkeletonPOrder(skeleton_order);
    
    DarcyControl->DivideSkeletonElements(0); //Insere material id do skeleton wrap
    bool Add_LagrangeAverageP = true;
    DarcyControl->SetLagrangeAveragePressure(Add_LagrangeAverageP);
    

    //std::vector<int> fracture_ids(1);
    //fracture_ids[0] = fmatFrac;
    //BreakH1Connectivity(DarcyControl->CMesh(),fracture_ids); //Open conect in the fracture
    //if (fsimData.GetNInterRefs()>0) {
    //    DarcyControl->SetCoarseAverageMultipliers(true);
    //}
    

    //Malha computacional
    DarcyControl->BuildComputationalMesh(0);
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo_1.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk1("MalhaGeo_1.vtk"); //Impressão da malha geométrica (formato vtk)
    DarcyControl->GMesh()->Print(fileg1);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk1,true);
#endif

    
#ifdef PZDEBUG
    {
        //Impressão da malha computacional da velocidade (formato txt)
//        std::ofstream filecv("MalhaC_v.txt");
//        std::ofstream filecpM("MalhaC_pM.txt");
//        std::ofstream filecgM("MalhaC_gM.txt");
//        cmesh_v->Print(filecv);
        
        
        
//        cmesh_p->Print(filecp);
//        cmesh_pM->Print(filecpM);
//        cmesh_gM->Print(filecgM);
        
        std::ofstream filecm("MalhaC_MHM.txt");
        DarcyControl->CMesh()->Print(filecm);
    }
#endif
    
    {

        if (Add_LagrangeAverageP) {
            TPZCompMesh *cmeshP = DarcyControl->GetMeshes()[0].operator->();
            std::ofstream filecpress("Malha_P_MHM.txt");
            cmeshP->Print(filecpress);
            
            TPZCompMesh *cmeshL = DarcyControl->GetMeshes()[1].operator->();
            std::ofstream filecL("Malha_L_MHM.txt");
            cmeshL->Print(filecL);
            
            TPZCompMesh *cmeshAp = DarcyControl->GetMeshes()[2].operator->();
            std::ofstream filecAp("MalhaC_pM.txt");
            cmeshAp->Print(filecAp);
        }
        
        
//        std::ofstream outp("Malha_P_MHM.vtk");
//        cmeshP->LoadReferences();
//        TPZVTKGeoMesh::PrintCMeshVTK(cmeshP, outp, false);
        
//
        
        TPZCompMesh *cmesh = DarcyControl->CMesh().operator->();
        std::ofstream outv("MalhaC_P_MHM.vtk");
        cmesh->LoadReferences();
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, outv, false);
        
//        TPZCompMesh *cmeshM = DarcyControl->CMesh().operator->();
//        std::ofstream out("MalhaC_MHM.vtk");
//        cmeshM->LoadReferences();
//        TPZVTKGeoMesh::PrintCMeshVTK(cmeshM, out, false);
        
    }
    
    
   
    std::cout << "MHM Hdiv Computational meshes created\n";
    std::cout << "Number of equations MHMDarcy " << DarcyControl->CMesh()->NEquations() << std::endl;
    std::string configuration;
    
    std::stringstream MHMDarcyPref;
    MHMDarcyPref << "MHMDarcy";
    
    SolveProblem(DarcyControl->CMesh(), DarcyControl->GetMeshes(), MHMDarcyPref.str());
    
    std::cout << "FINISHED!" << std::endl;
    
}

void MHMDarcyTest::SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix){
    
    bool shapetest = fsimData.GetShapeTest();
    //calculo solution
    bool shouldrenumber = false;
    TPZAnalysis an(cmesh,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(fsimData.GetNthreads());
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(fsimData.GetNthreads());
#endif
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();

    // Shape functions plot :
    if(shapetest){
        
        TPZVec<int64_t> equ_indexes(1);
        equ_indexes[0] = 38;
        //        for (int i=37; i<38; i++) {
        //            equ_indexes[i-37] = i;
        //        }
        std::string name_phi = "MHMDarcy_shape.vtk";
        TPZVec<std::string> var_name(2);
        var_name[0]="V";
        var_name[1]="P";
        
        //TPZBuildMultiphysicsMesh::ShowShape(f_mesh_vector,cmesh_m, an, name_phi, equ_indexes);
        an.ShowShape(name_phi, equ_indexes, 1, var_name);
        return;
        
    }
    
    std::cout << "Solving\n";
    an.Solve();
    

#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    if(1){
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif

    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmesh->Print(out);
    }

    
    std::cout << "Post Processing " << std::endl;
    std::string plotfile;
    std::stringstream sout_geo;
    std::stringstream sout;
    {
        sout << prefix << "Approx_";
        ConfigPrint(sout);
        plotfile = sout.str() + "_dim2.vtk";
    }
    plotfile = "DarcyMHMPlot.vtk";
    
    {
        sout_geo << prefix << "Geo_";
        ConfigPrint(sout_geo) << "_dim2.vtk";
    }
    
    std::ofstream plotfile3(sout_geo.str());
    TPZVTKGeoMesh::PrintGMeshVTK(cmesh.operator->()->Reference(), plotfile3, true);
    
    
    std::cout << "plotfiles " << " " << plotfile.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("P");
        scalnames.Push("P_exact");
    }
    
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 2;
    an.PostProcess(resolution,cmesh->Dimension());
    
    
#endif
    
    //Calculo do erro
    std::cout << "Comuting Error (need to check for MHM?) " << std::endl;
    TPZManVector<REAL,6> Errors;
    ofstream ErroOut("Error_Brinkman.txt", std::ofstream::app);
    an.SetExact(Sol_exact);
    an.SetThreadsForError(0);
    an.PostProcessError(Errors,false);
    
    ConfigPrint(ErroOut);
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - P = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - GradP = "<< Errors[1] << std::endl;
 //   ErroOut <<"Norma L2 - P = "<< Errors[4] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();
    
    
}

std::ostream &MHMDarcyTest::ConfigPrint(std::ostream &out)
{
    int int_order = fsimData.GetInternalOrder();
    int skeleton_order = fsimData.GetSkeletonOrder();
    TPZVec<int> n_s = fsimData.GetCoarseDivisions();
    TPZVec<REAL> h_s = fsimData.GetDomainSize();
    int nrefs = fsimData.GetNInterRefs();
    
    std::string elemName;
    
    if (feltype==EQuadrilateral) {
        elemName = " Quadrilateral elements : ";
        n_s[2] = 0;
    }else if(feltype==ETriangle){
        elemName = " Triangular elements : ";
        n_s[2] = 0;
    }else if(feltype==ETetraedro){
        elemName = " Tetrahedral elements : ";
    }else if(feltype==ECube){
        elemName = " Cubic elements : ";
    }
    
    out << elemName << n_s[0] <<" x "<< n_s[1] << " x " << n_s[2] << " - N refs : " << nrefs << " - Order Skel : " << skeleton_order << " - Order Intern : " << int_order<<"\n";
    return out;
}


void MHMDarcyTest::Rotate(TPZVec<REAL> &co, TPZVec<REAL> &co_r, bool rotate){
    
    if (rotate==true) {
        //rotação +
        co_r[0] = co[0]*cos(phi_r) - co[1]*sin(phi_r);
        co_r[1] = co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }else{
        
        co_r[0] = co[0]*cos(phi_r) + co[1]*sin(phi_r);
        co_r[1] = - co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }
    
    
}

void MHMDarcyTest::InsertLowerDimMaterial(TPZGeoMesh *gmesh){
    
    // Inserir elmentos fmatLambda and fmatLambdaBCs

            int64_t nel = gmesh->NElements();
            for (int64_t el = 0; el<nel; el++) {
                TPZGeoEl *gel = gmesh->Element(el);
                if(gel->HasSubElement()&&f_allrefine)
                {
                    continue;
                }
                if (gel->Dimension() != gmesh->Dimension()) {
                    continue;
                }
                int nsides = gel->NSides();
                for (int is = 0; is<nsides; is++) {
                    if (gel->SideDimension(is) != gmesh->Dimension() - 1) {
                        continue;
                    }
                    
                    TPZGeoElSide gelside(gel,is);
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    
                    if (neighbour == gelside && f_allrefine == false) {
                        continue;
                    }
                    

                    while (neighbour != gelside) {
                        if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
                            int neigh_matID = neighbour.Element()->MaterialId();
        
                            if(neigh_matID==fmatBCbott){
                                    TPZGeoElBC(gelside, fmatLambdaBC_bott);
                            }else if(neigh_matID==fmatBCtop){
                                    TPZGeoElBC(gelside, fmatLambdaBC_top);
                            }else if(neigh_matID==fmatBCleft){
                                    TPZGeoElBC(gelside, fmatLambdaBC_left);
                            }else if(neigh_matID==fmatBCright){
                                    TPZGeoElBC(gelside, fmatLambdaBC_right);
                            }else if(f_3Dmesh && neigh_matID==fmatBCbott_z){
                                    TPZGeoElBC(gelside, fmatLambdaBC_bott_z);
                            }else if(f_3Dmesh && neigh_matID==fmatBCtop_z){
                                    TPZGeoElBC(gelside, fmatLambdaBC_top_z);
                            }
        
                            break;
        
                        }
                        if(neighbour.Element()->HasSubElement()){
                            break;
                        }
                        
                        if (gel->LowestFather()->Index()!=neighbour.Element()->LowestFather()->Index()) {
                            break;
                        }
                        
                        neighbour = neighbour.Neighbour();
        
                    }
        
        
                    if (neighbour == gelside) {
                            TPZGeoElBC(gelside, fmatLambda);
                    }
                }
            }

}

void MHMDarcyTest::InsertFractureMaterial(TPZGeoMesh *gmesh){
    
    // Inserir elmentos fmatLambda and fmatLambdaBCs
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()&&f_allrefine)
        {
            continue;
        }
        if (gel->Dimension() != gmesh->Dimension()) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {

            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            
            if (neighbour == gelside) {
                continue;
            }
            
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == gmesh->Dimension()) {
                    TPZGeoElBC(gelside, fmatFrac);
                    break;
                    
                }
                neighbour = neighbour.Neighbour();
                
            }
        }
    }
    
}

bool MHMDarcyTest::IsSkellNeighbour(TPZGeoElSide neighbour){

    if (neighbour.Element()->Dimension() == f_mesh0->Dimension()) {
        int nskellneighs = f_skellNeighs.NElements();
    
        for (int iskell = 0; iskell < nskellneighs; iskell++) {
            TPZStack<TPZGeoElSide> sonSides;
            f_skellNeighs[iskell].GetAllSiblings(sonSides);
            for (int ison=0; ison<sonSides.NElements(); ison++) {
                if (neighbour == sonSides [ison]) {
                    return true;
                }
            }
        }
    
    }
    
    return false;
}


TPZGeoMesh *MHMDarcyTest::CreateGMesh(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{
    
    int dimmodel = 2;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    x0[0] = 0., x0[1] = -1.;
    x1[0] = 2., x1[1] = 1.;
    
//    x0[0] = 0., x0[1] = 0.;
//    x1[0] = 4., x1[1] = 2.;
    
    TPZGenGrid grid(n_div,x0,x1);
    
    //grid.SetDistortion(0.5);
    grid.SetRefpatternElements(true);
    if (feltype==ETriangle) {
        grid.SetElementType(ETriangle);
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    grid.Read(gmesh);
    grid.SetBC(gmesh, 4, fmatBCbott);
    grid.SetBC(gmesh, 5, fmatBCright);
    grid.SetBC(gmesh, 6, fmatBCtop);
    grid.SetBC(gmesh, 7, fmatBCleft);
    
    //Save the original mesh
    
    //SetAllRefine();
    
    TPZVec<REAL> centerCo(2,0.);
    centerCo[0]=1.;
    centerCo[1]=0.;
   // UniformRefine(1, gmesh, centerCo, true);

    //UniformRefine2(1, gmesh, n_div);
//    InsertLowerDimMaterial(gmesh);
    SetOriginalMesh(gmesh);
//    UniformRefine2(1, gmesh, n_div);
//    InsertLowerDimMaterial(gmesh);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();


        {
            std::ofstream Dummyfile("GeometricMesh2d.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
        }
    
    return gmesh;
    
}

TPZGeoMesh *MHMDarcyTest::CreateGMesh3D(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{
    
    int dimmodel = 2;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    x0[0] = 0., x0[1] = -1.;
    x1[0] = 2., x1[1] = 1.;
    
    x0[2] = 0.;
    x1[2] = 2.;
    
    TPZGenGrid grid(n_div,x0,x1);
    
    
    if (feltype == ETriangle|| feltype == EPrisma ) {
        grid.SetElementType(ETriangle);
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh;

    //grid.SetDistortion(0.2);
    
    if (feltype==EQuadrilateral||feltype==ECube||feltype==EPrisma||feltype==ETriangle) {
    
        grid.Read(gmesh);
        grid.SetBC(gmesh, 4, fmatBCbott);
        grid.SetBC(gmesh, 5, fmatBCright);
        grid.SetBC(gmesh, 6, fmatBCtop);
        grid.SetBC(gmesh, 7, fmatBCleft);
        
        
        REAL thickness = h_s[2]/n_div[2];
        TPZExtendGridDimension extend(gmesh,thickness);
        int numlayers = n_div[2];
        extend.SetElType(1);///???
        gmesh = extend.ExtendedMesh(numlayers,fmatBCbott_z,fmatBCtop_z);
        gmesh->SetDimension(3);
    
    } else if(feltype==ETetraedro){
        
        gmesh->SetDimension(3);
        int64_t id = 0;
        int nx = n_div[0]+1, ny = n_div[1]+1, nz = n_div[2]+1;
        TPZVec<REAL> coord(3,0.);
        int nnodes=nx*ny*nz;
        gmesh->NodeVec().Resize(nx*ny*nz);
        for (int k=0; k<nz; k++) {
            for(int i = 0; i < ny; i++){
                for(int j = 0; j < nx; j++){
                    id = i*nx + j+ k*nx*ny;
                    coord[0] = (j)*h_s[0]/(nx - 1);
                    coord[1] = (i)*h_s[1]/(ny - 1)-1;
                    coord[2] = (k)*h_s[2]/(nz - 1);
                    gmesh->NodeVec()[id].Initialize(coord, *gmesh);
                }
            }
        }
        
        
        TPZVec<int64_t> nodeindD1(4,0), nodeindD2(4,0), nodeindU1(4,0), nodeindU2(4,0), nodeindL1(4,0), nodeindL2(4,0);
        int64_t index=0;
        
        for(int kq=0; kq<n_div[2]; kq++){
            for(int iq = 0; iq < n_div[1]; iq++){
                for(int jq = 0; jq < n_div[0]; jq++){
                    
        
                    // Plano xy
                    nodeindD1[0] = (iq)*ny + (jq) + kq*nx*ny;
                    nodeindD1[1] = nodeindD1[0]+1;
                    nodeindD1[2] = nodeindD1[0]+nx;
                    nodeindD1[3] = nodeindD1[1] + (1)*nx*ny;
                    gmesh->CreateGeoElement(ETetraedro, nodeindD1, fmatID, index,1);
                    
                    index++;
                    
                    nodeindD2[0] = nodeindD1[1];
                    nodeindD2[1] = nodeindD1[2];
                    nodeindD2[2] = nodeindD1[1]+nx;
                    nodeindD2[3] = nodeindD1[1] + (1)*nx*ny;
                    gmesh->CreateGeoElement(ETetraedro, nodeindD2, fmatID, index,1);
                    
                    index++;
                    
                    nodeindU1[0] = nodeindD1[0] + (1)*nx*ny;
                    nodeindU1[1] = nodeindU1[0]+1;
                    nodeindU1[2] = nodeindU1[0]+nx;
                    nodeindU1[3] = nodeindD1[2];
                    gmesh->CreateGeoElement(ETetraedro, nodeindU1, fmatID, index,1);
                    
                    index++;
                    
                    nodeindU2[0] = nodeindU1[1];
                    nodeindU2[1] = nodeindU1[2];
                    nodeindU2[2] = nodeindU1[1]+nx;
                    nodeindU2[3] = nodeindD1[2];
                    gmesh->CreateGeoElement(ETetraedro, nodeindU2, fmatID, index,1);
                    
                    index++;
                    
                    // Plano xz
                    
                    nodeindL1[0] = nodeindD1[0];
                    nodeindL1[1] = nodeindD1[2];
                    nodeindL1[2] = nodeindL1[0]+nx*ny;
                    nodeindL1[3] = nodeindU1[1];
                    gmesh->CreateGeoElement(ETetraedro, nodeindL1, fmatID, index,1);
                    
                    index++;
                    
                    
                    nodeindL2[0] = nodeindD2[2];
                    nodeindL2[1] = nodeindD1[1]+nx*ny;
                    nodeindL2[2] = nodeindU2[2];
                    nodeindL2[3] = nodeindU2[3];
                    gmesh->CreateGeoElement(ETetraedro, nodeindL2, fmatID, index,1);
                    
                    index++;
                    
                }
            }
        }

        gmesh->BuildConnectivity();
        
        // Boundary Conditions
        const int numelements = gmesh->NElements();
        for(int el=0; el<numelements; el++)
        {
            TPZManVector <TPZGeoNode,4> Nodefinder(4);
            TPZManVector <REAL,3> nodecoord(3);
            TPZGeoEl *tetra = gmesh->ElementVec()[el];
            
            // na face x = 0
            TPZVec<int64_t> ncoordVec(0); int64_t sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[0]-x0[0])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCleft);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face x = 1
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[0]-x1[0])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCright);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face y = 0
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                
                
                if (fabs(nodecoord[1]-x0[1])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCbott);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face y = 1
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[1]-x1[1])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCtop);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face z = 0
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[2]-x0[2])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCbott_z);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face z = 1
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[2]-x1[2])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCtop_z);
            }
            
        }
        
    
    }
    //Save the original mesh
    grid.SetRefpatternElements(true);
    
    //SetAllRefine();
    
    TPZVec<REAL> centerCo(2,0.);
    centerCo[0]=1.;
    centerCo[1]=0.;
    // UniformRefine(1, gmesh, centerCo, true);
    
    //UniformRefine2(1, gmesh, n_div);
    
    //InsertLowerDimMaterial(gmesh);
    SetOriginalMesh(gmesh);
    
    //UniformRefine2(1, gmesh, n_div);
    //InsertLowerDimMaterial(gmesh);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    
    {
        std::ofstream Dummyfile("GeometricMesh3D.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    return gmesh;
    
}

TPZGeoMesh *MHMDarcyTest::CreateGMeshCurve()
{
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(2);
    
    int nodes = 6;
    REAL radius = 1.0;
    REAL innerradius = radius/2.0;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,6> TopolQQuadrilateral(6);
    TPZManVector<int64_t,8> TopolQuadrilateral(4);
    TPZManVector<int64_t,6> TopolQTriangle(6);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<int64_t,3> TopolArc(3);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    
    int64_t nodeindex = 0;
    
    for (int inode = 0; inode < 3 ; inode++) {
        // i node
        coord = ParametricCircle(radius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    for (int inode = 0; inode < 3 ; inode++) {
        // i node
        coord = ParametricCircle(innerradius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    //Ponto 1
    int64_t elementid = 0;
    
    TopolQuadrilateral[0] = 3;
    TopolQuadrilateral[1] = 0;
    TopolQuadrilateral[2] = 2;
    TopolQuadrilateral[3] = 5;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    elementid++;
    
    // outer arcs bc's
    
    TopolLine[0] = 3;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCbott,*geomesh);
    elementid++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCtop,*geomesh);
    elementid++;
    
    TopolArc[0] = 0;
    TopolArc[1] = 2;
    TopolArc[2] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCright,*geomesh);
    elementid++;
    
    TopolArc[0] = 5;
    TopolArc[1] = 3;
    TopolArc[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCleft,*geomesh);
    elementid++;
    
    
    
    geomesh->BuildConnectivity();
    
    int nref = 0;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    TPZCheckGeom check(geomesh);
    check.CheckUniqueId();
    geomesh->BuildConnectivity();
    
    
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    return geomesh;
    
}

TPZManVector<REAL,3>  MHMDarcyTest::ParametricCircle(REAL radius,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta);
    xcoor[1] = radius * sin(theta);
    xcoor[2] = 0.0 ;
    return xcoor;
}

void MHMDarcyTest::TetrahedralMeshCubo(TPZVec<int> &n_s){
    
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    GenerateNodes(gmesh,nelem);
//
//
//
//    return gmesh;
    
}

void MHMDarcyTest::BreakH1Connectivity(TPZCompMesh &cmesh, std::vector<int> fracture_ids)
{
    std::set<int> boundaries_ids;
    boundaries_ids.insert(fmatBCbott);
    boundaries_ids.insert(fmatBCleft);
    boundaries_ids.insert(fmatBCtop);
    boundaries_ids.insert(fmatBCright);
    
    for (unsigned int i_f = 0; i_f <  fracture_ids.size(); i_f++) {
        TPZFractureInsertion fracture(cmesh.Reference(),fracture_ids[i_f],boundaries_ids);
        fracture.ClassifyNeighboursofPivots();
        fracture.OpenFractureOnH1(&cmesh); // (ok)
        fracture.SetDiscontinuosFrac(&cmesh); // (ok)
        //fracture.SetInterfaces(&cmesh, fmatInterfaceLeft[i_f], fmatInterfaceRight[i_f]);
        //fractureInsert[i_f]=fracture;
    }
    cmesh.ComputeNodElCon();
    
}


void MHMDarcyTest::SubdomainRefine(int nrefine, TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices){
    
    int ncoarse = coarseindices.size();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int i_ref = 0; i_ref < nrefine; i_ref++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
         
            // BC elements
//            if(gel->MaterialId()<0){
//                    if(!gel->HasSubElement()) gel->Divide(filhos);
//            }
            
            if(gel->Dimension()!=gmesh->Dimension()){
                continue;
            }
            
            // Coarse volumetric elements
            TPZGeoEl * higher_el = gel->LowestFather();

            for (int i_coarse = 0; i_coarse<ncoarse; i_coarse++) {
                if(higher_el->Index()!=coarseindices[i_coarse]) continue;
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
        
    }
    
//    int64_t nel = coarseindices.size();
//    for (int64_t el=0; el<nel; el++) {
//        TPZGeoEl *gel = gmesh->Element(coarseindices[el]);
//
//        int nsubel = gel->NSubElements();
//
//        for (int iref = 0; iref<nrefine; iref++)
//        {
//            TPZManVector<TPZGeoEl *,10> gelsub;
//            gel->Divide(gelsub);
//        }
//    }
    
}




void MHMDarcyTest::UniformRefine4(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            higher_el->CenterPoint(side, centerMaster);
            higher_el->X(centerMaster,centerEuclid);
            
            
            if (fabs(centerCo[0]-centerEuclid[0]) > 1.e-9 &&  restriction == true) {
                continue;
            }
            if (fabs(centerCo[1]-centerEuclid[1]) > 1.e-9 && restriction == true) {
                continue;
            }
            
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        break;
             //           neighbour.Element()->Divide(filhos);
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}


void MHMDarcyTest::UniformRefine3(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        int count =0.;
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            if(higher_el->Dimension()!=2) continue;
            
            int intdiv = (higher_el->Index()/(n_div[0]*2))%2;
            if (intdiv==0) {
                count=0;
            }else{
                count=1;
            }

            count =0; //papapapa
            if((higher_el->Index()+count)%2!=0) continue;
            
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        
                        if(f_allrefine==false){
                            break;
                        }
                        neighbour.Element()->Divide(filhos);
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}


void MHMDarcyTest::UniformRefine2(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        int count =0.;
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            if(higher_el->Dimension()!=2) continue;
            
            int intdiv = (higher_el->Index()/n_div[0])%2;
            if (intdiv==0) {
                count=0;
            }else{
                count=1;
            }
            
            //if((higher_el->Index()+count)%2==0) continue;
          
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        
                        if (f_allrefine==false) {
                            break;
                        }
                        
                        neighbour.Element()->Divide(filhos);
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}



void MHMDarcyTest::UniformRefine(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            higher_el->CenterPoint(side, centerMaster);
            higher_el->X(centerMaster,centerEuclid);
            
            if (fabs(centerCo[0]-centerEuclid[0]) > 1.e-9 &&  restriction == true) {
                continue;
            }
            if (fabs(centerCo[1]-centerEuclid[1]) > 1.e-9 && restriction == true) {
                continue;
            }
            
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        if (f_allrefine==false) {
                            break;
                        }
                        neighbour.Element()->Divide(filhos);
                        
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

TPZCompEl *MHMDarcyTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}

void MHMDarcyTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
        sol.Resize(1);

        REAL x1 = x[0];
        REAL x2 = x[1];

//        STATE pressure = 1.-0.2*x1;
//        STATE pressure = 0.;
    STATE pressure= cos(x1)*sin(x2);
    
    sol[0]=pressure;



    
    
}

void MHMDarcyTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    //Applying rotation:
    TPZVec<REAL> x_in = x;
    TPZVec<REAL> x_rot(3,0.);
    
    f_InvT.Apply(x_in,x_rot);
    x[0] = x_rot[0];
    x[1] = x_rot[1];
    
    f.resize(1);
    REAL x1 = x[0];
    REAL x2 = x[1];

    
    STATE f_source = 2.*cos(x1)*sin(x2);
    
    f[0] = f_source;
    
    TPZVec<REAL> f_s(3,0), f_rot(3,0);
    
    // General form : : Artigo Botti, Di Pietro, Droniou
    
    //    REAL m_v= 1., m_u= 1.0;
    //
    //    REAL Cf=m_v/m_u;
    //
    //        f_1 = -sin(x1)*sin(x2)-exp(-Cf)*sin(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*sin(x1)*sin(x2)-m_u*(2.*exp(-Cf)*sin(x1)*sin(x2)-(1./m_v)*4.*(1-exp(-Cf))*sin(x1)*sin(x2));
    //
    //        f_2 = cos(x1)*cos(x2)-exp(-Cf)*cos(x1)*cos(x2)-(1./m_v)*(1.-exp(-Cf))*cos(x1)*cos(x2)-m_u*(2.*exp(-Cf)*cos(x1)*cos(x2)+(1./m_v)*4.*(1-exp(-Cf))*cos(x1)*cos(x2));
    //        STATE g_1 = (2./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
    //
    //        f[0] = f_1; // x direction
    //        f[1] = f_2; // y direction
    //
    //        f[2] = g_1; // g source
    
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou
    
    //    REAL e = exp(1.);
    //
    //    f_1 = (-8./e+ 4.)*sin(x1)*sin(x2);
    //    f_2 = (2./e- 4.)*cos(x1)*cos(x2);
    //    STATE g_1 = 2.*(1.-1./e)*cos(x1)*sin(x2);
    //
    //    f[0] = f_1; // x direction
    //    f[1] = f_2; // y direction
    //
    //    f[2] = g_1; // g source
    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
    
//    f_s[0] = -3.*sin(x1)*sin(x2);
//    f_s[1] = -1.*cos(x1)*cos(x2);
//
//    f_T.Apply(f_s, f_rot);
//    f_s = f_rot;
//
//
//    f[0] = f_s[0]; // x direction
//    f[1] = f_s[1]; // y direction
//    f[2] = f_s[2];
    
    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
    //        f_1 = 0.;
    //        f_2 = 0.;
    //
    //        f[0] = f_1; // x direction
    //        f[1] = f_2; // y direction
    //        f[2] = 2.*cos(x1)*sin(x2);
    
    
    // Darcy 3D : Artigo Botti, Di Pietro, Droniou
    
    
//        f_s[0] = 2.*cos(x1)*cos(x3) - 1.*(2. + cos(x3))*sin(x1)*sin(x2);
//        f_s[1] = cos(x1)*cos(x2)*(-2. + cos(x3));
//        f_s[2] = (2.*sin(x1) - cos(x1)*sin(x2))*sin(x3);
//
//        f_T.Apply(f_s, f_rot);
//        f_s = f_rot;
//
//
//        f[0] = f_s[0]; // x direction
//        f[1] = f_s[1]; // y direction
//        f[2] = f_s[2];
    
}

void MHMDarcyTest::ChangeExternalOrderConnects(TPZCompMesh *mesh, int addToOrder){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        int nshape2 = 0;
        
        if(cel->Dimension()== dim)
        {
            TPZConnect &conel = cel->Connect(ncon-1);
            corder = conel.Order();
            nshape = conel.NShape();
            
            int neworder = corder + addToOrder;//Aqui = +1
            int64_t cindex = cel->ConnectIndex(ncon-1);
            conel.SetOrder(neworder,cindex);
            
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            intel->SetPreferredOrder(neworder);
            nshape = intel->NConnectShapeF(ncon-1,neworder);
            
            if(dim==2 && addToOrder==1)
            {
                if(feltype==ETriangle){
                    nshape2 = (corder + 2)*(corder + 2)-1;
                }else{//Quadrilateral
                    nshape2 = 2*(corder + 1)*(corder + 2);
                }
                if(nshape2!=nshape)
                {
                    DebugStop();
                }
            }
            
            conel.SetNShape(nshape);
            mesh->Block().Set(conel.SequenceNumber(),nshape);
        }
    }
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
}

void MHMDarcyTest::InsertMaterialObjects(TPZMHMeshControl *control)
{

    TPZCompMesh &cmesh = control->CMesh();

    REAL visco = fsimData.GetViscosity();
    TPZMHMDarcyDFNMaterial *mat1 = new TPZMHMDarcyDFNMaterial(fmatID,fdim,1,visco,0,0);
    
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source, 9);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact,9);
    ((TPZDummyFunction<STATE>*)fp.operator->())->SetPolynomialOrder(9);
    ((TPZDummyFunction<STATE>*)solp.operator->())->SetPolynomialOrder(9);
    
    mat1->SetForcingFunction(fp); //Caso simples sem termo fonte
    mat1->SetForcingFunctionExact(solp);
    
    if (fsimData.GetShapeTest()) {
        mat1->SetForcingFunction(NULL);
        mat1->SetForcingFunctionExact(NULL);
    }
    
    //TPZMaterial * mat1(material);
    cmesh.InsertMaterialObject(mat1);
    
    int matSkeleton = 4;

    control->SwitchLagrangeMultiplierSign(false);
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    
    TPZBndCond * BCondD1 = mat1->CreateBC(mat1, fmatBCbott, fdirichlet, val1, val2);
    BCondD1->SetBCForcingFunction(0, solp);
    cmesh.InsertMaterialObject(BCondD1);
    //control->fMaterialBCIds.insert(fmatBCbott);
    
    TPZBndCond * BCondD2 = mat1->CreateBC(mat1, fmatBCtop, fdirichlet, val1, val2);
    BCondD2->SetBCForcingFunction(0, solp);
    cmesh.InsertMaterialObject(BCondD2);
    //control->fMaterialBCIds.insert(fmatBCtop);
    
    TPZBndCond * BCondD3 = mat1->CreateBC(mat1, fmatBCleft, fdirichlet, val1, val2);
    BCondD3->SetBCForcingFunction(0, solp);
    cmesh.InsertMaterialObject(BCondD3);
    //control->fMaterialBCIds.insert(fmatBCleft);
    
    TPZBndCond * BCondD4 = mat1->CreateBC(mat1, fmatBCright, fdirichlet, val1, val2);
    BCondD4->SetBCForcingFunction(0, solp);
    cmesh.InsertMaterialObject(BCondD4);
    //control->fMaterialBCIds.insert(fmatBCright);

    if (f_3Dmesh) {
        TPZBndCond * BCondD5 = mat1->CreateBC(mat1, fmatBCtop_z, fdirichlet, val1, val2);
        BCondD5->SetBCForcingFunction(0, solp);
        cmesh.InsertMaterialObject(BCondD5);
        
        TPZBndCond * BCondD6 = mat1->CreateBC(mat1, fmatBCbott_z, fdirichlet, val1, val2);
        BCondD6->SetBCForcingFunction(0, solp);
        cmesh.InsertMaterialObject(BCondD6);
    }
    
    
    //Skeleton::
    
    
//    TPZBndCond * bcFlux = mat1->CreateBC(mat1, matSkeleton, fneumann, val1, val2);
//    //bcFlux->SetBCForcingFunction(0, solp);
//    cmesh.InsertMaterialObject(bcFlux);
    
    
    
    // 2.1 - Material para tração tangencial 1D (Interior)
//    TPZNullMaterial *matLambda = new TPZNullMaterial(fmatLambda);
//    matLambda->SetDimension(fdim-1);
//    matLambda->SetNStateVariables(1);
//    control->CMesh()->InsertMaterialObject(matLambda);
 //   control->fMaterialIds.insert(fmatLambda);
    
    
    // 2.2 - Material for interfaces (Interior)
//    TPZMHMDarcyDFNMaterial *matInterfaceLeft = new TPZMHMDarcyDFNMaterial(control->fLagrangeMatIdLeft,fdim,1,visco,0,0);
//    matInterfaceLeft->SetMultiplier(1.);
//    cmesh.InsertMaterialObject(matInterfaceLeft);
 //   control->fLagrangeMatIdLeft=fmatInterfaceLeft;
    
//    TPZMHMDarcyDFNMaterial *matInterfaceRight = new TPZMHMDarcyDFNMaterial(control->fLagrangeMatIdRight,fdim,1,visco,0,0);
//    matInterfaceRight->SetMultiplier(-1.);
//    cmesh.InsertMaterialObject(matInterfaceRight);
 //   control->fLagrangeMatIdRight=fmatInterfaceRight;
    
    // 3.1 - Material para tração tangencial 1D nos contornos
//    TPZBndCond *matLambdaBC_bott = mat1->CreateBC(mat1, fmatLambdaBC_bott, fneumann, val1, val2);
//    matLambdaBC_bott->SetBCForcingFunction(0, solp);
//    cmesh.InsertMaterialObject(matLambdaBC_bott);
 //   control->fMaterialBCIds.insert(fmatLambdaBC_bott);
    
//    TPZBndCond *matLambdaBC_top = mat1->CreateBC(mat1, fmatLambdaBC_top, fneumann, val1, val2);
//    matLambdaBC_top->SetBCForcingFunction(0, solp);
//    cmesh.InsertMaterialObject(matLambdaBC_top);
 //   control->fMaterialBCIds.insert(fmatLambdaBC_top);
    
//    TPZBndCond *matLambdaBC_left = mat1->CreateBC(mat1, fmatLambdaBC_left, fneumann, val1, val2);
//    matLambdaBC_left->SetBCForcingFunction(0, solp);
//    cmesh.InsertMaterialObject(matLambdaBC_left);
 //   control->fMaterialBCIds.insert(fmatLambdaBC_left);
    
//    TPZBndCond *matLambdaBC_right = mat1->CreateBC(mat1, fmatLambdaBC_right, fneumann, val1, val2);
//    matLambdaBC_right->SetBCForcingFunction(0, solp);
//    cmesh.InsertMaterialObject(matLambdaBC_right);
 //   control->fMaterialBCIds.insert(fmatLambdaBC_right);
    
    if (f_3Dmesh) {
        TPZBndCond *matLambdaBC_top_z = mat1->CreateBC(mat1, fmatLambdaBC_top_z, fneumann, val1, val2);
        matLambdaBC_top_z->SetBCForcingFunction(0, solp);
        cmesh.InsertMaterialObject(matLambdaBC_top_z);
        
        TPZBndCond *matLambdaBC_bott_z = mat1->CreateBC(mat1, fmatLambdaBC_bott_z, fneumann, val1, val2);
        matLambdaBC_bott_z->SetBCForcingFunction(0, solp);
        cmesh.InsertMaterialObject(matLambdaBC_bott_z);
    }
    
}


void MHMDarcyTest::InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m){
    
    std::set<int> boundaries_ids;
    boundaries_ids.insert(fmatBCbott);
    boundaries_ids.insert(fmatBCleft);
    boundaries_ids.insert(fmatBCtop);
    boundaries_ids.insert(fmatBCright);
    if (f_3Dmesh) {
        boundaries_ids.insert(fmatBCbott_z);
        boundaries_ids.insert(fmatBCtop_z);
    }
    
    
    TPZInterfaceInsertion InterfaceInsertion(cmesh_m, fmatLambda, boundaries_ids, feltype);
    TPZManVector<int64_t,3> Interfaces(2,0);
    Interfaces[0] = fmatInterfaceLeft;
    Interfaces[1] = fmatInterfaceRight;
    InterfaceInsertion.SetInterfaceVectorId(Interfaces);
    
    if (f_allrefine) {
        InterfaceInsertion.AddMultiphysicsInterfacesLeftNRight(fmatLambda);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_bott,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_top,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_left,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_right,fmatInterfaceLeft);
    }else{
        InterfaceInsertion.AddMultiphysicsInterfacesLeftNRight2(fmatLambda);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_bott,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_top,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_left,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_right,fmatInterfaceLeft);
        if (f_3Dmesh) {
            InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_bott_z,fmatInterfaceLeft);
            InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_top_z,fmatInterfaceLeft);
        }
        
    }
    
    
    
}

void MHMDarcyTest::ComputeSkelNeighbours(){
    
    if (!f_mesh0) {
        DebugStop();
    }
    
    int64_t nel = f_mesh0->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = f_mesh0->Element(el);
//        if(gel->HasSubElement()&&f_allrefine)
//        {
//            continue;
//        }
        if (gel->MaterialId() != fmatLambda) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();

        if (neighbour == gelside) {
            continue;
        }
        
        while (neighbour != gelside) {
            int neigh_matID = neighbour.Element()->MaterialId();
            if (neighbour.Element()->Dimension() == f_mesh0->Dimension() && neigh_matID==fmatID) {
                f_skellNeighs.Push(neighbour);
            }
            if(neighbour.Element()->HasSubElement()){
                break;
            }
            neighbour = neighbour.Neighbour();
        }
    }
    
}


void MHMDarcyTest::GroupAndCondense(TPZMultiphysicsCompMesh *cmesh_m){
   
    
    //Criando apropamento de elementos
    
    int64_t ncompel = cmesh_m->ElementVec().NElements();
    int dim = cmesh_m->Reference()->Dimension();

    std::vector<int64_t> GroupIndex;
    TPZStack<TPZElementGroup *> elgroups;
    int count = 0;
    int64_t index =0;
    
    for(int64_t el = 0; el < ncompel; el++){
        
        TPZCompEl *cel = cmesh_m->Element(el);
        if (cel->Dimension()!=dim) {
            continue;
        }
        //GroupIndex[el] = cel->Index();
        count++;
        GroupIndex.resize(count);
        GroupIndex[count-1]=cel->Index();
        TPZElementGroup *GroupEl = new TPZElementGroup(*cmesh_m,index);
        elgroups.Push(GroupEl);
        elgroups[count-1]->AddElement(cel);
    }
    

    //Inserindo as respectivas interfaces e condições de contorno
    
    for(int64_t el = 0; el < ncompel; el++){
        TPZCompEl *cel = cmesh_m->Element(el);
        
        TPZMultiphysicsInterfaceElement *interel = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (interel) {
            TPZCompEl *Leftel = interel->LeftElement();
            
            if (Leftel->Dimension()!=dim) {
                continue;
            }
            int leftindex = Leftel->Index();
            
            for(int64_t iel = 0; iel < GroupIndex.size(); iel++){
                if (leftindex==GroupIndex[iel]) {
                    elgroups[iel]->AddElement(cel);
                }
            }
        }
        
        if (!cel) {
            continue;
        }
        
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension()==dim-1) {
            
            TPZBndCond *elBC = dynamic_cast<TPZBndCond *>(cel->Material());
            if (!elBC) {
                continue;
            }
            
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel, gel->NSides() - 1);
            
            gelside.EqualLevelCompElementList(celstack, 0, 0);
            
            for (auto &celstackindex : celstack) {
                if (celstackindex.Reference().Element()->Dimension() == dim) {
                    int bcindex = celstackindex.Element()->Index();
                    
                    for(int64_t iel = 0; iel < GroupIndex.size(); iel++){
                        if (bcindex==GroupIndex[iel]) {
                            elgroups[iel]->AddElement(cel);
                        }
                    }
                }
            }
        }
    }
    
    cmesh_m->ComputeNodElCon();
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation
//    for (int64_t ienv=0; ienv<nenvel; ienv++) {
//        TPZElementGroup *elgr = elgroups[ienv];
    int nenvel = elgroups.NElements();
    for (int64_t ienv=0; ienv<nenvel; ienv++) {
        TPZElementGroup *elgr = elgroups[ienv];
        
        int nc = elgroups[ienv]->GetElGroup()[0]->NConnects();
        elgroups[ienv]->GetElGroup()[0]->Connect(nc-1).IncrementElConnected();
        
        
//        for (int ic=0; ic<nc; ic++) {
//            TPZConnect &c = elgr->Connect(ic);
//            int connectpM = elgroups[ienv]->GetElGroup()[0]->NConnects();
//            int nc = elgr->NConnects();
//            TPZConnect &c = elgr->Connect(nc-1);
//            if (c.LagrangeMultiplier() > 0) {
//                c.IncrementElConnected();
//                break;
//            }
//        }
        new TPZCondensedCompEl(elgr);
    }

    
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
    
}


void MHMDarcyTest::GetElIndexCoarseMesh(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindex)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel=0;
    int dim = gmesh->Dimension();
    int eldim;
    int count =0.;
    
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        hassubel = gel->HasSubElement();
        eldim = gel->Dimension();
        if(!hassubel && eldim ==dim)
        {
            count++;
            coarseindex.resize(count);
            coarseindex[count-1] = gel->Index();
        }
    }
    
}

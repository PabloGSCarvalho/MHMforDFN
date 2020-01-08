/*
 *  HybridBrinkmanTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "HybridBrinkmanTest.h"
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
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGenSpecialGrid.h"
#include "tpzgeoelmapped.h"

using namespace std;

const REAL Pi=M_PI;

const REAL phi_r = 0.;

TPZTransform<REAL> HybridBrinkmanTest::f_T(3,3);

TPZTransform<REAL> HybridBrinkmanTest::f_InvT(3,3);

HybridBrinkmanTest::HybridBrinkmanTest()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    
    //Materiais das condições de contorno
    fmatBCbott=-1;
    fmatBCtop=-2;
    fmatBCleft=-3;
    fmatBCright=-4;
    fmatBChole=-7;

    fmatBCtop_z = -5; //3D
    fmatBCbott_z = -6; //3D normal negativa
    
    //Material do elemento de interface
    fmatLambda=4; // Multiplier material
    fmatLambdaBC=3;
    
    fmatLambdaBC_bott=11;
    fmatLambdaBC_top=12;
    fmatLambdaBC_left=13;
    fmatLambdaBC_right=14;
    fmatLambdaBC_hole=17;
    
    fmatLambdaBC_top_z=15;
    fmatLambdaBC_bott_z=16;
    
    fmatWrapBC_bott=21;
    fmatWrapBC_top=22;
    fmatWrapBC_left=23;
    fmatWrapBC_right=24;
    
    fmatInterfaceLeft=5;
    fmatInterfaceRight=6;
    fmatWrap = 7;
    
    //Materiais das condições de contorno (elementos de interface)
    fmatIntBCbott=-11;
    fmatIntBCtop=-12;
    fmatIntBCleft=-13;
    fmatIntBCright=-14;
    fmatIntBCtop_z=-15;
    fmatIntBCbott_z=-16;
    
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
    
    fSpaceV=0;
    
    fphi_r=0;
    
    f_inter_ref = 0.;
    
    f_is_hdivFull = false;
    
    f_hdivPlus = false;
    
    f_Holemesh = false;
    
    feltype = EQuadrilateral;
    
    f_mesh_vector.resize(2);
    
    f_T = TPZTransform<>(3,3);
    f_InvT = TPZTransform<>(3,3);
    
    f_HoleCoord.clear();
    f_ArcCentralNode.clear();
    
}

HybridBrinkmanTest::~HybridBrinkmanTest()
{
    
}

void HybridBrinkmanTest::Run(int Space, int pOrder, TPZVec<int> &n_s, TPZVec<REAL> &h_s, STATE visco)
{
    
    //Gerando malha geométrica:
    fSpaceV = Space;
    TPZGeoMesh *gmesh;
    
    if (f_3Dmesh) {
        gmesh = CreateGMesh3D(n_s, h_s);
    }else{
        gmesh = CreateGMesh(n_s, h_s);
        //gmesh = CreateGMeshRefPattern(n_s, h_s);
        //gmesh = CreateGMeshCurve();
        //gmesh = CreateGMeshCurveBlend();
        //gmesh = CreateGMeshCurveBlendSimple();
    }
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeoB.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeoB.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    int n_mais = 0;
    if (f_hdivPlus) {
        n_mais = 1;
    }
    
    TPZCompMesh *cmesh_v = this->CMesh_v(gmesh, Space, pOrder);
    TPZCompMesh *cmesh_p = this->CMesh_p(gmesh, Space, pOrder+n_mais);
    
    TPZCompMesh *cmesh_pM = this->CMesh_pM(gmesh, 0);
    TPZCompMesh *cmesh_gM = this->CMesh_gM(gmesh, 0);
    
    if (!f_mesh0) {
        DebugStop();
    }
    
//    TPZCompMesh *cmesh_pM_0 = this->CMesh_pM_0(f_mesh0, 0);
//    TPZCompMesh *cmesh_gM_0 = this->CMesh_gM_0(f_mesh0, 0);
    
    ChangeExternalOrderConnects(cmesh_v,n_mais);
    // ChangeExternalOrderConnects(cmesh_p,n_mais);
    
  
    f_mesh_vector[0]=cmesh_v;
    f_mesh_vector[1]=cmesh_p;
//    f_mesh_vector[2]=cmesh_pM;
//    f_mesh_vector[3]=cmesh_gM;
    
    TPZMultiphysicsCompMesh *cmesh_m = this->CMesh_m(gmesh, Space, pOrder, visco); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecv("MalhaC_v.txt");
        std::ofstream filecp("MalhaC_p.txt");
        std::ofstream filecpM("MalhaC_pM.txt");
        std::ofstream filecgM("MalhaC_gM.txt");
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
//        cmesh_pM->Print(filecpM);
//        cmesh_gM->Print(filecgM);
        
        std::ofstream filecm("MalhaC_m.txt");
        cmesh_m->Print(filecm);
    }
#endif
    
    
    cmesh_m->LoadReferences();
    InsertInterfaces(cmesh_m);
    
    
    //    AddMultiphysicsInterfacesLeftNRight(*cmesh_m,fmatLambda); // Rever isto aqui
    //    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCbott,fmatBCbott);
    //    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCtop,fmatBCtop);
    //    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCleft,fmatBCleft);
    //    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCright,fmatBCright);
    
    //    AddMultiphysicsInterfaces(*cmesh_m);
    
#ifdef PZDEBUG
    std::ofstream filecmbfCond("MalhaC_m_beforeCond.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecmbfCond);
#endif
    
    
    // Agrupar e condensar os elementos
    //GroupAndCondense(cmesh_m);

#ifdef PZDEBUG
    {
        std::ofstream fileg1("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
        std::ofstream filegvtk1("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
        gmesh->Print(fileg1);
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk1,true);
        
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
#endif
    
    
    //Resolvendo o Sistema:
    int numthreads = 3;
    
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    
    TPZSymetricSpStructMatrix struct_mat(cmesh_m);
    struct_mat.SetNumThreads(numthreads);
    an.SetStructuralMatrix(struct_mat);
    
    
    //TPZParSkylineStructMatrix matskl(cmesh_m, numthreads);
    
//        TPZSkylineStructMatrix matskl(cmesh_m); //OK para Hdiv
//        matskl.SetNumThreads(numthreads);
//        an.SetStructuralMatrix(matskl);
//    //
////        if (Space==1) {
//            TPZFStructMatrix matsklD(cmesh_m); //caso nao simetrico *** //OK para discont.
//            matsklD.SetNumThreads(numthreads);
//            an.SetStructuralMatrix(matsklD);
//        }
    
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    an.Assemble(); //Assembla a matriz de rigidez (e o vetor de carga) global
//    {
//        std::ofstream filestiff("stiffness_before.txt");
//        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
//
//    }
    std::cout << "Solving Matrix " << std::endl;
    an.Solve();
    
//        {
//            int eqGm= cmesh_gM->Solution().Rows();
//            TPZFMatrix<STATE> SolTriky(eqGm,1,1.0);
//            cmesh_gM->LoadSolution(SolTriky);
//            std::cout<<SolTriky<<std::endl;
//        }
    
#ifdef PZDEBUG
    if(0){
        std::ofstream filecv("MalhaC_v2.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p2.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
        
        std::ofstream filecm("MalhaC_m2.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
#endif
    
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    if(1){
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    if(0){
        TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
        
        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
#endif
    
    
    // Shape functions plot :
    
    if(0){
        int64_t target_index = 1;
        
        TPZVec<int64_t> equ_indexes(8);
        for (int i=8; i<16; i++) {
            equ_indexes[i-8] = i;
        }
        std::string name_phi = "Stokes_shape.vtk";
        TPZVec<std::string> var_name(2);
        var_name[0]="V";
        var_name[1]="P";
        
        //TPZBuildMultiphysicsMesh::ShowShape(f_mesh_vector,cmesh_m, an, name_phi, equ_indexes);
        an.ShowShape(name_phi, equ_indexes, 1, var_name);
    }
    
    
    
    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,6> Errors;
    ofstream ErroOut("Error_Brinkman.txt", std::ofstream::app);
    an.SetExact(Sol_exact);
    an.SetThreadsForError(3);
    an.PostProcessError(Errors,false);
    
    ErroOut <<"  //  Ordem = "<< pOrder << "  //  Tamanho da malha = "<< n_s[0] <<" x "<< n_s[1] << " x " << n_s[2] << std::endl;
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - V = "<< Errors[1] << std::endl;
    ErroOut <<"Semi-norma H1/Hdiv - V = "<< Errors[2] << std::endl;
    ErroOut <<"Norma L2 - P = "<< Errors[4] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();
    
    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("Brinkman.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    scalnames.Push("Div");
    
    
    int postProcessResolution = 3; //  keep low as possible

    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
}


void HybridBrinkmanTest::Rotate(TPZVec<REAL> &co, TPZVec<REAL> &co_r, bool rotate){
    
    if (rotate==true) {
        //rotação +
        co_r[0] = co[0]*cos(phi_r) - co[1]*sin(phi_r);
        co_r[1] = co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }else{
        
        co_r[0] = co[0]*cos(phi_r) + co[1]*sin(phi_r);
        co_r[1] = - co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }
    
    
}

void HybridBrinkmanTest::InsertLowerDimMaterial(TPZGeoMesh *gmesh){
    
    // Inserir elmentos fmatLambda and fmatLambdaBCs
    TPZManVector<int64_t,3> TopolArc(3);
    
            int64_t nel = gmesh->NElements();
            for (int64_t el = 0; el<nel; el++) {
                TPZGeoEl *gel = gmesh->Element(el);
                if(!gel){
                    continue;
                }
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
                    
                    if(neighbour.Element()->HasSubElement())
                    {
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
                                
                            }else if(neigh_matID==fmatBChole){
                                TPZGeoElBC(gelside, fmatLambdaBC_hole);
                                
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
                        
                        
                        if (f_skellNeighs.NElements()>0 && IsSkellNeighbour(neighbour)) {
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


bool HybridBrinkmanTest::IsSkellNeighbour(TPZGeoElSide neighbour){

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

TPZGeoMesh *HybridBrinkmanTest::CreateGMeshRefPattern(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    
    int nquadnods=n_div[0]+1;
    int nodes = nquadnods*nquadnods;
    gmesh->SetMaxNodeId(nodes-1);
    gmesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,6> TopolQTriangle(3);
    TPZManVector<int64_t,6> TopolQQuadrilateral(4);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<REAL,3> coord(3,0.),coordrot(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    int64_t nodeindex = 0;
    

    //Exterior coordinates
    for(int i = 0; i < nquadnods; i++){
        for(int j = 0; j < nquadnods; j++){
            coord[0] = (j)*h_s[0]/(nquadnods-1);
            coord[1] = (i)*h_s[1]/(nquadnods-1);
            gmesh->NodeVec()[nodeindex].SetCoord(coord);
            gmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }
    
    //Ponto 1
    int64_t elementid = 0;
    
    //Tringular elements with GeoBlend:
    
//    TopolQQuadrilateral[0] =0;
//    TopolQQuadrilateral[1] =n_div[0];
//    TopolQQuadrilateral[2] =n_div[1]*(n_div[0]+1);
//    TopolQQuadrilateral[3] =(n_div[0]+1)*(n_div[1]+1)-1;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQQuadrilateral, fmatID,*gmesh);
//    elementid++;
    
    
    
    //Conectividade dos elementos:
    
    for(int i = 0; i < n_div[1]; i++){
        for(int j = 0; j < n_div[0]; j++){
            TopolQQuadrilateral[0] = (i)*(n_div[1]+1) + (j);
            TopolQQuadrilateral[1] = TopolQQuadrilateral[0]+1;
            TopolQQuadrilateral[2] = TopolQQuadrilateral[1]+n_div[0]+1;
            TopolQQuadrilateral[3] = TopolQQuadrilateral[0]+n_div[0]+1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQQuadrilateral, fmatID,*gmesh);
            elementid++;
        }
    }
    
    gmesh->BuildConnectivity();
    
    
    if(f_Holemesh){
        TPZManVector<REAL,6> FirstCoord(3,0.), h_el(3,0.);
        h_el[0]=h_s[0]/n_div[0];
        h_el[1]=h_s[1]/n_div[1];
        
        int nreft = 1;
        TPZVec<TPZGeoEl *> sons1;
        for (int iref = 0; iref < nreft; iref++) {
            int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *gel = gmesh->ElementVec()[iel];
                if (gel->HasSubElement()) {
                    continue;
                }
                MElementType elType = gel->Type();
                if(gel->MaterialId()==1&&elType==EQuadrilateral){
                    
                    TPZAutoPointer<TPZRefPattern> RefpObst = CreateGMeshObstacle(FirstCoord,h_el);
                    gel->SetRefPattern(RefpObst);
                    gel->Divide(sons1);
                }
            }
        }
        
    }
    
    int nref1 = 0;
    TPZVec<TPZGeoEl *> sons2;
    for (int iref = 0; iref < nref1; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (!gel) {
                continue;
            }
            if (gel->HasSubElement()) {
                continue;
            }
            if (gel->MaterialId()==1){
                gel->Divide(sons2);
            }
        }
    }
    
    // BC and Hole elements
    int nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if (gel->HasSubElement()) {
            continue;
        }
        if (gel->MaterialId()==1){

            
            //Create Hole elements:
            TPZGeoElSide gelside(gel,7);
            TPZGeoElSide neighbour = gelside.Neighbour();
            
            bool near_hole = true;
            while (neighbour != gelside) {
                if(neighbour.Element()->Dimension()==2){
                    near_hole = false;
                }
                neighbour = neighbour.Neighbour();
            }

            if (near_hole) {
                TPZGeoElBC(gelside, fmatBChole);
            }
            

            TPZManVector<int64_t> TopolPlate(4);
            
            int64_t totalnodes = gel->NNodes();
            for (int i=0; i<totalnodes; i++){
                TopolPlate[i] = gel->NodeIndex(i);
            }
            
            //Set BC elements:
            TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
            TPZManVector <REAL,3> nodecoord(3);
            
            TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
            TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;
            TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
            TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;
            
            for (int64_t i = 0; i < totalnodes; i++)
            {
                Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
                Nodefinder[i].GetCoordinates(nodecoord);
               
                if (nodecoord[1] == 0.)
                {
                    sizeOfbottVec++;
                    ncoordzbottVec.Resize(sizeOfbottVec);
                    ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
                }
                if (nodecoord[1] == h_s[1])
                {
                    sizeOftopVec++;
                    ncoordztopVec.Resize(sizeOftopVec);
                    ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == 0.)
                {
                    sizeOfleftVec++;
                    ncoordzleftVec.Resize(sizeOfleftVec);
                    ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == h_s[0])
                {
                    sizeOfrightVec++;
                    ncoordzrightVec.Resize(sizeOfrightVec);
                    ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
                }
            }
            
            if (sizeOfbottVec == 2) {
                int sidesbott = gel->WhichSide(ncoordzbottVec);
                TPZGeoElSide platesidebott(gel, sidesbott);
                TPZGeoElBC(platesidebott,fmatBCbott);
            }
            if (sizeOftopVec == 2) {
                int sidestop = gel->WhichSide(ncoordztopVec);
                TPZGeoElSide platesidetop(gel, sidestop);
                TPZGeoElBC(platesidetop,fmatBCtop);
            }
            if (sizeOfleftVec == 2) {
                int sidesleft = gel->WhichSide(ncoordzleftVec);
                TPZGeoElSide platesideleft(gel, sidesleft);
                TPZGeoElBC(platesideleft,fmatBCleft);
            }
            if (sizeOfrightVec == 2) {
                int sidesright = gel->WhichSide(ncoordzrightVec);
                TPZGeoElSide platesideright(gel, sidesright);
                TPZGeoElBC(platesideright,fmatBCright);
            }
        }
    }
    
    
    TPZManVector<REAL,3> coordFather(3,0.);
    //Remove filiation:
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(gel->HasSubElement()) {
            int nsubel = gel->NSubElements();
            for (int isub=0; isub<nsubel; isub++) {
                //set center coord of father in subelements
                int sub_index = gel->SubElement(isub)->Index();
                TPZManVector<REAL> xcenter(3,0.);
                TPZManVector<REAL,3> xicenter(gel->Dimension(),0.);
                gel->CenterPoint(gel->NSides()-1, xicenter);
                gel->X(xicenter,xcenter);
                f_HoleCoord[sub_index] = xcenter;
                
                int lowest_index = gel->LowestFather()->Index();
                
                //Verifiy if father of father index was created
                if (f_HoleCoord.find(lowest_index)!=f_HoleCoord.end()) {
                    f_HoleCoord[sub_index] = f_HoleCoord[lowest_index];
                }

                
                gel->SetSubElement(isub, 0);
                gel->SetFatherIndex(-1);
            }
            gmesh->DeleteElement(gel);
        }
    }
    
    //Verify filiation removal
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            continue;
        }
        if (gel->HasSubElement()) {
            DebugStop();
        }
    }
    
     //Build Arcs in BC holes
    int nnodes = 0;
    TPZManVector<REAL,3> coord0(3,0.),coord1(3,0.),coordHole(3,0),coordNode(3,0);
    for(int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            continue;
        }
        if(gel->MaterialId()==fmatBChole){

            TPZManVector<int64_t> nodeindices;
            TPZManVector<int64_t,3> TopolArc(3);
            gel->GetNodeIndices(nodeindices);
            TopolArc[0]=nodeindices[0];
            TopolArc[1]=nodeindices[1];


            //Add new node:
            nnodes = gmesh->NNodes();
            gmesh->NodeVec().Resize(nnodes+1);
            nodeindex = nnodes;
            gmesh->SetNodeIdUsed(nodeindex);
            REAL radius = h_s[0]/(4*n_div[0]);
            
            //Find hole central coord
            int nsides = gel->NSides();

            TPZGeoElSide gelside(gel,nsides-1);
            TPZGeoElSide neighbour = gelside.Neighbour();

            while (neighbour != gelside) {
                if(neighbour.Element()->Dimension()==2){
                    int index_neigh = neighbour.Element()->Index();
                    coordHole = f_HoleCoord[index_neigh];
                    f_ArcCentralNode[index_neigh] = nodeindex;
                    break;
                }
                neighbour = neighbour.Neighbour();
            }

            gel->Node(0).GetCoordinates(coord0);
            gel->Node(1).GetCoordinates(coord1);

            REAL rel = (coord0[1]-coordHole[1])/radius;
            REAL theta0H = asin(rel);
            rel = (coord1[1]-coordHole[1])/radius;
            REAL theta1H = asin(rel);
            REAL theta = theta0H - (theta0H-theta1H)/2;

            if (fabs(coord1[0]-coordHole[0])<1.e-6) {
                coord1[0] = coordHole[0];
            }
            if (fabs(coord0[1]-coordHole[1])<1.e-6) {
                coord0[1] = coordHole[1];
            }

            if(coord0[0]<coordHole[0]||coord1[0]<coordHole[0]){
                theta = Pi - theta;
            }

            coordNode = ParametricCircle(radius, theta);
            coordNode[0] +=coordHole[0];
            coordNode[1] +=coordHole[1];
            gmesh->NodeVec()[nodeindex].SetCoord(coordNode);
            gmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);

            TopolArc[2]=nodeindex;

            elementid = gel->Id();

            delete gmesh->ElementVec()[elementid];
            gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*gmesh);
        }
    }

    gmesh->ResetConnectivities();

    //Seting Blend quadrilateral elements:

    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            continue;
        }
        TPZManVector<int64_t> nodeindices;
        MElementType elType = gel->Type();
        if(gel->MaterialId()==1&&elType==EQuadrilateral){
            if (gel->HasSubElement()) {
                DebugStop();
            }
            gel->GetNodeIndices(nodeindices);
            elementid = gel->Id();
            gmesh->ElementVec()[elementid]->SetFatherIndex(-1);
            delete gmesh->ElementVec()[elementid];
            gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad >> (elementid, nodeindices, fmatID,*gmesh);

         }
    }
    

    gmesh->BuildConnectivity();
    
//    InsertLowerDimMaterial(gmesh);
//    SetOriginalMesh(gmesh);
//    //UniformRefine2(1, gmesh, n_div);
//    //InsertLowerDimMaterial(gmesh);
//
//    gmesh->BuildConnectivity();
    
//    TPZCheckGeom check(gmesh);
//    check.CheckUniqueId();
    
    
    int nref2 = f_inter_ref;
    TPZVec<TPZGeoEl *> sons3;
    for (int iref = 0; iref < nref2; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (!gel) {
                continue;
            }
            if (gel->HasSubElement()) {
                continue;
            }
            //if(gel->MaterialId()==1){
                gel->Divide(sons3);
            //}
        }
    }

    gmesh->BuildConnectivity();
    InsertLowerDimMaterial(gmesh);
    SetOriginalMesh(gmesh);
    
    gmesh->BuildConnectivity();
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    
    {
        std::ofstream Dummyfile("GeometricMesh2d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    return gmesh;

    
}


TPZGeoMesh *HybridBrinkmanTest::CreateGMesh(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{
    
    int dimmodel = 2;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    x0[0] = 0., x0[1] = -1.;
    x1[0] = 2., x1[1] = 1.;
    
//    x0[0] = 0., x0[1] = 0.;
//    x1[0] = 4., x1[1] = 2.;
    
    TPZGenGrid grid(n_div,x0,x1);
    
    //grid.SetDistortion(0.2);
    
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
    //UniformRefine(1, gmesh, centerCo, true);
    //UniformRefine2(1, gmesh, n_div);


    InsertLowerDimMaterial(gmesh);
    SetOriginalMesh(gmesh);
    //UniformRefine2(1, gmesh, n_div);
    //InsertLowerDimMaterial(gmesh);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();


        {
            std::ofstream Dummyfile("GeometricMesh2d.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
        }
    
    return gmesh;
    
}

TPZGeoMesh *HybridBrinkmanTest::CreateGMesh3D(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{
    
    int dimmodel = 2;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    //    x0[0] = 0., x0[1] = -1.;
    //    x1[0] = 2., x1[1] = 1.;
    
    x0[0] = 0., x0[1] = 0., x0[2] = 0.;
    x1[0] = 2., x1[1] = 2., x1[2] = 2.;
    
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
                    coord[1] = (i)*h_s[1]/(ny - 1);
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
                    gmesh->CreateGeoElement(ETetraedro, nodeindD1, fmatID, index,0);
                    
                    index++;
                    
                    nodeindD2[0] = nodeindD1[1];
                    nodeindD2[1] = nodeindD1[2];
                    nodeindD2[2] = nodeindD1[1]+nx;
                    nodeindD2[3] = nodeindD1[1] + (1)*nx*ny;
                    gmesh->CreateGeoElement(ETetraedro, nodeindD2, fmatID, index,0);
                    
                    index++;
                    
                    nodeindU1[0] = nodeindD1[0] + (1)*nx*ny;
                    nodeindU1[1] = nodeindU1[0]+1;
                    nodeindU1[2] = nodeindU1[0]+nx;
                    nodeindU1[3] = nodeindD1[2];
                    gmesh->CreateGeoElement(ETetraedro, nodeindU1, fmatID, index,0);
                    
                    index++;
                    
                    nodeindU2[0] = nodeindU1[1];
                    nodeindU2[1] = nodeindU1[2];
                    nodeindU2[2] = nodeindU1[1]+nx;
                    nodeindU2[3] = nodeindD1[2];
                    gmesh->CreateGeoElement(ETetraedro, nodeindU2, fmatID, index,0);
                    
                    index++;
                    
                    // Plano xz
                    
                    nodeindL1[0] = nodeindD1[0];
                    nodeindL1[1] = nodeindD1[2];
                    nodeindL1[2] = nodeindL1[0]+nx*ny;
                    nodeindL1[3] = nodeindU1[1];
                    gmesh->CreateGeoElement(ETetraedro, nodeindL1, fmatID, index,0);
                    
                    index++;
                    
                    
                    nodeindL2[0] = nodeindD2[2];
                    nodeindL2[1] = nodeindD1[1]+nx*ny;
                    nodeindL2[2] = nodeindU2[2];
                    nodeindL2[3] = nodeindU2[3];
                    gmesh->CreateGeoElement(ETetraedro, nodeindL2, fmatID, index,0);
                    
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
    
    //SetAllRefine();
    
    TPZVec<REAL> centerCo(2,0.);
    centerCo[0]=1.;
    centerCo[1]=0.;
    // UniformRefine(1, gmesh, centerCo, true);
    
    //UniformRefine2(1, gmesh, n_div);
    
    InsertLowerDimMaterial(gmesh);
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

TPZGeoMesh *HybridBrinkmanTest::CreateGMeshCurve()
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
    TPZManVector<int64_t,6> TopolQTriangle(3);
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
    
 
//    char buf2[] =
//    "6     5  "
//    "-50       UnifTri3 "
//    " 1.    -1.     0. "
//    " 1.     0.     0. "
//    " 1.    1.     0. "
//    " -1     -1.     0. "
//    " -1.     0.     0. "
//    " -1.     1.     0. "
//    " 3     4     3     0     2     5 "
//    " 2     3     3     0     4 "
//    " 2     3     0     1     4 "
//    " 2     3     4     1     5 "
//    " 2     3     1     2     5 ";
//    " 2     3     3     0     1 "
//    " 2     3     3     1     4 "
//    " 2     3     4     1     2 "
//    " 2     3     4     2     5 ";
    
    char buf[] =
    "6     3  "
    "-50       UnifTri2 "
    " 1.    -1.     0. "
    " 1.     0.     0. "
    " 1.    1.     0. "
    " -1     -1.     0. "
    " -1.     0.     0. "
    " -1.     1.     0. "
    " 3     4     3     0     2     5 "
    " 2     3     3     0     5 "
    " 2     3     2     5     0 ";
    
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpattri = new TPZRefPattern(str);
    
    if(feltype==ETriangle){
        int nreft = 1;
        TPZVec<TPZGeoEl *> sons1;
        for (int iref = 0; iref < nreft; iref++) {
            int nel = geomesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *gel = geomesh->ElementVec()[iel];
                if (gel->HasSubElement()) {
                    continue;
                }
                MElementType elType = gel->Type();
                if(gel->MaterialId()==1&&elType==EQuadrilateral){
                    gel->SetRefPattern(refpattri);
                    gel->Divide(sons1);
                }
            }
        }
    }
    
    int nref = f_inter_ref;
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
    
    
    InsertLowerDimMaterial(geomesh);
    SetOriginalMesh(geomesh);

    TPZCheckGeom check(geomesh);
    check.CheckUniqueId();

    geomesh->BuildConnectivity();
    
//    int nref1 = 4;
//    TPZVec<TPZGeoEl *> sons2;
//    for (int iref = 0; iref < nref1; iref++) {
//        int nel = geomesh->NElements();
//        for (int iel = 0; iel < nel; iel++) {
//            TPZGeoEl *gel = geomesh->ElementVec()[iel];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            if(gel->MaterialId()==1){
//                gel->Divide(sons2);
//            }
//        }
//    }
    
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    return geomesh;
    
}

TPZAutoPointer<TPZRefPattern> HybridBrinkmanTest::CreateGMeshObstacle(TPZManVector<REAL,6> &FirstCoord, TPZManVector<REAL,6> &h_el)
{
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(2);
    
    int nodes = 25;
    REAL radius = h_el[0]/4.;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,6> TopolQQuadrilateral(6);
    TPZManVector<int64_t,8> TopolQuadrilateral(4);
    TPZManVector<int64_t,6> TopolQTriangle(3);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<int64_t,3> TopolArc(3);
    TPZManVector<REAL,3> coord(3,0.);
    TPZManVector<REAL,3> coordCentralnode(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    int nquadnods=3;
    REAL hx = h_el[0];
    REAL hy = h_el[1];
    
    int64_t nodeindex = 0;
    
    //Exterior coordinates
    for(int i = 0; i < nquadnods; i++){
        for(int j = 0; j < nquadnods; j++){
            coord[0] = (j)*hx/(nquadnods - 1);
            coord[1] = (i)*hy/(nquadnods - 1);
            geomesh->NodeVec()[nodeindex].SetCoord(coord);
            geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }
    
    //Circunference coordinates
    for (REAL inode = 0; inode < 16 ; inode++) {
        // i node
        coord = ParametricCircle(radius, inode * M_PI/8.0);
        coord[0] += hx/2.;
        coord[1] += hy/2.;
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    
    int64_t elementid = 0;
    
    TopolQuadrilateral[0] = 0;
    TopolQuadrilateral[1] = 2;
    TopolQuadrilateral[2] = 8;
    TopolQuadrilateral[3] = 6;
    
    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    elementid++;


    TopolQuadrilateral[0] = 9;
    TopolQuadrilateral[1] = 5;
    TopolQuadrilateral[2] = 8;
    TopolQuadrilateral[3] = 11;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son1 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son1->SetFather(father);
    son1->SetFatherIndex(father->Index());
    elementid++;

    TopolQuadrilateral[0] = 11;
    TopolQuadrilateral[1] = 8;
    TopolQuadrilateral[2] = 7;
    TopolQuadrilateral[3] = 13;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son2 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son2->SetFather(father);
    son2->SetFatherIndex(father->Index());
    elementid++;

    TopolQuadrilateral[0] = 13;
    TopolQuadrilateral[1] = 7;
    TopolQuadrilateral[2] = 6;
    TopolQuadrilateral[3] = 15;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son3 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son3->SetFather(father);
    son3->SetFatherIndex(father->Index());
    elementid++;

    TopolQuadrilateral[0] = 15;
    TopolQuadrilateral[1] = 6;
    TopolQuadrilateral[2] = 3;
    TopolQuadrilateral[3] = 17;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son4 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son4->SetFather(father);
    son4->SetFatherIndex(father->Index());
    elementid++;

    TopolQuadrilateral[0] = 17;
    TopolQuadrilateral[1] = 3;
    TopolQuadrilateral[2] = 0;
    TopolQuadrilateral[3] = 19;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son5 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son5->SetFather(father);
    son5->SetFatherIndex(father->Index());
    elementid++;

    TopolQuadrilateral[0] = 19;
    TopolQuadrilateral[1] = 0;
    TopolQuadrilateral[2] = 1;
    TopolQuadrilateral[3] = 21;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son6 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son6->SetFather(father);
    son6->SetFatherIndex(father->Index());
    elementid++;

    TopolQuadrilateral[0] = 21;
    TopolQuadrilateral[1] = 1;
    TopolQuadrilateral[2] = 2;
    TopolQuadrilateral[3] = 23;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son7 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son7->SetFather(father);
    son7->SetFatherIndex(father->Index());
    elementid++;

    TopolQuadrilateral[0] = 23;
    TopolQuadrilateral[1] = 2;
    TopolQuadrilateral[2] = 5;
    TopolQuadrilateral[3] = 9;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son8 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    son8->SetFather(father);
    son8->SetFatherIndex(father->Index());
    elementid++;
    
    //InsertLowerDimMaterial(geomesh);
    //SetOriginalMesh(geomesh);
    
    //Outer bc's
    
//    TopolLine[0] = 1;
//    TopolLine[1] = 0;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCbott,*geomesh);
//    elementid++;
//
//    TopolLine[0] = 1;
//    TopolLine[1] = 2;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCbott,*geomesh);
//    elementid++;
//
//    TopolLine[0] = 2;
//    TopolLine[1] = 5;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCright,*geomesh);
//    elementid++;
//    
//    TopolLine[0] = 5;
//    TopolLine[1] = 8;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCright,*geomesh);
//    elementid++;
//    
//    TopolLine[0] = 8;
//    TopolLine[1] = 7;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCtop,*geomesh);
//    elementid++;
//
//    TopolLine[0] = 7;
//    TopolLine[1] = 6;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCtop,*geomesh);
//    elementid++;
//
//    TopolLine[0] = 6;
//    TopolLine[1] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCleft,*geomesh);
//    elementid++;
//    
//    TopolLine[0] = 3;
//    TopolLine[1] = 0;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCleft,*geomesh);
//    elementid++;
//    
//    //Holes

    
    
//    TopolArc[0] = 9;
//    TopolArc[1] = 11;
//    TopolArc[2] = 10;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
//    
//    TopolArc[0] = 11;
//    TopolArc[1] = 13;
//    TopolArc[2] = 12;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
//
//    TopolArc[0] = 13;
//    TopolArc[1] = 15;
//    TopolArc[2] = 14;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
//    
//    TopolArc[0] = 15;
//    TopolArc[1] = 17;
//    TopolArc[2] = 16;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
//
//    TopolArc[0] = 17;
//    TopolArc[1] = 19;
//    TopolArc[2] = 18;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
//    
//    TopolArc[0] = 19;
//    TopolArc[1] = 21;
//    TopolArc[2] = 20;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
//    
//    TopolArc[0] = 21;
//    TopolArc[1] = 23;
//    TopolArc[2] = 22;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
//    
//    TopolArc[0] = 23;
//    TopolArc[1] = 9;
//    TopolArc[2] = 24;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*geomesh);
//    elementid++;
    
//    char buf[] =
//    "6     3  "
//    "-50       UnifTri2 "
//    " 1.    -1.     0. "
//    " 1.     0.     0. "
//    " 1.    1.     0. "
//    " -1     -1.     0. "
//    " -1.     0.     0. "
//    " -1.     1.     0. "
//    " 3     4     3     0     2     5 "
//    " 2     3     3     0     5 "
//    " 2     3     2     5     0 ";
//
//    std::istringstream str(buf);
//    TPZAutoPointer<TPZRefPattern> refpattri = new TPZRefPattern(str);
//
//    if(feltype==ETriangle){
//        int nreft = 1;
//        TPZVec<TPZGeoEl *> sons1;
//        for (int iref = 0; iref < nreft; iref++) {
//            int nel = geomesh->NElements();
//            for (int iel = 0; iel < nel; iel++) {
//                TPZGeoEl *gel = geomesh->ElementVec()[iel];
//                if (gel->HasSubElement()) {
//                    continue;
//                }
//                MElementType elType = gel->Type();
//                if(gel->MaterialId()==1&&elType==EQuadrilateral){
//                    gel->SetRefPattern(refpattri);
//                    gel->Divide(sons1);
//                }
//            }
//        }
//    }
    
//    int nref = 2;
//    TPZVec<TPZGeoEl *> sons;
//    for (int iref = 0; iref < nref; iref++) {
//        int nel = geomesh->NElements();
//        for (int iel = 0; iel < nel; iel++) {
//            TPZGeoEl *gel = geomesh->ElementVec()[iel];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            gel->Divide(sons);
//        }
//    }
    
//    geomesh->BuildConnectivity();
//    InsertLowerDimMaterial(geomesh);
//    SetOriginalMesh(geomesh);
    
    TPZCheckGeom check(geomesh);
    check.CheckUniqueId();
    
    geomesh->BuildConnectivity();
    
//        int nref1 = 0;
//        TPZVec<TPZGeoEl *> sons2;
//        for (int iref = 0; iref < nref1; iref++) {
//            int nel = geomesh->NElements();
//            for (int iel = 0; iel < nel; iel++) {
//                TPZGeoEl *gel = geomesh->ElementVec()[iel];
//                if (gel->HasSubElement()) {
//                    continue;
//                }
////                if(gel->MaterialId()==1){
//                    gel->Divide(sons2);
////                }
//            }
//        }
    
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    
    geomesh->SetMaxNodeId(geomesh->NNodes()-1);
    geomesh->SetMaxElementId(geomesh->NElements()-1);
    
    TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(*geomesh);
    TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
    if(!Found)
    {
        gRefDBase.InsertRefPattern(refPattern);
    }
    
    return refPattern;
    
}

TPZGeoMesh *HybridBrinkmanTest::CreateGMeshCurveBlendSimple()
{

    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(2);
    
    int nodes = 10;
    REAL radius = 1.0;
    REAL innerradius = radius/2.0;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,6> TopolQTriangle(3);
    TPZManVector<int64_t,6> TopolQQuadrilateral(4);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    
    int64_t nodeindex = 0;

    coord[0]=0.;
    coord[1]=0.;
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;

    coord[0]=1.;
    coord[1]=0.;
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;

    coord[0]=0.;
    coord[1]=1.;
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;

    //Ponto 1
    int64_t elementid = 0;

    //Tringular elements with GeoBlend:

    TopolQTriangle[0] =0;
    TopolQTriangle[1] =1;
    TopolQTriangle[2] =2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > (elementid,TopolQTriangle, fmatID,*geomesh);
    elementid++;

    // outer arcs bc's

    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCbott,*geomesh);
    elementid++;

    TopolLine[0] = 1;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCright,*geomesh);
    elementid++;

    TopolLine[0] = 2;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCleft,*geomesh);
    elementid++;

//    int64_t nodeindex = 0;
//
//    coord[0]=-1;
//    coord[1]=-1;
//    geomesh->NodeVec()[nodeindex].SetCoord(coord);
//    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//    nodeindex++;
//
//    coord[0]=1.;
//    coord[1]=-1;
//    geomesh->NodeVec()[nodeindex].SetCoord(coord);
//    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//    nodeindex++;
//
//    coord[0]=1.;
//    coord[1]=1.;
//    geomesh->NodeVec()[nodeindex].SetCoord(coord);
//    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//    nodeindex++;
//
//    coord[0]=-1.;
//    coord[1]=1.;
//    geomesh->NodeVec()[nodeindex].SetCoord(coord);
//    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//    nodeindex++;
//
//    //Ponto 1
//    int64_t elementid = 0;
//
//    //Tringular elements with GeoBlend:
//
//    TopolQQuadrilateral[0] =0;
//    TopolQQuadrilateral[1] =1;
//    TopolQQuadrilateral[2] =2;
//    TopolQQuadrilateral[3] =3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (elementid,TopolQQuadrilateral, fmatID,*geomesh);
//    elementid++;
//
//    // outer arcs bc's
//
//    TopolLine[0] = 0;
//    TopolLine[1] = 1;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCbott,*geomesh);
//    elementid++;
//
//    TopolLine[0] = 1;
//    TopolLine[1] = 2;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCright,*geomesh);
//    elementid++;
//
//    TopolLine[0] = 2;
//    TopolLine[1] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCtop,*geomesh);
//    elementid++;
//
//    TopolLine[0] = 3;
//    TopolLine[1] = 0;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCleft,*geomesh);
//    elementid++;
    
    
    geomesh->BuildConnectivity();
    
    //For GeoBlend
    //    int nreft = 1;
    //    TPZVec<TPZGeoEl *> sons1;
    //    for (int iref = 0; iref < nreft; iref++) {
    //        int nel = geomesh->NElements();
    //        for (int iel = 0; iel < nel; iel++) {
    //            TPZGeoEl *gel = geomesh->ElementVec()[iel];
    //            if (gel->HasSubElement()) {
    //                continue;
    //            }
    //            if(gel->MaterialId()==-3||gel->MaterialId()==-4){
    //                gel->Divide(sons1);
    //            }
    //        }
    //    }
    
    int nref = f_inter_ref;
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
    
    
    InsertLowerDimMaterial(geomesh);
    SetOriginalMesh(geomesh);
    
    TPZCheckGeom check(geomesh);
    check.CheckUniqueId();
    
    geomesh->BuildConnectivity();
    
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    return geomesh;

    
    
}


TPZGeoMesh *HybridBrinkmanTest::CreateGMeshCurveBlend()
{
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(2);
    
    int nodes = 10;
    REAL radius = 1.0;
    REAL innerradius = radius/2.0;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,6> TopolQQuadrilateral(6);
    TPZManVector<int64_t,8> TopolQuadrilateral(4);
    TPZManVector<int64_t,6> TopolQTriangle(3);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<int64_t,3> TopolArc(3);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    
    int64_t nodeindex = 0;
    
    for (int inode = 0; inode < 5 ; inode++) {
        // i node
        coord = ParametricCircle(radius, inode * M_PI/8.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    for (int inode = 0; inode < 5 ; inode++) {
        // i node
        coord = ParametricCircle(innerradius, inode * M_PI/8.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    //Ponto 1
    int64_t elementid = 0;
    
    //Tringular elements with GeoBlend:
    
    TopolQTriangle[0] =5;
    TopolQTriangle[1] =0;
    TopolQTriangle[2] =7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > (elementid,TopolQTriangle, fmatID,*geomesh);
    elementid++;
    
    TopolQTriangle[0] =0;
    TopolQTriangle[1] =2;
    TopolQTriangle[2] =7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > (elementid,TopolQTriangle, fmatID,*geomesh);
    elementid++;
    
    TopolQTriangle[0] =7;
    TopolQTriangle[1] =2;
    TopolQTriangle[2] =9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > (elementid,TopolQTriangle, fmatID,*geomesh);
    elementid++;
    
    TopolQTriangle[0] =2;
    TopolQTriangle[1] =4;
    TopolQTriangle[2] =9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > (elementid,TopolQTriangle, fmatID,*geomesh);
    elementid++;
    
    // outer arcs bc's
    
    TopolLine[0] = 5;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCbott,*geomesh);
    elementid++;
    
    TopolLine[0] = 4;
    TopolLine[1] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCtop,*geomesh);
    elementid++;
    
    TopolArc[0] = 0;
    TopolArc[1] = 2;
    TopolArc[2] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCright,*geomesh);
    elementid++;
    
    TopolArc[0] = 2;
    TopolArc[1] = 4;
    TopolArc[2] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCright,*geomesh);
    elementid++;

    TopolArc[0] = 5;
    TopolArc[1] = 7;
    TopolArc[2] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCleft,*geomesh);
    elementid++;
    
    TopolArc[0] = 7;
    TopolArc[1] = 9;
    TopolArc[2] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCleft,*geomesh);
    elementid++;
    
    geomesh->BuildConnectivity();
    
    //For GeoBlend
//    int nreft = 1;
//    TPZVec<TPZGeoEl *> sons1;
//    for (int iref = 0; iref < nreft; iref++) {
//        int nel = geomesh->NElements();
//        for (int iel = 0; iel < nel; iel++) {
//            TPZGeoEl *gel = geomesh->ElementVec()[iel];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            if(gel->MaterialId()==-3||gel->MaterialId()==-4){
//                gel->Divide(sons1);
//            }
//        }
//    }
    
    int nref = f_inter_ref;
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
    
    
    InsertLowerDimMaterial(geomesh);
    SetOriginalMesh(geomesh);
    
    TPZCheckGeom check(geomesh);
    check.CheckUniqueId();
    
    geomesh->BuildConnectivity();
    
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    return geomesh;
    
}


TPZManVector<REAL,3>  HybridBrinkmanTest::ParametricCircle(REAL radius,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta);
    xcoor[1] = radius * sin(theta);
    xcoor[2] = 0.0 ;
    return xcoor;
}

void HybridBrinkmanTest::TetrahedralMeshCubo(TPZVec<int> &n_s){
    
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    GenerateNodes(gmesh,nelem);
//
//
//
//    return gmesh;
    
}


void HybridBrinkmanTest::UniformRefine4(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction)
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


void HybridBrinkmanTest::UniformRefine3(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div)
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


void HybridBrinkmanTest::UniformRefine2(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div)
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



void HybridBrinkmanTest::UniformRefine(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction)
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

TPZCompEl *HybridBrinkmanTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}

void HybridBrinkmanTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
//    sol.resize(4);
//    dsol.Resize(3,3);
//
//    STATE xv = x[0];
//    STATE yv = x[1];
//    STATE theta = atan(yv/xv);
//    STATE r=sqrt(xv*xv+yv*yv);
//
//    STATE v_x =  -r*sin(theta);
//    STATE v_y =  r*cos(theta);
//    STATE v_norm =  sqrt(v_x*v_x+v_y*v_y);
//    STATE p =   3.;
//
//    //    STATE v_x =  xv;
//    //    STATE v_y =  0;
//    //    STATE p =   0.;
//
//
//    sol[0] = v_x; // x direction
//    sol[1] = v_y; // y direction
//    sol[2] = 0.;
//    sol[3] = p; //
//
//    dsol(0,1)= -(sqrt(1+(yv*yv)/(xv*xv))*xv)/r;
//    dsol(1,0)= (sqrt(1+(yv*yv)/(xv*xv))*xv)/r;
    
    
//    sol[0] = 1.; // x direction
//    sol[1] = 0.; // y direction
//    sol[3] = 0.;
//    dsol(0,0)= 0.;
//    dsol(1,1)= 0.;
    
//        dsol.Resize(3,3);
//        sol.Resize(4);
//
//        REAL x1 = x[0];
//        REAL x2 = x[1];
//
//        TPZVec<REAL> v_Dirichlet(3,0.);
//
////        v_Dirichlet[0] = -0.1*x2*x2+0.2*x2;
//        v_Dirichlet[0] = -1.+x2;
////        v_Dirichlet[0] = 1.;
//        v_Dirichlet[1] = 0.;
//        v_Dirichlet[2] = 0.;
////        STATE pressure = 1.-0.2*x1;
//        STATE pressure = 0.;
//
//        sol[0]=v_Dirichlet[0];
//        sol[1]=v_Dirichlet[1];
//        sol[2]=v_Dirichlet[2];
//        sol[3]=pressure;
//
//        // vx direction
//        dsol(0,0)= 0.;
////        dsol(0,1)= 0.2-0.2*x2;
//        dsol(0,1)= 1.;
//        //dsol(0,1)= 0.;
//        dsol(0,2)= 0.;
//
//        // vy direction
//        dsol(1,0)= 0.;
//        dsol(1,1)= 0.;
//        dsol(1,2)= 0.;
//
//        // vz direction
//        dsol(2,0)= 0.;
//        dsol(2,1)= 0.;
//        dsol(2,2)= 0.;
    
    // General form : : Artigo Botti, Di Pietro, Droniou
    
    //    dsol.Resize(3,3);
    //    sol.Resize(3);
    //
    //    REAL x1 = x[0];
    //    REAL x2 = x[1];
    //
    //    REAL m_v= 1., m_u= 1.0;
    //
    //    REAL Cf=m_v/m_u;
    //
    //    STATE v_1 = -exp(-Cf)*sin(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*sin(x1)*sin(x2);
    //    STATE v_2 = -exp(-Cf)*cos(x1)*cos(x2)-(1./m_v)*(1.-exp(-Cf))*cos(x1)*cos(x2);
    //    STATE pressure= cos(x1)*sin(x2);
    //
    //    sol[0]=v_1;
    //    sol[1]=v_2;
    //    sol[2]=pressure;
    //
    //    // vx direction
    //    dsol(0,0)= -exp(-Cf)*cos(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
    //    dsol(0,1)= exp(-Cf)*cos(x2)*sin(x1)+(1./m_v)*(1.-exp(-Cf))*cos(x2)*sin(x1);
    //
    //    // vy direction
    //    dsol(1,0)= -exp(-Cf)*cos(x2)*sin(x1)+(1./m_v)*(1.-exp(-Cf))*cos(x2)*sin(x1);
    //    dsol(1,1)= exp(-Cf)*cos(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
    //

    
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou
    
    //    dsol.Resize(3,3);
    //    sol.Resize(3);
    //
    //    REAL x1 = x[0];
    //    REAL x2 = x[1];
    //
    //    REAL e = exp(1.);
    //
    //    STATE v_1 = (1.-2./e)*sin(x1)*sin(x2);
    //    STATE v_2 = -1.*cos(x1)*cos(x2);
    //    STATE pressure= cos(x1)*sin(x2);
    //
    //    sol[0]=v_1;
    //    sol[1]=v_2;
    //    sol[2]=pressure;
    //
    //    // vx direction
    //    dsol(0,0)= (1.-2./e)*cos(x1)*sin(x2);
    //    dsol(0,1)= cos(x2)*sin(x1);
    //
    //    // vy direction
    //    dsol(1,0)= (1.-2./e)*cos(x2)*sin(x1);
    //    dsol(1,1)= cos(x1)*sin(x2);
    //

    
    // Stokes : : Artigo Botti, Di Pietro, Droniou
    
    dsol.Resize(3,3);
    sol.Resize(4);


    //Applying rotation:
    TPZVec<REAL> x_in = x;
    TPZVec<REAL> x_rot(3,0.);

    f_InvT.Apply(x_in,x_rot);
    x[0] = x_rot[0];
    x[1] = x_rot[1];

    REAL x1 = x[0];
    REAL x2 = x[1];

    REAL e = exp(1.);

    TPZVec<REAL> v_Dirichlet(3,0.), vbc_rot(3,0.);

    v_Dirichlet[0] = -1.*sin(x1)*sin(x2);
    v_Dirichlet[1] = -1.*cos(x1)*cos(x2);
    STATE pressure= cos(x1)*sin(x2);

    f_T.Apply(v_Dirichlet, vbc_rot);
    v_Dirichlet = vbc_rot;

    sol[0]=v_Dirichlet[0];
    sol[1]=v_Dirichlet[1];
    sol[2]=v_Dirichlet[2];
    sol[3]=pressure;


    // GradU * Rt
    TPZFMatrix<STATE> GradU(3,3,0.), GradURt(3,3,0.), RGradURt(3,3,0.);

    // vx direction
    GradU(0,0)= -1.*cos(x1)*sin(x2);
    GradU(0,1)= cos(x2)*sin(x1);

    // vy direction
    GradU(1,0)= -1.*cos(x2)*sin(x1);
    GradU(1,1)= cos(x1)*sin(x2);

    TPZFMatrix<STATE> R = f_T.Mult();
    TPZFMatrix<STATE> Rt(3,3,0.);
    R.Transpose(&Rt);

//    GradU.Print("GradU = ");
//    R.Print("R = ");
//    Rt.Print("Rt = ");

    GradU.Multiply(Rt,GradURt);
//    GradURt.Print("GradURt = ");

    R.Multiply(GradURt,RGradURt);
//    RGradURt.Print("RGradURt = ");

    // vx direction
    dsol(0,0)= RGradURt(0,0);
    dsol(0,1)= RGradURt(0,1);
    dsol(0,2)= RGradURt(0,2);

    // vy direction
    dsol(1,0)= RGradURt(1,0);
    dsol(1,1)= RGradURt(1,1);
    dsol(1,2)= RGradURt(1,2);

    // vz direction
    dsol(2,0)= RGradURt(2,0);
    dsol(2,1)= RGradURt(2,1);
    dsol(2,2)= RGradURt(2,2);
    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
    //        dsol.Resize(3,3);
    //        sol.Resize(3);
    //
    //        REAL x1 = x[0];
    //        REAL x2 = x[1];
    //
    //        STATE v_1 = sin(x1)*sin(x2);
    //        STATE v_2 = -1.*cos(x1)*cos(x2);
    //        STATE pressure= cos(x1)*sin(x2);
    //
    //        sol[0]=v_1;
    //        sol[1]=v_2;
    //        sol[2]=pressure;
    //
    //        // vx direction
    //        dsol(0,0)= cos(x1)*sin(x2);
    //        dsol(0,1)= cos(x2)*sin(x1);
    //
    //        // vy direction
    //        dsol(1,0)= cos(x2)*sin(x1);
    //        dsol(1,1)= cos(x1)*sin(x2);
    

    
    // Stokes 3D : Artigo Botti, Di Pietro, Droniou
    
//        dsol.Resize(3,3);
//        sol.Resize(4);
//
//        //Applying rotation:
//        TPZVec<REAL> x_in = x;
//        TPZVec<REAL> x_rot(3,0.);
//
//        f_InvT.Apply(x_in,x_rot);
//        x[0] = x_rot[0];
//        x[1] = x_rot[1];
//
//        REAL x1 = x[0];
//        REAL x2 = x[1];
//        REAL x3 = x[2];
//
//        TPZVec<REAL> v_Dirichlet(3,0.), vbc_rot(3,0.);
//
//        v_Dirichlet[0] = cos(x1)*cos(x3) -1.*sin(x1)*sin(x2);
//        v_Dirichlet[1] = -1.*cos(x1)*cos(x2);
//        v_Dirichlet[2] = sin(x1)*sin(x3);
//        STATE pressure= cos(x1)*sin(x2)*cos(x3);
//
//        f_T.Apply(v_Dirichlet, vbc_rot);
//        v_Dirichlet = vbc_rot;
//
//        sol[0]=v_Dirichlet[0];
//        sol[1]=v_Dirichlet[1];
//        sol[2]=v_Dirichlet[2];
//        sol[3]=pressure;
//
//
//        // GradU * Rt
//        TPZFMatrix<STATE> GradU(3,3,0.), GradURt(3,3,0.), RGradURt(3,3,0.);
//
//        // vx direction
//        GradU(0,0)= -cos(x3)*sin(x1)-1.*cos(x1)*sin(x2);
//        GradU(0,1)= cos(x2)*sin(x1);
//        GradU(0,2)= cos(x1)*sin(x3);
//
//        // vy direction
//        GradU(1,0)= -1.*cos(x2)*sin(x1);
//        GradU(1,1)= cos(x1)*sin(x2);
//        GradU(1,2)= 0.;
//
//        // vz direction
//        GradU(2,0)= -1.*cos(x1)*sin(x3);
//        GradU(2,1)= 0.;
//        GradU(2,2)= cos(x3)*sin(x1);
//
//
//        TPZFMatrix<STATE> R = f_T.Mult();
//        TPZFMatrix<STATE> Rt(3,3,0.);
//        R.Transpose(&Rt);
//
//    //    GradU.Print("GradU = ");
//    //    R.Print("R = ");
//    //    Rt.Print("Rt = ");
//
//        GradU.Multiply(Rt,GradURt);
//    //    GradURt.Print("GradURt = ");
//
//        R.Multiply(GradURt,RGradURt);
//    //    RGradURt.Print("RGradURt = ");
//
//        // vx direction
//        dsol(0,0)= RGradURt(0,0);
//        dsol(0,1)= RGradURt(0,1);
//        dsol(0,2)= RGradURt(0,2);
//
//        // vy direction
//        dsol(1,0)= RGradURt(1,0);
//        dsol(1,1)= RGradURt(1,1);
//        dsol(1,2)= RGradURt(1,2);
//
//        // vz direction
//        dsol(2,0)= RGradURt(2,0);
//        dsol(2,1)= RGradURt(2,1);
//        dsol(2,2)= RGradURt(2,2);
    
    
}

void HybridBrinkmanTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    //Applying rotation:
    TPZVec<REAL> x_in = x;
    TPZVec<REAL> x_rot(3,0.);
    
    f_InvT.Apply(x_in,x_rot);
    x[0] = x_rot[0];
    x[1] = x_rot[1];
    
    f.resize(3);
    REAL x1 = x[0];
    REAL x2 = x[1];
    REAL x3 = x[2];
    
    f[0] =0.;
    f[1] =0.;
    f[2] =0.;
    
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
    
    // Stokes : : Artigo Botti, Di Pietro, Droniou
    
    
    f_s[0] = -3.*sin(x1)*sin(x2);
    f_s[1] = -1.*cos(x1)*cos(x2);

    f_T.Apply(f_s, f_rot);
    f_s = f_rot;


    f[0] = f_s[0]; // x direction
    f[1] = f_s[1]; // y direction
    f[2] = f_s[2];
    
    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
    //        f_1 = 0.;
    //        f_2 = 0.;
    //
    //        f[0] = f_1; // x direction
    //        f[1] = f_2; // y direction
    //        f[2] = 2.*cos(x1)*sin(x2);
    
    
    // Stokes 3D : Artigo Botti, Di Pietro, Droniou
    
    
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

void HybridBrinkmanTest::ChangeExternalOrderConnects(TPZCompMesh *mesh, int addToOrder){
    
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


TPZCompMesh *HybridBrinkmanTest::CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    //BDM test papapaapapap
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    
    // 1 - Material volumétrico 2D
    TPZVecL2 *material = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material);
    
    if (Space==1) {
        cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
        //cmesh->ApproxSpace().CreateDisconnectedElements(true); //HDIV-Full:
        //Dimensões do material (para HDiv):
        //TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
        //material->SetMaterial(xkin, xcin, xfin);
        
    }else{
        DebugStop();
    }
    
    // 1 - Condições de contorno
    TPZFMatrix<STATE> val1(1,1,0.), val2(2,1,0.);
    {
        TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond0);
        
        TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond1);
        
        TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
        
        TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond3);
    }
    
    if (f_Holemesh){
        TPZMaterial * BCondHole = material->CreateBC(material, fmatBChole, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCondHole);
    }
    
    if (f_3Dmesh) {
        TPZMaterial * BCond4 = material->CreateBC(material, fmatBCtop_z, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond4);
        
        TPZMaterial * BCond5 = material->CreateBC(material, fmatBCbott_z, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond5);
    }
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
    
    
}


TPZCompMesh *HybridBrinkmanTest::CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    // 1 - Material volumétrico 2D
    TPZVecL2 *material = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material);
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xbin(1,1,0.), xfin(1,1,0.);
    //    material->SetMaterial(xkin, xcin, xfin);
    
    // 2 - Material para tração tangencial 1D
    
    TPZVecL2 *matLambda = new TPZVecL2(fmatLambda);
    //matLambda->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambda);
    
    // 3 - Material para tração tangencial nos contornos
    
    TPZVecL2 *matLambdaBC_bott = new TPZVecL2(fmatLambdaBC_bott);
    //matLambdaBC_bott->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_bott);
    
    TPZVecL2 *matLambdaBC_top = new TPZVecL2(fmatLambdaBC_top);
    //matLambdaBC_top->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_top);
    
    TPZVecL2 *matLambdaBC_left = new TPZVecL2(fmatLambdaBC_left);
    //matLambdaBC_left->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_left);
    
    TPZVecL2 *matLambdaBC_right = new TPZVecL2(fmatLambdaBC_right);
    //matLambdaBC_right->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_right);
    
    if (f_Holemesh){
        TPZVecL2 *matLambdaBC_hole = new TPZVecL2(fmatLambdaBC_hole);
        cmesh->InsertMaterialObject(matLambdaBC_hole);
    }

    
    if (f_3Dmesh) {

        TPZVecL2 *matLambdaBC_bott_z = new TPZVecL2(fmatLambdaBC_bott_z);
        //matLambdaBC_bott_z->SetMaterial(xkin, xcin, xbin, xfin);
        cmesh->InsertMaterialObject(matLambdaBC_bott_z);
        
        TPZVecL2 *matLambdaBC_top_z = new TPZVecL2(fmatLambdaBC_top_z);
        //matLambdaBC_top_z->SetMaterial(xkin, xcin, xbin, xfin);
        cmesh->InsertMaterialObject(matLambdaBC_top_z);

    }
    
//    //    Ponto de pressao:
//    //
//    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
//    ////
//    TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressao
//    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    //    TPZMaterial * BCPoint2 = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint2); //Insere material na malha
    
    //
    //    //    Ponto de pressao2:
    //    //
    //    TPZFMatrix<STATE> val5(1,1,0.), val6(1,1,0.);
    //    ////
    //    TPZMaterial * BCPoint2 = material->CreateBC(material, matPoint2, pointtype, val5, val6); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint2); //Insere material na malha
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(fmatID);
    
    // materialids.insert(fpointtype);
    cmesh->AutoBuild(materialids);
    
    gmesh->ResetReference();
    //  cmesh->LoadReferences();
    
    materialids.clear();
    materialids.insert(fmatLambda);
    materialids.insert(fmatLambdaBC_bott);
    materialids.insert(fmatLambdaBC_top);
    materialids.insert(fmatLambdaBC_left);
    materialids.insert(fmatLambdaBC_right);
    if (f_Holemesh){
        materialids.insert(fmatLambdaBC_hole);
    }
    
    if (f_3Dmesh) {
        materialids.insert(fmatLambdaBC_bott_z);
        materialids.insert(fmatLambdaBC_top_z);
    }
        
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(pOrder-1);
    cmesh->SetDimModel(fdim-1);
    cmesh->AutoBuild(materialids);
    
    cmesh->LoadReferences();
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *HybridBrinkmanTest::CMesh_pM(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    // 1 - Material volumétrico 2D
    TPZVecL2 *material_pM = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material_pM);
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }

    cmesh->AutoBuild();
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *HybridBrinkmanTest::CMesh_gM(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    // 1 - Material volumétrico 2D
    TPZVecL2 *material_pM = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material_pM);
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    cmesh->AutoBuild();
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


TPZCompMesh *HybridBrinkmanTest::CMesh_pM_0(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    // 1 - Material volumétrico 2D
    TPZVecL2 *material_pM = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material_pM);
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    cmesh->AutoBuild();
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *HybridBrinkmanTest::CMesh_gM_0(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    // 1 - Material volumétrico 2D
    TPZVecL2 *material_pM = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material_pM);
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    cmesh->AutoBuild();
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


TPZMultiphysicsCompMesh *HybridBrinkmanTest::CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco)
{
    
    //Criando malha computacional:
    
    TPZMultiphysicsCompMesh * cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    // Criando material:
    
    // 1 - Material volumétrico 2D
    TPZMHMBrinkmanMaterial *material = new TPZMHMBrinkmanMaterial(fmatID,fdim,Space,visco,0,0);
    int fexact_order = 5;
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source, fexact_order);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact,fexact_order);
    ((TPZDummyFunction<STATE>*)fp.operator->())->SetPolynomialOrder(fexact_order);
    ((TPZDummyFunction<STATE>*)solp.operator->())->SetPolynomialOrder(fexact_order);
    material->SetForcingFunction(fp); //Caso simples sem termo fonte
    material->SetForcingFunctionExact(solp);
    
    cmesh->InsertMaterialObject(material);
    
    // 1 - Condições de contorno:
    // Condições de contorno - Impõe v fortemente
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    val2(1,0) = 0.0;
    TPZBndCond * BC_bott = material->CreateBC(material, fmatBCbott, fneumann, val1, val2);
    BC_bott->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_bott);

    val2(1,0) = 0.0; // vx -> 0
    TPZBndCond * BC_top = material->CreateBC(material, fmatBCtop, fneumann, val1, val2);
    BC_top->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_top);

    //val2(0,0) = 1.0;
    TPZBndCond * BC_left = material->CreateBC(material, fmatBCleft, fneumann, val1, val2);
    BC_left->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_left);

    val2(0,0) = 0.0;
    TPZBndCond * BC_right = material->CreateBC(material, fmatBCright, fneumann, val1, val2);
    BC_right->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_right);

    if (f_Holemesh){
        val2(0,0) = 0.0;
        TPZBndCond * BC_hole = material->CreateBC(material, fmatBChole, fdirichlet, val1, val2);
        BC_hole->SetBCForcingFunction(0, solp);
        cmesh->InsertMaterialObject(BC_hole);
    }

    
    if (f_3Dmesh) {
        TPZBndCond * BC_bott_z = material->CreateBC(material, fmatBCbott_z, fdirichlet, val1, val2);
        BC_bott_z->SetBCForcingFunction(0, solp);
        cmesh->InsertMaterialObject(BC_bott_z);
        
        TPZBndCond * BC_top_z = material->CreateBC(material, fmatBCtop_z, fdirichlet, val1, val2);
        BC_top_z->SetBCForcingFunction(0, solp);
        cmesh->InsertMaterialObject(BC_top_z);
    }
    
    // 2.1 - Material para tração tangencial 1D (Interior)
    TPZNullMaterial *matLambda = new TPZNullMaterial(fmatLambda);
    matLambda->SetDimension(fdim-1);
    matLambda->SetNStateVariables(1);
    cmesh->InsertMaterialObject(matLambda);
    
    // 2.2 - Material for interfaces (Interior)
    TPZMHMBrinkmanMaterial *matInterfaceLeft = new TPZMHMBrinkmanMaterial(fmatInterfaceLeft,fdim,Space,visco,0,0);
    matInterfaceLeft->SetMultiplier(1.);
    cmesh->InsertMaterialObject(matInterfaceLeft);
    
    TPZMHMBrinkmanMaterial *matInterfaceRight = new TPZMHMBrinkmanMaterial(fmatInterfaceRight,fdim,Space,visco,0,0);
    matInterfaceRight->SetMultiplier(-1.);
    cmesh->InsertMaterialObject(matInterfaceRight);
    
    
    // 3.1 - Material para tração tangencial 1D nos contornos
    TPZBndCond *matLambdaBC_bott = material->CreateBC(material, fmatLambdaBC_bott, fneumann, val1, val2);
    matLambdaBC_bott->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(matLambdaBC_bott);
    
    TPZBndCond *matLambdaBC_top = material->CreateBC(material, fmatLambdaBC_top, fneumann, val1, val2);
    matLambdaBC_top->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(matLambdaBC_top);
  
    TPZBndCond *matLambdaBC_left = material->CreateBC(material, fmatLambdaBC_left, fneumann, val1, val2);
    matLambdaBC_left->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(matLambdaBC_left);
    
    TPZBndCond *matLambdaBC_right = material->CreateBC(material, fmatLambdaBC_right, fneumann, val1, val2);
    matLambdaBC_right->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(matLambdaBC_right);

    if (f_Holemesh){
        val2(0,0) = 0.0;
        TPZBndCond *matLambdaBC_hole = material->CreateBC(material, fmatLambdaBC_hole, fdirichlet, val1, val2);
        //    matLambdaBC_hole->SetBCForcingFunction(0, solp);
        cmesh->InsertMaterialObject(matLambdaBC_hole);
    }
    
    if (f_3Dmesh) {
        TPZBndCond *matLambdaBC_bott_z = material->CreateBC(material, fmatLambdaBC_bott_z, fneumann, val1, val2);
        matLambdaBC_bott_z->SetBCForcingFunction(0, solp);
        cmesh->InsertMaterialObject(matLambdaBC_bott_z);
        
        TPZBndCond *matLambdaBC_top_z = material->CreateBC(material, fmatLambdaBC_top_z, fneumann, val1, val2);
        matLambdaBC_top_z->SetBCForcingFunction(0, solp);
        cmesh->InsertMaterialObject(matLambdaBC_top_z);
    }
    
    //Ponto
//    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
//    val4(0,0)=1.;
//    TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressão
//    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    TPZManVector<int,5> active_approx_spaces(2,1);
    
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,f_mesh_vector);
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

void HybridBrinkmanTest::InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m){
    
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

void HybridBrinkmanTest::ComputeSkelNeighbours(){
    
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


void HybridBrinkmanTest::GroupAndCondense(TPZMultiphysicsCompMesh *cmesh_m){
   
    
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



//root -l MVA_code_2017.cpp\(\"A\"\)

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

//#include "TMVA/CrossValidation.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

#include "T3M_common.h"

using namespace TMVA;


void MVA_code_2017(TString categ){
    //Check on input argument
    if(!categ.Contains("A") && !categ.Contains("B") && !categ.Contains("C")){
        cout << "Please choose between any combination of 'A', 'B' and 'C'" << endl;
        return;
    }

    // Output file
    TFile *fout = new TFile("TMVA_"+TMVA_outputpath+categ+".root", "RECREATE");

    TString cat_name[] = {"A", "B", "C"};
    std::vector<TTree*> sigTree, bkgTree;

    // Get the signal and background trees from TFile source(s);
    //signal
    TString treeName = "FinalTree";
    //Ds
    TFile *f_sig_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f_sig_ds || !f_sig_ds->IsOpen()) {
         f_sig_ds = new TFile(inputpath_Ds);
    }
    sigTree.push_back( (TTree*)f_sig_ds->Get(treeName));
    //B0
    TFile *f_sig_b0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
    if (!f_sig_b0 || !f_sig_b0->IsOpen()) {
        f_sig_b0 = new TFile(inputpath_B0);
    }
    sigTree.push_back( (TTree*)f_sig_b0->Get(treeName));
    //Bp
    TFile *f_sig_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    if (!f_sig_bp || !f_sig_bp->IsOpen()) {
       f_sig_bp = new TFile(inputpath_Bp);
    }
    sigTree.push_back( (TTree*)f_sig_bp->Get(treeName));

    //background
    //Loop on run
    int n_bkg = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    for(int j = 0; j<n_bkg; j++){
        TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[j]);
        if (!f_bkg || !f_bkg->IsOpen()) {
            f_bkg = new TFile(inputpath_datarun[j]);
        }
        bkgTree.push_back((TTree*)f_bkg->Get(treeName));
    }

    //add TTreeFriend - muon quality variables
    TString friendTreeName[] = {"TreeMu1=TreeMu1", "TreeMu2=TreeMu2", "TreeMu3=TreeMu3"};
    for(int i=0; i<3; i++){
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_Ds<<endl;
        sigTree.at(0)->AddFriend(friendTreeName[i], inputpath_Ds);
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_B0<<endl;
        sigTree.at(1)->AddFriend(friendTreeName[i], inputpath_B0);
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_Bp<<endl;
        sigTree.at(2)->AddFriend(friendTreeName[i], inputpath_Bp);
    }
    for(int i=0; i<3; i++){
        for(int j=0; j<n_bkg; j++){
            cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_datarun[j]<<endl;
            bkgTree.at(j)->AddFriend(friendTreeName[i], inputpath_datarun[j]);
        }
    }

    //add TTreeFriend - MuonID evaluation
    TString muIdTreeName[] = {"MuonIDeval_Mu1=MuonIDeval_Mu1", "MuonIDeval_Mu2=MuonIDeval_Mu2", "MuonIDeval_Mu3=MuonIDeval_Mu3"};
    TString inputpath_Ds_muId = inputpath_Ds.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_B0_muId = inputpath_B0.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_Bp_muId = inputpath_Bp.ReplaceAll(".root", "_MuonID.root");
    for(int i=0; i<3; i++){
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_Ds_muId<<endl;
        sigTree.at(0)->AddFriend(muIdTreeName[i], inputpath_Ds_muId);
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_B0_muId<<endl;
        sigTree.at(1)->AddFriend(muIdTreeName[i], inputpath_B0_muId);
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_Bp_muId<<endl;
        sigTree.at(2)->AddFriend(muIdTreeName[i], inputpath_Bp_muId);
    }
    for(int j=0; j<n_bkg; j++){
        TString inputpath_data_muId = inputpath_datarun[j].ReplaceAll(".root", "_MuonID.root");
        for(int i=0; i<3; i++){
            cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_data_muId<<endl;
            bkgTree.at(j)->AddFriend(muIdTreeName[i], inputpath_data_muId);
        }
    }



    // Set the event weights per tree
    Double_t sigWeight1  = wNormDs/wNormDs; //1.0 Ds
    Double_t sigWeight2  = wNormB0/wNormDs; //B0
    Double_t sigWeight3  = wNormBp/wNormDs; //Bp
    cout<<"Ds sigWeight1: "<<sigWeight1<<" - B0 sigWeight2: "<<sigWeight2<<" - Bp sigWeight3: "<<sigWeight3<<endl;

    Double_t bkgWeight1 = 1.0;
    Double_t bkgWeight2 = 1.0;
    Double_t bkgWeight3 = 1.0;
    Double_t bkgWeight4 = 1.0;
    Double_t bkgWeight5 = 1.0;
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
    
    dataloader->AddSignalTree(sigTree.at(0), sigWeight1);
    dataloader->AddSignalTree(sigTree.at(1), sigWeight2);
    dataloader->AddSignalTree(sigTree.at(2), sigWeight3);
    dataloader->AddBackgroundTree(bkgTree.at(0), bkgWeight1);
    dataloader->AddBackgroundTree(bkgTree.at(1), bkgWeight2);
    dataloader->AddBackgroundTree(bkgTree.at(2), bkgWeight3);
    dataloader->AddBackgroundTree(bkgTree.at(3), bkgWeight4);
    dataloader->AddBackgroundTree(bkgTree.at(4), bkgWeight5);

    std::vector<TString> var_train_name;
    std::vector<TString> var_train_def;
    std::vector<TString> var_spec_name;
    std::vector<TString> var_spec_def;
    TString BDTinVar;
    if(categ.Contains("A")) BDTinVar = BDTinVar_A;
    if(categ.Contains("B")) BDTinVar = BDTinVar_B;
    if(categ.Contains("C")) BDTinVar = BDTinVar_C;

    readVarName(var_spec_name, var_spec_def, BDTspecVar);
    readVarName(var_train_name, var_train_def, BDTinVar);

    // Spectators declaration
    cout<<"Declaration of spectator variables - category "<<categ<<" from file:"<<BDTspecVar<<endl;
    for(int k = 0; k<var_spec_name.size(); k++){
        if(var_spec_name.at(k)=="" || var_spec_def.at(k)=="") continue;
        cout<<k<<" - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
        dataloader->AddSpectator(var_spec_def.at(k), var_spec_name.at(k), "", 'D');
    }

    // Variables declaration
    cout<<"Declaration of variables for training - category "<<categ<<" from file:"<<BDTinVar<<endl;
    for(int k = 0; k<var_train_name.size(); k++){
        if(var_train_name.at(k)=="" || var_train_def.at(k)=="") continue;
        cout<<k<<" - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
        dataloader->AddVariable(var_train_def.at(k), var_train_name.at(k), "", 'D');
    }


    //dataloader->SetSignalWeightExpression( "puFactor" );

    TCut cutS = "((tripletMass<2.0 && tripletMass>1.62))"; //Signal -> MC full range 
    TCut cutB = "(  ((tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)) )"; //Background -> data sidebands
    TCut preselCut = "";

    TCut reso_A = "tripletMassReso < 0.007";
    TCut reso_B = "tripletMassReso >= 0.007 && tripletMassReso <= 0.0105";
    TCut reso_C = "tripletMassReso > 0.0105";
    TCut reso_cat = "tripletMassReso < 0"; //always false    

    if(categ.Contains("A")) reso_cat = reso_cat || reso_A;
    if(categ.Contains("B")) reso_cat = reso_cat || reso_B;
    if(categ.Contains("C")) reso_cat = reso_cat || reso_C;

    TString prepareTrainTestOptions = ":SplitMode=Random"
                                      ":NormMode=NumEvents"
                                      ":nTrain_Signal=0"
                                      ":nTrain_Background=0";
    //no phi and omega vetos
    if(categ=="A") {
             prepareTrainTestOptions = prepareTrainTestOptions+
                                       ":nTest_Signal=9080" //20% total (45.4k) A
                                       ":nTest_Background=3900" //20% total (19.4k) A
                                       ":!V";
                        }   
    if(categ=="B") {
             prepareTrainTestOptions = prepareTrainTestOptions+
                                       ":nTest_Signal=16500" //20% total (82.5k) B
                                       ":nTest_Background=7700" //20% total (38.4k) B
                                       ":!V";
                        }   
    if(categ=="C") {
             prepareTrainTestOptions = prepareTrainTestOptions+
                                       ":nTest_Signal=5850" //20% total (21.9k) C
                                       ":nTest_Background=2600" //20% total (13k) C
                                       ":!V";
                        }   
//    //phi veto applied before BDT
//    if(categ=="A") {
//             prepareTrainTestOptions = prepareTrainTestOptions+
//                                       ":nTest_Signal=6724" //20% total (33622) A
//                                       ":nTest_Background=3944" //20% total (19722) A
//                                       ":!V";
//                        }   
//    if(categ=="B") {
//             prepareTrainTestOptions = prepareTrainTestOptions+
//                                       ":nTest_Signal=12759" //20% total (31897*2) A
//                                       ":nTest_Background=8305" //20% total (20763*2) A
//                                       ":!V";
//                        }   
//    if(categ=="C") {
//             prepareTrainTestOptions = prepareTrainTestOptions+
//                                       ":nTest_Signal=4712" //20% total (23560) A
//                                       ":nTest_Background=2890" //20% total (14448) A
//                                       ":!V";
//                        }   

    TString veto_phi = "";
    if (categ=="A") veto_phi = "!(dimu12>0.994 && dimu12<1.044) && !(dimu13>0.994 && dimu13<1.044) &&!(dimu23>0.994 && dimu23<1.044)";
    if (categ=="B") veto_phi = "!(dimu12>0.985 && dimu12<1.053) && !(dimu13>0.985 && dimu13<1.053) &&!(dimu23>0.985 && dimu23<1.053)";
    if (categ=="C") veto_phi = "!(dimu12>0.974 && dimu12<1.064) && !(dimu13>0.974 && dimu13<1.064) &&!(dimu23>0.974 && dimu23<1.064)";

    dataloader->PrepareTrainingAndTestTree(reso_cat&&preselCut&&cutS&&common_cut&&veto_phi, reso_cat&&preselCut&&cutB&&common_cut&&veto_phi, prepareTrainTestOptions);

     // Booking of MVA methods : BDT
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", 
          "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50" );
          //"!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=50" );

     // Booking of MVA methods : MLP
     //factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
 
     // Training the MVA methods
     factory->TrainAllMethods();
     
     // Testing the MVA methods
     factory->TestAllMethods();
     
     // Evaluating the MVA methods
     factory->EvaluateAllMethods();

    // Save the output
    fout->Close();
    
    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()){
        TMVAGui("TMVA_"+TMVA_outputpath+categ+".root");
    }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}

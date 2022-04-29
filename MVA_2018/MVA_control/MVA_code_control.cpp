//root -l MVA_code_control.cpp

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

#include "Control_common.h"

using namespace TMVA;


void MVA_code_control(){

    // Output file
    TFile *fout = new TFile("TMVA_"+TMVA_outputpath_control+".root", "RECREATE");

    std::vector<TTree*> sigTree, bkgTree;

    // Get the signal and background trees from TFile source(s);
    //signal
    TString treeName = "FinalTree_Control";
    //DsPhiPi MC
    TFile *f_sig_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_DsPhiPi);
    cout<<"Opening input file "<<inputpath_DsPhiPi<<endl;
    if (!f_sig_ds || !f_sig_ds->IsOpen()) {
         f_sig_ds = new TFile(inputpath_DsPhiPi);
    }
    sigTree.push_back( (TTree*)f_sig_ds->Get(treeName));

    //background
    //Loop on run
    int n_bkg = sizeof(inputpath_datarun_control)/sizeof(inputpath_datarun_control[0]);
    for(int j = 0; j<n_bkg; j++){
        TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun_control[j]);
        cout<<"Opening input file "<<inputpath_datarun_control[j]<<endl;
        if (!f_bkg || !f_bkg->IsOpen()) {
            f_bkg = new TFile(inputpath_datarun_control[j]);
        }
        bkgTree.push_back((TTree*)f_bkg->Get(treeName));
    }

    //add TTreeFriend - muon quality variables
    TString friendTreeName[] = {"TreeMu1=TreeMu1", "TreeMu2=TreeMu2"};
    for(int i=0; i<2; i++){
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_DsPhiPi<<endl;
        sigTree.at(0)->AddFriend(friendTreeName[i], inputpath_DsPhiPi);
    }
    for(int i=0; i<2; i++){
        for(int j=0; j<n_bkg; j++){
            cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_datarun_control[j]<<endl;
            bkgTree.at(j)->AddFriend(friendTreeName[i], inputpath_datarun_control[j]);
        }
    }

    //add TTreeFriend - MuonID evaluation
    TString muIdTreeName[] = {"MuonIDeval_Mu1=MuonIDeval_Mu1", "MuonIDeval_Mu2=MuonIDeval_Mu2"};
    TString inputpath_DsPhiPi_muId = inputpath_DsPhiPi.ReplaceAll(".root", "_MuonID.root");
    for(int i=0; i<2; i++){
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_DsPhiPi_muId<<endl;
        sigTree.at(0)->AddFriend(muIdTreeName[i], inputpath_DsPhiPi_muId);
    }
    for(int j=0; j<n_bkg; j++){
        TString inputpath_data_muId = inputpath_datarun_control[j].ReplaceAll(".root", "_MuonID.root");
        for(int i=0; i<2; i++){
            cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_data_muId<<endl;
            bkgTree.at(j)->AddFriend(muIdTreeName[i], inputpath_data_muId);
        }
    }


    // Set the event weights per tree
    Double_t sigWeight  = 1.; //1.0 DsPhiPi

    Double_t bkgWeight1 = 1.0;
    Double_t bkgWeight2 = 1.0;
    Double_t bkgWeight3 = 1.0;
    Double_t bkgWeight4 = 1.0;
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath_control);
    
    dataloader->AddSignalTree(sigTree.at(0), sigWeight);
    dataloader->AddBackgroundTree(bkgTree.at(0), bkgWeight1);
    dataloader->AddBackgroundTree(bkgTree.at(1), bkgWeight2);
    dataloader->AddBackgroundTree(bkgTree.at(2), bkgWeight3);
    dataloader->AddBackgroundTree(bkgTree.at(3), bkgWeight4);

    std::vector<TString> var_train_name;
    std::vector<TString> var_train_def;
    std::vector<TString> var_spec_name;
    std::vector<TString> var_spec_def;
    //TString BDTinVar_control; //taken from Control_common.h

    readVarName_control(var_spec_name, var_spec_def, BDTspecVar_control);
    readVarName_control(var_train_name, var_train_def, BDTinVar_control);

    // Spectators declaration
    cout<<"Declaration of spectator variables from file:"<<BDTspecVar_control<<endl;
    for(int k = 0; k<var_spec_name.size(); k++){
        if(var_spec_name.at(k)=="" || var_spec_def.at(k)=="") continue;
        cout<<k<<" - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
        dataloader->AddSpectator(var_spec_def.at(k), var_spec_name.at(k), "", 'D');
    }

    // Variables declaration
    cout<<"Declaration of variables for training from file:"<<BDTinVar_control<<endl;
    for(int k = 0; k<var_train_name.size(); k++){
        if(var_train_name.at(k)=="" || var_train_def.at(k)=="") continue;
        cout<<k<<" - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
        dataloader->AddVariable(var_train_def.at(k), var_train_name.at(k), "", 'D');
    }


    dataloader->SetSignalWeightExpression( "puFactor" );

    TCut cutS = "((tripletMass<2.02 && tripletMass>1.92))"; //Signal -> MC Ds peak 
    TCut cutB = "( (tripletMass<1.80 && tripletMass>1.72) )"; //Background -> data sidebands
    TCut preselCut = "";

    TString prepareTrainTestOptions = ":SplitMode=Random"
                                      ":NormMode=NumEvents"
                                      ":nTrain_Signal=0"
                                      ":nTest_Signal=0" //50/50 splitting used
                                      ":!V";

    //dataloader->PrepareTrainingAndTestTree(reso_cat&&preselCut&&cutS, reso_cat&&preselCut&&cutB, prepareTrainTestOptions);
    dataloader->PrepareTrainingAndTestTree(preselCut&&cutS&&common_cut_control, preselCut&&cutB&&common_cut_control, prepareTrainTestOptions);

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
        TMVAGui("TMVA_"+TMVA_outputpath_control+".root");
    }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}

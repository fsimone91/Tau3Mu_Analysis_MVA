//root -l MVA_code_2018.cpp\(\"A\"\)

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

void MVA_code_2018(TString categ){
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
    std::vector<Double_t> bkgWeight; 
    for(int j = 0; j<n_bkg; j++){
        TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[j]);
        if (!f_bkg || !f_bkg->IsOpen()) {
            f_bkg = new TFile(inputpath_datarun[j]);
        }
        bkgTree.push_back((TTree*)f_bkg->Get(treeName));
        bkgWeight.push_back(1.0);
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

    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
    
    dataloader->AddSignalTree(sigTree.at(0), sigWeight1);
    dataloader->AddSignalTree(sigTree.at(1), sigWeight2);
    dataloader->AddSignalTree(sigTree.at(2), sigWeight3);
    for(int j=0; j<n_bkg; j++){
        dataloader->AddBackgroundTree(bkgTree.at(j), bkgWeight.at(j));
    }

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

    dataloader->SetSignalWeightExpression( "puFactor" );

    TCut cutS = "tripletMass<2.0 && tripletMass>1.62"; //Signal -> MC full range 
    //TCut cutB = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)"; //Background -> data sidebands
    TCut cutB = "";
    if(categ.Contains("A")) cutB = "(tripletMass<=1.753 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.801)";
    if(categ.Contains("B")) cutB = "(tripletMass<=1.739 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.815)";
    if(categ.Contains("C")) cutB = "(tripletMass<=1.727 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.827)";
    TCut preselCut = "";

    TCut reso_A = "tripletMassReso < 0.007";
    TCut reso_B = "tripletMassReso >= 0.007 && tripletMassReso <= 0.0105";
    TCut reso_C = "tripletMassReso > 0.0105";
    TCut reso_cat = "tripletMassReso < 0"; //always false    

    TCut phasespace = "";
    if(categ.Contains("highpt")) phasespace = "Ptmu1 >= 7";
    if(categ.Contains("lowpt"))  phasespace = "Ptmu1 < 7";

    if(categ.Contains("A")) reso_cat = reso_cat || reso_A;
    if(categ.Contains("B")) reso_cat = reso_cat || reso_B;
    if(categ.Contains("C")) reso_cat = reso_cat || reso_C;

    TString prepareTrainTestOptions = ":SplitMode=Random"
                                      ":NormMode=NumEvents";
    if(categ=="A") { 
             prepareTrainTestOptions = prepareTrainTestOptions+
                                       //":nTrain_Signal=0"
                                       ":nTest_Signal=17300" //tot 86.6k
                                       //":nTest_Signal=14000" //tot 70.3k
                                       ":nTest_Background=6400" //tot 32.3k
                                       ":!V";
                        }
    else if(categ=="B"){
             prepareTrainTestOptions = prepareTrainTestOptions+
                                       //":nTrain_Signal=0"
                                       ":nTest_Signal=34000" //tot 171k
                                       //":nTest_Signal=28400" //tot 142k
                                       ":nTest_Background=15000" //tot 75.7k
                                       ":!V";
                        }
    else if(categ=="C"){
             prepareTrainTestOptions = prepareTrainTestOptions+
                                       //":nTrain_Signal=0"
                                       ":nTest_Signal=12800" //tot 64k
                                       //":nTest_Signal=11000" //tot 55.3k
                                       ":nTest_Background=6800" //tot 34.3k
                                       ":!V";
                        }
    TString veto_phi = "",  phi_low = "", phi_high = "";
        if(categ=="A") { phi_low = "0.998"; phi_high = "1.040"; }
        if(categ=="B") { phi_low = "0.991"; phi_high = "1.047"; }
        if(categ=="C") { phi_low = "0.984"; phi_high = "1.054"; }
    veto_phi = "!(dimu12_ref>"+phi_low+" && dimu12_ref<"+phi_high+") && "
               "!(dimu13_ref>"+phi_low+" && dimu13_ref<"+phi_high+") && "
               "!(dimu23_ref>"+phi_low+" && dimu23_ref<"+phi_high+")";

    dataloader->PrepareTrainingAndTestTree(reso_cat&&preselCut&&phasespace&&cutS&&common_cut&&veto_phi, reso_cat&&preselCut&&phasespace&&cutB&&common_cut&&veto_phi, prepareTrainTestOptions);

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

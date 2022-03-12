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


void MVA_code_run2(TString categ){
    //Check on input argument
    if(!categ.Contains("A") && !categ.Contains("B") && !categ.Contains("C")){
        cout << "Please choose between any combination of 'A', 'B' and 'C'" << endl;
        return;
    }

    // Output file
    TFile *fout = new TFile("TMVA_"+TMVA_outputpath+categ+".root", "RECREATE");

    TString cat_name[] = {"A", "B", "C"};
    std::vector<TTree*> sigTree_2017, bkgTree_2017;
    std::vector<TTree*> sigTree_2018, bkgTree_2018;
    TString treeName = "FinalTree";
    TString friendTreeName[] = {"TreeMu1=TreeMu1", "TreeMu2=TreeMu2", "TreeMu3=TreeMu3"};
    TString muIdTreeName[] = {"MuonIDeval_Mu1=MuonIDeval_Mu1", "MuonIDeval_Mu2=MuonIDeval_Mu2", "MuonIDeval_Mu3=MuonIDeval_Mu3"};

    //////////////////////// 2017 /////////////////////////////
    // Get the signal and background trees from TFile source(s);
    //signal
    //Ds
    TFile *f_sig_2017_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds_2017);
    if (!f_sig_2017_ds || !f_sig_2017_ds->IsOpen()) {
         f_sig_2017_ds = new TFile(inputpath_Ds_2017);
    }
    sigTree_2017.push_back( (TTree*)f_sig_2017_ds->Get(treeName));
    //B0
    TFile *f_sig_2017_b0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0_2017);
    if (!f_sig_2017_b0 || !f_sig_2017_b0->IsOpen()) {
        f_sig_2017_b0 = new TFile(inputpath_B0_2017);
    }
    sigTree_2017.push_back( (TTree*)f_sig_2017_b0->Get(treeName));
    //Bp
    TFile *f_sig_2017_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp_2017);
    if (!f_sig_2017_bp || !f_sig_2017_bp->IsOpen()) {
       f_sig_2017_bp = new TFile(inputpath_Bp_2017);
    }
    sigTree_2017.push_back( (TTree*)f_sig_2017_bp->Get(treeName));

    //background
    //Loop on run
    int n_bkg_2017 = sizeof(inputpath_datarun_2017)/sizeof(inputpath_datarun_2017[0]);
    std::vector<Double_t> bkgWeight_2017; 
    for(int j = 0; j<n_bkg_2017; j++){
        TFile *f_bkg_2017 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun_2017[j]);
        if (!f_bkg_2017 || !f_bkg_2017->IsOpen()) {
            f_bkg_2017 = new TFile(inputpath_datarun_2017[j]);
        }
        bkgTree_2017.push_back((TTree*)f_bkg_2017->Get(treeName));
        bkgWeight_2017.push_back(1.0);
    }

    //add TTreeFriend - muon quality variables
    for(int i=0; i<3; i++){
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_Ds_2017<<endl;
        sigTree_2017.at(0)->AddFriend(friendTreeName[i], inputpath_Ds_2017);
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_B0_2017<<endl;
        sigTree_2017.at(1)->AddFriend(friendTreeName[i], inputpath_B0_2017);
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_Bp_2017<<endl;
        sigTree_2017.at(2)->AddFriend(friendTreeName[i], inputpath_Bp_2017);
    }
    for(int i=0; i<3; i++){
        for(int j=0; j<n_bkg_2017; j++){
            cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_datarun_2017[j]<<endl;
            bkgTree_2017.at(j)->AddFriend(friendTreeName[i], inputpath_datarun_2017[j]);
        }
    }

    //add TTreeFriend - MuonID evaluation
    TString inputpath_Ds_muId_2017 = inputpath_Ds_2017.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_B0_muId_2017 = inputpath_B0_2017.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_Bp_muId_2017 = inputpath_Bp_2017.ReplaceAll(".root", "_MuonID.root");
    for(int i=0; i<3; i++){
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_Ds_muId_2017<<endl;
        sigTree_2017.at(0)->AddFriend(muIdTreeName[i], inputpath_Ds_muId_2017);
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_B0_muId_2017<<endl;
        sigTree_2017.at(1)->AddFriend(muIdTreeName[i], inputpath_B0_muId_2017);
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_Bp_muId_2017<<endl;
        sigTree_2017.at(2)->AddFriend(muIdTreeName[i], inputpath_Bp_muId_2017);
    }
    for(int j=0; j<n_bkg_2017; j++){
        TString inputpath_data_muId = inputpath_datarun_2017[j].ReplaceAll(".root", "_MuonID.root");
        for(int i=0; i<3; i++){
            cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_data_muId<<endl;
            bkgTree_2017.at(j)->AddFriend(muIdTreeName[i], inputpath_data_muId);
        }
    }


    // Set the event weights per tree
    Double_t sigWeight1_2017  = wNormDs_2017/wNormDs_2017; //1.0 Ds
    Double_t sigWeight2_2017  = wNormB0_2017/wNormDs_2017; //B0
    Double_t sigWeight3_2017  = wNormBp_2017/wNormDs_2017; //Bp
    cout<<"Ds sigWeight1: "<<sigWeight1_2017<<" - B0 sigWeight2: "<<sigWeight2_2017<<" - Bp sigWeight3: "<<sigWeight3_2017<<endl;

    Factory *factory = new Factory("TMVA_new", fout, "V:!Silent:Color:DrawProgressBar");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
    
    dataloader->AddSignalTree(sigTree_2017.at(0), sigWeight1_2017);
    dataloader->AddSignalTree(sigTree_2017.at(1), sigWeight2_2017);
    dataloader->AddSignalTree(sigTree_2017.at(2), sigWeight3_2017);
    for(int j=0; j<n_bkg_2017; j++){
        dataloader->AddBackgroundTree(bkgTree_2017.at(j), bkgWeight_2017.at(j));
    }
    ///////////////////////////////////////////////////////////////

    //////////////////////// 2018 /////////////////////////////
    // Get the signal and background trees from TFile source(s);
    //signal
    //Ds
    TFile *f_sig_2018_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds_2018);
    if (!f_sig_2018_ds || !f_sig_2018_ds->IsOpen()) {
         f_sig_2018_ds = new TFile(inputpath_Ds_2018);
    }
    sigTree_2018.push_back( (TTree*)f_sig_2018_ds->Get(treeName));
    //B0
    TFile *f_sig_2018_b0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0_2018);
    if (!f_sig_2018_b0 || !f_sig_2018_b0->IsOpen()) {
        f_sig_2018_b0 = new TFile(inputpath_B0_2018);
    }
    sigTree_2018.push_back( (TTree*)f_sig_2018_b0->Get(treeName));
    //Bp
    TFile *f_sig_2018_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp_2018);
    if (!f_sig_2018_bp || !f_sig_2018_bp->IsOpen()) {
       f_sig_2018_bp = new TFile(inputpath_Bp_2018);
    }
    sigTree_2018.push_back( (TTree*)f_sig_2018_bp->Get(treeName));

    //background
    //Loop on run
    int n_bkg_2018 = sizeof(inputpath_datarun_2018)/sizeof(inputpath_datarun_2018[0]);
    std::vector<Double_t> bkgWeight_2018; 
    for(int j = 0; j<n_bkg_2018; j++){
        TFile *f_bkg_2018 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun_2018[j]);
        if (!f_bkg_2018 || !f_bkg_2018->IsOpen()) {
            f_bkg_2018 = new TFile(inputpath_datarun_2018[j]);
        }
        bkgTree_2018.push_back((TTree*)f_bkg_2018->Get(treeName));
        bkgWeight_2018.push_back(1.0);
    }

    //add TTreeFriend - muon quality variables
    for(int i=0; i<3; i++){
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_Ds_2018<<endl;
        sigTree_2018.at(0)->AddFriend(friendTreeName[i], inputpath_Ds_2018);
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_B0_2018<<endl;
        sigTree_2018.at(1)->AddFriend(friendTreeName[i], inputpath_B0_2018);
        cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_Bp_2018<<endl;
        sigTree_2018.at(2)->AddFriend(friendTreeName[i], inputpath_Bp_2018);
    }
    for(int i=0; i<3; i++){
        for(int j=0; j<n_bkg_2018; j++){
            cout<<"Adding FriendTree "<<friendTreeName[i]<<" from file "<<inputpath_datarun_2018[j]<<endl;
            bkgTree_2018.at(j)->AddFriend(friendTreeName[i], inputpath_datarun_2018[j]);
        }
    }

    //add TTreeFriend - MuonID evaluation
    TString inputpath_Ds_muId_2018 = inputpath_Ds_2018.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_B0_muId_2018 = inputpath_B0_2018.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_Bp_muId_2018 = inputpath_Bp_2018.ReplaceAll(".root", "_MuonID.root");
    for(int i=0; i<3; i++){
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_Ds_muId_2018<<endl;
        sigTree_2018.at(0)->AddFriend(muIdTreeName[i], inputpath_Ds_muId_2018);
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_B0_muId_2018<<endl;
        sigTree_2018.at(1)->AddFriend(muIdTreeName[i], inputpath_B0_muId_2018);
        cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_Bp_muId_2018<<endl;
        sigTree_2018.at(2)->AddFriend(muIdTreeName[i], inputpath_Bp_muId_2018);
    }
    for(int j=0; j<n_bkg_2018; j++){
        TString inputpath_data_muId = inputpath_datarun_2018[j].ReplaceAll(".root", "_MuonID.root");
        for(int i=0; i<3; i++){
            cout<<"Adding FriendTree "<<muIdTreeName[i]<<" from file "<<inputpath_data_muId<<endl;
            bkgTree_2018.at(j)->AddFriend(muIdTreeName[i], inputpath_data_muId);
        }
    }


    // Set the event weights per tree
    Double_t sigWeight1_2018  = wNormDs_2018/wNormDs_2018; //1.0 Ds
    Double_t sigWeight2_2018  = wNormB0_2018/wNormDs_2018; //B0
    Double_t sigWeight3_2018  = wNormBp_2018/wNormDs_2018; //Bp
    cout<<"Ds sigWeight1: "<<sigWeight1_2018<<" - B0 sigWeight2: "<<sigWeight2_2018<<" - Bp sigWeight3: "<<sigWeight3_2018<<endl;

    dataloader->AddSignalTree(sigTree_2018.at(0), sigWeight1_2018);
    dataloader->AddSignalTree(sigTree_2018.at(1), sigWeight2_2018);
    dataloader->AddSignalTree(sigTree_2018.at(2), sigWeight3_2018);
    for(int j=0; j<n_bkg_2018; j++){
        dataloader->AddBackgroundTree(bkgTree_2018.at(j), bkgWeight_2018.at(j));
    }
    ///////////////////////////////////////////////////////////////

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

    dataloader->AddSpectator("evt:=int(evt)%4096", 'I');

    // Variables declaration
    cout<<"Declaration of variables for training - category "<<categ<<" from file:"<<BDTinVar<<endl;
    for(int k = 0; k<var_train_name.size(); k++){
        if(var_train_name.at(k)=="" || var_train_def.at(k)=="") continue;
        cout<<k<<" - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
        dataloader->AddVariable(var_train_def.at(k), var_train_name.at(k), "", 'D');
    }

    dataloader->SetSignalWeightExpression( "puFactor" );

    TCut cutS = "tripletMass<2.0 && tripletMass>1.62"; //Signal -> MC full 
    //TCut cutS = "tripletMass<1.80 && tripletMass>1.754"; //Signal -> MC peak 
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

    TString prepareTrainTestOptions = "";
    if(doCV) prepareTrainTestOptions = 
                                       "nTest_Signal=1" //with CV the test set is unused. For this reason, we assign 1 event only to it (0 cannot be used, would split test/training 50/50//
                                       ":nTest_Background=1" //same comment as above
                                       ":NormMode=NumEvents"
                                       ":SplitMode=Random"
                                       ":V";
    else     prepareTrainTestOptions = ":nTrain_Signal=0"
                                       ":nTrain_Background=0"
                                       ":SplitMode=Random"
                                       ":NormMode=NumEvents"
                                       ":!V";
    
    TString veto_phi = "",  phi_low = "", phi_high = "";
        if(categ=="A") { phi_low = "0.99"; phi_high = "1.050"; }
        if(categ=="B") { phi_low = "0.98"; phi_high = "1.060"; }
        if(categ=="C") { phi_low = "0.965"; phi_high = "1.075"; }
    veto_phi = "!(dimu12_ref>"+phi_low+" && dimu12_ref<"+phi_high+") && "
               "!(dimu13_ref>"+phi_low+" && dimu13_ref<"+phi_high+") && "
               "!(dimu23_ref>"+phi_low+" && dimu23_ref<"+phi_high+")";

    TString common_cut = "(fv_nC>0 && !(TMath::IsNaN(fv_d3Dsig)) && (MuonIDeval_Mu3.year==2017?(bs_sv_d2Dsig>=3.75):(bs_sv_d2Dsig>=2.0) ))"; 
    dataloader->PrepareTrainingAndTestTree(reso_cat&&preselCut&&phasespace&&cutS&&common_cut&&veto_phi, reso_cat&&preselCut&&phasespace&&cutB&&common_cut&&veto_phi, prepareTrainTestOptions);

    //UInt_t numFolds = 3;
    TString analysisType = "Classification";
    TString splitExpr = "int(fabs([evt]))%int([NumFolds])"; //"int([evt])\%int([numFolds])";
    TString outputEnsembling = "None"; //"Avg";
    TString cvOptions = Form("!V"
                           ":!Silent"
                           //":VerboseLevel=Debug"
                           ":DrawProgressBar=True"
                           ":ModelPersistence"
                           ":AnalysisType=%s"
                           ":SplitType=Deterministic"
                           ":NumFolds=%i"
                           ":SplitExpr=%s"
                           ":OutputEnsembling=%s"
                           ":FoldFileOutput=True",
                           analysisType.Data(), 
                           numFolds,
                           splitExpr.Data(),
                           outputEnsembling.Data()
                           );
    
    TMVA::CrossValidation cv{"TMVACrossValidation", dataloader, fout, cvOptions};

    if(doCV){
        // Booking of MVA methods : BDTG
        cv.BookMethod( TMVA::Types::kBDT, "BDTG",
             "!H:V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:DoBoostMonitor=True"
             ":AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50");
       // cv.BookMethod(TMVA::Types::kBDT, "BDTG",
       //          "!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad"
       //          ":SkipNormalization=True"
       //          ":NegWeightTreatment=Pray:Shrinkage=0.10:nCuts=20"
       //          ":MaxDepth=2");
        cv.Evaluate();

        size_t iMethod = 0;
        for (auto && result : cv.GetResults()) {
            std::cout << "Summary for method " << cv.GetMethods()[iMethod++].GetValue<TString>("MethodName")
                    << std::endl;
            for (UInt_t iFold = 0; iFold<cv.GetNumFolds(); ++iFold) {
                std::cout << "\tFold " << iFold << ": "
                       << "ROC int: " << result.GetROCValues()[iFold]
                       << ", "
                       << "BkgEff@SigEff=0.3: " << result.GetEff30Values()[iFold]
                       << std::endl;
            }
        }
    }
    else {
        // Booking of MVA methods : BDT
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", 
             "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50" );
        // Training the MVA methods
        factory->TrainAllMethods();
        // Testing the MVA methods
        factory->TestAllMethods();
        // Evaluating the MVA methods
        factory->EvaluateAllMethods();
    } 

    // Save the output
    fout->Close();
    std::cout << "==> Wrote root file: " << fout->GetName() << std::endl;

    if (!gROOT->IsBatch()) {
       // gui
       cout<<"TMVA::TMVAGui("<<fout->GetName()<<")"<<endl;
       TMVA::TMVAGui(fout->GetName());
    }
       
    delete factory;
    delete dataloader;
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();

    return 0;
}

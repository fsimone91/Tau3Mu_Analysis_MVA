#include <stdio.h>
#include "Control_common.h"
#include "/lustrehome/fsimone/MuonID_study/MuonMVA_2017/MuonID_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

//root -l -b evalMuonID_addbranch.cpp

void upd(TString filename, TString treename){

    TFile *f = new TFile(filename,"update"); 
    TTree *tin = (TTree*)f->Get(treename); 

    //Create TTreeReader
    cout<<"Accessing to input tree "<<treename<<endl;
    TTreeReader treeReader(tin);
    
    Double_t muonID_value;
    TBranch *b_muid = tin->Branch("mu_MVAid",&muonID_value,"mu_MVAid/D"); 

    //create the TMVA::Reader
    TMVA::Tools::Instance();
    TMVA::Reader *MVAreader_barrel = new TMVA::Reader("!Color:!Silent:!V");
    TMVA::Reader *MVAreader_endcap = new TMVA::Reader("!Color:!Silent:!V");

    //number of spectators
    size_t n_spec = var_MuonID_spec_names.size();
    std::vector<Float_t> var_spec;
    for(int j = 0; j<n_spec; j++) var_spec.push_back(0);
    // Spectators declaration
    for(int k = 0; k<n_spec; k++){
        MVAreader_barrel->TMVA::Reader::AddSpectator(var_MuonID_spec_names.at(k), &var_spec.at(k));
        MVAreader_endcap->TMVA::Reader::AddSpectator(var_MuonID_spec_names.at(k), &var_spec.at(k));
        cout<<var_MuonID_spec_names.at(k)<<endl;
    }
    //number of variables used for training
    size_t n_train = var_MuonID_train_def.size();
    std::vector<Float_t> var_train;
    for(int j = 0; j<n_train; j++) var_train.push_back(0);
    // Variables declaration
    for(int k = 0; k<n_train; k++){
        MVAreader_barrel->TMVA::Reader::AddVariable(var_MuonID_train_def.at(k), &var_train.at(k));
        MVAreader_endcap->TMVA::Reader::AddVariable(var_MuonID_train_def.at(k), &var_train.at(k));
        cout<<var_MuonID_train_def.at(k)<<endl;
    }
    TString category = "barrel";
    //Book the MVA method used for MuonID
    TString weight_path_barrel = "/lustrehome/fsimone/MuonID_study/MuonMVA_2017/"+TMVA_MuonID_inputpath+category+TMVA_MuonID_weightfilename;
    Bool_t weightfileExists = (gSystem->AccessPathName(weight_path_barrel) == kFALSE);
    if (weightfileExists) {
       MVAreader_barrel->TMVA::Reader::BookMVA(method_MuonID, weight_path_barrel);
       cout<<"Using weights in "<<weight_path_barrel<<endl;
    } else {
       std::cout << "Weightfile " <<weight_path_barrel<<" for method " << method_MuonID << " not found."
                    " Did you run TMVACrossValidation with a specified splitExpr?" << std::endl;
       exit(0);
    }
    category = "endcap";
    //Book the MVA method used for MuonID
    TString weight_path_endcap = "/lustrehome/fsimone/MuonID_study/MuonMVA_2017/"+TMVA_MuonID_inputpath+category+TMVA_MuonID_weightfilename;
    weightfileExists = (gSystem->AccessPathName(weight_path_endcap) == kFALSE);
    if (weightfileExists) {
       MVAreader_endcap->TMVA::Reader::BookMVA(method_MuonID, weight_path_endcap);
       cout<<"Using weights in "<<weight_path_endcap<<endl;
    } else {
       std::cout << "Weightfile " <<weight_path_endcap<<" for method " << method_MuonID << " not found."
                    " Did you run TMVACrossValidation with a specified splitExpr?" << std::endl;
       exit(0);
    }

    //Read branches
    std::vector<TTreeReaderValue<double>> reader_spec;
    for(int k = 0; k<n_spec; k++){
        reader_spec.emplace_back(treeReader, var_MuonID_spec_names.at(k));
    }
    std::vector<TTreeReaderValue<double>> reader_train;
    for(int k = 0; k<n_train; k++){
        reader_train.emplace_back(treeReader, var_MuonID_train_names.at(k));
    }
    TTreeReaderValue<double> reader_mueta(treeReader, "mu_eta");

    //prepare branch in output tree
    //Loop on input Tree
    Double_t score = 0;
    while (treeReader.Next()) {
        for(int k = 0; k<n_spec; k++){
            var_spec.at(k)= *reader_spec.at(k);
        }
        for(int k = 0; k<n_train; k++){
            var_train.at(k) = *reader_train.at(k);
        }
        
        //Evaluate method(s) and fill histogram or MiniTree
        if(abs(*reader_mueta) < 1.2) score = MVAreader_barrel->EvaluateMVA(method_MuonID);
        else score  = MVAreader_endcap->EvaluateMVA(method_MuonID);
        muonID_value = score;
        b_muid->Fill();
    }
    delete MVAreader_barrel;
    delete MVAreader_endcap;
    f->cd();
    tin->Print();
    tin->Write(); 
    f->Close();
}


void evalMuonID_addbranch(){

    //loop input files
    size_t n_bkg = sizeof(inputpath_datarun_control)/sizeof(inputpath_datarun_control[0]);
    //Bkg samples
    //for(auto i=0; i<n_bkg; i++){
    //    upd(inputpath_datarun_control[i], "TreeMu1");
    //    upd(inputpath_datarun_control[i], "TreeMu2");
    //    std::cout<<"Updated input file: "<<inputpath_datarun_control[i]<<std::endl;
    //}
    //MC Tau3Mu
    upd(inputpath_DsPhiPi, "TreeMu1");
    upd(inputpath_DsPhiPi, "TreeMu2");
    std::cout<<"Updated input file: "<<inputpath_DsPhiPi<<std::endl;
}

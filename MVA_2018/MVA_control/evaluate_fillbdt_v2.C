#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "Control_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TRandom.h"
#include "TRandom3.h"

void fill_BDT_score(TChain *t, TH1F* hBDTdecision, Int_t isMC, TTree *tout, TString* outvar_name, Double_t* outvar_value){

    //Create TTreeReader
    cout<<"Accessing to input tree"<<endl;
    TTreeReader treeReader(t);

    //create the TMVA::Reader
    TMVA::Tools::Instance();
    TMVA::Reader *MVAreader = new TMVA::Reader("!Color:!Silent:!V");

    //Variables
    std::vector<TString> var_train_name;
    std::vector<TString> var_train_def;
    std::vector<TString> var_spec_name;
    std::vector<TString> var_spec_def;
    //TString BDTinVar; //taken from Control_common.h

    readVarName_control(var_spec_name, var_spec_def, BDTspecVar_control);
    readVarName_control(var_train_name, var_train_def, BDTinVar_control);

    // Spectators declaration
    cout<<"Declaration of spectator variables from file:"<<BDTspecVar_control<<endl;
    for(int k = 0; k<var_spec_name.size(); k++){
        cout<<k<<" - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
    }

    // Variables declaration
    cout<<"Declaration of variables for training from file:"<<BDTinVar_control<<endl;
    for(int k = 0; k<var_train_name.size(); k++){
        cout<<k<<" - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
    }

    //number of spectators
    size_t n_spec = var_spec_name.size();
    std::vector<Float_t> var_spec;
    for(int j = 0; j<n_spec; j++) var_spec.push_back(0);
    // Spectators declaration
    for(int k = 0; k<n_spec; k++){
        MVAreader->TMVA::Reader::AddSpectator(var_spec_def.at(k), &var_spec.at(k));
        cout<<k<<" added to Reader - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
    }
    //number of variables used for training
    size_t n_train = var_train_name.size();
    std::vector<Float_t> var_train;
    for(int j = 0; j<n_train; j++) var_train.push_back(0);
    // Variables declaration
    for(int k = 0; k<n_train; k++){
        MVAreader->TMVA::Reader::AddVariable(var_train_def.at(k), &var_train.at(k));
        cout<<k<<" added to Reader - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
    }

    //Book the MVA method(s)
    TString weight_path = TMVA_inputpath_control+TMVA_weightfilename_control;
    Bool_t weightfileExists = (gSystem->AccessPathName(weight_path) == kFALSE);
    if (weightfileExists) {
       MVAreader->TMVA::Reader::BookMVA(method_control, weight_path);
       cout<<"Using weights in "<<weight_path<<endl;
    } else {
       std::cout << "Weightfile " <<weight_path<<" for method " << method_control << " not found."
                    " Did you run TMVACrossValidation with a specified splitExpr?" << std::endl;
       exit(0);
    }

    TTreeReaderValue<double> reader_triplMass(treeReader, "tripletMass");

    //Branches for output tree
    int nvar = 30;
    std::vector<TTreeReaderValue<double>> reader_out;
    for(auto v = 0; v<nvar-11; v++){
        reader_out.emplace_back(treeReader, outvar_name[v]);
        cout<<"   reading from TChain branch "<<outvar_name[v]<<endl;
    }

    cout<<"booleans must be treated separately"<<endl;
    std::vector<TTreeReaderValue<bool>> reader_out_bool;
    for(auto v = nvar-11; v<nvar-7; v++){
        reader_out_bool.emplace_back(treeReader, outvar_name[v]);
        cout<<"   reading from TChain branch "<<outvar_name[v]<<endl;
    }
    cout<<"\n------------------------"<<endl;

    TTreeReaderValue<double> reader_dxy1(treeReader, "dxy1");
    TTreeReaderValue<double> reader_dxyErr1(treeReader, "dxyErr1");
    TTreeReaderValue<double> reader_dxy2(treeReader, "dxy2");
    TTreeReaderValue<double> reader_dxyErr2(treeReader, "dxyErr2");
    TTreeReaderValue<bool> reader_doubleMu4(treeReader, "l1double_DoubleMu4_fired");
    TTreeReaderValue<bool> reader_doubleMu0(treeReader, "l1double_DoubleMu0_fired");

    TTreeReaderValue<float> reader_muid1(treeReader, "MuonIDeval_Mu1.MuonID");
    TTreeReaderValue<float> reader_muid2(treeReader, "MuonIDeval_Mu2.MuonID");

    //nBranches for MVA
    std::vector<TTreeReaderValue<double>> reader_spec;
    for(int k = 0; k<n_spec; k++){
        cout<<"   reading from TChain branch "<<var_spec_name.at(k)<<endl;
        reader_spec.emplace_back(treeReader, var_spec_name.at(k));
    }
    std::vector<TTreeReaderValue<double>> reader_train;
    for(int k = 0; k<n_train; k++){
        if(var_train_name.at(k)=="d0sig_min") continue;
        if(var_train_name.at(k)=="d0sig_max") continue;
        if(var_train_name.at(k)=="d0_min") continue;
        if(var_train_name.at(k)=="MuonIDeval_Mu1.MuonID") continue;
        if(var_train_name.at(k)=="MuonIDeval_Mu2.MuonID") continue;
        cout<<"   reading from TChain branch "<<var_train_name.at(k)<<endl;
        reader_train.emplace_back(treeReader, var_train_name.at(k));
    }

    //Loop on input Tree
    while (treeReader.Next()) {

        for(int k = 0; k<n_spec; k++){
            var_spec.at(k)= *reader_spec.at(k);
        }
        for(int k = 0; k<n_train; k++){
            if(var_train_name.at(k)=="d0sig_min"){ 
                if((*reader_dxyErr1)!=0 && (*reader_dxyErr2)!=0) var_train.at(k) = std::fmin ( abs(*reader_dxy1)/(*reader_dxyErr1) , abs(*reader_dxy2)/(*reader_dxyErr2) );
                else var_train.at(k) = 0;
            }
            else if(var_train_name.at(k)=="d0sig_max"){ 
                if((*reader_dxyErr1)!=0 && (*reader_dxyErr2)!=0) var_train.at(k) = std::fmax ( abs(*reader_dxy1)/(*reader_dxyErr1), abs(*reader_dxy2)/(*reader_dxyErr2) );
                else var_train.at(k) = 0;
            }
            else if(var_train_name.at(k)=="d0_min"){ 
                if((*reader_dxyErr1)!=0 && (*reader_dxyErr2)!=0) var_train.at(k) = std::fmin ( abs(*reader_dxy1), abs(*reader_dxy2) );
                else var_train.at(k) = 0;
            }
            else if(var_train_name.at(k)=="MuonIDeval_Mu1.MuonID") var_train.at(k) = *reader_muid1;
            else if(var_train_name.at(k)=="MuonIDeval_Mu2.MuonID") var_train.at(k) = *reader_muid2;
            else var_train.at(k) = *reader_train.at(k);
        }
        //Evaluate method(s) and fill histogram or MiniTree
        auto BDTscore = MVAreader->EvaluateMVA(method_control);
        auto mass = *reader_triplMass;
        Int_t isSB = 0;
        if( ( mass > 1.70 && mass < 1.80 ) ) isSB = 1; 

        if((!isMC && isSB) || isMC) hBDTdecision->Fill(BDTscore);
        // for(int k = 0; k<n_spec; k++){
        //     cout<<"   spect. "<<var_spec.at(k)<<endl;
        // }
        // for(int k = 0; k<n_train; k++){
        //     cout<<"   train "<<var_train.at(k)<<endl;
        // }
        // cout<<"reso "<<reso<<" bdt "<<MVAreader->EvaluateMVA(method)<<endl;        

        for(auto v = 0; v<nvar-11; v++){
            outvar_value[v] = *reader_out.at(v);
            //cout<<"     to output tree -> "<<outvar_name[v]<<"="<<outvar_value[v]<<endl;
        }
        for(auto v = nvar-11; v<nvar-7; v++){
            outvar_value[v] = *reader_out_bool.at(v-(nvar-11));
            //cout<<"     to output tree -> "<<outvar_name[v]<<"="<<outvar_value[v]<<endl;
        }
        for(auto v = nvar-7; v<nvar; v++){
            if(outvar_name[v] == "MVA_glbmu1")  outvar_value[v] = *reader_muid1;
            if(outvar_name[v] == "MVA_glbmu2")  outvar_value[v] = *reader_muid2;
            if(outvar_name[v] == "bdt")  outvar_value[v] = BDTscore;
            if(outvar_name[v] == "isMC") outvar_value[v] = isMC;
            if(outvar_name[v] == "isSB") outvar_value[v] = isSB;
            if(outvar_name[v] == "weight" && isMC == 0) outvar_value[v] = 1.; //data
            if(outvar_name[v] == "weight" && isMC == 1) outvar_value[v] = 1.; //Ds
            if(outvar_name[v] == "category") outvar_value[v] = 0; //no categorisation in DsPhiPi
            //cout<<"   to output tree -> "<<outvar_name[v]<<"="<<outvar_value[v]<<endl;
        }
        //reject 31.5% of MC events exclusively triggered by DM4
        bool dM4_excl = !(*reader_doubleMu0) && (*reader_doubleMu4);
        if(isMC && dM4_excl){
            TRandom3 rand(0); Double_t x_rand = rand.Rndm( ); //uniformly distributed random number in 0..1
            if(x_rand < 0.315) { continue; }
        }

        tout->Fill();

    }
    delete MVAreader;
}


void evaluate_fillbdt_v2() 
{
    TString run[] = {"A", "B", "C", "D"};
    size_t nrun = sizeof(run)/sizeof(run[0]);
    //open input files
    //data
    TChain *tdata = new TChain("FinalTree_Control");
    for(auto i=0; i<nrun; i++){
        tdata->Add(inputpath_datarun_control[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun_control[i]<<std::endl;
    }
    //MC DsPhiPi
    TChain *tmc1 = new TChain("FinalTree_Control");
    tmc1->Add(inputpath_DsPhiPi); 
    std::cout<<"Opened input file: "<<inputpath_DsPhiPi<<std::endl;

    //data
    TChain *tdata_mu1 = new TChain("TreeMu1");
    TChain *tdata_mu2 = new TChain("TreeMu2");
    //TChain *tdata_mu3 = new TChain("TreeMu3");
    for(auto i=0; i<nrun; i++){
        tdata_mu1->Add(inputpath_datarun_control[i]); 
        tdata_mu2->Add(inputpath_datarun_control[i]); 
        //tdata_mu3->Add(inputpath_datarun_control[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun_control[i]<<std::endl;
    }
    //MC Ds
    TChain *tmc1_mu1 = new TChain("TreeMu1");
    TChain *tmc1_mu2 = new TChain("TreeMu2");
    //TChain *tmc1_mu3 = new TChain("TreeMu3");
    tmc1_mu1->Add(inputpath_DsPhiPi); 
    tmc1_mu2->Add(inputpath_DsPhiPi); 
    //tmc1_mu3->Add(inputpath_DsPhiPi); 
    std::cout<<"Opened input file: "<<inputpath_DsPhiPi<<std::endl;

    //data
    TChain *tdata_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tdata_muid2 = new TChain("MuonIDeval_Mu2");
    for(auto i=0; i<nrun; i++){
        TString inputpath_data_muId = inputpath_datarun_control[i].ReplaceAll(".root", "_MuonID.root");
        tdata_muid1->Add(inputpath_datarun_control[i]);
        tdata_muid2->Add(inputpath_datarun_control[i]);
        std::cout<<"Opened input file: "<<inputpath_data_muId<<std::endl;
    }
    //MC Ds
    TChain *tmc1_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc1_muid2 = new TChain("MuonIDeval_Mu2");
    TString inputpath_Ds_muId = inputpath_DsPhiPi.ReplaceAll(".root", "_MuonID.root");
    tmc1_muid1->Add(inputpath_Ds_muId);
    tmc1_muid2->Add(inputpath_Ds_muId);
    std::cout<<"Opened input file: "<<inputpath_Ds_muId<<std::endl;


    tdata->AddFriend(tdata_mu1);
    tdata->AddFriend(tdata_mu2);
    //tdata->AddFriend(tdata_mu3);

    tmc1->AddFriend(tmc1_mu1);
    tmc1->AddFriend(tmc1_mu2);
    //tmc1->AddFriend(tmc1_mu3);

    tdata->AddFriend(tdata_muid1);
    tdata->AddFriend(tdata_muid2);
    //tdata->AddFriend(tdata_muid3);

    tmc1->AddFriend(tmc1_muid1);
    tmc1->AddFriend(tmc1_muid2);
    //tmc1->AddFriend(tmc1_muid3);

    //Book output MiniTree
    TString fout_tree_path = TMVA_inputpath_control+"outputTree.root";
    cout<<"Output file for final tree: "<<fout_tree_path<<endl;
    TFile *fout_tree = new TFile(fout_tree_path, "recreate");
    fout_tree->cd();
    TTree *tout = new TTree("outputTree","outputTree");

    //variables for output tree taken from input tree
    TString outvar_name [] = { 
                                "evt",
                                "run",
                                "lumi",
                                "puFactor",
                                "tripletMass",
                                "tripletMassReso",
                                "phiMass",
                                //"dimu12",
                                //"dimu13",
                                //"dimu23",
                                "Ptmu1",
                                "Etamu1",
                                "Ptmu2",
                                "Etamu2",
                                "Ptmu3",
                                "Etamu3",
                                "pv_sv_dxy",
                                "pv_sv_dxy_sig",
                                "fv_d3D",
                                "fv_d3Dsig",
                                "bs_sv_d2D",
                                "bs_sv_d2Dsig",

                                "l1triple_fired",
                                "l1double_fired",
                                "l1double_DoubleMu0_fired",
                                "l1double_DoubleMu4_fired",

                                "MVA_glbmu1",
                                "MVA_glbmu2",

                                "isMC",
                                "isSB",
                                "weight",
                                "category",
                                "bdt"
                                };

    int nvar = 30;
    Double_t outvar_value[30] = {0.};

    //Set branches output tree
    for(auto v = 0; v<nvar; v++){
        tout->Branch(outvar_name[v], &outvar_value[v]);
    }

    //Make sure TChain points to first event
    tdata->LoadTree(0);
    tmc1->LoadTree(0);
  
    //Book output histogram or MiniTree
    TString fout_path = TMVA_inputpath_control+"/BDTdecision_v2.root";
    TFile *fout = new TFile(fout_path, "recreate");
    fout->cd();

    TH1F * hBDTdecision_data = new TH1F ("BDTdecision_data_obs","BDTdecision_data_obs", 240, -0.6, 0.6);
    TH1F * hBDTdecision_Ds = new TH1F ("BDTdecision_signalDs","BDTdecision_signalDs", 240, -0.6, 0.6);

    //Read branches, loop on tree, perform evaluate and fill histo
    fill_BDT_score(tdata, hBDTdecision_data, 0, tout, outvar_name, outvar_value); 
    fill_BDT_score(tmc1,  hBDTdecision_Ds,   1, tout, outvar_name, outvar_value); 

    //rescale signal
    Double_t Lumi = 0;
    for(int i=0; i<nrun; i++) Lumi = Lumi + Lumi_data_control[i]; 
    hBDTdecision_Ds->Scale(Lumi*xsection_mc*BR/N_MC);

    //add signals to common histo
    TH1F * hBDTdecisionSignal = new TH1F ("BDTdecision_signal", "BDTdecision_signal", 240, -0.6, 0.6);
    hBDTdecisionSignal->Add(hBDTdecision_Ds);

    //Write and close the file
    fout->cd();
    hBDTdecision_data->Write();
    hBDTdecisionSignal->Write();
    hBDTdecision_Ds->Write();

    fout->Close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    
    fout_tree->Write();
    fout_tree->Close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_tree_path<<"\n\n"<<endl;
    return 0;
}

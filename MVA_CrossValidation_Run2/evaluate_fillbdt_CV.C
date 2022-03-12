#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

void fill_BDT_score(TChain *t, TH1F* hBDTdecision, TString category, Int_t isMC, TTree *tout, TString* outvar_name, Double_t* outvar_value, TString year){

    Double_t wNormDs = 0, Ds_correction = 0, Dplus_correction = 0, wNormB0 = 0, wNormBp = 0, Bs_correction = 0, wNormW = 0;
    if(year == "2017"){
        wNormDs = wNormDs_2017;
        wNormB0 = wNormB0_2017;
        wNormBp = wNormBp_2017;
        wNormW = wNormW_2017;
        
        Ds_correction = Ds_correction_2017;
        Dplus_correction = Dplus_correction_2017;
        Bs_correction = Bs_correction_2017;
    } else if(year == "2018"){
        wNormDs = wNormDs_2018;
        wNormB0 = wNormB0_2018;
        wNormBp = wNormBp_2018;
        wNormW = wNormW_2018;
        
        Ds_correction = Ds_correction_2018;
        Dplus_correction = Dplus_correction_2018;
        Bs_correction = Bs_correction_2018;
    }

    Int_t count = 0;
    cout<<"category "<<category<<" sample "<<isMC<<endl;
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
    TString BDTinVar;
    if(category.Contains("A")) BDTinVar = BDTinVar_A;
    if(category.Contains("B")) BDTinVar = BDTinVar_B;
    if(category.Contains("C")) BDTinVar = BDTinVar_C;

    readVarName(var_spec_name, var_spec_def, BDTspecVar);
    readVarName(var_train_name, var_train_def, BDTinVar);

    // Spectators declaration
    cout<<"Declaration of spectator variables - category "<<category<<" from file:"<<BDTspecVar<<endl;
    for(int k = 0; k<var_spec_name.size(); k++){
        cout<<k<<" - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
    }

    // Variables declaration
    cout<<"Declaration of variables for training - category "<<category<<" from file:"<<BDTinVar<<endl;
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
    Int_t evt=0;
    MVAreader->TMVA::Reader::AddSpectator("evt:=int(evt)%4096", &evt);

    //number of variables used for training
    size_t n_train = var_train_name.size();
    std::vector<Float_t> var_train;
    for(int j = 0; j<n_train; j++) var_train.push_back(0);
    // Variables declaration
    cout<<" -> variables for training"<<endl;
    for(int k = 0; k<n_train; k++){
        MVAreader->TMVA::Reader::AddVariable(var_train_def.at(k), &var_train.at(k));
        cout<<k<<" added to Reader - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
    }

    //Book the MVA method(s)
    TString weight_path = TMVA_inputpath+category+TMVA_weightfilename;
    Bool_t weightfileExists = (gSystem->AccessPathName(weight_path) == kFALSE);
    if (weightfileExists) {
       MVAreader->TMVA::Reader::BookMVA(method, weight_path);
       cout<<"Using weights in "<<weight_path<<endl;
    } else {
       std::cout << "Weightfile " <<weight_path<<" for method " << method << " not found."
                    " Did you run TMVACrossValidation with a specified splitExpr?" << std::endl;
       exit(0);
    }

    TTreeReaderValue<double> reader_triplMassReso(treeReader, "tripletMassReso");
    TTreeReaderValue<double> reader_triplMass(treeReader, "tripletMass");

    //Branches for preselections
    TTreeReaderValue<bool> reader_doubleMu4(treeReader, "l1double_DoubleMu4_fired");
    TTreeReaderValue<bool> reader_doubleMu0(treeReader, "l1double_DoubleMu0_fired");
    TTreeReaderValue<bool> reader_tripleMu(treeReader, "l1triple_fired");

    TTreeReaderValue<double> reader_mu1_pt(treeReader, "TreeMu1.mu_pt");
    TTreeReaderValue<double> reader_mu2_pt(treeReader, "TreeMu2.mu_pt");
    TTreeReaderValue<double> reader_mu3_pt(treeReader, "TreeMu3.mu_pt");
    TTreeReaderValue<double> reader_mu1_eta(treeReader, "TreeMu1.mu_eta");
    TTreeReaderValue<double> reader_mu2_eta(treeReader, "TreeMu2.mu_eta");
    TTreeReaderValue<double> reader_mu3_eta(treeReader, "TreeMu3.mu_eta");

    TTreeReaderValue<double> reader_mu1_tLWM(treeReader, "TreeMu1.mu_trackerLayersWithMeasurement");
    TTreeReaderValue<double> reader_mu2_tLWM(treeReader, "TreeMu2.mu_trackerLayersWithMeasurement");
    TTreeReaderValue<double> reader_mu3_tLWM(treeReader, "TreeMu3.mu_trackerLayersWithMeasurement");
    TTreeReaderValue<double> reader_mu1_validFra(treeReader, "TreeMu1.mu_innerTrack_validFraction");
    TTreeReaderValue<double> reader_mu2_validFra(treeReader, "TreeMu2.mu_innerTrack_validFraction");
    TTreeReaderValue<double> reader_mu3_validFra(treeReader, "TreeMu3.mu_innerTrack_validFraction");
    TTreeReaderValue<double> reader_mu1_iTnC(treeReader, "TreeMu1.mu_innerTrack_normalizedChi2");
    TTreeReaderValue<double> reader_mu2_iTnC(treeReader, "TreeMu2.mu_innerTrack_normalizedChi2");
    TTreeReaderValue<double> reader_mu3_iTnC(treeReader, "TreeMu3.mu_innerTrack_normalizedChi2");
    TTreeReaderValue<double> reader_mu1_tKink(treeReader, "TreeMu1.mu_combinedQuality_trkKink");
    TTreeReaderValue<double> reader_mu2_tKink(treeReader, "TreeMu2.mu_combinedQuality_trkKink");
    TTreeReaderValue<double> reader_mu3_tKink(treeReader, "TreeMu3.mu_combinedQuality_trkKink");
    TTreeReaderValue<double> reader_fv_nC(treeReader, "fv_nC");

    TTreeReaderValue<double> reader_dxy1(treeReader, "dxy1");
    TTreeReaderValue<double> reader_dxyErr1(treeReader, "dxyErr1");
    TTreeReaderValue<double> reader_dxy2(treeReader, "dxy2");
    TTreeReaderValue<double> reader_dxyErr2(treeReader, "dxyErr2");
    TTreeReaderValue<double> reader_dxy3(treeReader, "dxy3");
    TTreeReaderValue<double> reader_dxyErr3(treeReader, "dxyErr3");

    TTreeReaderValue<double> reader_SVx(treeReader, "SVx");
    TTreeReaderValue<double> reader_SVy(treeReader, "SVy");
    TTreeReaderValue<double> reader_SVz(treeReader, "SVz");
    TTreeReaderValue<double> reader_PVx(treeReader, "RefVx1");
    TTreeReaderValue<double> reader_PVy(treeReader, "RefVy1");
    TTreeReaderValue<double> reader_PVz(treeReader, "RefVz1");

    TTreeReaderValue<double> reader_dimu12(treeReader, "dimu12_ref");
    TTreeReaderValue<double> reader_dimu13(treeReader, "dimu13_ref");
    TTreeReaderValue<double> reader_dimu23(treeReader, "dimu23_ref");

    TTreeReaderValue<double> reader_fv_d3Dsig(treeReader, "fv_d3Dsig");
    TTreeReaderValue<double> reader_bs_sv_d2Dsig(treeReader, "bs_sv_d2Dsig");

    TTreeReaderValue<float> reader_muid1(treeReader, "MuonIDeval_Mu1.MuonID");
    TTreeReaderValue<float> reader_muid2(treeReader, "MuonIDeval_Mu2.MuonID");
    TTreeReaderValue<float> reader_muid3(treeReader, "MuonIDeval_Mu3.MuonID");
    TTreeReaderValue<int> reader_year(treeReader, "MuonIDeval_Mu3.year");

    TTreeReaderValue<double> reader_evt(treeReader, "evt");

    //Branches for output tree
    int nvar = 36;
    std::vector<TTreeReaderValue<double>> reader_out;
    for(auto v = 0; v<nvar-11; v++){
        reader_out.emplace_back(treeReader, outvar_name[v]);
        cout<<"   reading from TChain branch "<<outvar_name[v]<<endl;
    }

    //Branches for MVA
    std::vector<TTreeReaderValue<double>> reader_spec;
    for(int k = 0; k<n_spec; k++){
        cout<<"   reading from TChain branch "<<var_spec_name.at(k)<<endl;
        reader_spec.emplace_back(treeReader, var_spec_name.at(k));
    }
    std::vector<TTreeReaderValue<double>> reader_train;
    cout<<" -> variables for training"<<endl;
    for(int k = 0; k<n_train; k++){
        if(var_train_name.at(k)=="d0_min") continue;
        if(var_train_name.at(k)=="d0sig_min") continue;
        if(var_train_name.at(k)=="d0sig_max") continue;
        if(var_train_name.at(k)=="MuonIDeval_Mu1.MuonID") continue;
        if(var_train_name.at(k)=="MuonIDeval_Mu2.MuonID") continue;
        if(var_train_name.at(k)=="MuonIDeval_Mu3.MuonID") continue;
        if(var_train_name.at(k)=="MuonIDeval_Mu3.year")   continue;
        if(var_train_name.at(k)=="l1double_DoubleMu0_fired") continue;
        if(var_train_name.at(k)=="l1double_DoubleMu4_fired") continue;
        if(var_train_name.at(k)=="l1triple_fired") continue;
        cout<<"   reading from TChain branch "<<var_train_name.at(k)<<endl;
        reader_train.emplace_back(treeReader, var_train_name.at(k));
    }

    double phi_low, phi_high;
    //Loop on input Tree
    while (treeReader.Next()) {

        bool passedPresel = false;
        if(category=="A") { phi_low = 0.998; phi_high = 1.040; }
        if(category=="B") { phi_low = 0.991; phi_high = 1.047; }
        if(category=="C") { phi_low = 0.984; phi_high = 1.054; }
        //if(category=="A") { phi_low = 0.994; phi_high = 1.044; }
        //if(category=="B") { phi_low = 0.985; phi_high = 1.053; }
        //if(category=="C") { phi_low = 0.974; phi_high = 1.064; }
        if (
            ( (*reader_dimu12<phi_high && *reader_dimu12>phi_low) ||
              (*reader_dimu13<phi_high && *reader_dimu13>phi_low) ||
              (*reader_dimu23<phi_high && *reader_dimu23>phi_low) ) 
                                      ) continue;
        //excluding pT<3.5 GeV and |eta|<1.2
        if (
            ( (*reader_mu1_pt<3.5 && abs(*reader_mu1_eta)<1.2) ||
              (*reader_mu2_pt<3.5 && abs(*reader_mu2_eta)<1.2) ||
              (*reader_mu3_pt<3.5 && abs(*reader_mu3_eta)<1.2) ) 
                                      ) continue;

        if (TMath::IsNaN(*reader_fv_d3Dsig) ) {cout<<"warning TMath::IsNaN(*reader_fv_d3Dsig) in mass bin "<<*reader_triplMass<<endl; continue;}
        if (( *reader_fv_nC>0 && *reader_fv_nC<100) &&
            ( *reader_mu1_tLWM>8 ) &&
            ( *reader_mu2_tLWM>8 ) &&
            ( *reader_mu3_tLWM>8 ) &&
            ( *reader_mu1_validFra>0.5 ) &&
            ( *reader_mu2_validFra>0.5 ) &&
            ( *reader_mu3_validFra>0.5 ) &&
            ( *reader_mu1_iTnC<40 )&&
            ( *reader_mu2_iTnC<40 )&&
            ( *reader_mu3_iTnC<40 )&&
            ( *reader_mu1_tKink<900 ) &&
            ( *reader_mu2_tKink<900 ) &&
            ( *reader_mu3_tKink<900 ) ) passedPresel = true;
        if(!passedPresel) continue;

        for(int k = 0; k<n_spec; k++){
            var_spec.at(k)= *reader_spec.at(k);
        }
        evt = int(*reader_evt)%4096;
        for(int k = 0; k<n_train; k++){
            if(var_train_name.at(k)=="d0sig_min"){ 
                if((*reader_dxyErr1)!=0 && (*reader_dxyErr2)!=0  && (*reader_dxyErr3)!=0) var_train.at(k) = std::fmin ( abs(*reader_dxy1)/(*reader_dxyErr1) , std::fmin ( abs(*reader_dxy2)/(*reader_dxyErr2), abs(*reader_dxy3)/(*reader_dxyErr3) ) );
                else var_train.at(k) = 0;
            }
            else if(var_train_name.at(k)=="d0sig_max"){ 
                if((*reader_dxyErr1)!=0 && (*reader_dxyErr2)!=0  && (*reader_dxyErr3)!=0) var_train.at(k) = std::fmax ( abs(*reader_dxy1)/(*reader_dxyErr1) , std::fmax ( abs(*reader_dxy2)/(*reader_dxyErr2), abs(*reader_dxy3)/(*reader_dxyErr3) ) );
                else var_train.at(k) = 0;
            }
            else if(var_train_name.at(k)=="d0_min"){ 
                if((*reader_dxyErr1)!=0 && (*reader_dxyErr2)!=0  && (*reader_dxyErr3)!=0) var_train.at(k) = std::fmin ( abs(*reader_dxy1) , std::fmin ( abs(*reader_dxy2), abs(*reader_dxy3) ));
                else var_train.at(k) = 0;
            }
            else if(var_train_name.at(k)=="MuonIDeval_Mu1.MuonID") var_train.at(k) = *reader_muid1;
            else if(var_train_name.at(k)=="MuonIDeval_Mu2.MuonID") var_train.at(k) = *reader_muid2;
            else if(var_train_name.at(k)=="MuonIDeval_Mu3.MuonID") var_train.at(k) = *reader_muid3;
            else if(var_train_name.at(k)=="MuonIDeval_Mu3.year")   var_train.at(k) = *reader_year;
            else var_train.at(k) = *reader_train.at(k);
        }
        //Evaluate method(s) and fill histogram or MiniTree
        auto BDTscore = MVAreader->EvaluateMVA(method);
        auto reso = *reader_triplMassReso;
        auto mass = *reader_triplMass;
        Int_t isSB = 0;
        bool dM4_excl = !(*reader_tripleMu) && !(*reader_doubleMu0) && (*reader_doubleMu4);
        if( category == "A" && ( ( mass >= 1.62 && mass <= 1.746 ) || ( mass >= 1.794 && mass <= 2.0 ) )) isSB = 1;
        if( category == "B" && ( ( mass >= 1.62 && mass <= 1.732 ) || ( mass >= 1.808 && mass <= 2.0 ) )) isSB = 1;
        if( category == "C" && ( ( mass >= 1.62 && mass <= 1.720 ) || ( mass >= 1.820 && mass <= 2.0 ) )) isSB = 1;

        Double_t Ds_SF = Ds_correction; //0.9 by default
        if (dM4_excl and year=="2018") Ds_SF = Ds_correction_dM4_2018; //0.7 for events excl. triggered by dM4

        if(((category == "A" && reso < 0.007 ) ||
            (category == "B" && reso>= 0.007 && reso <= 0.0105 ) ||
            (category == "C" && reso > 0.0105 ))){
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
            for(auto v = nvar-11; v<nvar; v++){
                if(outvar_name[v] == "l1double_DoubleMu4_fired") outvar_value[v] = *reader_doubleMu4;
                if(outvar_name[v] == "l1double_DoubleMu0_fired") outvar_value[v] = *reader_doubleMu0;
                if(outvar_name[v] == "l1triple_fired") outvar_value[v] = *reader_tripleMu;
                if(outvar_name[v] == "MuonIDeval_Mu1.MuonID") outvar_value[v] = *reader_muid1;
                if(outvar_name[v] == "MuonIDeval_Mu2.MuonID") outvar_value[v] = *reader_muid2;
                if(outvar_name[v] == "MuonIDeval_Mu3.MuonID") outvar_value[v] = *reader_muid3;
                if(outvar_name[v] == "MuonIDeval_Mu3.year")   outvar_value[v] = *reader_year;
                if(outvar_name[v] == "bdt")  outvar_value[v] = BDTscore;
                if(outvar_name[v] == "isMC") outvar_value[v] = isMC;
                if(outvar_name[v] == "isSB") outvar_value[v] = isSB;
                if(outvar_name[v] == "weight" && isMC == 0) outvar_value[v] = 1.; //data
                if(outvar_name[v] == "weight" && isMC == 1) outvar_value[v] = wNormDs * Ds_SF * Dplus_correction; //Ds
                if(outvar_name[v] == "weight" && isMC == 2) outvar_value[v] = wNormB0 * Ds_SF * Bs_correction; //B0
                if(outvar_name[v] == "weight" && isMC == 3) outvar_value[v] = wNormBp * Ds_SF * Bs_correction; //Bp
                if(outvar_name[v] == "weight" && isMC == 4) outvar_value[v] = wNormW; //W, no corrections here
                if(outvar_name[v] == "category" && category == "A") outvar_value[v] = 0;
                if(outvar_name[v] == "category" && category == "B") outvar_value[v] = 1;
                if(outvar_name[v] == "category" && category == "C") outvar_value[v] = 2;
                //if(isMC){
                //     if(outvar_name[v] == "rew_Da_1") outvar_value[v] = reader_rew_Da[0];
                //     if(outvar_name[v] == "rew_Da_2") outvar_value[v] = reader_rew_Da[1];
                //     if(outvar_name[v] == "rew_Da_3") outvar_value[v] = reader_rew_Da[2];
                //     if(outvar_name[v] == "rew_Da_4") outvar_value[v] = reader_rew_Da[3];
                //     if(outvar_name[v] == "rew_Da_5") outvar_value[v] = reader_rew_Da[4];
                //}
                //else{
                //    if(outvar_name[v] == "rew_Da_1" || outvar_name[v] == "rew_Da_2" || outvar_name[v] == "rew_Da_3" || outvar_name[v] == "rew_Da_4" || outvar_name[v] == "rew_Da_5" ) outvar_value[v] = 1.0;
                //}
                //cout<<"   to output tree -> "<<outvar_name[v]<<"="<<outvar_value[v]<<endl;
            }
            count++;
            tout->Fill();
        }

    }
    cout<<"counted "<<count<<endl;
    delete MVAreader;
}


void evaluate_fillbdt_CV(TString year) 
{
    int nrun = 0;
    if(year == "2017") nrun = sizeof(inputpath_datarun_2017)/sizeof(inputpath_datarun_2017[0]);
    if(year == "2018") nrun = sizeof(inputpath_datarun_2018)/sizeof(inputpath_datarun_2018[0]);

    TString inputpath_datarun[nrun];
    TString inputpath_Ds = "";
    TString inputpath_B0 = "";
    TString inputpath_Bp = "";
    TString inputpath_W  = "";

    if(year == "2017"){
        for(int i=0; i<nrun; i++) inputpath_datarun[i] = inputpath_datarun_2017[i];
        inputpath_Ds      = inputpath_Ds_2017;
        inputpath_B0      = inputpath_B0_2017;
        inputpath_Bp      = inputpath_Bp_2017;
        inputpath_W       = inputpath_W_2017;
    } else if(year == "2018"){
        for(int i=0; i<nrun; i++) inputpath_datarun[i] = inputpath_datarun_2018[i];
        inputpath_Ds      = inputpath_Ds_2018;
        inputpath_B0      = inputpath_B0_2018;
        inputpath_Bp      = inputpath_Bp_2018;
        inputpath_W       = inputpath_W_2018;
    }
        
    Double_t wNormDs = 0, Ds_correction = 0, Dplus_correction = 0, wNormB0 = 0, wNormBp = 0, Bs_correction = 0, wNormW = 0;
    if(year == "2017"){
        wNormDs = wNormDs_2017;
        wNormB0 = wNormB0_2017;
        wNormBp = wNormBp_2017;
        wNormW = wNormW_2017;
        
        Ds_correction = Ds_correction_2017;
        Dplus_correction = Dplus_correction_2017;
        Bs_correction = Bs_correction_2017;
    } else if(year == "2018"){
        wNormDs = wNormDs_2018;
        wNormB0 = wNormB0_2018;
        wNormBp = wNormBp_2018;
        wNormW = wNormW_2018;
        
        Ds_correction = Ds_correction_2018;
        Dplus_correction = Dplus_correction_2018;
        Bs_correction = Bs_correction_2018;
    }
    int ncat = sizeof(cat_label)/sizeof(cat_label[0]);
    //open input files
    //data
    TChain *tdata = new TChain("FinalTree");
    for(auto i=0; i<nrun; i++){
        tdata->Add(inputpath_datarun[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC Ds
    TChain *tmc1 = new TChain("FinalTree");
    tmc1->Add(inputpath_Ds); 
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    //MC B0
    TChain *tmc2 = new TChain("FinalTree");
    tmc2->Add(inputpath_B0); 
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    //MC Bp
    TChain *tmc3 = new TChain("FinalTree");
    tmc3->Add(inputpath_Bp); 
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;
    //MC W
    TChain *tmc4 = new TChain("FinalTree");
    tmc4->Add(inputpath_W); 
    std::cout<<"Opened input file: "<<inputpath_W<<std::endl;

    //data
    TChain *tdata_mu1 = new TChain("TreeMu1");
    TChain *tdata_mu2 = new TChain("TreeMu2");
    TChain *tdata_mu3 = new TChain("TreeMu3");
    for(auto i=0; i<nrun; i++){
        tdata_mu1->Add(inputpath_datarun[i]); 
        tdata_mu2->Add(inputpath_datarun[i]); 
        tdata_mu3->Add(inputpath_datarun[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC Ds
    TChain *tmc1_mu1 = new TChain("TreeMu1");
    TChain *tmc1_mu2 = new TChain("TreeMu2");
    TChain *tmc1_mu3 = new TChain("TreeMu3");
    tmc1_mu1->Add(inputpath_Ds); 
    tmc1_mu2->Add(inputpath_Ds); 
    tmc1_mu3->Add(inputpath_Ds); 
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    //MC B0
    TChain *tmc2_mu1 = new TChain("TreeMu1");
    TChain *tmc2_mu2 = new TChain("TreeMu2");
    TChain *tmc2_mu3 = new TChain("TreeMu3");
    tmc2_mu1->Add(inputpath_B0); 
    tmc2_mu2->Add(inputpath_B0); 
    tmc2_mu3->Add(inputpath_B0); 
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    //MC Bp
    TChain *tmc3_mu1 = new TChain("TreeMu1");
    TChain *tmc3_mu2 = new TChain("TreeMu2");
    TChain *tmc3_mu3 = new TChain("TreeMu3");
    tmc3_mu1->Add(inputpath_Bp); 
    tmc3_mu2->Add(inputpath_Bp); 
    tmc3_mu3->Add(inputpath_Bp); 
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;
    //MC W
    TChain *tmc4_mu1 = new TChain("TreeMu1");
    TChain *tmc4_mu2 = new TChain("TreeMu2");
    TChain *tmc4_mu3 = new TChain("TreeMu3");
    tmc4_mu1->Add(inputpath_W); 
    tmc4_mu2->Add(inputpath_W); 
    tmc4_mu3->Add(inputpath_W); 
    std::cout<<"Opened input file: "<<inputpath_W<<std::endl;


    //data
    TChain *tdata_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tdata_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tdata_muid3 = new TChain("MuonIDeval_Mu3");
    for(auto i=0; i<nrun; i++){
        TString inputpath_data_muId = inputpath_datarun[i].ReplaceAll(".root", "_MuonID.root");
        tdata_muid1->Add(inputpath_datarun[i]);
        tdata_muid2->Add(inputpath_datarun[i]);
        tdata_muid3->Add(inputpath_datarun[i]);
        std::cout<<"Opened input file: "<<inputpath_data_muId<<std::endl;
    }
    //MC Ds
    TChain *tmc1_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc1_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc1_muid3 = new TChain("MuonIDeval_Mu3");
    TString inputpath_Ds_muId = inputpath_Ds.ReplaceAll(".root", "_MuonID.root");
    tmc1_muid1->Add(inputpath_Ds_muId);
    tmc1_muid2->Add(inputpath_Ds_muId);
    tmc1_muid3->Add(inputpath_Ds_muId);
    std::cout<<"Opened input file: "<<inputpath_Ds_muId<<std::endl;
    //MC B0
    TChain *tmc2_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc2_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc2_muid3 = new TChain("MuonIDeval_Mu3");
    TString inputpath_B0_muId = inputpath_B0.ReplaceAll(".root", "_MuonID.root");
    tmc2_muid1->Add(inputpath_B0_muId);
    tmc2_muid2->Add(inputpath_B0_muId);
    tmc2_muid3->Add(inputpath_B0_muId);
    std::cout<<"Opened input file: "<<inputpath_B0_muId<<std::endl;
    //MC Bp
    TChain *tmc3_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc3_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc3_muid3 = new TChain("MuonIDeval_Mu3");
    TString inputpath_Bp_muId = inputpath_Bp.ReplaceAll(".root", "_MuonID.root");
    tmc3_muid1->Add(inputpath_Bp_muId);
    tmc3_muid2->Add(inputpath_Bp_muId);
    tmc3_muid3->Add(inputpath_Bp_muId);
    std::cout<<"Opened input file: "<<inputpath_Bp_muId<<std::endl;
    //MC W
    TChain *tmc4_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc4_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc4_muid3 = new TChain("MuonIDeval_Mu3");
    TString inputpath_W_muId = inputpath_W.ReplaceAll(".root", "_MuonID.root");
    tmc4_muid1->Add(inputpath_W_muId);
    tmc4_muid2->Add(inputpath_W_muId);
    tmc4_muid3->Add(inputpath_W_muId);
    std::cout<<"Opened input file: "<<inputpath_W_muId<<std::endl;

    tdata->AddFriend(tdata_mu1);
    tdata->AddFriend(tdata_mu2);
    tdata->AddFriend(tdata_mu3);

    tmc1->AddFriend(tmc1_mu1);
    tmc1->AddFriend(tmc1_mu2);
    tmc1->AddFriend(tmc1_mu3);

    tmc2->AddFriend(tmc2_mu1);
    tmc2->AddFriend(tmc2_mu2);
    tmc2->AddFriend(tmc2_mu3);

    tmc3->AddFriend(tmc3_mu1);
    tmc3->AddFriend(tmc3_mu2);
    tmc3->AddFriend(tmc3_mu3);

    tmc4->AddFriend(tmc4_mu1);
    tmc4->AddFriend(tmc4_mu2);
    tmc4->AddFriend(tmc4_mu3);

    tdata->AddFriend(tdata_muid1);
    tdata->AddFriend(tdata_muid2);
    tdata->AddFriend(tdata_muid3);

    tmc1->AddFriend(tmc1_muid1);
    tmc1->AddFriend(tmc1_muid2);
    tmc1->AddFriend(tmc1_muid3);

    tmc2->AddFriend(tmc2_muid1);
    tmc2->AddFriend(tmc2_muid2);
    tmc2->AddFriend(tmc2_muid3);

    tmc3->AddFriend(tmc3_muid1);
    tmc3->AddFriend(tmc3_muid2);
    tmc3->AddFriend(tmc3_muid3);
        
    tmc4->AddFriend(tmc4_muid1);
    tmc4->AddFriend(tmc4_muid2);
    tmc4->AddFriend(tmc4_muid3);

    cout<<"entries in input MC:"<<endl;
    cout<<"MC1 "<<tmc1->GetEntries()<<endl;
    cout<<"MC2 "<<tmc2->GetEntries()<<endl;
    cout<<"MC3 "<<tmc3->GetEntries()<<endl;
    cout<<"MC4 "<<tmc4->GetEntries()<<endl;
    
    //Book output MiniTree
    TString fout_tree_path = TMVA_inputpath+"_"+year+"_outputTree.root";
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
                                "Ptmu1",
                                "Etamu1",
                                "Ptmu2",
                                "Etamu2",
                                "Ptmu3",
                                "Etamu3",   //12

                                "fv_nC",
                                "fv_dphi3D",
                                "fv_d3Dsig",
                                "fv_d3D",
                                "Pt_tripl",
                                "cLP",
                                "tKink",
                                "d0sig",
                                "mindca_iso",
                                "trkRel",
                                "segmComp", //23

                                "dimu_OS1",
                                "dimu_OS2",

                                "MuonIDeval_Mu1.MuonID",
                                "MuonIDeval_Mu2.MuonID",
                                "MuonIDeval_Mu3.MuonID",
                                "l1double_DoubleMu0_fired",
                                "l1double_DoubleMu4_fired",
                                "l1triple_fired", //31

                                "isMC",
                                "isSB",
                                "weight",
                                "category",
                                "bdt"
                                };

    int nvar = 36;
    Double_t outvar_value[36] = {0.};

    //Set branches output tree
    for(auto v = 0; v<nvar; v++){
        //TString outvar_name_mod = outvar_name[v];
        //outvar_name_mod = outvar_name_mod.ReplaceAll(".", "_");
        //tout->Branch(outvar_name_mod, &outvar_value[v]);
        tout->Branch(outvar_name[v], &outvar_value[v]);
    }

    //Loop on categories A, B, C
    for(auto i = 0; i<3; i++){
        TString category = cat_label[i];
        cout<<"Category "<<category<<endl;

        //Make sure TChain points to first event
        tdata->LoadTree(0);
        tmc1->LoadTree(0);
        tmc2->LoadTree(0);
        tmc3->LoadTree(0);
        tmc4->LoadTree(0);
  
        //Book output histogram or MiniTree
        TString fout_path = TMVA_inputpath+category+"/BDTdecision_"+category+"_v2.root";
        TFile *fout = new TFile(fout_path, "recreate");
        fout->cd();

        TH1F * hBDTdecision_data = new TH1F ("BDTdecision_data_obs"+category,"BDTdecision_data_obs"+category, 240, -0.6, 0.6);
        TH1F * hBDTdecision_Ds = new TH1F ("BDTdecision_signalDs"+category,"BDTdecision_signalDs"+category, 240, -0.6, 0.6);
        TH1F * hBDTdecision_B0 = new TH1F ("BDTdecision_signalB0"+category,"BDTdecision_signalB0"+category, 240, -0.6, 0.6);
        TH1F * hBDTdecision_Bp = new TH1F ("BDTdecision_signalBp"+category,"BDTdecision_signalBp"+category, 240, -0.6, 0.6);
        TH1F * hBDTdecision_W  = new TH1F ("BDTdecision_signalW"+category,"BDTdecision_signalW"+category, 240, -0.6, 0.6);

        //Read branches, loop on tree, perform evaluate and fill histo
        fill_BDT_score(tdata, hBDTdecision_data, category, 0, tout, outvar_name, outvar_value, year); 
        fill_BDT_score(tmc1,  hBDTdecision_Ds,   category, 1, tout, outvar_name, outvar_value, year); 
        fill_BDT_score(tmc2,  hBDTdecision_B0,   category, 2, tout, outvar_name, outvar_value, year); 
        fill_BDT_score(tmc3,  hBDTdecision_Bp,   category, 3, tout, outvar_name, outvar_value, year); 
        fill_BDT_score(tmc4,  hBDTdecision_W,    category, 4, tout, outvar_name, outvar_value, year); 

        //rescale signal
        hBDTdecision_Ds->Scale(wNormDs * Ds_correction * Dplus_correction);
        hBDTdecision_B0->Scale(wNormB0 * Ds_correction * Bs_correction);
        hBDTdecision_Bp->Scale(wNormBp * Ds_correction * Bs_correction);
        hBDTdecision_W->Scale(wNormW);

        //add signals to common histo
        TH1F * hBDTdecisionSignal = new TH1F ("BDTdecision_signal"+category, "BDTdecision_signal"+category, 240, -0.6, 0.6);
        hBDTdecisionSignal->Add(hBDTdecision_Ds);
        hBDTdecisionSignal->Add(hBDTdecision_B0);
        hBDTdecisionSignal->Add(hBDTdecision_Bp);
        hBDTdecisionSignal->Add(hBDTdecision_W);

        //Write and close the file
        fout->cd();
        hBDTdecision_data->Write();
        hBDTdecisionSignal->Write();
        hBDTdecision_Ds->Write();
        hBDTdecision_B0->Write();
        hBDTdecision_Bp->Write();
        hBDTdecision_W->Write();

        fout->Close();
        cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    }
    fout_tree->Write();
    fout_tree->Close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_tree_path<<"\n\n"<<endl;
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}

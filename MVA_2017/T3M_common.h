#include <iostream>

using namespace std;
    TString workdir = "/lustrehome/fsimone/MVA_2017/";

//XGBoost 
    TString XGBoost_tag = "";
//TMVA Training options
    TString cat_label[] = {"A", "B", "C"};

    TString TMVA_outputpath = "dataset_2017_16dec_vetoByCatPreBDT_cutOnDispl_addSample_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = false;
    TString method = "BDT";
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    //TString TMVA_inputpath = "dataset_2017_14oct_vetoPreBDT_";  //name to load TMVA results for evaluation
    //TString TMVA_inputpath = "dataset_2017_19nov_vetoByCatPreBDT_";  //name to load TMVA results for evaluation
    TString TMVA_inputpath = "dataset_2017_16dec_vetoByCatPreBDT_cutOnDispl_addSample_";  //name to load TMVA results for evaluation

//data rootfiles

    TString inputpath_datarun[] = {
           //updated bs_sv_d2D renamed
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201216_1208/AnalysedTree_data_2017B_tau3mu_16dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201216_1209/AnalysedTree_data_2017C_tau3mu_16dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201216_1210/AnalysedTree_data_2017D_tau3mu_16dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201216_1211/AnalysedTree_data_2017E_tau3mu_16dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201216_1212/AnalysedTree_data_2017F_tau3mu_16dec.root",
           };

    //no vetoes + PU reweighting
    TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201215_1527/AnalysedTree_MC_2017Ds_tau3mu_15dec_mergedOffJian.root"; //merged official Ds+Jian's sample
    //TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1512/AnalysedTree_MC_2017Ds_tau3mu_25nov.root";
    TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1533/AnalysedTree_MC_2017B0_tau3mu_25nov.root";
    TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1533/AnalysedTree_MC_2017Bp_tau3mu_25nov.root";

//Coefficients for signal normalisation
    Double_t Dplus_correction = 1.05; // to be applied to D signal  
    Double_t Bs_correction = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction = 0.775; //2017 Ds correction factor updated on 28 November

    Double_t Lumi_data[] = {4.79, 9.63, 4.24, 9.3, 10.04};
    Double_t wNormDs = 7.203E-04; // Ds 3.67+2.65 initial events
    //Double_t wNormDs = 1.242E-03; // Ds 3.67M initial events
    Double_t wNormB0 = 4.160E-04; // B0 3.00M initial events
    Double_t wNormBp = 6.203E-04; // Bp 2.01M initial events
    Double_t sig_norm = 0;


//TMVA settings
    TString common_cut = "(fv_nC>0 && bs_sv_d2Dsig>=3.75)";

  //  TString common_cut = "( "
  //                       " (fv_nC>0) &&"
  //                       "!(dimu12>0.9975 && dimu12<1.0415) && "
  //                       "!(dimu13>0.9975 && dimu13<1.0415) && "
  //                       "!(dimu23>0.9975 && dimu23<1.0415) "
  //                       " )";  

    TString BDTinVar_A = "BDTinputVar_v3_muId.txt";
    TString BDTinVar_B = "BDTinputVar_v3_muId.txt";
    TString BDTinVar_C = "BDTinputVar_v3_muId.txt";
    TString BDTspecVar = "BDTspecVar_2017.txt";

//Utility to read variables from text file
void readVarName(std::vector<TString> &var_name, std::vector<TString> &var_def,  TString filename ){
     TString var[2];
     ifstream inputVar;
     inputVar.open(filename);
     while (!inputVar.fail() && !inputVar.eof()){
         inputVar >> var[0] >> var[1];
         if(! (var[0]=="" || var[1]=="" )){
             var_name.push_back(var[0]);
             var_def.push_back(var[1]);
         }
     }
     inputVar.close();
}

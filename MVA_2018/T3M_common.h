#include <iostream>
#include <fstream>

using namespace std;
TString workdir = "/lustrehome/fsimone/MVA_2018/"; 

//XGBoost 
    TString XGBoost_tag = ""; //9june2020_optimisedvar_5fold";


//TMVA Training options

    //TString cat_label[] = {"Alowpt", "Ahighpt", "Blowpt", "Bhighpt", "C"};
    TString cat_label[] = {"A", "B", "C"};

    TString TMVA_outputpath = "dataset_2018_17dec_addSample_"; //dataset_2018_6may_opt_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = false;
    TString method = "BDT";
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
 
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath = "dataset_2018_17dec_addSample_"; //dataset_2018_6may_opt_";  //name to load TMVA results for evaluation

//data rootfiles

    //no omega no phi vetos
    TString inputpath_datarun[] = {
    //no omega no phi vetos //added dimu info in finalTree //bs_sv_d2Dsig>2
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201217_1432/AnalysedTree_data_2018A_tau3mu_17dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201217_1433/AnalysedTree_data_2018B_tau3mu_17dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201217_1434/AnalysedTree_data_2018C_tau3mu_17dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201217_1435/AnalysedTree_data_2018D_tau3mu_17dec.root",
           };
    //additional signal samples
    //no omega no phi vetos //added dimu info in finalTree //bs_sv_d2Dsig>2 //additional Ds and Bp samples
    TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201215_1529/AnalysedTree_MC_2018Ds_tau3mu_15dec_addSample_merged.root";
    TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201217_1437/AnalysedTree_MC_2018B0_tau3mu_17dec.root";
    TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201218_1737/AnalysedTree_MC_2018Bp_tau3mu_18dec_addSample_merged.root";

//Coefficients for signal normalisation
    Double_t Lumi_data[] = {13.98, 7.06, 6.90, 31.75};
    Double_t Dplus_correction = 1.05; // to be applied to D signal  
    Double_t Bs_correction = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction = 0.93; //2018 updated on 22 march

    Double_t wNormDs = 7.58E-04; //additional sample 3.4M + new 4.6E+06 MC DsTau
    //Double_t wNormDs = 1.323E-03; //new 4.6E+06 MC DsTau
  //        Double_t wNormDs = 7.15E-03; //old 2.13E+06 MC DsTau
    Double_t wNormB0 = 4.784E-04; //new 3.5E+06 MC B0Tau
  //        Double_t wNormB0 = 3.03E-03;  //old 0.655E+06 MC B0Tau
    Double_t wNormBp = 5.37E-04; //additional sample 2.23M + usual 1.33E+06 MC BpTau
    //Double_t wNormBp = 1.437E-03; //usual 1.33E+06 MC BpTau
    Double_t sig_norm = 0;
//(wNormDs+wNormB0+wNormBp) * Dplus_correction / 3.; //average normalization factor for the three signal samples


//TMVA settings
    //TString common_cut = "( "
    //                     " (fv_nC>0) &&"
    //                     "!(dimu12>0.9975 && dimu12<1.0415) && "
    //                     "!(dimu13>0.9975 && dimu13<1.0415) && "
    //                     "!(dimu23>0.9975 && dimu23<1.0415) "
    //                     " )";  
    TString common_cut = "( fv_nC>0 )";

    TString BDTinVar_A = "BDTinputVar_dataset_2018_4Nov_A.txt"; 
    TString BDTinVar_B = "BDTinputVar_dataset_2018_4Nov_B.txt"; 
    TString BDTinVar_C = "BDTinputVar_dataset_2018_4Nov_C.txt"; 
    TString BDTspecVar = "BDTspecVar_12nov.txt";

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

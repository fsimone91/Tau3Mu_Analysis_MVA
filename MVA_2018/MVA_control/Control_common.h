#include <iostream>

using namespace std;
    TString workdir_control = "/lustrehome/fsimone/MVA_2018/MVA_control/";
 
//TMVA Training options
    TString TMVA_outputpath_control = "MVA_control_2018_16nov_chi2cut"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    TString method_control = "BDT";
    TString TMVA_weightfilename_control = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath_control = "MVA_control_2018_16nov_chi2cut";  //name to load TMVA results for evaluation

//data rootfiles

    TString inputpath_datarun_control[] = {
           //DsPhiPi 2018 Data - no cut on SV chi2 - isGlobal+isPF - prescaled HLT //only DoubleMu0
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2155/AnalysedTree_data_2018A_control_14oct.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2156/AnalysedTree_data_2018B_control_14oct.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2157/AnalysedTree_data_2018C_control_14oct.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2158/AnalysedTree_data_2018D_control_14oct.root"
           ////DsPhiPi 2018 Data - no cut on SV chi2 - isGlobal+isPF - prescaled HLT
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1033/AnalysedTree_data_2018A_control_5ott.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1033/AnalysedTree_data_2018B_control_5ott.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1558/AnalysedTree_data_2018C_control_5ott.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1033/AnalysedTree_data_2018D_control_5ott.root"
           };

    //DsPhiPi 2018 MC - no cut on SV chi2 - isGlobal+isPF - prescaled HLT //only DoubleMu0
    TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2154/AnalysedTree_MC_2018DsPhiPi_control_14oct.root";
    ////DsPhiPi 2018 MC - no cut on SV chi2 - isGlobal+isPF
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1032/AnalysedTree_MC_2018DsPhiPi_control_5ott.root";


//Coefficients for signal normalisation
// Norm = Lumi_data[i]*xsection_mc*BR/N_MC;
    TString run_name_control[] = {"2018A", "2018B", "2018C", "2018D"};
    Double_t Lumi_data_control[] = {0.7025, 0.3522, 0.3458, 1.58803}; //recorded lumi by HLT_DoubleMu3_Trk_*
    //Double_t Lumi_data_control[] = {13.98, 7.06, 6.90, 31.75};
    Double_t xsection_mc = 2.32e10; //Ds Production Cross section
    int N_MC = 2050462;  //rereco Total number of events in MC sample
    Double_t BR = 1.29e-5;  //Branching ratio Ds to Phi Pi


//TMVA settings
    //common preselections
    TCut common_cut_control = "(fv_nC > 0 && fv_nC < 10)";
    //input variables
    TString BDTinVar_control = "../BDTinputVar_2016.txt";
    TString BDTspecVar_control = "../BDTspecVar_2018_noevt.txt";
//Utility to read variables from text file
void readVarName_control(std::vector<TString> &var_name, std::vector<TString> &var_def,  TString filename ){
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

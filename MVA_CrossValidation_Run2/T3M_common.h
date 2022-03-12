#include <iostream>

using namespace std;
    TString workdir = "/lustrehome/fsimone/MVA_CrossValidation_Run2/";

//XGBoost 
    TString XGBoost_tag = "";
//TMVA Training options
    TString cat_label[] = {"A", "B", "C"};

    TString TMVA_outputpath = "training_run2_CV_ntrees500depth3beta0p5cuts50_5folds_";
    //change it to perform 5-fold Cross Validation
    bool doCV = true;
    UInt_t numFolds = 5;
    TString method = "BDTG";
    TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath = "training_run2_CV_ntrees500depth3beta0p5cuts50_5folds_";

//data rootfiles

    TString inputpath_datarun_2017[] = {
           //UL datasets
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1549/AnalysedTree_data_2017B_tau3mu_28july.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1550/AnalysedTree_data_2017C_tau3mu_28july.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1551/AnalysedTree_data_2017D_tau3mu_28july.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1552/AnalysedTree_data_2017E_tau3mu_28july.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1553/AnalysedTree_data_2017F_tau3mu_28july.root",
           };
    TString inputpath_datarun_2018[] = {
           //UL - Added nPV
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1615/AnalysedTree_data_2017F_tau3mu_14dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1544/AnalysedTree_data_2018A_tau3mu_14dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1545/AnalysedTree_data_2018B_tau3mu_14dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1546/AnalysedTree_data_2018C_tau3mu_14dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1547/AnalysedTree_data_2018D_tau3mu_14dec.root",
           };
    //UL MC
    TString inputpath_Ds_2017 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1602/AnalysedTree_MC_2017Ds_tau3mu_30july.root";
    TString inputpath_B0_2017 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1603/AnalysedTree_MC_2017B0_tau3mu_30july.root";
    TString inputpath_Bp_2017 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1604/AnalysedTree_MC_2017Bp_tau3mu_30july.root";
    TString inputpath_W_2017 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_2305/AnalysedTree_MC_2017W_tau3mu_28july.root";

    TString inputpath_Ds_2018 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1539/AnalysedTree_MC_2018Ds_tau3mu_14dec.root";
    TString inputpath_B0_2018 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1540/AnalysedTree_MC_2018B0_tau3mu_14dec.root";
    TString inputpath_Bp_2018 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20211214_1541/AnalysedTree_MC_2018Bp_tau3mu_14dec.root";
    TString inputpath_W_2018  = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1831/AnalysedTree_MC_2018W_tau3mu_30july.root";

    //Coefficients for signal normalisation
    Double_t Dplus_correction_2017 = 1.05; // to be applied to D signal  
    Double_t Bs_correction_2017 = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction_2017 = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction_2017 = 0.89; //UL 2017 Ds correction factor updated Feb 2021
    Double_t Lumi_data_2017[] = {4.79, 9.63, 4.24, 9.3, 9.89};

    Double_t wNormDs_2017 = 5.344E-04; // Ds 7.23M initial events
    Double_t wNormB0_2017 = 4.991E-04; //5.256E-04; // B0 1.95M initial events
    Double_t wNormBp_2017 = 5.096E-04; //5.259E-04; // Bp 1.87M initial events
    Double_t wNormW_2017  = 1.611E-04; //W 0.5M initial events initial events

    Double_t Dplus_correction_2018 = 1.05; // to be applied to D signal  
    Double_t Bs_correction_2018 = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction_2018 = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction_2018 = 0.91; //2018 updated feb 2021
    Double_t Ds_correction_dM4_2018 = 0.81; //2018 events excl triggered by DoubleMu4
    Double_t Lumi_data_2018[] = {3.67, 13.98, 7.06, 6.90, 31.75};

    Double_t wNormDs_2018 = 5.70E-04; //using 63.36 tot lumi 5.37E-04; //Ds UL
    Double_t wNormB0_2018 = 5.42E-04; //using 63.36 tot lumi 5.1093E-04; //5.36E-04; //B0 UL
    Double_t wNormBp_2018 = 5.58E-04; //using 63.36 tot lumi 5.2549E-04; //5.40E-04; //Bp UL
    Double_t wNormW_2018 = 3.05E-04;

    TString common_cut_2017 = "(fv_nC>0 && bs_sv_d2Dsig>=3.75 && !(TMath::IsNaN(fv_d3Dsig)))";
    TString common_cut_2018 = "(fv_nC>0 && !(TMath::IsNaN(fv_d3Dsig)) )";
     //&& !(Ptmu1<3.5 && Etamu1<1.2) && !(Ptmu2<3.5 && Etamu2<1.2) && !(Ptmu3<3.5 && Etamu3<1.2))"; //&& MuonIDeval_Mu1.MuonID>0 && MuonIDeval_Mu2.MuonID>0 && MuonIDeval_Mu3.MuonID>0)";

    TString BDTspecVar = "BDTspecVar_2017.txt";
    TString BDTinVar_A = "BDTinputVar_dataset_Run2.txt"; 
    TString BDTinVar_B = "BDTinputVar_dataset_Run2_dropCorr1.txt"; 
    TString BDTinVar_C = "BDTinputVar_dataset_Run2_dropCorr1.txt"; 

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

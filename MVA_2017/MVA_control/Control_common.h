#include <iostream>

using namespace std;
    TString workdir_control = "/lustrehome/fsimone/MVA_2017/MVA_control/";

//TMVA Training options
    TString TMVA_outputpath_control = "MVA_control_2017_baseline_PU_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    TString method_control = "BDT";
    TString TMVA_weightfilename_control = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath_control = "MVA_control_2017_baseline_PU_";  //name to load TMVA results for evaluation

//data rootfiles

    TString inputpath_datarun_control[] = {
 //          //DsPhiPi 2017 Data - no cut on SV chi2 - isGlobal+isPF
 //          "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1603/AnalysedTree_data_2017B_control_25nov.root",
 //          "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1604/AnalysedTree_data_2017C_control_25nov.root",
 //          "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1605/AnalysedTree_data_2017D_control_25nov.root",
 //          "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1606/AnalysedTree_data_2017E_control_25nov.root",
 //          "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1607/AnalysedTree_data_2017F_control_25nov.root"
 //          };

 //   //DsPhiPi 2017 MC - no cut on SV chi2
 //   TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201125_1549/AnalysedTree_MC_2017DsPhiPi_control_25nov.root";

//           //DsPhiPi 2017 Data - no cut on SV chi2 - isGlobal+isPF
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200928_1812/AnalysedTree_data_2017B_control_28set.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200928_1813/AnalysedTree_data_2017C_control_28set.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200928_1814/AnalysedTree_data_2017D_control_28set.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200928_1815/AnalysedTree_data_2017E_control_28set.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200928_1816/AnalysedTree_data_2017F_control_28set.root"
//           };
//
//    //DsPhiPi 2017 MC - no cut on SV chi2
//    TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200928_1817/AnalysedTree_MC_2017DsPhiPi_control_28set.root";

           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201201_1850/AnalysedTree_data_2017B_control_1dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201201_1851/AnalysedTree_data_2017C_control_1dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201201_1852/AnalysedTree_data_2017D_control_1dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201201_1853/AnalysedTree_data_2017E_control_1dec.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201201_1854/AnalysedTree_data_2017F_control_1dec.root"
           };

    //DsPhiPi 2017 MC - no cut on SV chi2
    TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201201_1849/AnalysedTree_MC_2017DsPhiPi_control_1dec.root";

//Coefficients for signal normalisation
// Norm = Lumi_data[i]*xsection_mc*BR/N_MC;
    TString run_name_control[] = {"2017B", "2017C", "2017D", "2017E", "2017F"};
    Double_t Lumi_data_control[] = {4.79, 9.63, 4.24, 9.3, 10.04};
    Double_t xsection_mc = 1.17e10; //Ds Production Cross section
    //int N_MC = 414357;  //UL --> Total number of events in MC sample
    int N_MC = 1832740;  //rereco --> Total number of events in MC sample
    Double_t BR = 1.29e-5;  //Branching ratio Ds to Phi Pi


//TMVA settings
    //common preselections
    TCut common_cut_control = "";
    //input variables
    TString BDTinVar_control = "../BDTinputVar_baseline.txt";
    TString BDTspecVar_control = "../BDTspecVar_2017_noevt.txt";
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

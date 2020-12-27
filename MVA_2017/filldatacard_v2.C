#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "BDT_optimal_cut_v2.h"


void filldatacard_v2() 
{
    TString cat_label[] = {"A", "B", "C"};

    //Open output MiniTree
    TString fout_tree_path = TMVA_inputpath+"outputTree.root";
    TChain *ttree = new TChain("outputTree");
    ttree->Add(fout_tree_path); 
    cout<<"Reading final tree from file: "<<fout_tree_path<<endl;

    //Book output mass histograms
    TString fout_path = "datacardT3Mu_"+TMVA_inputpath+"v2.root";
    TFile *fout = new TFile(fout_path,"recreate");

    //Loop on categories A, B, C
    for(auto i = 0; i<3; i++){
        TString category = cat_label[i];
        cout<<"Category "<<category<<endl;

        //Make sure TChain points to firts event
        ttree->LoadTree(0);

        //set optimal categorisation
        BDTcut BDTcutvalues = Get_BDT_cut_v2(category, false);   

        //Book output mass histograms
        fout->cd();
        TH1F * hTriplMass_data1;// = new TH1F ("data_obs"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_data2;// = new TH1F ("data_obs"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);
        TH1F * hSignal1        ;// = new TH1F ("signal"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hSignal2        ;// = new TH1F ("signal"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);        

        TString c = ""; 
        if (category=="A") c = "0";
        if (category=="B") c = "1";
        if (category=="C") c = "2";
        //TString veto_phi = "";
        //if (category=="A") veto_phi = "!(dimu12>0.994 && dimu12<1.044) && !(dimu13>0.994 && dimu13<1.044) &&!(dimu23>0.994 && dimu23<1.044)";
        //if (category=="B") veto_phi = "!(dimu12>0.985 && dimu12<1.053) && !(dimu13>0.985 && dimu13<1.053) &&!(dimu23>0.985 && dimu23<1.053)";
        //if (category=="C") veto_phi = "!(dimu12>0.974 && dimu12<1.064) && !(dimu13>0.974 && dimu13<1.064) &&!(dimu23>0.974 && dimu23<1.064)";
        //TString veto = " && "+veto_phi; 
        TString veto = "";
        
        //bdt score distribution
        TString binning = "(42, 1.600000, 2.020000)"; 
        //data
        TString BDTa = to_string(BDTcutvalues.a);
        TString BDTb = to_string(BDTcutvalues.b);
        ttree->Draw("tripletMass>>hTriplMass_data1"+binning, "weight*(isMC==0  && category =="+c+" && bdt >= "+BDTa+veto+")");
        hTriplMass_data1 = (TH1F *)gDirectory->Get("hTriplMass_data1");
        ttree->Draw("tripletMass>>hTriplMass_data2"+binning, "weight*(isMC==0  && category =="+c+" && bdt < "+BDTa+" && bdt >= "+BDTb+veto+")");
        hTriplMass_data2 = (TH1F *)gDirectory->Get("hTriplMass_data2");
        //signal
        ttree->Draw("tripletMass>>hSignal1"+binning, "weight*(isMC>0  && category =="+c+" && bdt >= "+BDTa+veto+")");
        hSignal1 = (TH1F *)gDirectory->Get("hSignal1");
        ttree->Draw("tripletMass>>hSignal2"+binning, "weight*(isMC>0  && category =="+c+" && bdt < "+BDTa+" && bdt >= "+BDTb+veto+")");
        hSignal2 = (TH1F *)gDirectory->Get("hSignal2");


        hTriplMass_data1->SetDirectory(0);
        hTriplMass_data2->SetDirectory(0);
        hSignal1->SetDirectory(0);
        hSignal2->SetDirectory(0);

        //Write and close the file
        fout->cd();
        hTriplMass_data1->Write("data_obs"+category+"1");
        hTriplMass_data2->Write("data_obs"+category+"2");
        hTriplMass_data1->Write("background"+category+"1");
        hTriplMass_data2->Write("background"+category+"2");
        hSignal1->Write("signal"+category+"1");
        hSignal2->Write("signal"+category+"2");

    }
    fout->Close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    return 0;
}

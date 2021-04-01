#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "BDT_optimal_cut_v3.h"


void filldatacard_v3() 
{
    TString cat_label[] = {"A", "B", "C"};

    //Open output MiniTree
    TString fout_tree_path = TMVA_inputpath+"outputTree.root";
    TChain *ttree = new TChain("outputTree");
    ttree->Add(fout_tree_path); 
    cout<<"Reading final tree from file: "<<fout_tree_path<<endl;

    //Book output mass histograms
    TString fout_path = "datacardT3Mu_"+TMVA_inputpath+"v3.root";
    TFile *fout = new TFile(fout_path,"recreate");

    //Loop on categories A, B, C
    for(auto i = 0; i<3; i++){
        TString category = cat_label[i];
        cout<<"Category "<<category<<endl;

        //Make sure TChain points to firts event
        ttree->LoadTree(0);

        //set optimal categorisation
        BDTcut3d BDTcutvalues = Get_BDT_cut_3D(category);   

        //Book output mass histograms
        fout->cd();
        TH1F * hTriplMass_data1;// = new TH1F ("data_obs"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_data2;// = new TH1F ("data_obs"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);
        TH1F * hTriplMass_data3;// = new TH1F ("data_obs"+category+"3","Triplet mass "+category+"3",42, 1.600000, 2.020000);
        TH1F * hSignal1        ;// = new TH1F ("signal"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hSignal2        ;// = new TH1F ("signal"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);        
        TH1F * hSignal3        ;// = new TH1F ("signal"+category+"3","Triplet mass "+category+"3",42, 1.600000, 2.020000);        

        TString c = ""; 
        if (category=="A") c = "0";
        if (category=="B") c = "1";
        if (category=="C") c = "2";

        //bdt score distribution
        TString binning = "(42, 1.600000, 2.020000)"; 
        //data
        TString BDTa = to_string(BDTcutvalues.a);
        TString BDTb = to_string(BDTcutvalues.b);
        TString BDTc = to_string(BDTcutvalues.c);
        ttree->Draw("tripletMass>>hTriplMass_data1"+binning, "weight*(isMC==0  && category =="+c+" && bdt >= "+BDTa+")");
        hTriplMass_data1 = (TH1F *)gDirectory->Get("hTriplMass_data1");
        ttree->Draw("tripletMass>>hTriplMass_data2"+binning, "weight*(isMC==0  && category =="+c+" && bdt < "+BDTa+" && bdt >= "+BDTb+")");
        hTriplMass_data2 = (TH1F *)gDirectory->Get("hTriplMass_data2");
        ttree->Draw("tripletMass>>hTriplMass_data3"+binning, "weight*(isMC==0  && category =="+c+" && bdt < "+BDTb+" && bdt >= "+BDTc+")");
        hTriplMass_data3 = (TH1F *)gDirectory->Get("hTriplMass_data3");
        //signal
        ttree->Draw("tripletMass>>hSignal1"+binning, "weight*(isMC>0  && category =="+c+" && bdt >= "+BDTa+")");
        hSignal1 = (TH1F *)gDirectory->Get("hSignal1");
        ttree->Draw("tripletMass>>hSignal2"+binning, "weight*(isMC>0  && category =="+c+" && bdt < "+BDTa+" && bdt >= "+BDTb+")");
        hSignal2 = (TH1F *)gDirectory->Get("hSignal2");
        ttree->Draw("tripletMass>>hSignal3"+binning, "weight*(isMC>0  && category =="+c+" && bdt < "+BDTb+" && bdt >= "+BDTc+")");
        hSignal3 = (TH1F *)gDirectory->Get("hSignal3");


        hTriplMass_data1->SetDirectory(0);
        hTriplMass_data2->SetDirectory(0);
        hTriplMass_data3->SetDirectory(0);
        hSignal1->SetDirectory(0);
        hSignal2->SetDirectory(0);
        hSignal3->SetDirectory(0);

        //Write and close the file
        fout->cd();
        hTriplMass_data1->Write("data_obs"+category+"1");
        hTriplMass_data2->Write("data_obs"+category+"2");
        hTriplMass_data3->Write("data_obs"+category+"3");
        hTriplMass_data1->Write("background"+category+"1");
        hTriplMass_data2->Write("background"+category+"2");
        hTriplMass_data3->Write("background"+category+"3");
        hSignal1->Write("signal"+category+"1");
        hSignal2->Write("signal"+category+"2");
        hSignal3->Write("signal"+category+"3");

    }
    fout->Close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    return 0;
}

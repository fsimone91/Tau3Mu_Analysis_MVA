//root -l -b plot_TMVA_bdtcorr.C\(\"C\"\)

#include "TH1F.h"
#include <cmath>
#include <string> 
#include "../T3M_common.h"

void plot_TMVA_bdtcorr(TString category) 
{
    TString TMVA_filename = TMVA_inputpath+category;
    TString bdt_range_signal[] = {
                                  "tripletMass<1.75 && tripletMass>1.74",
                                  "tripletMass<1.76 && tripletMass>1.75",
                                  "tripletMass<1.77 && tripletMass>1.76",
                                  "tripletMass<1.78 && tripletMass>1.77",
                                  "tripletMass<1.79 && tripletMass>1.78",
                                  "tripletMass<1.80 && tripletMass>1.79",
                                  "tripletMass<1.81 && tripletMass>1.80",
                                  "tripletMass<1.82 && tripletMass>1.81"
                                 };
    TString bdt_range_bkg[] = {
                                  "tripletMass<1.70 && tripletMass>1.65",
                                  "tripletMass<1.75 && tripletMass>1.70",
                                  "tripletMass<1.85 && tripletMass>1.80",
                                  "tripletMass<1.90 && tripletMass>1.85"
                                 };

    TString file_name = workdir+TMVA_inputpath+"outputTree.root";

    //open input files
    TChain *t = new TChain("outputTree");
    t->Add(file_name); 
    std::cout<<"Opened input file: "<<file_name<<std::endl;

    TString c = ""; 
    if (category=="A") c = "0";
    if (category=="B") c = "1";
    if (category=="C") c = "2";
    TString signal = "weight*(isMC>0  && category =="+c+")";
    TString bkg    = "weight*(isMC==0 && category =="+c+")";
    
    int n_signal = sizeof(bdt_range_signal)/sizeof(bdt_range_signal[0]);
    int n_bkg = sizeof(bdt_range_bkg)/sizeof(bdt_range_bkg[0]);
    TH1F *hTrain_signal[n_signal];
    TH1F *hTrain_bkg[n_bkg];

    TString binning;
    binning = "(80, -0.4, 0.4)";
    TString varname = "bdt";

    gStyle->SetOptTitle(0);
 
    //SIGNAL
    TCanvas *c1 = new TCanvas("c1","Signal - cat "+category,150,10,990,660);
    THStack hs_signal("hs_signal","hs_signal");
    TLegend*leg_signal = new TLegend(0.1,0.7,0.45,0.9);

    for(int i = 0; i<n_signal; i++){
        TString range = bdt_range_signal[i];

        TString s = std::to_string(i);
        cout<<"bdt range "<<range<<endl;

        binning = "(40, -0.4, 0.4)";
        t->Draw(varname+">>hTrain_signal"+s+binning, range+" && "+signal);
        hTrain_signal[i] = (TH1F *)gDirectory->Get("hTrain_signal"+s); 
        hTrain_signal[i]->SetLineColor(1+i);
        hTrain_signal[i]->Scale(1/hTrain_signal[i]->GetEntries());
        hs_signal.Add(hTrain_signal[i]);
        leg_signal->AddEntry(hTrain_signal[i],range,"l");

    }
    hs_signal.Draw("hist nostack");
    leg_signal->Draw();
    hs_signal.GetXaxis()->SetTitle("BDT score");

    gPad->Modified();
    c1->Update();
    c1->SaveAs(workdir+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"correlation_signal.png");

    //BACKGROUND
    TCanvas *c2 = new TCanvas("c2","Bkg - cat "+category,150,10,990,660);
    gStyle->SetOptTitle(0);
    THStack hs_bkg("hs_bkg","hs_bkg");
    TLegend*leg_bkg = new TLegend(0.55,0.7,0.9,0.9);

    for(int i = 0; i<n_bkg; i++){
        TString range = bdt_range_bkg[i];

        TString s = std::to_string(i);
        cout<<"bdt range "<<range<<endl;

        binning = "(40, -0.4, 0.4)";
        t->Draw(varname+">>hTrain_bkg"+s+binning, range+" && "+bkg);
        hTrain_bkg[i] = (TH1F *)gDirectory->Get("hTrain_bkg"+s); 
        hTrain_bkg[i]->GetXaxis()->SetTitle("BDT score");
        hTrain_bkg[i]->SetLineColor(1+i);
        hTrain_bkg[i]->Scale(1/hTrain_bkg[i]->GetEntries());
        hs_bkg.Add(hTrain_bkg[i]);
        leg_bkg->AddEntry(hTrain_bkg[i],range,"l");

    }
    hs_bkg.Draw("hist nostack");
    leg_bkg->Draw();
    hs_bkg.GetXaxis()->SetTitle("BDT score");

    gPad->Modified();
    c2->Update();
    c2->SaveAs(workdir+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"correlation_bkg.png");
}

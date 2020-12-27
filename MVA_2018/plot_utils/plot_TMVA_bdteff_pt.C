//root -l -b plot_TMVA_bdtcorr.C\(\"C\"\)

#include "TH1F.h"
#include <cmath>
#include <string> 
#include "../T3M_common.h"
#include "../BDT_optimal_cut_v2.h"

void plot_TMVA_bdteff_pt(TString category) 
{
    TString TMVA_filename = TMVA_inputpath+category;

    TString file_name = workdir+TMVA_inputpath+"outputTree.root";

    //open input files
    TChain *t = new TChain("outputTree");
    t->Add(file_name); 
    std::cout<<"Opened input file: "<<file_name<<std::endl;

    TString c = ""; 
    if (category=="A") c = "0";
    if (category=="B") c = "1";
    if (category=="C") c = "2";
    TString sample = "DsTau3Mu"; 
    TString s = ""; 
    if (sample=="DsTau3Mu") s = "1";
    if (sample=="B0Tau3Mu") s = "2";
    if (sample=="BpTau3Mu") s = "3";
    
    //BDTcut BDTcutvalues = Get_BDT_cut_v2(category, false);
    //TString bdt_cut = std::to_string(BDTcutvalues.b);
    TString bdt_cut[7] = {"-0.20", "-0.10", "0.0", "0.10", "0.20", "0.30", "0.40"};
    TH1F *hratio_sig[7];
    TH1F *hratio_bkg[7];
    int n_cut = 7; 
    THStack hs_ratio_sig("hs_ratio_sig","hs_ratio_sig");

    //denominators
    TH1F *hden_sig;
    TH1F *hden_bkg;
    TH1F *hden_sig_plot;
    TH1F *hden_bkg_plot;
    TString signal = "weight*(isMC=="+s+"  && category =="+c+")";
    TString bkg    = "weight*(isMC==0 && category =="+c+" && isSB == 1)";

    TString binning_sig, binning_bkg;
    binning_bkg = "(80, 0.5, 40.5)";
    binning_sig = "(80, 0.5, 40.5)";
    TString varname = "Ptmu3";

    t->Draw(varname+">>hden_sig"+binning_sig, signal);
    hden_sig = (TH1F *)gDirectory->Get("hden_sig"); 
    hden_sig_plot = (TH1F*)hden_sig->Clone("hden_sig_plot");
    t->Draw(varname+">>hden_bkg"+binning_bkg, bkg);
    hden_bkg = (TH1F *)gDirectory->Get("hden_bkg"); 
    hden_bkg_plot = (TH1F*)hden_bkg->Clone("hden_bkg_plot");


    for(int i=0; i<n_cut; i++){
       TString signal_num = "weight*(isMC=="+s+"  && category =="+c+" && bdt>"+bdt_cut[i]+")";
       TH1F *hnum_sig;
   
       //SIGNAL
       t->Draw(varname+">>hnum_sig"+binning_sig, signal_num);
       hnum_sig = (TH1F *)gDirectory->Get("hnum_sig"); 
   
       // Define the ratio plot
       hratio_sig[i] = (TH1F*)hnum_sig->Clone("hratio_sig");
       hratio_sig[i]->Sumw2();
       hratio_sig[i]->Divide(hden_sig);
       hratio_sig[i]->SetStats(0);
       // Ratio plot settings
       gStyle->SetLineWidth(2);
       hratio_sig[i]->SetTitle(""); // Remove the ratio title
       hratio_sig[i]->GetYaxis()->SetTitle("bdt efficiency");
       hratio_sig[i]->GetXaxis()->SetTitle("#mu_{3} p_{T} (GeV)");
       hratio_sig[i]->SetLineColor(1+i);
       hratio_sig[i]->GetYaxis()->SetTitleSize(20);
       hratio_sig[i]->GetYaxis()->SetTitleFont(43);
       hratio_sig[i]->GetYaxis()->SetTitleOffset(1.25);
       hratio_sig[i]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       hratio_sig[i]->GetYaxis()->SetLabelSize(15);
   
       // X axis ratio plot settings
       hratio_sig[i]->GetXaxis()->SetTitle(varname);
       hratio_sig[i]->GetXaxis()->SetTitleSize(20);
       hratio_sig[i]->GetXaxis()->SetTitleFont(43);
       hratio_sig[i]->GetXaxis()->SetTitleOffset(3);
       hratio_sig[i]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       hratio_sig[i]->GetXaxis()->SetLabelSize(15);
    } 
    //SIGNAL
    TCanvas *c_sig = new TCanvas("c_sig","Signal "+sample+" - cat "+category,150,10,800,800);
    TLegend*leg_signal = new TLegend(0.1,0.3,0.3,0.5);
    hden_sig_plot->GetYaxis()->SetTitle("BDT cut efficiency");
    hden_sig_plot->GetXaxis()->SetTitle("#mu_{3} p_{T} (GeV)");
    hden_sig_plot->SetTitle("Signal "+sample+" - cat "+category);
    cout<<"hden_sig_plot->GetMaximum() "<<hden_sig_plot->GetMaximum()<<" hden_sig_plot->GetMean(1) "<<hden_sig_plot->GetMean(1)<<" hden_sig_plot->GetNbinsX() "<<hden_sig_plot->GetNbinsX()<<endl;
    double mean_sig = hden_sig_plot->GetMean() / hden_sig_plot->GetNbinsX(); //avg bin height
    hden_sig_plot->Scale(0.5 / hden_sig_plot->GetMaximum());
    //hden_sig_plot->Scale(1.0 / hden_sig_plot->Integral());
    hden_sig_plot->GetYaxis()->SetRangeUser(0.,1.);
    hden_sig_plot->GetYaxis()->SetRange(0.,1.);
    hden_sig_plot->SetMaximum(1.);
    hden_sig_plot->SetStats(0);
    hden_sig_plot->Draw("hist");
    leg_signal->AddEntry(hden_sig_plot, sample+" #mu_{3} p_{T} (GeV)", "f");

    for(int i=0; i<n_cut; i++){
       hratio_sig[i]->Draw("same EL");
       leg_signal->AddEntry(hratio_sig[i], "Eff. BTD>"+bdt_cut[i], "epl");
    } 
   gPad->Modified();
   leg_signal->Draw();
   c_sig->SaveAs(workdir+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"bdteff_signal"+sample+".png");
 
    //BACKGROUND
    for(int i=0; i<n_cut; i++){
       TString bkg_num    = "weight*(isMC==0 && category =="+c+" && isSB == 1 && bdt>"+bdt_cut[i]+")";
       TH1F *hnum_bkg;
       t->Draw(varname+">>hnum_bkg"+binning_bkg, bkg_num);
       hnum_bkg = (TH1F *)gDirectory->Get("hnum_bkg"); 
   
       // Define the ratio plot
       hratio_bkg[i] = (TH1F*)hnum_bkg->Clone("hratio_bkg");
       hratio_bkg[i]->Sumw2();
       hratio_bkg[i]->Divide(hden_bkg);
       hratio_bkg[i]->SetStats(0);
       // Ratio plot settings
       gStyle->SetLineWidth(2);
       hratio_bkg[i]->SetTitle(""); // Remove the ratio title
       hratio_bkg[i]->GetYaxis()->SetTitle("bdt efficiency");
       hratio_bkg[i]->GetXaxis()->SetTitle("#mu_{3} p_{T} (GeV)");
       hratio_bkg[i]->SetLineColor(1+i);
       hratio_bkg[i]->GetYaxis()->SetTitleSize(20);
       hratio_bkg[i]->GetYaxis()->SetTitleFont(43);
       hratio_bkg[i]->GetYaxis()->SetTitleOffset(1.25);
       hratio_bkg[i]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       hratio_bkg[i]->GetYaxis()->SetLabelSize(15);
   
       // X axis ratio plot settings
       hratio_bkg[i]->GetXaxis()->SetTitle(varname);
       hratio_bkg[i]->GetXaxis()->SetTitleSize(20);
       hratio_bkg[i]->GetXaxis()->SetTitleFont(43);
       hratio_bkg[i]->GetXaxis()->SetTitleOffset(3);
       hratio_bkg[i]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       hratio_bkg[i]->GetXaxis()->SetLabelSize(15);
    }  
    TCanvas *c_bkg = new TCanvas("c_bkg","Bkg - cat "+category,150,10,800,800);
    TLegend*leg_bkg = new TLegend(0.1,0.8,0.3,0.95);
    hden_bkg_plot->GetYaxis()->SetTitle("BDT cut efficiency");
    hden_bkg_plot->GetXaxis()->SetTitle("#mu_{3} p_{T} (GeV)");
    hden_bkg_plot->SetTitle("Data 2018 - cat "+category);
    hden_bkg_plot->Scale(0.5 / hden_bkg_plot->GetMaximum());
    //hden_bkg_plot->Scale(3.0 / hden_bkg->Integral());
    hden_bkg_plot->GetYaxis()->SetRangeUser(0.,1.);
    hden_bkg_plot->GetYaxis()->SetRange(0.,1.);
    hden_bkg_plot->SetMaximum(1.);
    hden_bkg_plot->SetStats(0);
    hden_bkg_plot->Draw("hist");
    leg_bkg->AddEntry(hden_bkg_plot, "Data #mu_{3} p_{T} (GeV)", "f");
    for(int i=0; i<n_cut; i++){
       hratio_bkg[i]->Draw("ep same");
       leg_bkg->AddEntry(hratio_bkg[i], "Eff. BTD>"+bdt_cut[i], "epl");
    } 
   gPad->Modified();
   leg_bkg->Draw();
   c_bkg->SaveAs(workdir+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"bdteff_bkg.png");
}

//root -l -b plot_TMVA_bdtcorr.C\(\"C\"\)

#include "TH1F.h"
#include <cmath>
#include <string> 
#include "../T3M_common.h"
#include "../BDT_optimal_cut_v2.h"

void plot_TMVA_bdteff(TString category) 
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
    
    BDTcut BDTcutvalues = Get_BDT_cut_v2(category, false);
    TString bdt_cut = std::to_string(BDTcutvalues.b);

    TString signal = "weight*(isMC>0  && category =="+c+")";
    TString bkg    = "weight*(isMC==0 && category =="+c+" && isSB == 1)";
    TString signal_num = "weight*(isMC>0  && category =="+c+" && bdt>"+bdt_cut+")";
    TString bkg_num    = "weight*(isMC==0 && category =="+c+" && isSB == 1 && bdt>"+bdt_cut+")";
    TH1F *hnum_sig;
    TH1F *hnum_bkg;
    TH1F *hden_sig;
    TH1F *hden_bkg;

    TString binning;
    binning = "(76, 1.62, 2.0)";
    TString varname = "tripletMass";

    gStyle->SetOptTitle(0);
    Double_t norm = 1.38;
 
    //SIGNAL
    TCanvas *c1 = new TCanvas("c1","Signal - cat "+category,150,10,990,660);
    TLegend*leg_signal = new TLegend(0.1,0.7,0.45,0.9);

    t->Draw(varname+">>hnum_sig"+binning, signal_num);
    hnum_sig = (TH1F *)gDirectory->Get("hnum_sig"); 
    t->Draw(varname+">>hden_sig"+binning, signal);
    hden_sig = (TH1F *)gDirectory->Get("hden_sig"); 

    // Define the ratio plot
    TH1F *hratio_sig = (TH1F*)hnum_sig->Clone("hratio_sig");
    hratio_sig->Sumw2();
    hratio_sig->Divide(hden_sig);
    hratio_sig->SetStats(0);
    // Ratio plot settings
    gStyle->SetLineWidth(2);
    hratio_sig->SetTitle(""); // Remove the ratio title
    hratio_sig->GetYaxis()->SetTitle("bdt efficiency");
    hratio_sig->GetXaxis()->SetTitle("m(3#mu) (GeV)");
    hratio_sig->SetLineColor(kBlack);
    hratio_sig->GetYaxis()->SetTitleSize(20);
    hratio_sig->GetYaxis()->SetTitleFont(43);
    hratio_sig->GetYaxis()->SetTitleOffset(1.25);
    hratio_sig->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio_sig->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    hratio_sig->GetXaxis()->SetTitle(varname);
    hratio_sig->GetXaxis()->SetTitleSize(20);
    hratio_sig->GetXaxis()->SetTitleFont(43);
    hratio_sig->GetXaxis()->SetTitleOffset(3);
    hratio_sig->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio_sig->GetXaxis()->SetLabelSize(15);

    hden_sig->GetYaxis()->SetRangeUser(0.,1.);
    hden_sig->GetYaxis()->SetRange(0.,1.);
    hden_sig->Scale(norm / hden_sig->Integral());
    hden_sig->SetMaximum(1);
    hden_sig->Draw("hist");
    hratio_sig->Draw("ep same");
    gPad->Modified();
    c1->Update();
    c1->SaveAs(workdir+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"bdteff_signal.png");


    //BACKGROUND
    TCanvas *c2 = new TCanvas("c2","Bkg - cat "+category,150,10,990,660);
    TLegend*leg_bkg = new TLegend(0.1,0.7,0.45,0.9);

    t->Draw(varname+">>hnum_bkg"+binning, bkg_num);
    hnum_bkg = (TH1F *)gDirectory->Get("hnum_bkg"); 
    t->Draw(varname+">>hden_bkg"+binning, bkg);
    hden_bkg = (TH1F *)gDirectory->Get("hden_bkg"); 

    // Define the ratio plot
    TH1F *hratio_bkg = (TH1F*)hnum_bkg->Clone("hratio_bkg");
    hratio_bkg->Sumw2();
    hratio_bkg->Divide(hden_bkg);
    hratio_bkg->SetStats(0);
    // Ratio plot settings
    gStyle->SetLineWidth(2);
    hratio_bkg->SetTitle(""); // Remove the ratio title
    hratio_bkg->GetYaxis()->SetTitle("bdt efficiency");
    hratio_bkg->GetXaxis()->SetTitle("m(3#mu) (GeV)");
    hratio_bkg->SetLineColor(kBlack);
    hratio_bkg->GetYaxis()->SetTitleSize(20);
    hratio_bkg->GetYaxis()->SetTitleFont(43);
    hratio_bkg->GetYaxis()->SetTitleOffset(1.25);
    hratio_bkg->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio_bkg->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    hratio_bkg->GetXaxis()->SetTitle(varname);
    hratio_bkg->GetXaxis()->SetTitleSize(20);
    hratio_bkg->GetXaxis()->SetTitleFont(43);
    hratio_bkg->GetXaxis()->SetTitleOffset(3);
    hratio_bkg->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio_bkg->GetXaxis()->SetLabelSize(15);

    hden_bkg->GetYaxis()->SetRangeUser(0.,1.);
    hden_bkg->GetYaxis()->SetRange(0.,1.);
    hden_bkg->Scale(norm / hden_bkg->Integral());
    hden_bkg->SetMaximum(0.2);
    hden_bkg->Draw("hist");
    hratio_bkg->Draw("ep same");
    gPad->Modified();
    c2->Update();
    c2->SaveAs(workdir+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"bdteff_bkg.png");

}

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1F.h"
#include "../../T3M_common.h"
#include "../Control_common.h"

using namespace RooFit;

void Ds_syst() 
{
    //reading Ds yield measured from control channel
    vector<Double_t> Dsyield;
    vector<Double_t> Dsyield_err;
    Double_t temp1;
    Double_t temp2;

    ifstream file_dsyield("dsphipi_yield.txt");
    while ((file_dsyield >> temp1) && (file_dsyield >> temp2))
    {
        Dsyield.push_back(temp1);
        Dsyield_err.push_back(temp2);
        //Dsyield_err.push_back(sqrt(temp1));
    }

    //reading SB entries in data from the different L1 seeds
    int nrun = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    TH1F *hdata;

    TCut cutSB = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)"; // data sidebands
    TCut L1_double  = "l1double_fired";
    TCut L1_triple  = "l1triple_fired";
    TCut L1_double0 = "l1double_DoubleMu0_fired";
    TCut L1_double4 = "l1double_DoubleMu4_fired";

    std::vector<Double_t> sideByield_OR, sideByield_double, sideByield_triple, sideByield_double0, sideByield_double4;
    for(auto i=0; i<nrun; i++){
        TChain *tdata = new TChain("FinalTree");
        tdata->Add(inputpath_datarun[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
        //sideband OR
        tdata->Draw("tripletMass>>hdata", cutSB);
        hdata = (TH1F *)gDirectory->Get("hdata");
        sideByield_OR.push_back(hdata->GetEntries());
        //sideband double
        tdata->Draw("tripletMass>>hdata", cutSB&&L1_double);
        hdata = (TH1F *)gDirectory->Get("hdata");
        sideByield_double.push_back(hdata->GetEntries());
        //sideband triple
        tdata->Draw("tripletMass>>hdata", cutSB&&L1_triple);
        hdata = (TH1F *)gDirectory->Get("hdata");
        sideByield_triple.push_back(hdata->GetEntries());
        //sideband double0
        tdata->Draw("tripletMass>>hdata", cutSB&&L1_double0);
        hdata = (TH1F *)gDirectory->Get("hdata");
        sideByield_double0.push_back(hdata->GetEntries());
        //sideband double4
        tdata->Draw("tripletMass>>hdata", cutSB&&L1_double4);
        hdata = (TH1F *)gDirectory->Get("hdata");
        sideByield_double4.push_back(hdata->GetEntries());
    }

    TCut cutPeak = "(tripletMass<2.01 && tripletMass>1.93)"; // DsPhiPi MC peak
    TH1F *hMC;
    TChain *tMC = new TChain("FinalTree_Control");
    tMC->Add(inputpath_DsPhiPi); 
    std::cout<<"Opened input file: "<< inputpath_DsPhiPi <<std::endl;
    //entries DsPhiPiMC
    tMC->Draw("tripletMass>>hMC", cutPeak);
    hMC = (TH1F *)gDirectory->Get("hMC");
    Double_t Dsyield_MC = hMC->GetEntries();

    std::vector<TString> era = {'B','C','D','E','F'};

    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    Double_t Ds_tot = 0;
    Double_t Ds_tot_MC = 0;
    TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);
    TH1F *hDs    = new TH1F("h","Ds",nrun,0,nrun);
    TH1F *hDs_MC = new TH1F("h","Ds_MC",nrun,0,nrun);
    for (Int_t i=0;i<nrun;i++) {

      hDs->Fill(i,Dsyield[i]/Lumi_data[i]);
      hDs->GetXaxis()->SetBinLabel(i+1, era[i]);
      Double_t err = pow(pow(Dsyield_err[i]/Lumi_data[i],2) + pow(0.023*Lumi_data[i]*Dsyield[i]/pow(Lumi_data[i],2),2),0.5);
      hDs->SetBinError(i+1, err);

      cout<<"MC "<<Lumi_data[i]*xsection_mc*BR/N_MC * Dsyield_MC<<endl;
      hDs_MC->Fill(i,xsection_mc*BR/N_MC * Dsyield_MC);
      hDs_MC->GetXaxis()->SetBinLabel(i+1, era[i]);
      Double_t err_MC = pow( ( xsection_mc*BR/N_MC * Dsyield_MC ),0.5);
      hDs_MC->SetBinError(i+1, err_MC);

      Ds_tot+=Dsyield[i];
      Ds_tot_MC+=Lumi_data[i]*xsection_mc*BR/N_MC * Dsyield_MC;
    }
    cout<<"Ds yield scale factor = "<<Ds_tot/Ds_tot_MC<<endl;

    hDs->SetLineColor(kRed);
    hDs-> SetMarkerStyle(22);
    hDs-> SetMarkerSize(1.5);
    hDs-> SetMarkerColor(kRed);
    hDs->Draw();
    //hDs->GetYaxis()->SetRangeUser(0,4000);
    hDs->GetXaxis()->SetTitle("2017 era");
    hDs->GetXaxis()->SetTitleSize(0.045);
    hDs->GetXaxis()->SetLabelSize(0.05);
    hDs->GetYaxis()->SetTitle("Events/fb^{-1}");
    hDs->GetYaxis()->SetTitleSize(0.045);
    hDs->GetYaxis()->SetTitleOffset(1);
    hDs->GetYaxis()->SetLabelSize(0.04);
    hDs->SetStats(0);

    hDs_MC->SetLineColor(kBlue);
    //hDs_MC-> SetMarkerStyle(22);
    //hDs_MC-> SetMarkerSize(1.5);
    //hDs_MC-> SetMarkerColor(kBlue);
    hDs_MC->Draw("same");

    TLegend*leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(hDs,"Ds yield data","l");
    leg->AddEntry(hDs_MC,"Ds expected","l");
    leg->Draw();
    c1->Update();

    TH1F *h3mu = new TH1F("h","3mu",nrun,0,nrun);
    for (Int_t i=0;i<nrun;i++) {
      cout<<sideByield_double[i]<<endl;
      cout<<Lumi_data[i]<<endl;
      h3mu->Fill(i,sideByield_double[i]/Lumi_data[i]);
      h3mu->GetXaxis()->SetBinLabel(i+1, era[i]);
      Double_t err = pow(sideByield_double[i]/pow(Lumi_data[i],2) + pow(0.023*Lumi_data[i]*sideByield_double[i]/pow(Lumi_data[i],2),2),0.5);
      h3mu->SetBinError(i+1, err);
    }
    h3mu->SetLineColor(kBlue);
    h3mu-> SetMarkerStyle(22);
    h3mu-> SetMarkerSize(1.5);
    h3mu-> SetMarkerColor(kBlue);
 //   h3mu->Draw("same");

 //   TLegend*leg = new TLegend(0.1,0.7,0.48,0.9);
 //   leg->AddEntry(hDs,"Ds yield","l");
 //   leg->AddEntry(h3mu,"3mu(SB)_l1DoubleMu","l");
 //   leg->Draw();
 //   c1->Update();

    // Create ratio histograms data dsphipi vs data sideband
    TCanvas *c2 = new TCanvas("c2","c2",200,10,700,500);
    TH1F *h3mu_ratio = (TH1F*)h3mu->Clone("h3mu_ratio");
    TH1F *hDs2 = (TH1F*)hDs->Clone("hDs2");
    h3mu_ratio->Divide(hDs2);
    h3mu_ratio->Draw();
    //Fitting the ratio with horizontal line
    TF1 *func = new TF1("fit","[0]",0,nrun);
    func->SetParNames("Ratio_avg");
    h3mu_ratio->Fit("fit");

    h3mu_ratio->SetLineColor(kBlack);
    h3mu_ratio-> SetMarkerStyle(22);
    h3mu_ratio-> SetMarkerSize(1.5);
    h3mu_ratio-> SetMarkerColor(kBlack);
    h3mu_ratio->GetXaxis()->SetTitle("2017 era");
    h3mu_ratio->GetXaxis()->SetTitleSize(0.045);
    h3mu_ratio->GetXaxis()->SetLabelSize(0.05);
    h3mu_ratio->GetYaxis()->SetTitle("Ds/3#muSB yield ratio");
    h3mu_ratio->GetYaxis()->SetTitleSize(0.045);
    h3mu_ratio->GetYaxis()->SetTitleOffset(1);
    h3mu_ratio->GetYaxis()->SetLabelSize(0.04);
    TF1 *f = h3mu_ratio->GetFunction("fit");
    Double_t ndf = nrun;
    Double_t scale = pow(f->GetChisquare()/(ndf-1),0.5);
    cout<<"Scale factor "<<scale<<endl;   
    cout<<"Relative error after rescaling "<<f->GetParError(0)/f->GetParameter(0)*scale<<endl;
    h3mu_ratio->Draw();
    c2->Update();

    // Create ratio histograms data dsphipi vs yield expected from MC
    TCanvas *c4 = new TCanvas("c4","c4",200,10,700,500);
    TH1F *hDs_ratio = (TH1F*)hDs->Clone("hDs_ratio");
    hDs_ratio->Divide(hDs_MC);
    hDs_ratio->Draw();
    Double_t ave_ratio = Ds_tot/Ds_tot_MC;
    //Drawing line corresponding to average ratio
    TLine *line = new TLine(0,ave_ratio,nrun,ave_ratio);
    line->SetLineColor(kRed);

    hDs_ratio->SetLineColor(kBlack);
    hDs_ratio-> SetMarkerStyle(22);
    hDs_ratio-> SetMarkerSize(1.5);
    hDs_ratio-> SetMarkerColor(kBlack);
    hDs_ratio->GetXaxis()->SetTitle("2017 era");
    hDs_ratio->GetXaxis()->SetTitleSize(0.045);
    hDs_ratio->GetXaxis()->SetLabelSize(0.05);
    hDs_ratio->GetYaxis()->SetTitle("Ds yield data/MC ratio");
    hDs_ratio->GetYaxis()->SetTitleSize(0.045);
    hDs_ratio->GetYaxis()->SetTitleOffset(1);
    hDs_ratio->GetYaxis()->SetLabelSize(0.04);
    hDs_ratio->Draw();
    line->Draw("same");
    c4->Update();

    TCanvas *c3 = new TCanvas("c3","c3",200,10,700,500);

    TH1F *h3mu_or = new TH1F("h","3mu_or",nrun,0,nrun);
    for (Int_t i=0;i<nrun;i++) {
      h3mu_or->Fill(i,sideByield_OR[i]/Lumi_data[i]);
      h3mu_or->GetXaxis()->SetBinLabel(i+1, era[i]);
      Double_t err = pow(sideByield_OR[i]/pow(Lumi_data[i],2) + pow(0.023*Lumi_data[i]*sideByield_OR[i]/pow(Lumi_data[i],2),2),0.5);
      h3mu_or->SetBinError(i+1, err);
    }
    h3mu_or->SetLineColor(kBlack);
    h3mu_or-> SetMarkerStyle(22);
    h3mu_or-> SetMarkerSize(1.5);
    h3mu_or-> SetMarkerColor(kBlack);
    h3mu_or->Draw();
    //h3mu_or->GetYaxis()->SetRangeUser(0,2000);
    h3mu_or->GetXaxis()->SetTitle("2017 era");
    h3mu_or->GetXaxis()->SetTitleSize(0.045);
    h3mu_or->GetXaxis()->SetLabelSize(0.05);
    h3mu_or->GetYaxis()->SetTitle("Events/fb^{-1}");
    h3mu_or->GetYaxis()->SetTitleSize(0.045);
    h3mu_or->GetYaxis()->SetTitleOffset(1);
    h3mu_or->GetYaxis()->SetLabelSize(0.04);
    h3mu_or->SetStats(0);

    TH1F *h3mu_double = new TH1F("h","3mu_double",nrun,0,nrun);
    for (Int_t i=0;i<nrun;i++) {
      h3mu_double->Fill(i,sideByield_double[i]/Lumi_data[i]);
      h3mu_double->GetXaxis()->SetBinLabel(i+1, era[i]);
      Double_t err = pow(sideByield_double[i]/pow(Lumi_data[i],2) + pow(0.023*Lumi_data[i]*sideByield_double[i]/pow(Lumi_data[i],2),2),0.5);
      h3mu_double->SetBinError(i+1, err);
    }
    h3mu_double->SetLineColor(kBlue);
    h3mu_double-> SetMarkerStyle(22);
    h3mu_double-> SetMarkerSize(1.5);
    h3mu_double-> SetMarkerColor(kBlue);
    h3mu_double->Draw("same");

    TH1F *h3mu_triple = new TH1F("h","3mu_triple",nrun,0,nrun);
    for (Int_t i=0;i<nrun;i++) {
      h3mu_triple->Fill(i,sideByield_triple[i]/Lumi_data[i]);
      h3mu_triple->GetXaxis()->SetBinLabel(i+1, era[i]);
      Double_t err = pow(sideByield_triple[i]/pow(Lumi_data[i],2) + pow(0.023*Lumi_data[i]*sideByield_triple[i]/pow(Lumi_data[i],2),2),0.5);
      h3mu_triple->SetBinError(i+1, err);
    }
    h3mu_triple->SetLineColor(kRed);
    h3mu_triple-> SetMarkerStyle(22);
    h3mu_triple-> SetMarkerSize(1.5);
    h3mu_triple-> SetMarkerColor(kRed);
    h3mu_triple->Draw("same");

    TLegend*leg3 = new TLegend(0.1,0.7,0.48,0.9);
    leg3->AddEntry(h3mu_or,"3mu_l1OR","l");
    leg3->AddEntry(h3mu_double,"3mu_l1DoubleMu","l");
    leg3->AddEntry(h3mu_triple,"3mu_l1TripleMu","l");
    leg3->Draw();
    c3->Update();
}

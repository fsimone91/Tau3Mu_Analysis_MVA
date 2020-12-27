#include "TH1F.h"
#include <cmath>
#include <string> 

double TH1_integral(TH1F *h, float xmin, float xmax){
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/ axis->GetBinWidth(bmax);

    return integral;
}

void plot_TMVA_inputvariables(TString category, bool apply_weights) 
{
    category = "_"+category;
    TString TMVA_filename = "MuonMVA_29april"+category;

    std::vector<TString> var = {
                                      "abs_mu_eta_",
                                      "mu_pt",
                                      "mu_phi",
                                      //Muon reconstruction
                                      "mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum",
                                      "mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition",
                                      "mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", //cut >0 ?
                                      "mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", //cut >0 ?
                                      "mu_combinedQuality_globalDeltaEtaPhi",
                                      "log_mu_combinedQuality_trkKink_",
                                      "log_mu_combinedQuality_glbKink_",
                                      "mu_combinedQuality_glbTrackProbability",

                                      //collection of hits in the HitPattern
                                      "mu_Numberofvalidtrackerhits", //Valid Tracker Hits
                                      "mu_Numberofvalidpixelhits",
                                      "mu_trackerLayersWithMeasurement",
                                      "mu_GLhitPattern_numberOfValidMuonHits", //cut >0 ?
                                      "mu_validMuonHitComb", //Hits in DT, CSC, RPC //cut >0 ?

                                      //muon track reconstruction
                                      "mu_numberOfMatchedStations",
                                      "mu_segmentCompatibility", //cut >0.05 ?
                                      "mu_timeAtIpInOut", //to be studied
                                      "mu_timeAtIpInOutErr", //to be studied

                                      //general track properties
                                      "mu_GLnormChi2>50?50:mu_GLnormChi2", 
                                      "mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2",
                                      "mu_outerTrack_normalizedChi2",
                                      "mu_innerTrack_validFraction", //Inner Valid Fraction
                                      "mu_QInnerOuter"
                                      //"", //dxyRef
                                      //"", //dzRef

                                      //custom variables track multiplicity
                                      };


    cout<<"Opening data file"<<endl;
    TFile *f_tmva = new TFile("TMVA_"+TMVA_filename+".root","READ");
    cout<<"Opened TMVA file"<<endl;
    
    TTree *tTrain = (TTree*)f_tmva->Get(TMVA_filename+"/TrainTree");
    //TTree *tTest = (TTree*)f_tmva->Get(TMVA_filename+"/TestTree");

    int n = var.size();

    TH1F *hTrain_signal[n];
    TH1F *hTrain_bkg[n];

    TString binning;
 
   for(int i = 0; i<n; i++){
        TString varname = var[i];

        TString s = std::to_string(i);
        cout<<"Input variable "<<varname<<endl;

        TString binning = "";
        TString w = "";
        if(apply_weights) w = "weight*";
        if(varname=="mu_pt") binning = binning+"(90, 0, 30)"; //binning AN
        tTrain->Draw(varname+">>hTrain_signal"+s+binning, w+"(classID==0)");
        tTrain->Draw(varname+">>hTrain_bkg"+s+binning,    w+"(classID==1)");

        hTrain_signal[i] = (TH1F *)gDirectory->Get("hTrain_signal"+s); 
        hTrain_bkg[i] = (TH1F *)gDirectory->Get("hTrain_bkg"+s); 

        //Data
        TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(111111);
        hTrain_signal[i]->SetLineColor(kBlue);
        hTrain_signal[i]->SetLineWidth(3);
        hTrain_signal[i]->SetFillStyle(3004);
        hTrain_signal[i]->SetFillColor(kBlue);
        hTrain_bkg[i]->SetLineColor(kRed);
        hTrain_bkg[i]->SetLineWidth(3);
        hTrain_bkg[i]->SetFillStyle(3005);
        hTrain_bkg[i]->SetFillColor(kRed);
       
        Double_t norm = 1; 
        double X_min = std::min(hTrain_signal[i]->GetXaxis()->GetXmin(), hTrain_bkg[i]->GetXaxis()->GetXmin());
        double X_max = std::max(hTrain_signal[i]->GetXaxis()->GetXmax(), hTrain_bkg[i]->GetXaxis()->GetXmax());
        hTrain_signal[i]->Scale(norm / TH1_integral(hTrain_signal[i], X_min, X_max));
        hTrain_bkg[i]->Scale(norm / TH1_integral(hTrain_bkg[i], X_min, X_max));

        hTrain_signal[i]->GetXaxis()->SetTitle(varname);

        hTrain_signal[i]->GetXaxis()->SetRange(1, hTrain_signal[i]->GetNbinsX() + 1); // will draw with the overflow bin
        hTrain_bkg[i]->GetXaxis()->SetRange(1, hTrain_bkg[i]->GetNbinsX() + 1);    // will draw with the overflow bin

        ////Adjust Y range
        //if(varname=="mu_phi")      hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.022);
        //if(varname=="cLP")       hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.275); 
        //if(varname=="segmComp")  hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 3.6);
        //if(varname=="tKink")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.11);
        //if(varname=="fv_nC")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.55);
        //if(varname=="fv_dphi3D") hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 50);
        //if(varname=="fv_d3Dsig") hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.3);
        //if(varname=="d0sig")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.8);
        //if(varname=="mindca_iso") hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 31);
        //if(varname=="trkRel")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 1.5);

        THStack hs("hs",varname);
        hs.Add(hTrain_signal[i]);
        hs.Add(hTrain_bkg[i]);
        hs.Draw("hist nostack");

        //hTrain_signal[i]->Draw("hist");
        //hTrain_bkg[i]->Draw("same hist");

        //position of legend = top right
        Double_t x1 = 0.45;
        Double_t x2 = 0.9;
        Double_t y1 = 0.7;
        Double_t y2 = 0.9;

        //bottom centre
        if(varname=="abs_mu_eta_" || varname=="mu_phi"){
            x1=0.3; x2=0.75;
            y1=0.1; y2=0.3;
        }
 
        TLegend*leg = new TLegend(x1, y1, x2, y2);
        leg->AddEntry(hTrain_signal[i],varname+"_signal","f");
        leg->AddEntry(hTrain_bkg[i],varname+"_bkg","f");
        leg->Draw();

        //c1->SetLogy();
        c1->Update();
        TString out_filename = "plots/TMVA_train_"+varname+category+".png";
        if(apply_weights) out_filename = "plots/TMVA_train_weighted_"+varname+category+".png";
        c1->SaveAs(out_filename);

    }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}

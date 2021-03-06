#include "TH1F.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "T3M_common.h"
#include "BDT_optimal_cut.h"

using namespace std;

void BDT_optimal_cut() 
{
    int ncat = sizeof(cat_label) / sizeof(cat_label[0]);
    
    for(int k=0; k<ncat; k++)
    {
  /*
            //full path of root files generated by TMVA
            TString file_name = "TMVA_"+TMVA_inputpath+cat_label[k]+".root";
	    TFile *f = new TFile(file_name,"READ");

	    TH1F *h_test_signal;
	    TH1F *h_test_bkg;
	    TH1F *h_test_signal2;
	    TH1F *h_test_bkg2;
	    TH1F *h_train_signal2;
	    TH1F *h_train_bkg2;
            //high granularity BDT scores
	    h_test_signal = (TH1F*)f->Get(TMVA_inputpath+cat_label[k]+"/Method_"+method+"/"+method+"/MVA_"+method+"_S_high");
	    h_test_bkg = (TH1F*)f->Get(TMVA_inputpath+cat_label[k]+"/Method_"+method+"/"+method+"/MVA_"+method+"_B_high");
            //Normal binning BDT scores
	    h_test_signal2 = (TH1F*)f->Get(TMVA_inputpath+cat_label[k]+"/Method_"+method+"/"+method+"/MVA_"+method+"_S");
	    h_test_bkg2 = (TH1F*)f->Get(TMVA_inputpath+cat_label[k]+"/Method_"+method+"/"+method+"/MVA_"+method+"_B");
	    h_train_signal2 = (TH1F*)f->Get(TMVA_inputpath+cat_label[k]+"/Method_"+method+"/"+method+"/MVA_"+method+"_Train_S");
	    h_train_bkg2 = (TH1F*)f->Get(TMVA_inputpath+cat_label[k]+"/Method_"+method+"/"+method+"/MVA_"+method+"_Train_B");
            h_test_signal->SetDirectory(0);
            h_test_bkg->SetDirectory(0);
            h_test_signal2->SetDirectory(0);
            h_test_bkg2->SetDirectory(0);
            h_train_signal2->SetDirectory(0);
            h_train_bkg2->SetDirectory(0);
            f->Close();
            //KolmogorovTest
            Double_t ks_signal = h_test_signal2->KolmogorovTest(h_train_signal2);
            Double_t ks_bkg = h_test_bkg2->KolmogorovTest(h_train_bkg2);
            cout<<"+++++++++ "<<cat_label[k]<<" KS test on signal "<<ks_signal<<" and bkg "<<ks_bkg<<endl;

*/
            //full path of root files containing BDT decision distribution
            TString file_name = TMVA_inputpath+cat_label[k]+"/BDTdecision_"+cat_label[k]+".root"; 
            TFile *f = new TFile(file_name,"READ");
            cout<<"Opened file "<<file_name<<endl;
            TH1F *h_test_signal;
            TH1F *h_test_bkg;
            TH1F *h_test_signal2;
            TH1F *h_test_bkg2;
            h_test_signal = (TH1F*)f->Get("BDTdecision_signal"+cat_label[k]);
            h_test_bkg = (TH1F*)f->Get("BDTdecision_data_obs"+cat_label[k]);
            h_test_signal2 = (TH1F*)f->Get("BDTdecision_signal"+cat_label[k]);
            h_test_bkg2 = (TH1F*)f->Get("BDTdecision_data_obs"+cat_label[k]);
            h_test_signal->SetDirectory(0);
            h_test_bkg->SetDirectory(0);
            h_test_signal2->SetDirectory(0);
            h_test_bkg2->SetDirectory(0);
            f->Close();

	    //Make up on plots
	    h_test_signal->GetXaxis()->SetRangeUser(-0.5,0.5);
	    h_test_bkg->SetLineColor(kBlack);
	    h_test_signal->SetLineColor(kRed);

            double X_min = std::min(h_test_signal->GetXaxis()->GetXmin(), h_test_signal->GetXaxis()->GetXmin());
            double X_max = std::max(h_test_signal->GetXaxis()->GetXmax(), h_test_signal->GetXaxis()->GetXmax());
	    
            //Compute cut and make colz plot
            BDTcut cut_value = Get_BDT_cut(cat_label[k], h_test_signal, h_test_bkg, true);
            //Draw BDT working point on ROC curve
            plot_ROC(cat_label[k], h_test_signal, h_test_bkg);

	    TCanvas *c1 = new TCanvas("c1","c1",150,10,800,800);
	    h_test_bkg->Draw("HISTE");
	    h_test_signal->Draw("same HISTE");
	    c1->Update();
            TLine l;
            l.DrawLine(cut_value.a,0,cut_value.a,0.1);
            l.DrawLine(cut_value.b,0,cut_value.b,0.1);

	    TLegend*leg = new TLegend(0.1,0.75,0.45,0.9);
	    leg->AddEntry(h_test_signal,cat_label[k]+"_signal","f");
	    leg->AddEntry(h_test_bkg,cat_label[k]+"_bkg","f");
	    leg->Draw();
	    c1->Update();
            c1->SaveAs(TMVA_inputpath+cat_label[k]+"/"+method+"_"+cat_label[k]+"normalizedSignal.png");

            //Drawing BDT score from scratch without signal normalization
	    TCanvas *c2 = new TCanvas("c2","c2",150,10,800,800);
	    h_test_signal2->GetXaxis()->SetRangeUser(-0.4,0.4);
	    h_test_bkg2->GetXaxis()->SetRangeUser(-0.4,0.4);
	    //h_test_signal2->GetXaxis()->SetRangeUser(X_min,X_max);
	    h_test_bkg2->SetLineColor(kBlack);
	    h_test_signal2->SetLineColor(kRed);
	    h_test_signal2->Scale(1/h_test_signal2->Integral());
	    h_test_bkg2->Scale(1/h_test_bkg2->Integral());
	    h_test_signal2->Rebin(4);
	    h_test_bkg2->Rebin(4);
            double Y_max = std::max(h_test_signal2->GetMaximum(), h_test_bkg2->GetMaximum());
	    h_test_bkg2->Draw("HISTE");
	    //h_test_bkg2->GetYaxis()->SetRangeUser(0, Y_max);
	    h_test_bkg2->GetYaxis()->SetRangeUser(0, 0.12);
	    h_test_bkg2->GetXaxis()->SetRangeUser(-0.4,0.4);
	    h_test_signal2->Draw("same HISTE");
	    c2->Update();
            l.DrawLine(cut_value.a,0,cut_value.a,0.12);
            l.DrawLine(cut_value.b,0,cut_value.b,0.12);
            TLatex ta(cut_value.a,0,"a");
            ta.Draw();
            TLatex tb(cut_value.b,0,"b");
            tb.Draw();

	    TLegend*leg2 = new TLegend(0.1,0.76,0.45,0.9);
	    leg2->AddEntry(h_test_signal2,cat_label[k]+"_signal","f");
	    leg2->AddEntry(h_test_bkg2,cat_label[k]+"_bkg","f");
	    leg2->Draw();
	    c2->Update();
            c2->SaveAs(TMVA_inputpath+cat_label[k]+"/"+method+"_"+cat_label[k]+".png");

     }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}

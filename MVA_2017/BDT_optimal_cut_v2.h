#ifndef BDT_OPTIMAL_CUT_H_
#define BDT_OPTIMAL_CUT_H_

#include "TH1F.h"
#include "TGraph2D.h"
#include "TLine.h"
#include "TCanvas.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <array>

using namespace std;

struct BDTcut { 
  float a; 
  float b; 
};

double TH1_integral(TH1F *h, float xmin, float xmax){
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/ axis->GetBinWidth(bmax);

    return integral;
}

double log_significance(double S, double B){
    double significance = 0;
    significance = sqrt(2*( (S+B) * log( 1+S/B ) - S));
    //cout<<"log sign is "<<significance<<" while S/sqrt(S + B) gives "<< S/sqrt(S + B) <<endl;
    return significance;
}

BDTcut Get_BDT_cut(TString categ, TH1F *h_test_signal, TH1F *h_test_bkg, bool make_plot) 
{
    float a,b;
    double N_s_1, N_b_1;
    double N_s_2, N_b_2;
    double S1, S2, S;
    std::vector<double> S1_list, S2_list, S_list, a_list, b_list;

    //Bkg is scaled depending on categ to fit signal mass window
    double bkg_scale = 1;
    if(categ.Contains("A")) bkg_scale = 4. * 12./(380.);
    if(categ.Contains("B")) bkg_scale = 4. * 19./(380.);
    if(categ.Contains("C")) bkg_scale = 4. * 25./(380.);
    cout<<"bkg_scale = "<<bkg_scale<<endl;

    double X_min = std::min(h_test_signal->GetXaxis()->GetXmin(), h_test_signal->GetXaxis()->GetXmin());
    double X_max = std::max(h_test_signal->GetXaxis()->GetXmax(), h_test_signal->GetXaxis()->GetXmax());
    //cout<<"X_min = "<<X_min<<" X_max = "<<X_max<<endl;
    //Loop on both cuts in [X_min;X_max]
    Int_t dim = 0;
    //Increase N to increase (a,b) scan granularity!
    Int_t N = 350; double step = (X_max - X_min)/N;
    for(int i=0; i<N; i++){
	a = X_min + i * step;
	for(int j=0; j<N; j++){
	    b = X_min + j * step;
            if(a<b) continue;
            //computing areas in range [a;X_max]
	    N_s_1 = TH1_integral(h_test_signal,a,X_max);
	    N_b_1 = TH1_integral(h_test_bkg,a,X_max);
            //skip iteration if integral in the tails is < 0.001% of total (sensitive to fluctuations!)
            if(N_s_1 < TH1_integral(h_test_signal,X_min,X_max)*0.00001) continue;
            if(N_b_1 < TH1_integral(h_test_bkg,X_min,X_max)*0.00001) continue;
            if(N_b_1 < 2) continue;
            //computing areas in range [b;a]
	    N_s_2 = TH1_integral(h_test_signal,b,a);
            N_b_2 = TH1_integral(h_test_bkg,b,a);
            //skip iteration if integral in cat 2 is < 0.001% of total (sensitive to fluctuations!)
            if(N_s_2 < TH1_integral(h_test_signal,X_min,X_max)*0.00001) continue;
            if(N_b_2 < TH1_integral(h_test_bkg,X_min,X_max)*0.00001) continue;
            //scale n(b) depending on categ
            N_b_1 = N_b_1*bkg_scale;
            N_b_2 = N_b_2*bkg_scale;
            if(a<b) continue;
 	    if ( (N_b_1)>0 && (N_b_2)>0 ) {
                double S = log_significance(N_s_1, N_b_1);
                S1 = log_significance(N_s_1, N_b_1);
                S2 = log_significance(N_s_2, N_b_2);
                //S1 = N_s_1 / sqrt(N_s_1 + N_b_1);
		//S2 = N_s_2 / sqrt(N_s_2 + N_b_2);
	        //Combined significance
	        S = sqrt(S1*S1 + S2*S2);
	        a_list.push_back(a);
	        b_list.push_back(b);
	        S_list.push_back(S);
                //cout<<"(a;b)=("<<a<<" ;"<<b<<") - S="<<S<<endl;
   	        dim++;
	    }
	}
    }
    //Taking absolute maximum of the combined significance
    double S_max = *max_element(S_list.begin(), S_list.end());
    int S_maxIndex = std::max_element(S_list.begin(),S_list.end()) - S_list.begin();
    float a_max = a_list.at(S_maxIndex);
    float b_max = b_list.at(S_maxIndex);
    cout<<"Category "<<categ<<endl;
    cout<<"b cut: "<<b_max<<", a cut: "<<a_max<<endl;
    cout<<"S value: "<<S_max<<endl;

    //Computing cut efficiency on signal
    Double_t N_S_12 = TH1_integral(h_test_signal,b_max,X_max);
    Double_t N_S_tot = TH1_integral(h_test_signal,X_min,X_max);
    cout<<"Signal events kept by BDT "<<N_S_12<<" over "<<N_S_tot<<" ratio: "<<N_S_12/N_S_tot<<endl;
    //Computing cut efficiency on backgroup
    Double_t N_B_12 = TH1_integral(h_test_bkg,b_max,X_max);
    Double_t N_B_tot = TH1_integral(h_test_bkg,X_min,X_max);
    cout<<"Background events kept by BDT "<<N_B_12<<" over "<<N_B_tot<<" ratio: "<<N_B_12/N_B_tot<<endl;
    
    if(make_plot){
        TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
        TGraph2D *g2d = new TGraph2D(dim, &a_list[0], &b_list[0], &S_list[0]);
        g2d->GetXaxis()->SetRangeUser(X_min, X_max);
        g2d->GetYaxis()->SetRangeUser(X_min, X_max);
        g2d->SetTitle(categ+";a;b");
        g2d->GetZaxis()->SetRangeUser(0.3, 1);
        g2d->Draw("colz");    
        g2d->GetZaxis()->SetRangeUser(0.3, 1);
        g2d->GetXaxis()->SetRangeUser(X_min, X_max);
        g2d->GetYaxis()->SetRangeUser(X_min, X_max);
        g2d->SetMinimum(0.25);
        c3->Update();
        g2d->GetXaxis()->SetRangeUser(X_min, X_max);
        g2d->GetYaxis()->SetRangeUser(X_min, X_max);
        TLine l;
        l.DrawLine(a_max,X_min,a_max,X_max);
        l.DrawLine(X_min,b_max,X_max,b_max);

	c3->Update();
	c3->SaveAs(TMVA_inputpath+categ+"/"+method+"_2Dmap_"+categ+"_v2.png");
    }

    S_list.clear();
    S1_list.clear();
    S2_list.clear();
    a_list.clear();
    b_list.clear();
    return {a_max,b_max};
}

BDTcut Get_BDT_cut_v2(TString categ, bool make_plot) 
{
    TString file_name = workdir+TMVA_inputpath+"outputTree.root";
    TH1F *h_test_signal;
    TH1F *h_test_bkg;
    TH1F *h_test_signal2;
    TH1F *h_test_bkg2;

    //open input files
    TChain *t = new TChain("outputTree");
    t->Add(file_name); 
    std::cout<<"Opened input file: "<<file_name<<std::endl;

    TString c = ""; 
    if (categ=="A") c = "0";
    if (categ=="B") c = "1";
    if (categ=="C") c = "2";
    TString signal = "weight*(isMC>0  && category =="+c+")";
    TString bkg    = "weight*(isMC==0 && category =="+c+" && isSB==1)";
    
    //bdt score distribution
    TString binning = "(240, -0.6, 0.6)"; 
    t->Draw("bdt>>h_test_bkg"+binning, bkg);
    h_test_bkg = (TH1F *)gDirectory->Get("h_test_bkg");
    h_test_bkg2 = (TH1F *)gDirectory->Get("h_test_bkg");
    t->Draw("bdt>>h_test_signal"+binning, signal);
    h_test_signal = (TH1F *)gDirectory->Get("h_test_signal");
    h_test_signal2 = (TH1F *)gDirectory->Get("h_test_signal");

    h_test_signal->SetDirectory(0);
    h_test_bkg->SetDirectory(0);
    h_test_signal2->SetDirectory(0);
    h_test_bkg2->SetDirectory(0);

//    h_test_signal->GetXaxis()->SetRangeUser(0,1.0);

    BDTcut BDTcutvalues = Get_BDT_cut(categ, h_test_signal, h_test_bkg, true); 
    return BDTcutvalues;
}


TGraph* plot_ROC(TString categ, TH1F *h_signal, TH1F *h_bkg){

    float bdt_min = std::min(h_signal->GetXaxis()->GetXmin(), h_bkg->GetXaxis()->GetXmin());
    float bdt_max = std::max(h_signal->GetXaxis()->GetXmax(), h_bkg->GetXaxis()->GetXmax());
    int n = 200; //number of points
    float step = abs(bdt_max-bdt_min)/(n*1.0);

    Double_t x[n], y[n];
    double n_true_positive = 0;
    double n_false_positive = 0;
    double n_true_negative_rej = 0;
    double n_false_negative_rej = 0;
    float bdt_wp = 0;

    for(int i = 0; i<n; i++){
        bdt_wp = bdt_min + step*i;

        n_true_positive  = TH1_integral(h_signal, bdt_wp, bdt_max);
        n_false_positive = TH1_integral(h_bkg, bdt_wp, bdt_max);
        n_true_negative_rej  = TH1_integral(h_bkg,    bdt_min, bdt_wp);
        n_false_negative_rej = TH1_integral(h_signal, bdt_min, bdt_wp);

        x[i] = n_true_positive/(n_true_positive+n_false_negative_rej);
        y[i] = n_true_negative_rej/(n_false_positive+n_true_negative_rej);
    }

    TGraph* gr = new TGraph(n,x,y);
    gr->SetTitle("ROC;Signal Efficiency;Background rejection");
    return gr;
}

array<double, 4> get_wp_XGBoost(TChain *t, BDTcut BDTcutvalues){

   TString bdt_branch = "bdt";
   TString target_branch = "target";
   TTreeReader reader (t);

   TTreeReaderValue<double> reader_bdt (reader, bdt_branch);
   TTreeReaderValue<long long> reader_classid (reader, target_branch);

   double n_b_true_positive = 0;
   double n_b_false_positive = 0;
   double n_b_true_negative_rej = 0;
   double n_b_false_negative_rej = 0;
   double n_a_true_positive = 0;
   double n_a_false_positive = 0;
   double n_a_true_negative_rej = 0;
   double n_a_false_negative_rej = 0;

   while (reader.Next()) {
       if ((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 1) {
           n_b_true_positive += 1.;
       } else if((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 0){
           n_b_false_positive += 1.;
       }
       if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 0) {
           n_b_true_negative_rej += 1.;
       } else if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 1) {
           n_b_false_negative_rej += 1.;
       }

       if ((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 1) {
           n_a_true_positive += 1.;
       } else if((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 0){
           n_a_false_positive += 1.;
       }
       if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 0) {
           n_a_true_negative_rej += 1.;
       } else if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 1) {
           n_a_false_negative_rej += 1.;
       }
    }

    array<double, 4> c;
    c[0] = n_b_true_positive/(n_b_true_positive+n_b_false_negative_rej); //x[0]
    c[1] = n_b_true_negative_rej/(n_b_false_positive+n_b_true_negative_rej); //y[0]
    c[2] = n_a_true_positive/(n_a_true_positive+n_a_false_negative_rej); //x[1]
    c[3] = n_a_true_negative_rej/(n_a_false_positive+n_a_true_negative_rej); //y[1]
    return c;
}

array<double, 4> get_wp_TMVA(TChain *t, BDTcut BDTcutvalues){

   TString bdt_branch = "BDT";
   TString target_branch = "classID";
   TTreeReader reader (t);
   
   TTreeReaderValue<float> reader_bdt (reader, bdt_branch);
   TTreeReaderValue<int> reader_classid (reader, target_branch);

   double n_b_true_positive = 0;
   double n_b_false_positive = 0;
   double n_b_true_negative_rej = 0;
   double n_b_false_negative_rej = 0;
   double n_a_true_positive = 0;
   double n_a_false_positive = 0;
   double n_a_true_negative_rej = 0;
   double n_a_false_negative_rej = 0;

   while (reader.Next()) {
       if ((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 0) {
           n_b_true_positive += 1.;
       } else if((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 1){
           n_b_false_positive += 1.;
       }
       if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 1) {
           n_b_true_negative_rej += 1.;
       } else if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 0) {
           n_b_false_negative_rej += 1.;
       }

       if ((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 0) {
           n_a_true_positive += 1.;
       } else if((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 1){
           n_a_false_positive += 1.;
       }
       if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 1) {
           n_a_true_negative_rej += 1.;
       } else if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 0) {
           n_a_false_negative_rej += 1.;
       }
    }

    array<double, 4> c;
    c[0] = n_b_true_positive/(n_b_true_positive+n_b_false_negative_rej); //x[0]
    c[1] = n_b_true_negative_rej/(n_b_false_positive+n_b_true_negative_rej); //y[0]
    c[2] = n_a_true_positive/(n_a_true_positive+n_a_false_negative_rej); //x[1]
    c[3] = n_a_true_negative_rej/(n_a_false_positive+n_a_true_negative_rej); //y[1]
    return c;
}

void plot_wp_ROC(TString categ, TH1F *h_test_signal, TH1F *h_test_bkg, bool isXGBoost){
 
     TGraph* ROC =  plot_ROC(categ, h_test_signal, h_test_bkg);
 
     TCanvas *c1 = new TCanvas("c1","ROC "+categ,150,10,990,660);
 //    ROC->GetYaxis()->SetRangeUser(0.5,1.05);
     ROC->Draw();
 
     BDTcut BDTcutvalues = Get_BDT_cut(categ, h_test_signal, h_test_bkg, false);
 
     TChain *tTrain = new TChain("tree");
 
     if(isXGBoost) {
         TString file_name_signal = "Tau3Mu_BDT_XGBoost/ntuples/signal_"+XGBoost_tag+"_"+categ+".root";
         TString file_name_bkg = "Tau3Mu_BDT_XGBoost/ntuples/background_"+XGBoost_tag+"_"+categ+".root";
         tTrain->Add(file_name_bkg); //data 
         tTrain->Add(file_name_signal); //signal
 
     } else { tTrain->Add("TMVA_"+TMVA_inputpath+categ+".root/"+TMVA_inputpath+categ+"/TrainTree"); }

     tTrain->LoadTree(-1);
     int nentries = tTrain->GetEntries();
     cout<<"nentries "<<nentries<<endl;

     array<double, 4> coordinates;
     if(isXGBoost) coordinates = get_wp_XGBoost(tTrain, BDTcutvalues);
     else coordinates = get_wp_TMVA(tTrain, BDTcutvalues);

     Double_t x[2], y[2];
     x[0] = coordinates[0];
     y[0] = coordinates[1];
     x[1] = coordinates[2];
     y[1] = coordinates[3];

     TGraph* gr = new TGraph(2,x,y);
     TLatex *l = new TLatex(0.5, 0.5, "label");
     l->SetTextSize(0.025);
     l->SetTextFont(42);
     l->SetTextAlign(21);
     l->SetTextColor(kBlue);
     l->DrawLatex(x[0],y[0]+0.01,Form("bdt.b %4.3f",BDTcutvalues.b));
     l->DrawLatex(x[0],y[0]-0.02,Form("(%4.2f, %4.2f)",x[0], y[0]));
     l->DrawLatex(x[1],y[1]+0.01,Form("bdt.a %4.3f",BDTcutvalues.a));
     l->DrawLatex(x[1],y[1]-0.02,Form("(%4.2f, %4.2f)",x[1], y[1]));
     gr->Draw("same *p");

     c1->Update();
     c1->SaveAs(TMVA_inputpath+categ+"/"+method+"_ROC_"+categ+".png");
     delete tTrain;
}
#endif // BDT_OPTIMAL_CUT_H ///:~

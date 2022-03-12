#include "TH1F.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "../T3M_common.h"
#include "plot_ROC_test_train.h"

using namespace std;

void plot_ROC_test_train() 
{

    int ncat = sizeof(cat_label) / sizeof(cat_label[0]);

    for(int k=0; k<ncat; k++)
    { //if(k!=2) continue;//only cat C
        for(int j=0; j<numFolds; j++)
        {
            TString category = cat_label[k];
            TString fold = std::to_string(j+1);

            TFile *f_tmva = new TFile("../"+TMVA_inputpath+category+"/BDTG_fold"+fold+".root","READ");
            cout<<"Opened TMVA file"<<endl;

            TTree *tTrain = (TTree*)f_tmva->Get(TMVA_inputpath+category+"/TrainTree");
            TTree *tTest = (TTree*)f_tmva->Get(TMVA_inputpath+category+"/TestTree");

            TString binning;
            binning = "(80, -0.4, 0.4)";
            TString varname = "BDTG_fold"+fold;

            TH1F *h_test_signal;
            TH1F *h_test_bkg;
            TH1F *h_train_signal;
            TH1F *h_train_bkg;

            //take out histograms 
            //Bkg
            tTrain->Draw(varname+">>h_train_bkg"+binning, "classID==1");
            h_train_bkg = (TH1F *)gDirectory->Get("h_train_bkg"); 

            tTest->Draw(varname+">>h_test_bkg"+binning, "classID==1");
            h_test_bkg = (TH1F *)gDirectory->Get("h_test_bkg"); 
            //Signal
            tTrain->Draw(varname+">>h_train_signal"+binning, "weight*(classID==0)");
            h_train_signal = (TH1F *)gDirectory->Get("h_train_signal"); 

            tTest->Draw(varname+">>h_test_signal"+binning, "weight*(classID==0)");
            h_test_signal = (TH1F *)gDirectory->Get("h_test_signal"); 

            h_test_signal->SetDirectory(0);
            h_test_bkg->SetDirectory(0);
            h_train_signal->SetDirectory(0);
            h_train_bkg->SetDirectory(0);

            //take out ROC curve for training and test sets
            TGraph* ROC_test =  plot_ROC(category, h_test_signal, h_test_bkg);
            TGraph* ROC_train =  plot_ROC(category, h_train_signal, h_train_bkg);

            ROC_test->SetLineColor(kBlue);
            ROC_test->SetLineWidth(2);
            ROC_train->SetLineColor(kBlack);
            ROC_train->SetLineWidth(2);

            ROC_train->GetXaxis()->SetRangeUser(0,0.6);            
            ROC_train->GetYaxis()->SetRangeUser(0.9,1.05);            
            TCanvas *c1 = new TCanvas("c1","c1",150,10,800,800);
            ROC_train->Draw();
            ROC_test->Draw("same");

            TLegend*leg = new TLegend(0.12,0.12,0.45,0.45);
            leg->AddEntry(ROC_train, "Train - fold "+fold,"l");
            leg->AddEntry(ROC_test, "Test - fold "+fold,"l");
            leg->Draw();
            c1->Update();
            c1->SaveAs("../"+TMVA_inputpath+cat_label[k]+"/"+method+"_"+cat_label[k]+"_fold"+fold+"_ROC.png");

         }
     }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}

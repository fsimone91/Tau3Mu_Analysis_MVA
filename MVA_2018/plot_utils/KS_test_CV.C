#include "TH1F.h"
#include "TGraph2D.h"
#include "TLine.h"
#include "TCanvas.h"

#include "TMVA/Tools.h"

#include "../T3M_common.h"

using namespace TMVA;

void KS_test_parser(TString category, TString TMVA_inputpath, TString fold){

    cout<<"Opening data file"<<endl;
    TFile *f_tmva = new TFile("../"+TMVA_inputpath+category+"/BDTG_fold"+fold+".root","READ");
    cout<<"Opened TMVA file"<<endl;

    TH1F *h_MVA_BDT_Test_S;
    TH1F *h_MVA_BDT_Train_S;
    TH1F *h_MVA_BDT_Test_B;
    TH1F *h_MVA_BDT_Train_B;

    TTree *tTrain = (TTree*)f_tmva->Get(TMVA_inputpath+category+"/TrainTree");
    TTree *tTest = (TTree*)f_tmva->Get(TMVA_inputpath+category+"/TestTree");

    TString binning;
    binning = "(40, 0.0, 0.4)";
    TString varname = "BDTG_fold"+fold;

    TString signal_label_all = "MC #tau#rightarrow3#mu";
    //take out histograms 
    //Bkg
    tTrain->Draw(varname+">>hTrain_bkg"+binning, "classID==1");
    h_MVA_BDT_Train_B = (TH1F *)gDirectory->Get("hTrain_bkg"); 

    tTest->Draw(varname+">>hTest_bkg"+binning, "classID==1");
    h_MVA_BDT_Test_B = (TH1F *)gDirectory->Get("hTest_bkg"); 

    //Signal
    tTrain->Draw(varname+">>hTrain_signal"+binning, "weight*(classID==0)");
    h_MVA_BDT_Train_S = (TH1F *)gDirectory->Get("hTrain_signal"); 
    tTest->Draw(varname+">>hTest_signal"+binning, "weight*(classID==0)");
    h_MVA_BDT_Test_S = (TH1F *)gDirectory->Get("hTest_signal"); 
    
    //SIGNAL
    TCanvas *c3 = new TCanvas("c3","2018 "+signal_label_all+" MC - cat "+category,150,10,800,800);
    gStyle->SetOptStat(0);
    TLegend*leg_signal = new TLegend(0.55,0.7,0.9,0.9);

    //Train
    h_MVA_BDT_Train_S->SetLineColor(kBlue);
    h_MVA_BDT_Train_S->SetLineWidth(2);
    h_MVA_BDT_Train_S->Scale(1/h_MVA_BDT_Train_S->GetEntries());
    h_MVA_BDT_Train_S->GetXaxis()->SetTitle("BDT score");
    //Test
    h_MVA_BDT_Test_S->SetLineColor(kBlue);
    h_MVA_BDT_Test_S->SetLineWidth(2);
    h_MVA_BDT_Test_S->Scale(1/h_MVA_BDT_Test_S->GetEntries());
    h_MVA_BDT_Test_S->GetXaxis()->SetTitle("BDT score");

    h_MVA_BDT_Train_S->SetTitle("2018 "+signal_label_all+" MC - cat "+category);
    h_MVA_BDT_Train_S->Draw("hist");
    h_MVA_BDT_Test_S->Draw("lep same");

    leg_signal->AddEntry(h_MVA_BDT_Train_S,signal_label_all+" - fold "+fold,"f");
    leg_signal->AddEntry(h_MVA_BDT_Test_S, signal_label_all+" - fold "+fold,"lep");
    leg_signal->Draw();
    c3->Update();

    //KS test
    cout<<"Signal KS test: "<<h_MVA_BDT_Test_S->KolmogorovTest(h_MVA_BDT_Train_S)<<endl;
    Double_t KS_value_S = h_MVA_BDT_Test_S->KolmogorovTest(h_MVA_BDT_Train_S);
    std::stringstream stream_ks_s;
    stream_ks_s << std::fixed << std::setprecision(4) << KS_value_S;
    std::string string_ks_s = stream_ks_s.str();
    TString KS_signal = "KS test: "+string_ks_s;
    TLatex* text_KS_S = new TLatex(0.55,0.65, KS_signal);
    text_KS_S->SetTextSize(0.040);
    text_KS_S->SetNDC(kTRUE);
    text_KS_S->Draw("same");

    c3->SaveAs("../"+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"_KStest_signalStack.png");

    //BACKGROUND
    TCanvas *c2 = new TCanvas("c2","2018 Bkg - cat "+category,150,10,800,800);
    gStyle->SetOptStat(0);
    TLegend*leg_bkg = new TLegend(0.6,0.75,0.9,0.9);

    //Train
    h_MVA_BDT_Train_B->SetLineColor(kRed);
    h_MVA_BDT_Train_B->SetLineWidth(2);
    h_MVA_BDT_Train_B->Scale(1/h_MVA_BDT_Train_B->GetEntries());
    h_MVA_BDT_Train_B->GetXaxis()->SetTitle("BDT score");
    c2->cd();
    h_MVA_BDT_Train_B->Draw("hist");
    leg_bkg->AddEntry(h_MVA_BDT_Train_B,"Bkg. train - fold "+fold,"f");
    //Test
    h_MVA_BDT_Test_B->SetLineColor(kRed);
    h_MVA_BDT_Test_B->SetLineWidth(2);
    h_MVA_BDT_Test_B->Scale(1/h_MVA_BDT_Test_B->GetEntries());
    h_MVA_BDT_Test_B->GetXaxis()->SetTitle("BDT score");
    c2->cd();
    h_MVA_BDT_Test_B->Draw("lep same");
    leg_bkg->AddEntry(h_MVA_BDT_Test_B,"Bkg. test - fold "+fold,"lep");
    leg_bkg->Draw();
    c2->Update();

    //KS test
    cout<<"Background KS test: "<<h_MVA_BDT_Test_B->KolmogorovTest(h_MVA_BDT_Train_B)<<endl;
    Double_t KS_value_B = h_MVA_BDT_Test_B->KolmogorovTest(h_MVA_BDT_Train_B);
    std::stringstream stream_ks_b;
    stream_ks_b << std::fixed << std::setprecision(4) << KS_value_B;
    std::string string_ks_b = stream_ks_b.str();
    TString KS_bkg = "KS test: "+string_ks_b;
    TLatex* text_KS_B = new TLatex(0.6,0.65, KS_bkg);
    text_KS_B->SetTextSize(0.045);
    text_KS_B->SetNDC(kTRUE);
    text_KS_B->Draw("same");

    c2->SaveAs("../"+TMVA_inputpath+category+"/"+TMVA_inputpath+category+"_"+varname+"_KStest_bkg.png");

    f_tmva->Close();
}

void KS_test_CV(){
    KS_test_parser("C",TMVA_inputpath, "1");   
    KS_test_parser("C",TMVA_inputpath, "2");   
    KS_test_parser("C",TMVA_inputpath, "3");   
}


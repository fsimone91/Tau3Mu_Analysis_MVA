#include "TH1F.h"
#include <cmath>
#include <string> 
#include "../T3M_common.h"

//macro modified to compare different MC (DsTau3Mu with shifted tau masses)
void plot_BDT_inputvariables_MCcomp() 
{
    TString inputpath_Ds_1p77 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201211_1617/AnalysedTree_MC_2018Ds_tau3mu_11Dec_1p77.root";
    TString inputpath_Ds_1p67 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201211_1621/AnalysedTree_MC_2018Ds_tau3mu_11Dec_1p67.root";
    TString inputpath_Ds_1p87 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201211_1620/AnalysedTree_MC_2018Ds_tau3mu_11Dec_1p87.root";

    TString cat_label[] = {"A", "B", "C"};
    TString var[] = {
                     "Pmu3","cLP","tKink","segmComp","fv_nC","fv_dphi3D","fv_d3Dsig","d0sig","mindca_iso","trkRel","nMatchesMu3",
                     "fv_d3D",
                   //  "abs(dxy1/dxyErr1)", "abs(dxy2/dxyErr2)", "abs(dxy3/dxyErr3)",
                   //  "TreeMu1.mu_SoftMVA",
                   //  "TreeMu2.mu_SoftMVA",
                   //  "TreeMu3.mu_SoftMVA",
                     "MuonIDeval_Mu1.MuonID",
                     "MuonIDeval_Mu2.MuonID",
                     "MuonIDeval_Mu3.MuonID"
                     };

    //MC Ds - 
    //tau mass 1.77
    TChain *tmc_1p77 = new TChain("FinalTree");
    tmc_1p77->Add(inputpath_Ds_1p77); 
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    //tau mass 1.67
    TChain *tmc_1p67 = new TChain("FinalTree");
    tmc_1p67->Add(inputpath_Ds_1p67); 
    std::cout<<"Opened input file: "<<inputpath_Ds_1p67<<std::endl;
    //tau mass 1.87
    TChain *tmc_1p87 = new TChain("FinalTree");
    tmc_1p87->Add(inputpath_Ds_1p87); 
    std::cout<<"Opened input file: "<<inputpath_Ds_1p87<<std::endl;

    //MC Ds, B0, Bp
    TChain *tmc_1p77_mu1 = new TChain("TreeMu1");
    TChain *tmc_1p67_mu1 = new TChain("TreeMu1");
    TChain *tmc_1p87_mu1 = new TChain("TreeMu1");
    tmc_1p77_mu1->Add(inputpath_Ds_1p77);
    tmc_1p67_mu1->Add(inputpath_Ds_1p67);
    tmc_1p87_mu1->Add(inputpath_Ds_1p87);
    TChain *tmc_1p77_mu2 = new TChain("TreeMu2");
    TChain *tmc_1p67_mu2 = new TChain("TreeMu2");
    TChain *tmc_1p87_mu2 = new TChain("TreeMu2");
    tmc_1p77_mu2->Add(inputpath_Ds_1p77);
    tmc_1p67_mu2->Add(inputpath_Ds_1p67);
    tmc_1p87_mu2->Add(inputpath_Ds_1p87);
    TChain *tmc_1p77_mu3 = new TChain("TreeMu3");
    TChain *tmc_1p67_mu3 = new TChain("TreeMu3");
    TChain *tmc_1p87_mu3 = new TChain("TreeMu3");
    tmc_1p77_mu3->Add(inputpath_Ds_1p77);
    tmc_1p67_mu3->Add(inputpath_Ds_1p67);
    tmc_1p87_mu3->Add(inputpath_Ds_1p87);

    //MC Ds, B0, Bp
    TString inputpath_Ds_1p77_muId = inputpath_Ds_1p77.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_Ds_1p67_muId = inputpath_Ds_1p67.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_Ds_1p87_muId = inputpath_Ds_1p87.ReplaceAll(".root", "_MuonID.root");
    TChain *tmc_1p77_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc_1p67_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc_1p87_muid1 = new TChain("MuonIDeval_Mu1");
    tmc_1p77_muid1->Add(inputpath_Ds_1p77_muId);
    tmc_1p67_muid1->Add(inputpath_Ds_1p67_muId);
    tmc_1p87_muid1->Add(inputpath_Ds_1p87_muId);
    TChain *tmc_1p77_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc_1p67_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc_1p87_muid2 = new TChain("MuonIDeval_Mu2");
    tmc_1p77_muid2->Add(inputpath_Ds_1p77_muId);
    tmc_1p67_muid2->Add(inputpath_Ds_1p67_muId);
    tmc_1p87_muid2->Add(inputpath_Ds_1p87_muId);
    TChain *tmc_1p77_muid3 = new TChain("MuonIDeval_Mu3");
    TChain *tmc_1p67_muid3 = new TChain("MuonIDeval_Mu3");
    TChain *tmc_1p87_muid3 = new TChain("MuonIDeval_Mu3");
    tmc_1p77_muid3->Add(inputpath_Ds_1p77_muId);
    tmc_1p67_muid3->Add(inputpath_Ds_1p67_muId);
    tmc_1p87_muid3->Add(inputpath_Ds_1p87_muId);
    std::cout<<"Opened input file: "<<inputpath_Ds_1p77_muId<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_Ds_1p67_muId<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_Ds_1p87_muId<<std::endl;

    tmc_1p77->AddFriend(tmc_1p77_mu1);
    tmc_1p77->AddFriend(tmc_1p77_mu2);
    tmc_1p77->AddFriend(tmc_1p77_mu3);
    tmc_1p77->AddFriend(tmc_1p77_muid1);
    tmc_1p77->AddFriend(tmc_1p77_muid2);
    tmc_1p77->AddFriend(tmc_1p77_muid3);

    tmc_1p67->AddFriend(tmc_1p67_mu1);
    tmc_1p67->AddFriend(tmc_1p67_mu2);
    tmc_1p67->AddFriend(tmc_1p67_mu3);
    tmc_1p67->AddFriend(tmc_1p67_muid1);
    tmc_1p67->AddFriend(tmc_1p67_muid2);
    tmc_1p67->AddFriend(tmc_1p67_muid3);

    tmc_1p87->AddFriend(tmc_1p87_mu1);
    tmc_1p87->AddFriend(tmc_1p87_mu2);
    tmc_1p87->AddFriend(tmc_1p87_mu3);
    tmc_1p87->AddFriend(tmc_1p87_muid1);
    tmc_1p87->AddFriend(tmc_1p87_muid2);
    tmc_1p87->AddFriend(tmc_1p87_muid3);

    int n = sizeof(var)/sizeof(var[0]);
    TH1F *hmc_1p77[3];
    TH1F *hmc_1p77_merged;

    TH1F *hmc_1p67[3];
    TH1F *hmc_1p67_merged;

    TH1F *hmc_1p87[3];
    TH1F *hmc_1p87_merged;

    TString binning;

    TCut cutS = "tripletMass<2.0 && tripletMass>1.62"; //Signal -> MC full range 
    TCut cutB = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)"; //Background -> data sidebands

    TCut reso_A = "tripletMassReso < 0.007";
    TCut reso_B = "tripletMassReso >= 0.007 && tripletMassReso <= 0.0105";
    TCut reso_C = "tripletMassReso > 0.0105";
    TCut reso_cat = "tripletMassReso < 0"; //always false 

    //Loop on variables
    for(int i = 0; i<n; i++){
        TString varname = var[i];
        cout<<"Input variable "<<varname<<endl;
        binning = "";
        if(varname=="Pmu3") binning = "(50,0,50)";
        if(varname=="cLP") binning = "(100,0,100)";
        if(varname=="segmComp") binning = "(100,-0.1,1.1)";
        if(varname=="tKink") binning = "(150,0,150)";
        if(varname=="fv_nC") binning = "(60,-0.5,15.5)";
        if(varname=="fv_dphi3D") binning = "(100,-0.1,0.3)";
        if(varname=="fv_d3Dsig") binning = "(200,-0.1,200)";
        if(varname=="d0sig") binning = "(60,-0.1,30)";
        if(varname=="mindca_iso") binning = "(100,-0.1,1)";
        if(varname=="trkRel") binning = "(100,-0.5,10)";
        if(varname=="nMatchesMu3") binning = "(10,0,10)";
        if(varname=="fv_d3D") binning = "(48,0,6)";
        if(varname.Contains("mu_SoftMVA")) binning = "(200,-1,1)";
        if(varname.Contains("MuonID"))     binning = "(200,-1,1)";
 
        //Loop on categories A, B, C
        for(auto j = 0; j<3; j++){
            TString s = std::to_string(j);
            TString category = cat_label[j];
            cout<<"Category "<<category<<endl;
            if(category == "A") reso_cat = reso_cat || reso_A;
            if(category == "B") reso_cat = reso_cat || reso_B;
            if(category == "C") reso_cat = reso_cat || reso_C;

            tmc_1p77->Draw(varname+">>hmc_1p77"+s+binning, cutS&&reso_cat);
            hmc_1p77[j] = (TH1F *)gDirectory->Get("hmc_1p77"+s);

            tmc_1p67->Draw(varname+">>hmc_1p67"+s+binning, cutS&&reso_cat);
            hmc_1p67[j] = (TH1F *)gDirectory->Get("hmc_1p67"+s);

            tmc_1p87->Draw(varname+">>hmc_1p87"+s+binning, cutS&&reso_cat);
            hmc_1p87[j] = (TH1F *)gDirectory->Get("hmc_1p87"+s);

            reso_cat = "tripletMassReso < 0";
        }
        tmc_1p77->Draw(varname+">>hmc_1p77_merged"+binning, cutS);
        hmc_1p77_merged = (TH1F *)gDirectory->Get("hmc_1p77_merged");
 
        tmc_1p67->Draw(varname+">>hmc_1p67_merged"+binning, cutS);
        hmc_1p67_merged = (TH1F *)gDirectory->Get("hmc_1p67_merged");

        tmc_1p87->Draw(varname+">>hmc_1p87_merged"+binning, cutS);
        hmc_1p87_merged = (TH1F *)gDirectory->Get("hmc_1p87_merged");

        //Loop on categories A, B, C
        for(auto j = 0; j<3; j++){
           TString s = std::to_string(j);
           TString category = cat_label[j];
           cout<<"Category "<<category<<endl;
           //plot categories on same canvas - MC
           TCanvas *c2 = new TCanvas("c2","c2",150,10,990,660);
           gStyle->SetOptTitle(0);
           gStyle->SetOptStat(0);
           hmc_1p77[j]->SetLineColor(kBlue);
           hmc_1p67[j]->SetLineColor(kRed);
           hmc_1p87[j]->SetLineColor(kBlack);
    
           Double_t min_mc = std::min( std::min(hmc_1p77[j]->GetXaxis()->GetXmin(),hmc_1p67[j]->GetXaxis()->GetXmin()), hmc_1p87[j]->GetXaxis()->GetXmin() );
           min_mc = min_mc - min_mc*0.1;
           Double_t max_mc = std::max( std::max(hmc_1p77[j]->GetXaxis()->GetXmax(),hmc_1p67[j]->GetXaxis()->GetXmax()), hmc_1p87[j]->GetXaxis()->GetXmax() );
           max_mc = max_mc + max_mc*0.1;
    
           hmc_1p77[j]->GetXaxis()->SetRangeUser(min_mc, max_mc);
           hmc_1p77[j]->GetXaxis()->SetTitle(varname);

           hmc_1p77[j]->Scale(1.0/( hmc_1p77[j]->Integral() ));
           hmc_1p67[j]->Scale(1.0/( hmc_1p67[j]->Integral() ));
           hmc_1p87[j]->Scale(1.0/( hmc_1p87[j]->Integral() ));
    
           hmc_1p77[j]->Draw("hist");
           hmc_1p67[j]->Draw("same hist");
           hmc_1p87[j]->Draw("same hist");
    
           TLegend*leg2 = new TLegend(0.4,0.7,0.7,0.9);
           leg2->AddEntry(hmc_1p77[j],"2018 "+varname+"_DsTau3Mu_1p77_"+category,"f");
           leg2->AddEntry(hmc_1p67[j],"2018 "+varname+"_DsTau3Mu_1p67_"+category,"f");
           leg2->AddEntry(hmc_1p87[j],"2018 "+varname+"_DsTau3Mu_1p87_"+category,"f");
           leg2->Draw();
    
           //c2->SetLogy();
           c2->Update();
           c2->SaveAs("../plots/"+TMVA_inputpath+"_MC_"+varname+"_"+category+".png");
        }
        //TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
        //gStyle->SetOptTitle(0);
        //gStyle->SetOptStat(0);
        //hmc_merged->SetLineColor(kBlue);
        //hmc_merged->SetLineWidth(3);
        //hmc_merged->SetFillStyle(3004);
        //hmc_merged->SetFillColor(kBlue);
        //hdata_merged->SetLineColor(kRed);
        //hdata_merged->SetLineWidth(3);
        //hdata_merged->SetFillStyle(3005);
        //hdata_merged->SetFillColor(kRed);

        //Double_t min_compare = std::min(min_data, min_mc);
        //Double_t max_compare = std::max(max_data, max_mc);

        //hdata_merged->GetXaxis()->SetRangeUser(min_compare, max_compare);
        //hdata_merged->GetXaxis()->SetTitle(varname);
        //hdata_merged->Scale(1.0/(hdata_merged->Integral()));

        //hmc_merged->GetXaxis()->SetRangeUser(min_compare, max_compare);
        //hmc_merged->GetXaxis()->SetTitle(varname);
        //hmc_merged->Scale(1.0/(hmc_merged->Integral()));

        //c3->Update();
        //hmc_merged->Draw();
        //hdata_merged->Draw("same");

        //TLegend*leg3 = new TLegend(0.2,0.7,0.45,0.9);
        //leg3->AddEntry(hmc_merged,"2018 "+varname+" mc","f");
        //leg3->AddEntry(hdata_merged,"2018 "+varname+" data","f");
        //leg3->Draw();

        ////c3->SetLogy();
        //c3->Update();
        //c3->SaveAs("../plots/"+TMVA_inputpath+"_compare_"+varname+".png");

    }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}


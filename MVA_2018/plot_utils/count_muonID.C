#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

double log_significance(double S, double B){
    double significance = 0;
    significance = sqrt(2*( (S+B) * log( 1+S/B ) - S));
    //cout<<"log sign is "<<significance<<" while S/sqrt(S + B) gives "<< S/sqrt(S + B) <<endl;
    return significance;
}
void count_muonID() 
{
    TString run[] = {"A", "B", "C", "D"};
    TString category[] = {"A", "B", "C"};

    //open input files 
    //data eraA
    TFile *f_data0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[0]);
    if (!f_data0 || !f_data0->IsOpen()) f_data0 = new TFile(inputpath_datarun[0]);
    std::cout<<"Opened input file: "<<inputpath_datarun[0]<<std::endl;
    //data eraB
    TFile *f_data1 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[1]);
    if (!f_data1 || !f_data1->IsOpen()) f_data1 = new TFile(inputpath_datarun[1]);
    std::cout<<"Opened input file: "<<inputpath_datarun[1]<<std::endl;
    //data eraC
    TFile *f_data2 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[2]);
    if (!f_data2 || !f_data2->IsOpen()) f_data2 = new TFile(inputpath_datarun[2]);
    std::cout<<"Opened input file: "<<inputpath_datarun[2]<<std::endl;
    //data eraD
    TFile *f_data3 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[3]);
    if (!f_data3 || !f_data3->IsOpen()) f_data3 = new TFile(inputpath_datarun[3]);
    std::cout<<"Opened input file: "<<inputpath_datarun[3]<<std::endl;
    //MC Ds
    TFile *f_mc1 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f_mc1 || !f_mc1->IsOpen()) f_mc1 = new TFile(inputpath_Ds);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    //MC B0
    TFile *f_mc2 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
    if (!f_mc2 || !f_mc2->IsOpen()) f_mc2 = new TFile(inputpath_B0);
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    //MC Bp
    TFile *f_mc3 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    if (!f_mc3 || !f_mc3->IsOpen()) f_mc3 = new TFile(inputpath_Bp);
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;

    //Access TTrees
    TTree *tdata0 = (TTree*)f_data0->Get("FinalTree");
    TTree *tdata1 = (TTree*)f_data1->Get("FinalTree");
    TTree *tdata2 = (TTree*)f_data2->Get("FinalTree");
    TTree *tdata3 = (TTree*)f_data3->Get("FinalTree");
    
    TTree *tmcDs = (TTree*)f_mc1->Get("FinalTree");
    TTree *tmcB0 = (TTree*)f_mc2->Get("FinalTree");
    TTree *tmcBp = (TTree*)f_mc3->Get("FinalTree");

    TCut full = "tripletMass<2.0 && tripletMass>1.62";
    TCut sideband = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)";
    TCut massreso[] = {
                       "tripletMassReso < 0.007",
                       "tripletMassReso >= 0.007 && tripletMassReso <= 0.0105",
                       "tripletMassReso > 0.0105"
    };
    TCut muonID[] = {
                     "isPF3",
                     "isGlb3",
                     "isTracker3",
                     "isLoose3",
                     "isSoft3",
                     "isMedium3",
                     "isPF3 && isGlb3",
                     "isPF3 && isTracker3",
                     "isPF3 && isLoose3",
                     "isPF3 && isSoft3",
                     "isPF3 && isMedium3",
    };
    size_t n = sizeof(muonID)/sizeof(muonID[0]);
    float sensitivity[n];

    TString varname = "tripletMass";

    //count reparately for each category
    for(auto j = 0; j<3; j++){
        cout<<"category "<<category[j]<<endl;
        cout<<"MuonID\tsignificance\tbkg\tsignal"<<endl; 
        //for each muonID i count the events, sum up the entries, apply normalisation and compute the sensitivity before BDT
        for(auto i = 0; i<n; i++){

            tdata0->Draw(varname+">>hdata0", sideband && muonID[i] && massreso[j]);
            tdata1->Draw(varname+">>hdata1", sideband && muonID[i] && massreso[j]);
            tdata2->Draw(varname+">>hdata2", sideband && muonID[i] && massreso[j]);
            tdata3->Draw(varname+">>hdata3", sideband && muonID[i] && massreso[j]);

            tmcDs->Draw(varname+">>hmcDs", full && muonID[i] && massreso[j]);
            tmcB0->Draw(varname+">>hmcB0", full && muonID[i] && massreso[j]);
            tmcBp->Draw(varname+">>hmcBp", full && muonID[i] && massreso[j]);

            TH1F *hdata0 = (TH1F*)gDirectory->Get("hdata0");
            TH1F *hdata1 = (TH1F*)gDirectory->Get("hdata1");
            TH1F *hdata2 = (TH1F*)gDirectory->Get("hdata2");
            TH1F *hdata3 = (TH1F*)gDirectory->Get("hdata3");

            TH1F *hmcDs = (TH1F*)gDirectory->Get("hmcDs");
            TH1F *hmcB0 = (TH1F*)gDirectory->Get("hmcB0");
            TH1F *hmcBp = (TH1F*)gDirectory->Get("hmcBp");

            //cout<<hdata0->GetEntries()<<"\t"<<hdata1->GetEntries()<<"\t"<<hdata2->GetEntries()<<"\t"<<hdata3->GetEntries()<<endl; 
            //cout<<hmcDs->GetEntries()<<"\t"<<hmcB0->GetEntries()<<"\t"<<hmcBp->GetEntries()<<endl; 
            float entries_data = (hdata0->GetEntries()) + (hdata1->GetEntries()) + (hdata2->GetEntries()) + (hdata3->GetEntries());
            float entries_signal = (hmcDs->GetEntries()*wNormDs) + (hmcB0->GetEntries()*wNormB0) + (hmcBp->GetEntries()*wNormBp);

            //normalisation data //380MeV total //330MeV sideband //50MeV signal region
            entries_data = entries_data / 330.0 * 50.0;

            sensitivity[i] = entries_signal / sqrt(entries_data); 
            double signif = log_significance(entries_signal, entries_data); 
            cout<<muonID[i]<<"\t"<<signif<<"\t"<<entries_data<<"\t"<<entries_signal<<endl; 

        }
    cout<<"__________________________________________"<<endl;
    }
    cout<<"all categories"<<endl;
    cout<<"muonID\tsignificance\tbkg\tsignal"<<endl; 
    //for each muonID i count the events, sum up the entries, apply normalisation and compute the sensitivity before BDT
    for(auto i = 0; i<n; i++){

        tdata0->Draw(varname+">>hdata0", sideband && muonID[i]);
        tdata1->Draw(varname+">>hdata1", sideband && muonID[i]);
        tdata2->Draw(varname+">>hdata2", sideband && muonID[i]);
        tdata3->Draw(varname+">>hdata3", sideband && muonID[i]);

        tmcDs->Draw(varname+">>hmcDs", full && muonID[i]);
        tmcB0->Draw(varname+">>hmcB0", full && muonID[i]);
        tmcBp->Draw(varname+">>hmcBp", full && muonID[i]);

        TH1F *hdata0 = (TH1F*)gDirectory->Get("hdata0");
        TH1F *hdata1 = (TH1F*)gDirectory->Get("hdata1");
        TH1F *hdata2 = (TH1F*)gDirectory->Get("hdata2");
        TH1F *hdata3 = (TH1F*)gDirectory->Get("hdata3");

        TH1F *hmcDs = (TH1F*)gDirectory->Get("hmcDs");
        TH1F *hmcB0 = (TH1F*)gDirectory->Get("hmcB0");
        TH1F *hmcBp = (TH1F*)gDirectory->Get("hmcBp");

        //cout<<hdata0->GetEntries()<<"\t"<<hdata1->GetEntries()<<"\t"<<hdata2->GetEntries()<<"\t"<<hdata3->GetEntries()<<endl; 
        //cout<<hmcDs->GetEntries()<<"\t"<<hmcB0->GetEntries()<<"\t"<<hmcBp->GetEntries()<<endl; 
        float entries_data = (hdata0->GetEntries()) + (hdata1->GetEntries()) + (hdata2->GetEntries()) + (hdata3->GetEntries());
        float entries_signal = (hmcDs->GetEntries()*wNormDs) + (hmcB0->GetEntries()*wNormB0) + (hmcBp->GetEntries()*wNormBp);

        //normalisation data //380MeV total //330MeV sideband //50MeV signal region
        entries_data = entries_data / 330.0 * 50.0;

        sensitivity[i] = entries_signal / sqrt(entries_data); 
        double signif = log_significance(entries_signal, entries_data); 
        cout<<muonID[i]<<"\t"<<signif<<"\t"<<entries_data<<"\t"<<entries_signal<<endl; 
        cout<<muonID[i]<<"\t"<<sensitivity[i]<<"\t"<<entries_data<<"\t"<<entries_signal<<endl; 

    }
}

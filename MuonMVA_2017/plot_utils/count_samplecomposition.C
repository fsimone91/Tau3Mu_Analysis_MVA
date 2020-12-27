#include "TH1F.h"
#include <cmath>
#include <string> 
#include "../MuonID_common.h"

void count_samplecomposition() 
{
    //number of signal samples
    size_t nbkg = sizeof(inputpath_MuonID_bkg)/sizeof(inputpath_MuonID_bkg[0]);
    //open input files
    //bkg
    TChain *tbkg = new TChain("FinalTree");
    for(auto i=0; i<nbkg; i++){
        tbkg->Add(inputpath_MuonID_bkg[i]);
        std::cout<<"Opened input file: "<<inputpath_MuonID_bkg[i]<<std::endl;
    }
    //signal sgn 
    TChain *tsgn = new TChain("FinalTree");
    tsgn->Add(inputpath_MuonID_Ds);
    tsgn->Add(inputpath_MuonID_B0);
    tsgn->Add(inputpath_MuonID_Bp);
    std::cout<<"Opened input file: "<<inputpath_MuonID_Ds<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_MuonID_B0<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_MuonID_Bp<<std::endl;

    //numbero of input variables
    TString varname = "mu_pt";

    TH1F *hbkg;
    TH1F *hsgn;

    TString binning = "";
    Int_t n_bkg = 0, n_sgn = 0;
    TString title_bkg, title_sgn;

    //signal muons
    TCut cutS = "abs(mu_simPdgId)==13 && abs(mu_simMotherPdgId)==15"; 
    //bkg muons = pions, kaon, muons frmo decay in flight
    TCut cutB_pi = "( abs(mu_simPdgId) == 211 )"; //pion
    TCut cutB_k = "( abs(mu_simPdgId) == 321 )"; //kaon
    TCut cutB_mu_from_pi = "( abs(mu_simPdgId) == 13 && abs(mu_simMotherPdgId) == 211 )"; //true mu from pion decay
    TCut cutB_mu_from_k = "( abs(mu_simPdgId) == 13 && abs(mu_simMotherPdgId) == 321 )"; //true mu from kaon decay
    //running on global muons with associated simInfo
    TCut preselCutS = "mu_simType != 0 && mu_isGlobal == 1 && mu_pt > 2 && abs(mu_eta)<2.4";
    TCut preselCutB = "ptetaWeight*(mu_simType != 0 && mu_isGlobal == 1 && mu_pt > 2 && abs(mu_eta)<2.4)";

    TCut eta = "mu_eta<1.2";

    cout<<"BARREL:"<<endl;
    //BKG
    cout<<"BACKGROUND:"<<endl;
    //preselections 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //pions 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_pi);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //kaons 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_k);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //mu from pions 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_mu_from_pi);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //mu from kaons
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_mu_from_k);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;

    //SIGNAL
    cout<<"\nSIGNAL:"<<endl;
    //preselections 
    tsgn->Draw(varname+">>hsgn"+binning, preselCutS&&eta);
    hsgn = (TH1F *)gDirectory->Get("hsgn");
    n_sgn = hsgn->GetEntries();
    title_sgn = hsgn->GetTitle();
    cout<<title_sgn<<" | "<<n_sgn<<endl;
    //tau3mu 
    tsgn->Draw(varname+">>hsgn"+binning, preselCutS&&eta&&cutS);
    hsgn = (TH1F *)gDirectory->Get("hsgn");
    n_sgn = hsgn->GetEntries();
    title_sgn = hsgn->GetTitle();
    cout<<title_sgn<<" | "<<n_sgn<<endl;

    eta = "mu_eta>1.2";

    cout<<"ENDCAP:"<<endl;
    //BKG
    cout<<"BACKGROUND:"<<endl;
    //preselections 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //pions 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_pi);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //kaons 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_k);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //mu from pions 
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_mu_from_pi);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;
    //mu from kaons
    tbkg->Draw(varname+">>hbkg"+binning, preselCutB&&eta&&cutB_mu_from_k);
    hbkg = (TH1F *)gDirectory->Get("hbkg");
    n_bkg = hbkg->GetEntries();
    title_bkg = hbkg->GetTitle();
    cout<<title_bkg<<" | "<<n_bkg<<endl;

    //SIGNAL
    cout<<"\nSIGNAL:"<<endl;
    //preselections 
    tsgn->Draw(varname+">>hsgn"+binning, preselCutS&&eta);
    hsgn = (TH1F *)gDirectory->Get("hsgn");
    n_sgn = hsgn->GetEntries();
    title_sgn = hsgn->GetTitle();
    cout<<title_sgn<<" | "<<n_sgn<<endl;
    //tau3mu 
    tsgn->Draw(varname+">>hsgn"+binning, preselCutS&&eta&&cutS);
    hsgn = (TH1F *)gDirectory->Get("hsgn");
    n_sgn = hsgn->GetEntries();
    title_sgn = hsgn->GetTitle();
    cout<<title_sgn<<" | "<<n_sgn<<endl;

    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}


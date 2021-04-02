#include <stdio.h>
#include "../MuonID_common.h"

//root -l -b Add_branch_weights.cpp

void upd(TString filename) { 

    bool isT3M = false;
    if(filename.Contains("2018Ds") || filename.Contains("2018B0") || filename.Contains("2018Bp")) isT3M = true;
    if(filename.Contains("2017Ds") || filename.Contains("2017B0") || filename.Contains("2017Bp")) isT3M = true;

    TFile *f = new TFile(filename,"update"); 
    TTree *T = (TTree*)f->Get("FinalTree"); 

    Double_t pt,eta; 
    Double_t myweight_pt, myweight_eta;
    Double_t myweight_pteta;
    //TBranch *b_wpt = T->Branch("myweight_pt",&myweight_pt,"myweight_pt/D"); 
    //TBranch *b_weta = T->Branch("myweight_eta",&myweight_eta,"myweight_eta/D"); 
    TBranch *b_wpteta = T->Branch("myweight_pteta",&myweight_pteta,"myweight_pteta/D"); 
    T->SetBranchAddress("mu_pt",&pt); 
    T->SetBranchAddress("mu_eta",&eta);

    TFile *f_weights_barrel = TFile::Open("/lustrehome/fsimone/MuonID_study/MuonMVA_2017/Phase_space_reweighing_tools/PT_reweighting_muonid_barrel.root", "read");
    TH2D *h_weights_barrel = (TH2D*)f_weights_barrel->Get("hweight_pteta");
    TFile *f_weights_endcap = TFile::Open("/lustrehome/fsimone/MuonID_study/MuonMVA_2017/Phase_space_reweighing_tools/PT_reweighting_muonid_endcap.root", "read");
    TH2D *h_weights_endcap = (TH2D*)f_weights_endcap->Get("hweight_pteta");
    Int_t bin_pt, bin_eta;

    Long64_t nentries = T->GetEntries(); 
    for (Long64_t i=0; i<nentries; i++) { 
        T->GetEntry(i);
        //compute weights
        if(!isT3M){
            if(abs(eta)<1.2) {
                bin_pt =  h_weights_barrel->GetXaxis()->FindBin(pt); 
                bin_eta = h_weights_barrel->GetYaxis()->FindBin(eta); 
                myweight_pteta = h_weights_barrel->GetBinContent(bin_pt, bin_eta);
            }
            else {
                bin_pt =  h_weights_endcap->GetXaxis()->FindBin(pt); 
                bin_eta = h_weights_endcap->GetYaxis()->FindBin(eta); 
                myweight_pteta = h_weights_endcap->GetBinContent(bin_pt, bin_eta);
            }
            cout<<"mu_pt "<<pt<<" mu_eta "<<eta<<" bin_pt "<<bin_pt<<" bin_eta "<<bin_eta<<" myweight_pteta "<<myweight_pteta<<endl;
       }else{myweight_pteta = 1;}

       //b_wpt->Fill();
       //b_weta->Fill();
       b_wpteta->Fill();
     } 
     f_weights_barrel->Close();
     f_weights_endcap->Close();

     f->cd();
     T->Print();
     T->Write(); 
     f->Close();
}

void Add_branch_weights(){

    //loop input files
    size_t n_bkg = sizeof(inputpath_MuonID_bkg)/sizeof(inputpath_MuonID_bkg[0]);
    //Bkg samples
    for(auto i=0; i<n_bkg; i++){
        upd(inputpath_MuonID_bkg[i]);
        std::cout<<"Updated input file: "<<inputpath_MuonID_bkg[i]<<std::endl;
    }
    //MC Tau3Mu
    upd(inputpath_MuonID_Ds);
    std::cout<<"Updated input file: "<<inputpath_MuonID_Ds<<std::endl;
    upd(inputpath_MuonID_B0);
    std::cout<<"Updated input file: "<<inputpath_MuonID_B0<<std::endl;
    upd(inputpath_MuonID_Bp);
    std::cout<<"Updated input file: "<<inputpath_MuonID_Bp<<std::endl;
}

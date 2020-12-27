#include <stdio.h>

//root -l PileUp_Reweighting.cpp\(\"2018Ds\"\)

void PileUp_Reweighting(TString datasetName, TString option){
    
    if (option=="PUint") cout<<"Distributio of number of PileUp interactions."<<endl;
    else if (option=="nPV") cout<<"Distribution of number of primary vertices."<<endl;
    else cout<<"Choose between \"PUint\" or \"nPV\" "<<endl;

    TTree *tree;
    TString filename1, filename2;

    //2017
    if(datasetName == "2017Ds") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_1333/AnalysedTree_MC_2017Ds_tau3mu_14oct.root"; //MC DsTau3Mu
    if(datasetName == "2017B0") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_1341/AnalysedTree_MC_2017B0_tau3mu_14oct.root"; //MC B0Tau3Mu
    if(datasetName == "2017Bp") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_1342/AnalysedTree_MC_2017Bp_tau3mu_14oct.root"; //MC B0Tau3Mu
    if(datasetName == "2017DsPhiPi") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200928_1817/AnalysedTree_MC_2017DsPhiPi_control_28set.root"; // MC_DsPhiPi


    //2018
    //if(datasetName == "2018Ds") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1141/AnalysedTree_MC_2018Ds_tau3mu_pileup_21march.root";
    //if(datasetName == "2018B0") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1145/AnalysedTree_MC_2018B0_tau3mu_pileup_21march.root";
    //if(datasetName == "2018Bp") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1146/AnalysedTree_MC_2018Bp_tau3mu_pileup_21march.root";
    //if(datasetName == "2018DsPhiPi") filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1146/AnalysedTree_MC_2018Bp_tau3mu_pileup_21march.root";

    TString dataSample = datasetName;
    dataSample.ReplaceAll("DsPhiPi","");
    dataSample.ReplaceAll("Ds","");
    dataSample.ReplaceAll("B0","");
    dataSample.ReplaceAll("Bp","");

    //distribution from data if using nPU interactions
    if(option=="PUint"){
        filename2 = "MyDataPileupHistogram_";
        filename2 += dataSample;   filename2 += ".root"; // Data
    }

    //distribution from data if using nPrimary Vertices
    if(option=="nPV") filename2 = "/lustrehome/fsimone/MVA_2018/AnalysedTree_data_"+dataSample+"_merged_tau3mu.root";

    //open files
    TFile *f1 = new TFile(filename1);
    TFile *f2 = new TFile(filename2);
    cout<<"Opened input file:\ndata = "<<filename2<<"\nMC = "<<filename1<<endl; 
    // Creation of the output file
    TString foutName = "PileUp_ReweightingStudy_";
    foutName += datasetName; foutName += ".root";
    TFile *fout = new TFile(foutName, "RECREATE");
    fout->cd();
   
    Int_t NBINS = 100; 
    
    //TH1 *hPileUp_MC, *hPileUp_data_temp;
    TH1 *hPileUp_MC, *hPileUp_data;
    double weight = 0, scale1 = 0, scale2 = 0;
    // MC pile-up histo
    if(option=="PUint") hPileUp_MC = (TH1*)f1->Get("BeforeCuts/hNPileUp");
    if(option=="nPV")   hPileUp_MC = (TH1*)f1->Get("BeforeCuts/hNPrimaryVertices");
    if(hPileUp_MC->GetEntries() != 0)
        scale1 = 1/(hPileUp_MC->Integral());
    else scale1 = 1;
    cout << "Scale factor MC : " << scale1 << endl;
    hPileUp_MC->Scale(scale1);

    // Data pile-up histo
    if(option=="PUint") hPileUp_data = (TH1*)f2->Get("pileup");
    //TH1 *hPileUp_data = new TH1D("Pileup distributions", "Pileup distributions", NBINS, -0.5, NBINS-0.5);
    //for(int m=1; m<hPileUp_MC->GetXaxis()->GetXmax(); m++){
    //    hPileUp_data->SetBinContent(m+1, hPileUp_data_temp->GetBinContent(m));
    //}

    if(option=="nPV")   hPileUp_data = (TH1*)f2->Get("BeforeCuts/hNPrimaryVertices");
        
    if(hPileUp_data->GetEntries() != 0)
        scale2 = 1/(hPileUp_data->Integral());
    else scale2 = 1;
    cout << "Scale factor data : " << scale2 << endl;
    hPileUp_data->Scale(scale2);
    
    int nBins1 = hPileUp_MC->GetNbinsX();
    cout << "MC pile-up histo has " << nBins1 << " bins." << endl;
    int nBins2 = hPileUp_data->GetNbinsX();
    cout << "Data pile-up histo has " << nBins2 << " bins." << endl;
    double Xmin = hPileUp_MC->GetXaxis()->GetXmin();
    cout << "Xmin : " << Xmin << endl;
    double Xmax = hPileUp_MC->GetXaxis()->GetXmax();
    cout << "Xmax : " << Xmax << endl;
   
    int nbin = std::min(nBins1, nBins2); 
    std::vector<double> val1, val2;
    for(int k=0; k<nbin; k++){
        val1.push_back(hPileUp_MC->GetBinContent(k));
        val2.push_back(hPileUp_data->GetBinContent(k));
    }
    
    // Canvas with PileUp distributions normalized to 1
    TString canvasName = "PileUp_distr_"; canvasName += datasetName;
    TCanvas *canv = new TCanvas(canvasName, canvasName, 0, 0, 1200, 1000);
    hPileUp_data->SetStats(0);
    hPileUp_data->SetMarkerStyle(22);
    hPileUp_data->SetMarkerColor(1);
    hPileUp_data->GetXaxis()->SetTitle("nPileUp");
    hPileUp_data->GetYaxis()->SetTitle("Entries");
    hPileUp_data->Draw();

    hPileUp_MC->SetLineColor(kBlue+1);
    hPileUp_MC->SetLineWidth(2);
    hPileUp_MC->SetFillColorAlpha(kBlue, 0.35);
    hPileUp_MC->Draw("bar same");
    
    TLegend *leg = new TLegend(0.80,0.80,0.99,0.94);
    TString legLabel[2];
    legLabel[0] = "Pileup_MC_"; legLabel[0] += datasetName;
    legLabel[1] = "Pileup_data_"; legLabel[1] += dataSample;
    leg->AddEntry(hPileUp_MC, legLabel[0], "f");
    leg->AddEntry(hPileUp_data, legLabel[1], "p");
    leg->Draw();
    
    canv->Write();
    canv->Close();
    
    //Pile-up Reweighting
    TH1D *hweight = new TH1D("PileUp_Reweighting", "PileUp_Reweighting", NBINS, -0.5, NBINS-0.5);
    for(int i=0; i<NBINS; i++){
        cout << " BIN n. " << i << endl;
        cout << " val_MC : " << val1[i] << endl;
        cout << " val_data : " << val2[i] << endl;
        if(val1[i] == 0 || val2[i] < 1e-7) weight = 0;
        else weight = val2[i]/val1[i];
        cout << "weight : " << weight << endl << endl;
        hweight->SetBinContent(i, weight);
    }
    
     // Write and close the file
    fout->Write();
    fout->Close();
    f1->Close();
    f2->Close();
    
}

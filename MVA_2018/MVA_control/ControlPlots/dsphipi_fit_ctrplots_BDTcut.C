#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include "../Control_common.h"

using namespace RooFit;

void dsphipi_fit_ctrplots_BDTcut() 
{
    bool doPUrew = false;

    //open root files where to store all plots!!
    TFile *fout = new TFile("plots_BDTcut/dsphipi_all_output_"+TMVA_inputpath_control+"_BDT_cut.root", "RECREATE");
    fout->cd();

    //set here list of bdt cuts
    TString bdt_cutlist[] = {"", "bdt>-0.02","bdt>0", "bdt>0.02", "bdt>0.03", "bdt>0.04", "bdt>0.06", "bdt>0.07", "bdt>0.075", "bdt>0.08", "bdt>0.10"};
    int ncut = sizeof(bdt_cutlist)/sizeof(bdt_cutlist[0]);
    Double_t Dsyield_data[ncut];
    Double_t Dsyield_data_err[ncut];
    Double_t Dsyield_MC[ncut];

    //TFile *fin = new TFile("../MVA_control_2018_30july_dphi3DoutputTree.root", "READ");
    TFile *fin = new TFile("../"+TMVA_inputpath_control+"outputTree.root", "READ");
    cout<<"opened input file ../"+TMVA_inputpath_control+"outputTree.root"<<endl;
    TTree *tin = (TTree*)fin->Get("outputTree");

    int nrun = sizeof(inputpath_datarun_control)/sizeof(inputpath_datarun_control[0]);

    TChain *tdata_all = new TChain("FinalTree_Control");

    for(auto i=0; i<nrun; i++){
        TFile *f = new TFile(inputpath_datarun_control[i],"READ");
        tdata_all->Add(inputpath_datarun_control[i]);
        f->Close();
    }
    Double_t lumi = 0;
    for(int k=0; k<nrun; k++) lumi = lumi+ Lumi_data_control[k];
    TString run_lable = "2018";

    TString common_cut = " && bs_sv_d2Dsig>2.0 && Ptmu3 > 1.2 && " 
                         "((Ptmu1>3.5 && Etamu1<1.2) || (Ptmu1>2.0 && Etamu1>=1.2 && Etamu1<=2.4)) && "
                         "((Ptmu2>3.5 && Etamu2<1.2) || (Ptmu2>2.0 && Etamu2>=1.2 && Etamu2<=2.4)) && "
                         "abs(phiMass-1.02)<0.045 && "
                         "!(l1double_DoubleMu4_fired && !l1double_DoubleMu0_fired)";
                         //"(l1double_DoubleMu4_fired || l1double_DoubleMu0_fired)";

    TString invmass_peak_MC = "";
    TString invmass_peak_data = "";
    TString invmass_all_data = "";

    if(doPUrew){
    invmass_all_data  = "puFactor*(tripletMass<2.02 && tripletMass>1.62"+common_cut+" && isMC==0";
    invmass_peak_MC = "puFactor*(tripletMass<2.01 && tripletMass>1.93"+common_cut+" && isMC==1";
    invmass_peak_data = "puFactor*(tripletMass<2.01 && tripletMass>1.93"+common_cut+" && isMC==0";
    } else {
    invmass_all_data  = "(tripletMass<2.02 && tripletMass>1.62"+common_cut+" && isMC==0";
    invmass_peak_MC = "(tripletMass<2.01 && tripletMass>1.93"+common_cut+" && isMC==1";
    invmass_peak_data = "(tripletMass<2.01 && tripletMass>1.93"+common_cut+" && isMC==0";
    }
    TString binning_mass = "(58, 1.72, 2.01)";

    Double_t events_SB = 0; //will save value for BDT shape scaling only for ncut==0
    for(int i = 0; i<ncut; i++){

        TH1F *h_tripletmass_mc;
        TH1F *h_tripletmass;
        TH1F *h_tripletmass_sign;

        TString bdt_cut = bdt_cutlist[i];
        TString bdt_cut_label = bdt_cutlist[i];
        if(bdt_cut!="") bdt_cut = +"&& "+bdt_cut+")";
        else            bdt_cut = ")";

        bdt_cut_label = bdt_cut_label.ReplaceAll(".", "p");
        bdt_cut_label = bdt_cut_label.ReplaceAll(">", "_");
        bdt_cut_label = bdt_cut_label.ReplaceAll("<", "_");

        cout<<"tripletMass>>h_tripletmass"+binning_mass+", "+invmass_all_data+bdt_cut<<endl;
        tin->Draw("tripletMass>>h_tripletmass"+binning_mass, invmass_all_data+bdt_cut);
        tin->Draw("tripletMass>>h_tripletmass_sign"+binning_mass, invmass_peak_data+bdt_cut);
        tin->Draw("tripletMass>>h_tripletmass_mc"+binning_mass, invmass_peak_MC+bdt_cut);

        h_tripletmass     = (TH1F *)gDirectory->Get("h_tripletmass");
        h_tripletmass_sign = (TH1F *)gDirectory->Get("h_tripletmass_sign");
        h_tripletmass_mc = (TH1F *)gDirectory->Get("h_tripletmass_mc");

        cout<<"Events passing selections "<<bdt_cutlist[i]<<" = "<<h_tripletmass->GetEntries()<<endl;
        //drawing triplet mass in MC for region mass selection and integral
        TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
        h_tripletmass_mc->Draw();
        auto f1  = new TF1("f1","gaus",1.93,2.01);
        h_tripletmass_mc->Fit("f1", "R");
        f1->Draw("same");
        
        Double_t n_mc_peak = f1->Integral(1.93, 2.01) / h_tripletmass_mc->Integral(h_tripletmass_mc->FindFixBin(1.93),h_tripletmass_mc->FindFixBin(2.01),"width") * h_tripletmass_mc->Integral(h_tripletmass_mc->FindFixBin(1.93),h_tripletmass_mc->FindFixBin(2.01));
        cout<<"n_mc_peak "<<n_mc_peak<<endl;
        c3->Update();
        fout->WriteObject(c3,"2mu1trk_invmass_mc");

        TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
        h_tripletmass->Draw();
        c1->Update();
        h_tripletmass_sign->Draw("same");
        c1->Update();

        // Declare observable x
        TCanvas *c5 = new TCanvas("c5","c5",150,10,800,800);
        c5->SetLeftMargin(0.15);
        c5->Update();
        RooRealVar x("x","2mu+1trk inv. mass (GeV)",1.65,2.01);

        // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
        RooDataHist dh("dh","dh",x,Import(*h_tripletmass)); 

        // Make plot of binned dataset showing Poisson error bars (RooFit default)
        RooPlot* frame = x.frame(Title(" "));
        dh.plotOn(frame);
        //dh.statOn(frame,Layout(0.55,0.99,0.8)) ;

        //set ranges
        x.setRange("R1",1.83,1.89); //first peak D+(1.87GeV)
        x.setRange("R2",1.93,2.01); //second peak Ds(1.97)
        x.setRange("R3",1.72,1.84); //background    
        x.setRange("R4",1.90,1.925); //background    
        x.setRange("R5",1.99,2.01); //background    
        x.setRange("R6",1.72,2.01); //full range    

        RooRealVar mass2("mass2","Central value of Gaussian",1.965,1.94,2.0);
        RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.06);
        Double_t mass1_low = 1.83;
        if(i==0) mass1_low = 1.825; 
        else  mass1_low = 1.83; 
        RooRealVar mass1("mass1","Central value of Gaussian",1.87,mass1_low,1.89);
        RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0003,0.1);


        // Fit a Crystal ball p.d.f to the data
        RooRealVar alpha2("alpha2","alpha value CB",1,-30,30);
        RooRealVar n2("n2", "n2", 2, 0, 20);
        RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
        signal_CB2.fitTo(dh, Range("R2"));
        
        
        // Fit a Crystal ball p.d.f to the data
        RooRealVar alpha1("alpha1","alpha value CB",1,-30,30);
        RooRealVar n1("n1", "n1", 2, 0, 10);
        RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
        signal_CB1.fitTo(dh, Range("R1"));
       

        // Fit an exponential to the background
        Double_t a_high = 0.;
        //if(i==0) a_high = 1.; 
        RooRealVar a("a", "a", -5, -20, a_high);
        RooExponential bg_exp("bg_exp", "bg_exp", x, a);
        bg_exp.fitTo(dh, Range("R3,R4,R5"));
        //bg_exp.fitTo(dh, Range("R3,R4"));
        //bg_exp.fitTo(dh, Range("R3"));


        // Combine the models
        RooRealVar nsig2("nsig2","#signal2 events",80000,500,5000000);
        RooRealVar nsig1("nsig1","#signal1 events",40000,500,5000000);
        RooRealVar nbkg("nbkg","#background events",400000,500,10000000);
        RooAddPdf model("model","g+a",RooArgList(signal_CB1,signal_CB2,bg_exp),RooArgList(nsig1,nsig2,nbkg));
        RooFitResult * r = model.fitTo(dh, Save(true));
        r->Print();
        //model.paramOn(frame,Layout(0.12, 0.4, 0.9));

        //plot
        model.plotOn(frame);
        model.plotOn(frame, Components(bg_exp), LineColor(kGreen), LineStyle(kDashed));
        model.plotOn(frame, Components(RooArgSet(signal_CB2, signal_CB1)), LineColor(kRed), LineStyle(kDashed) );
        frame->Draw();
        //Chi2
        cout<<"frame->chiSquare(3) "<<frame->chiSquare(3)<<endl;
        RooChi2Var chi2("chi2","chi2",model,dh);
        cout << "chi2.getVal() " << chi2.getVal() << endl ;

        //Compute integrals
        x.setRange("signal",1.93,2.01);
        x.setRange("sideband",1.7,1.8);

        //fraction of total events in 1.93,2.01 (n_signal_region_events/n_total_events)
        RooAbsReal* fsigregion_model = model.createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fs = fsigregion_model->getVal();
        Double_t fs_err = fsigregion_model->getPropagatedError(*r);
        //fraction of total events in 1.70,1.80 (n_sideband_region_events/n_total_events)
        RooAbsReal* fsidebandregion_model = model.createIntegral(x,NormSet(x),Range("sideband")); 

        //fraction of background events in 1.93,2.01 
        RooAbsReal* fsigregion_bkg = bg_exp.createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fb = fsigregion_bkg->getVal();
        Double_t fb_err = fsigregion_bkg->getPropagatedError(*r);
        //fraction of background events in 1.70, 1.80 
        RooAbsReal* fsidebandregion_bkg = bg_exp.createIntegral(x,NormSet(x),Range("sideband")); 


        Double_t nsigevents = fs * (nsig2.getVal()+nsig1.getVal()+nbkg.getVal()) - fb*nbkg.getVal(); 
        Double_t nsig_err = pow( pow(fs_err,2) * pow(nsig2.getVal()+nsig1.getVal()+nbkg.getVal(),2)  + ( pow(nsig2.getPropagatedError(*r),2)+pow(nsig1.getPropagatedError(*r),2)+pow(nbkg.getPropagatedError(*r),2)) * pow(fs,2) + pow(fb_err,2) * pow(nbkg.getVal(),2) + pow(nbkg.getPropagatedError(*r),2)*pow(fb,2) , 0.5);

        Double_t fsig = nsigevents/(fsigregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal()));

        //cout<<"fsigregion_model "<<  fsigregion_model->getVal()<<endl;
        //cout<<"fsigregion_bkg "<<  fsigregion_bkg->getVal()<<endl;

        //cout<<"fsidebandregion_model "<<  fsidebandregion_model->getVal()<<endl;
        //cout<<"fsidebandregion_bkg "<<  fsidebandregion_bkg->getVal()<<endl;

        cout<<"n background events in 1.70,1.80 "<< fsidebandregion_bkg->getVal()*nbkg.getVal() <<endl;
        cout<<"n total events in 1.70,1.80 "<< fsidebandregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal())<<endl;

        cout<<"+++++++++++++++++++++++++"<<endl;
        cout<<"Ds yield "+run_lable<<" lumi "<<lumi<<endl;
        cout<<"n signal events in 1.93,2.01 "<< nsigevents<<endl;
        cout<<"error on signal events in 1.93,2.01 "<< nsig_err<<endl;
        cout<<"fraction of signal events in 1.93,2.01 "<< fsig <<endl;
        cout<<"n background events in 1.93,2.01 "<< fsigregion_bkg->getVal()*nbkg.getVal() <<endl;
        Double_t ntotalevents = fsigregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal());
        cout<<"n total events in 1.93,2.01 "<< ntotalevents <<endl;

        if(i==0) events_SB = fsigregion_bkg->getVal()*nbkg.getVal();
        Dsyield_data[i] = nsigevents;
        Dsyield_data_err[i] = nsig_err;
        Dsyield_MC[i] = h_tripletmass_mc->GetEntries();

        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << lumi;
        std::string strLumi = stream.str();
        TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");
        TLatex* bdtcut = new TLatex(0.20,0.71, bdt_cutlist[i]);
        bdtcut->SetNDC(kTRUE);
        bdtcut->Draw("same");
        TLatex* text = new TLatex(0.20,0.91, "\\text{data }"+run_lable+"\n\\text{ }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
        text->SetNDC(kTRUE);
        text->Draw("same");
        c5->Update();
        c5->SaveAs("plots_BDTcut/DsPhiPi_invmass_"+TMVA_inputpath_control+"_"+run_lable+"_BDTcut_"+bdt_cut_label+".png");

        if(i==0) c5->Print(TMVA_inputpath_control+"_"+run_lable+"_BDTcut.gif");
        c5->Print(TMVA_inputpath_control+"_"+run_lable+"_BDTcut.gif+50");
        //if end of loop, save last gif
        if(i==(ncut-1)) c5->Print(TMVA_inputpath_control+"_"+run_lable+"_BDTcut.gif++");

        TString dir_lable = "nocut";
        if(i>0) dir_lable = bdt_cutlist[i];
        TDirectory *dirRun = fout->mkdir(dir_lable);
        dirRun->cd();
        dirRun->WriteObject(c5,"2mu1trk_invmass_data");

        fout->cd();
    }

    fout->cd();
    cout<<"Dsyield_data | Dsyield_data_err | Dsyield_MC | bdt_cut "<<endl;
    for(int i = 0; i<ncut; i++){
        cout<<Dsyield_data[i]<<" | "<<Dsyield_data_err[i]<<" | "<<Dsyield_MC[i]<<" | "<<bdt_cutlist[i]<<endl;
    }

    TCanvas *c4 = new TCanvas("c4","c4",150,10,990,660);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1F *hMC_num = new TH1F("hMC_num", "hMC_num", ncut, -0.5,(ncut-0.5));
    TH1F *hMC_den = new TH1F("hMC_den", "hMC_den", ncut, -0.5,(ncut-0.5));
    for(int i = 0; i<ncut; i++){
        hMC_num->SetBinContent(i+1, Dsyield_MC[i]);
        hMC_num->SetBinError(i+1, sqrt(Dsyield_MC[i]));
        hMC_den->SetBinContent(i+1, Dsyield_MC[0]);
        hMC_den->SetBinError(i+1, sqrt(Dsyield_MC[0]));
        hMC_num->GetXaxis()->SetBinLabel(i+1, bdt_cutlist[i]);
    }
    // Define the ratio plot
    TH1F *hMC_eff = (TH1F*)hMC_num->Clone("hMC_eff");
    hMC_eff->Sumw2();
    hMC_eff->Divide(hMC_den);
    hMC_eff->GetYaxis()->SetTitle("BDT cut efficiency");
    hMC_eff->GetYaxis()->SetRangeUser(0, 1);
    hMC_eff->Draw("ep text0");

    TH1F *hdata_num = new TH1F("hdata_num", "hdata_num", ncut, -0.5,(ncut-0.5));
    TH1F *hdata_den = new TH1F("hdata_den", "hdata_den", ncut, -0.5,(ncut-0.5));
    for(int i = 0; i<ncut; i++){
        hdata_num->SetBinContent(i+1, Dsyield_data[i]);
        hdata_num->SetBinError(i+1, Dsyield_data_err[i]);
        hdata_den->SetBinContent(i+1, Dsyield_data[0]);
        hdata_den->SetBinError(i+1, Dsyield_data_err[0]);
    }
    // Define the ratio plot
    TH1F *hdata_eff = (TH1F*)hdata_num->Clone("hdata_eff");
    hdata_eff->Sumw2();
    hdata_eff->Divide(hdata_den);
    hdata_eff->SetLineColor(kRed+1);
    hdata_eff->SetLineWidth(2);
    hdata_eff->Draw("ep text0 same");
    c4->Update();

    TLegend* leg = new TLegend(0.6, 0.7, .4, .80);
    leg->AddEntry(hMC_eff, "BDT cut eff. - 2018 MC", "lep");
    leg->AddEntry(hdata_eff, "BDT cut eff. - 2018 data", "lep");
    leg->Draw("same");

    c4->SaveAs("plots_BDTcut/BDT_cut_"+TMVA_inputpath_control+".png");
    fout->cd();
    fout->WriteObject(c4,"BDT_cut_"+TMVA_inputpath_control);

    //Drawing BDT score 
    TH1F *h_signal;
    TH1F *h_data;
    TH1F *h_data_bkg;

    TString invmass_SB_data = "(tripletMass<1.80 && tripletMass>1.70"+common_cut+" && isMC==0";
 
    //bdt score distribution
    TString binning = "(240, -0.6, 0.6)"; 
    tin->Draw("bdt>>h_data"+binning, invmass_peak_data+")");
    h_data = (TH1F *)gDirectory->Get("h_data");
    tin->Draw("bdt>>h_data_bkg"+binning, invmass_SB_data+")");
    h_data_bkg = (TH1F *)gDirectory->Get("h_data_bkg");
    tin->Draw("bdt>>h_signal"+binning, invmass_peak_MC+")");
    h_signal = (TH1F *)gDirectory->Get("h_signal");

    h_signal->SetDirectory(0);
    h_data->SetDirectory(0);
    h_data_bkg->SetDirectory(0);

    //scaling the SB distribution to number of background events in 1.93,2.01
    Double_t normSB = h_data_bkg->GetEntries();
    h_data_bkg->Scale(events_SB/normSB);
    
    h_data->Add(h_data_bkg,-1); //subtract h2 from h1 : h1->Add(h2,-1)

    TCanvas *c2 = new TCanvas("c2","c2",150,10,800,800);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //h_signal->GetXaxis()->SetRangeUser(X_min,X_max);
    //h_data->GetXaxis()->SetRangeUser(X_min,X_max);
    h_data->SetLineColor(kBlack);
    h_signal->SetLineColor(kRed);
    h_signal->SetLineWidth(2);
    h_data->SetLineWidth(2);
    h_data->GetXaxis()->SetTitle("BDT score");
    h_signal->Scale(1/h_signal->Integral());
    h_data->Scale(1/h_data->Integral());
    h_signal->Rebin(4);
    h_data->Rebin(4);

    double Y_max;
    Y_max = 0.12;

    h_data->Draw("HISTE");
    h_data->GetYaxis()->SetRangeUser(0, Y_max);
    //h_data->GetXaxis()->SetRangeUser(X_min,X_max);
    h_signal->Draw("same HISTE");
    c2->Update();

    TLegend*leg2 = new TLegend(0.12,0.75,0.47,0.9);
    leg2->AddEntry(h_signal,"2018 control - signal","f");
    leg2->AddEntry(h_data,   "2018 control - data (SB subtracted)","f");
    leg2->Draw();
    c2->Update();
    c2->SaveAs("plots_BDTcut/BDT_shape_"+TMVA_inputpath_control+".png");
        
    fout->cd();
    fout->WriteObject(c2,"BDT_shape_"+TMVA_inputpath_control);

    fout->Write();
    fout->Close();
 
}

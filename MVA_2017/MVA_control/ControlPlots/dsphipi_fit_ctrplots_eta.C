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

void dsphipi_fit_ctrplots_eta() 
{
    bool doPUrew = false;

    gROOT->SetBatch(1);
    //open root files where to store all plots!!
    TFile *fout = new TFile("plots_eta/dsphipi_all_output_"+TMVA_inputpath_control+"_eta.root", "RECREATE");
    fout->cd();

    TString mfeta = "std::fmax( std::fmax(Etamu1, Etamu2), Etamu3)";
    //set here list of most forward mu eta cuts
    TString eta_cutlist[] = {mfeta+">=0 && "+mfeta+"<0.6", mfeta+">=0.6 && "+mfeta+"<1.2", mfeta+">=1.2 && "+mfeta+"<2.4"};
    int ncut = sizeof(eta_cutlist)/sizeof(eta_cutlist[0]);
    Double_t Dsyield_data[ncut];
    Double_t Dsyield_data_err[ncut];
    Double_t Dsyield_MC[ncut];

    Double_t Dsmass_data[ncut];
    Double_t Dsreso_data[ncut];
    Double_t DsresoMeV_data[ncut];
    Double_t Dsmass_MC[ncut];
    Double_t Dsreso_MC[ncut];
    Double_t DsresoMeV_MC[ncut];

    TFile *fin = new TFile("../"+TMVA_inputpath_control+"outputTree.root", "READ");
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
    TString run_lable = "2017";

//    TString common_cut = 
//                         //" && bs_sv_d2Dsig>2.00 && Ptmu3 > 1.2 && "
//                         " && bs_sv_d2Dsig>3.75 && Ptmu3 > 1.2 && "
//                         "((Ptmu1>3.5 && Etamu1<1.2) || (Ptmu1>2.0 && Etamu1>=1.2 && Etamu1<=2.4)) && "
//                         "((Ptmu2>3.5 && Etamu2<1.2) || (Ptmu2>2.0 && Etamu2>=1.2 && Etamu2<=2.4)) && "
//                         //"(abs(phiMass-1.02)<0.06) &&"
//                         "l1double_DoubleMu0_fired";
    TString common_cut = " && bs_sv_d2Dsig>3.75 && Ptmu3 > 2 && l1double_DoubleMu0_fired";

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
    TString binning_mass = "(72, 1.65, 2.01)";

    for(int i = 0; i<ncut; i++){

        TH1F *h_tripletmass_mc;
        TH1F *h_tripletmass;
        TH1F *h_tripletmass_sign;

        TString eta_cut = eta_cutlist[i];
        if(eta_cut!="") eta_cut = +"&& "+eta_cut+")";
        else            eta_cut = ")";

        cout<<"tripletMass>>h_tripletmass"+binning_mass+", "+invmass_all_data+eta_cut<<endl;
        tin->Draw("tripletMass>>h_tripletmass"+binning_mass, invmass_all_data+eta_cut);
        tin->Draw("tripletMass>>h_tripletmass_sign"+binning_mass, invmass_peak_data+eta_cut);
        tin->Draw("tripletMass>>h_tripletmass_mc"+binning_mass, invmass_peak_MC+eta_cut);

        h_tripletmass     = (TH1F *)gDirectory->Get("h_tripletmass");
        h_tripletmass_sign = (TH1F *)gDirectory->Get("h_tripletmass_sign");
        h_tripletmass_mc = (TH1F *)gDirectory->Get("h_tripletmass_mc");

        cout<<"Events passing selections "<<eta_cutlist[i]<<" = "<<h_tripletmass->GetEntries()<<endl;
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
        TCanvas *c5 = new TCanvas("c5","c5",150,10,990,660);
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
        x.setRange("R3",1.65,1.84); //background    
        x.setRange("R4",1.90,1.925); //background    
        x.setRange("R5",1.99,2.01); //background    
        x.setRange("R6",1.65,2.01); //full range    

        RooRealVar mass2("mass2","Central value of Gaussian",1.965,1.94,2.0);
        RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.06);
        Double_t mass1_low = 1.83;
        if(i==0 || i==1) mass1_low = 1.825; 
        else  mass1_low = 1.83; 
        RooRealVar mass1("mass1","Central value of Gaussian",1.87,mass1_low,1.89);
        RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0003,0.3);


        // Fit a Crystal ball p.d.f to the data
        RooRealVar alpha2("alpha2","alpha value CB",1,-20,20);
        RooRealVar n2("n2", "n2", 2, 0, 20);
        RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
        signal_CB2.fitTo(dh, Range("R2"));
        
        
        // Fit a Crystal ball p.d.f to the data
        RooRealVar alpha1("alpha1","alpha value CB",1,-20,20);
        RooRealVar n1("n1", "n1", 2, 0, 5);
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
        RooRealVar nsig2("nsig2","#signal2 events",80000,50,5000000);
        RooRealVar nsig1("nsig1","#signal1 events",40000,50,5000000);
        RooRealVar nbkg("nbkg","#background events",400000,50,10000000);
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

        Dsyield_data[i] = nsigevents;
        Dsyield_data_err[i] = nsig_err;
        Dsyield_MC[i] = h_tripletmass_mc->GetEntries();

        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << lumi;
        std::string strLumi = stream.str();
        TLatex* cmslabel = new TLatex(0.15,0.81, "#bf{CMS Preliminary}");
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");
        //TLatex* text = new TLatex(0.5,0.91, "\n\\text{data }"+run_lable+"\n\\text{ }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
        //text->SetNDC(kTRUE);
        //text->Draw("same");
        c5->Update();
        c5->SaveAs("plots_eta/DsPhiPi_invmass_"+TMVA_inputpath_control+"_"+run_lable+"_eta_"+eta_cutlist[i]+".png");

        TString dir_lable = eta_cutlist[i];
        TDirectory *dirRun = fout->mkdir(dir_lable);
        dirRun->cd();
        dirRun->WriteObject(c5,"2mu1trk_invmass_data");

        fout->cd();
        cout<<"gaus "<<f1->GetParameter(1)<<" | "<<f1->GetParameter(2)<<endl;
        cout<<"CB  "<<mass2.getVal()<<" | "<<sigma2.getVal()<<endl;
        Dsmass_data[i] = mass2.getVal();
        Dsreso_data[i] = sigma2.getVal();
        DsresoMeV_data[i] = sigma2.getVal()*1000;
        Dsmass_MC[i] = f1->GetParameter(1);
        Dsreso_MC[i] = f1->GetParameter(2);
        DsresoMeV_MC[i] = f1->GetParameter(2)*1000;

    }

    fout->cd();
    cout<<"Dsyield_data | Dsyield_data_err | Dsyield_MC | eta_cut "<<endl;
    for(int i = 0; i<ncut; i++){
        cout<<Dsyield_data[i]<<" | "<<Dsyield_data_err[i]<<" | "<<Dsyield_MC[i]<<" | "<<eta_cutlist[i]<<endl;
    }
    gROOT->SetBatch(0);
    TCanvas *c4 = new TCanvas("c4","c4",150,10,800,800);
    //auto g_mass = new TGraphErrors(ncut, Dsmass_data[i], Dsmass_MC[i], Dsreso_data[i], Dsreso_MC[i]);
    gStyle->SetOptFit();
    gStyle->SetOptTitle(0);

    auto g_reso = new TGraph(ncut, DsresoMeV_data, DsresoMeV_MC);
    g_reso->GetXaxis()->SetTitle("Mass Resolution (MeV), 2017 data");
    g_reso->GetYaxis()->SetTitle("Mass Resolution (MeV), MC");
    g_reso->SetMarkerStyle(21);
    g_reso->SetMarkerSize(1.5);
    g_reso->Draw("AP");
    TF1 *flin = new TF1("flin", "[0]*x + [1]", gPad->GetUxmin(), gPad->GetUxmax());
    g_reso->Fit(flin);
    flin->Draw("same");
    g_reso->Draw("P same");
    c4->Update();

      TLegend *legmc = new TLegend(0.48,0.60,0.86,0.74);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->SetTextFont(42);
      //legmc->SetTextSize(0.05);
      legmc->AddEntry(g_reso,"D_{s} mass resolution from m(#mu#mu#pi)","PE");
      legmc->AddEntry(flin,"linear fit (y=1.05x-0.45)","L");
      legmc->Draw();

      //add lumi TDR style
      TLatex latex;
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(kBlack);
      float l  = c4->GetLeftMargin();
      float t  = c4->GetTopMargin();
      float ri = c4->GetRightMargin();
      float bo = c4->GetBottomMargin();
      float lumiTextSize     = 0.50;;
      float lumiTextOffset   = 0.2; 
      TString lumiText = "2017, 38 fb^{-1} (13 TeV)";
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(lumiTextSize*t);
      latex.DrawLatex(1-ri,1-t+lumiTextOffset*t,lumiText);  
      //add CMS text
      float cmsTextFont   = 61;
      float cmsTextSize   = 0.65;
      float extraTextFont = 52;
      float relPosX    = 0.045;
      float relPosY    = 0.035;
      float posX_ =  l + 0.05*(1-l-ri);
      float posY_ = 0.95 -t - relPosY*(1-t-bo);
      latex.SetTextFont(cmsTextFont);
      latex.SetTextSize(cmsTextSize*t);
      latex.SetTextAlign(12);
      latex.DrawLatex(posX_, posY_, "CMS");
      //add Preliminary  text
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(cmsTextSize*t);
      latex.DrawLatex(posX_, posY_ - 1.2*cmsTextSize*t, "Preliminary");

    c4->SaveAs("plots_eta/MassScale_"+TMVA_inputpath_control+".png");
    c4->SaveAs("plots_eta/MassScale_"+TMVA_inputpath_control+".root");

    Double_t eta_bin[] = {0.3, 0.9, 1.8};
    Double_t eta_err[] = {0.3, 0.3, 0.6};
    TCanvas *c2 = new TCanvas("c2","c2",150,10,800,800);
    auto g_mass_data = new TGraphErrors(ncut, eta_bin, Dsmass_data, eta_err, Dsreso_data);
    auto g_mass_MC   = new TGraphErrors(ncut, eta_bin, Dsmass_MC,   eta_err, Dsreso_MC);
    //auto g_reso = new TGraph(ncut, DsresoMeV_data, DsresoMeV_MC);
    g_mass_data->GetXaxis()->SetTitle("|#eta| of most forward decay product");
    g_mass_data->GetYaxis()->SetTitle("Reconstructed D_{s} mass (GeV)");
    g_mass_data->SetMarkerColor(2);
    g_mass_data->SetLineColor(2);
    g_mass_data->SetMarkerStyle(53);
    g_mass_data-> SetMarkerStyle(20);
    g_mass_data-> SetMarkerSize(1.5);
    g_mass_MC->SetMarkerColor(4);
    g_mass_MC->SetLineColor(4);
    g_mass_MC->SetMarkerStyle(54);
    g_mass_MC->SetMarkerStyle(21);
    g_mass_MC->SetMarkerSize(1.5);
    g_mass_data->Draw("AP");
    g_mass_MC->Draw("P");
    c2->Update();
    TLegend* leg = new TLegend(0.6, 0.7, .4, .80);
    leg->AddEntry(g_mass_data, "2017 data", "lep");
    leg->AddEntry(g_mass_MC,   "MC", "lep");
    leg->Draw();
    c2->SaveAs("plots_eta/ResoScale_"+TMVA_inputpath_control+".png");

    fout->Write();
    fout->Close();
 
}

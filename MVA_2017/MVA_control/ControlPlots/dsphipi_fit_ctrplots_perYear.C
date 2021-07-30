#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include "../Control_common.h"

using namespace RooFit;

void dsphipi_fit_ctrplots_perYear() 
{
    bool doPUrew = true;

    //open root files where to store all plots!!
    TFile *fout = new TFile("dsphipi_all_output_"+TMVA_inputpath_control+"perYear.root", "RECREATE");
    fout->cd();

    TString filename_text = "dsphipi_yield_"+TMVA_inputpath_control+"perYear.txt";
    //open text file where to write Ds yields
    ofstream fout_yield(filename_text);

    TFile *f_mc = new TFile(inputpath_DsPhiPi, "READ");
    TChain *tmc = new TChain("FinalTree_Control");
    tmc->Add(inputpath_DsPhiPi); 
    //MC Ds
    TChain *tmc_mu1 = new TChain("TreeMu1");
    TChain *tmc_mu2 = new TChain("TreeMu2");
    tmc_mu1->Add(inputpath_DsPhiPi); 
    tmc_mu2->Add(inputpath_DsPhiPi); 
    tmc->AddFriend(tmc_mu1);
    tmc->AddFriend(tmc_mu2);
    f_mc->Close();

    //adding +1 to nrun: last will be full 2017 fit
    const int nrun = 1 + sizeof(inputpath_datarun_control)/sizeof(inputpath_datarun_control[0]);

    TChain *tdata_all = new TChain("FinalTree_Control");
    TChain *tdata_all_mu1 = new TChain("TreeMu1");
    TChain *tdata_all_mu2 = new TChain("TreeMu2");

    for(auto i=0; i<nrun-1; i++){
        TFile *f = new TFile(inputpath_datarun_control[i],"READ");
        tdata_all->Add(inputpath_datarun_control[i]);
        tdata_all_mu1->Add(inputpath_datarun_control[i]); 
        tdata_all_mu2->Add(inputpath_datarun_control[i]); 
        f->Close();
    }
    tdata_all->AddFriend(tdata_all_mu1);
    tdata_all->AddFriend(tdata_all_mu2);

    Double_t lumi[nrun] = {0.0};
    TString run_lable[nrun];

    for(int k=0; k<nrun-1; k++){
        lumi[k] = Lumi_data_control[k];
        lumi[nrun-1] += Lumi_data_control[k];
        run_lable[k] = run_name_control[k];
    }
    run_lable[nrun-1] = "2017";
 
    TString common_cut = " && bs_sv_d2Dsig>3.75 && Ptmu3 > 1.2 && " 
                         "((Ptmu1>3.5 && Etamu1<1.2) || (Ptmu1>2.0 && Etamu1>=1.2 && Etamu1<=2.4)) && "
                         "((Ptmu2>3.5 && Etamu2<1.2) || (Ptmu2>2.0 && Etamu2>=1.2 && Etamu2<=2.4)) && "
                         "(abs(phiMass-1.02)<0.045) &&"
                         "l1double_DoubleMu0_fired";

    Double_t Dsyield_data[nrun];
    Double_t Dsyield_data_err[nrun];
    Double_t Dsyield_MC[nrun];
    Double_t Dsyield_MC_err[nrun];
    Double_t Bkgyield_data[nrun];
    Double_t Bkgyield_data_err[nrun];

    TH1F *h_tripletmass[nrun];
    TH1F *h_tripletmass_mc;
    Double_t events_SB = 0;
    TString binning_mass = "(62, 1.72, 2.01)";
    Int_t n_bins = 65;

    TString invmass_SB   = "";
    TString invmass_peak = "";
    TString invmass_all = "";
    if(doPUrew){
        invmass_all  = "puFactor*(tripletMass<2.02 && tripletMass>1.62"+common_cut+")";
        invmass_SB   = "puFactor*(tripletMass<1.80 && tripletMass>1.73"+common_cut+")";
        invmass_peak = "puFactor*(tripletMass<2.01 && tripletMass>1.93"+common_cut+")";
    } else {
        invmass_all  = "(tripletMass<2.02 && tripletMass>1.62"+common_cut+")";
        invmass_SB   = "(tripletMass<1.80 && tripletMass>1.73"+common_cut+")";
        invmass_peak = "(tripletMass<2.01 && tripletMass>1.93"+common_cut+")";
    }

    for(int i=0; i<nrun-1; i++){
        TString s = std::to_string(i);

        TFile *f = new TFile(inputpath_datarun_control[i],"READ");
        TTree *tdata = (TTree*)f->Get("FinalTree_Control");
        TTree *tdata_mu1 = (TTree*)f->Get("TreeMu1");
        TTree *tdata_mu2 = (TTree*)f->Get("TreeMu2");
        tdata->AddFriend(tdata_mu1);
        tdata->AddFriend(tdata_mu2);
        tdata->Draw("tripletMass>>h_tripletmass"+s+binning_mass, invmass_all);
        h_tripletmass[i] = (TH1F *)gDirectory->Get("h_tripletmass"+s);
    }
     
    //Fill histogram for full 2017
    tdata_all->Draw("tripletMass>>h_tripletmass"+std::to_string(nrun-1)+binning_mass, invmass_all);
    TString s = std::to_string(nrun-1);
    h_tripletmass[nrun-1] = (TH1F *)gDirectory->Get("h_tripletmass"+s);

    //Fill histogram for MC
    tmc->Draw("tripletMass>>h_tripletmass_mc"+binning_mass, invmass_peak);
    h_tripletmass_mc = (TH1F *)gDirectory->Get("h_tripletmass_mc");

    std::map<std::string, TH1 *> hmap;
    RooCategory c("c", "c");
    for(int i = 0; i<nrun; i++){
        std::stringstream category;
        TString s = std::to_string(i);
        category << "cat" <<i;
        //define category
        c.defineType(category.str().c_str());
        //map category - TH1
        hmap[category.str()] = h_tripletmass[i];
    }
    RooRealVar x("x","2mu+1trk inv. mass (GeV)",1.72,2.01);
    RooDataHist dh("dh", "dh", x, Index(c), Import(hmap));

    // Construct a simultaneous pdf using category "c" as index
    RooSimultaneous simPdf("simPdf", "simultaneous pdf", c);

    //set ranges
    x.setRange("R1",1.83,1.89); //first peak D+(1.87GeV)
    x.setRange("R2",1.93,2.01); //second peak Ds(1.97)
    x.setRange("R3",1.72,1.82); //background    
    x.setRange("R4",1.90,1.925); //background    
    x.setRange("R5",1.99,2.01); //background    
    x.setRange("R6",1.65,2.01); //full range    

    RooRealVar mass2("mass2","Central value of Gaussian",1.965,1.94,2.0);
    RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.06);
    RooRealVar mass1("mass1","Central value of Gaussian",1.87,1.85,1.89);
    RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0003,0.3);

    // Fit a Crystal ball p.d.f to the data
    RooRealVar alpha2("alpha2","alpha value CB",1.5,-20,20);
    RooRealVar n2("n2", "n2", 2, 0, 20);
    RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
    
    // Fit a Crystal ball p.d.f to the data
    RooRealVar alpha1("alpha1","alpha value CB",1.5,-20,20);
    RooRealVar n1("n1", "n1", 2, 0, 10);
    RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
       
    // Fit an exponential to the background
    RooRealVar a("a", "a", -5, -20, 0.0);
    RooExponential bg_exp("bg_exp", "bg_exp", x, a);
    bg_exp.fitTo(dh, Range("R3"));

    // Combine the models
    // relative contributions of background and signal will depend on category
    RooRealVar *nsig2[nrun];
    RooRealVar *nsig1[nrun];
    RooRealVar *nbkg[nrun];
    RooAddPdf *model[nrun];

    // Associate model with the categories. In our case, model is the same for all categories
    for(int i = 0; i<nrun; i++){
        TString category = "cat"+std::to_string(i);
        Int_t entries = h_tripletmass[i]->GetEntries(); //reference for normalisation variables
        nsig2[i] = new RooRealVar("nsig2_"+category,"#signal2 events",int(entries/4),5,int(entries/2)); 
        nsig1[i] = new RooRealVar("nsig1_"+category,"#signal1 events",int(entries/4),5,int(entries/2)); 
        nbkg[i] = new RooRealVar("nbkg_"+category,"#bkacground events",int(entries/2),5,entries);
        model[i] = new RooAddPdf("model_"+category,"g+a",RooArgList(signal_CB1,signal_CB2,bg_exp),RooArgList(*nsig1[i],*nsig2[i],*nbkg[i]));
        simPdf.addPdf(*model[i], category);
    }
    // P e r f o r m   a   s i m u l t a n e o u s   f i t
    // ---------------------------------------------------
    // Perform simultaneous fit of model to data and model_ctl to data_ctl
    RooFitResult * r = simPdf.fitTo(dh, Save(true));
    cout<<"Fit results"<<endl;
    r->floatParsFinal().Print("s");

    TCanvas *c5 = new TCanvas("c5","c5",150,10,800,800);
    c5->SetLeftMargin(0.15);
    c5->Update();

    x.setRange("signal",1.93,2.01);
    x.setRange("sideband",1.72,1.8);

    //for each category: do plot and store integrals:
    for(int i = 0; i<nrun; i++){
        TString category = "cat"+std::to_string(i);
        // Make plot of binned dataset showing Poisson error bars (RooFit default)
        //RooPlot* frame = x.frame(Title(category));
        RooPlot* frame = x.frame(Title(" "));
        dh.plotOn(frame, Cut("c==c::"+category), DataError(RooAbsData::SumW2), Name("data_"+category));
        simPdf.plotOn(frame, Slice(c, category), ProjWData(c, dh), Name("model_"+category), Range("chi2"));
        simPdf.plotOn(frame, Slice(c, category), Components(bg_exp), LineColor(kGreen), LineStyle(kDashed), ProjWData(c, dh));
    //    simPdf.plotOn(frame, Slice(c, category), Components(signal_CB2, signal_CB1), LineColor(kRed), LineStyle(kDashed), ProjWData(c, dh));
        //simPdf.paramOn(frame,Layout(0.12, 0.4, 0.9));
        frame->Draw();

        //add Lumi and Chi2 to plot
        std::stringstream stream_lumi;
        stream_lumi << std::fixed << std::setprecision(1) << lumi[i];
        std::string strLumi = stream_lumi.str();
        TLatex* text_lumi = new TLatex(0.10,0.91, "\n\\text{data }"+run_lable[i]+"\n\\text{    }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
        text_lumi->SetTextSize(0.05);
        text_lumi->SetNDC(kTRUE);
        text_lumi->Draw("same");

        cout<<"frame->chiSquare() "<<frame->chiSquare("model_"+category, "data_"+category, r->floatParsFinal().getSize())<<endl;
        Double_t Chi2 = frame->chiSquare("model_"+category, "data_"+category, r->floatParsFinal().getSize());
        std::stringstream stream_chi2;
        stream_chi2 << std::fixed << std::setprecision(2) << Chi2;
        std::string strChi2 = stream_chi2.str();
        TString chi2tstring = "\\chi^{2}\\text{/NDOF} = "+strChi2;
        TLatex* text_chi2 = new TLatex(0.20,0.74, chi2tstring);
        text_chi2->SetTextSize(0.04);
        text_chi2->SetNDC(kTRUE);
        text_chi2->Draw("same");

        TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");



        fout->WriteObject(c5,category);
        c5->SaveAs("dsphipi_fit_perYear_"+category+".png");

        //fraction of total events in 1.93,2.01 (n_signal_region_events/n_total_events)
        RooAbsReal* fsigregion_model = model[i]->createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fs = fsigregion_model->getVal();
        Double_t fs_err = fsigregion_model->getPropagatedError(*r);
        //fraction of total events in 1.70,1.80 (n_sideband_region_events/n_total_events)
        RooAbsReal* fsidebandregion_model = model[i]->createIntegral(x,NormSet(x),Range("sideband")); 

        //fraction of background events in 1.93,2.01
        RooAbsReal* fsigregion_bkg = ((RooAbsPdf*)model[i]->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fb = fsigregion_bkg->getVal();
        Double_t fb_err = fsigregion_bkg->getPropagatedError(*r);
        //fraction of background events in 1.70, 1.80 
        RooAbsReal* fsidebandregion_bkg = ((RooAbsPdf*)model[i]->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("sideband")); 

        Dsyield_data[i] = fs * (nsig2[i]->getVal()+nsig1[i]->getVal()+nbkg[i]->getVal()) - fb*nbkg[i]->getVal();;
        Dsyield_data_err[i] = pow( pow(fs_err,2) * pow(nsig2[i]->getVal()+nsig1[i]->getVal()+nbkg[i]->getVal(),2)  + ( pow(nsig2[i]->getPropagatedError(*r),2)+pow(nsig1[i]->getPropagatedError(*r),2)+pow(nbkg[i]->getPropagatedError(*r),2)) * pow(fs,2) + pow(fb_err,2) * pow(nbkg[i]->getVal(),2) + pow(nbkg[i]->getPropagatedError(*r),2)*pow(fb,2) , 0.5);;

        Bkgyield_data[i] = fb*nbkg[i]->getVal();
        Bkgyield_data_err[i] = sqrt( pow(fb_err, 2.0) * pow(nbkg[i]->getVal(), 2.0) + pow(fb, 2.0) * pow(nbkg[i]->getPropagatedError(*r), 2.0));
    }

    //drawing triplet mass in MC for region mass selection and integral
    TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
    h_tripletmass_mc->Draw();
    auto f1  = new TF1("f1","gaus",1.93,2.01);
    h_tripletmass_mc->Fit("f1", "R");
    f1->Draw("same");

    for(int i = 0; i<nrun; i++){
        Double_t n_mc_peak = f1->Integral(1.93, 2.01) / h_tripletmass_mc->Integral(h_tripletmass_mc->FindFixBin(1.93),h_tripletmass_mc->FindFixBin(2.01),"width") * h_tripletmass_mc->Integral(h_tripletmass_mc->FindFixBin(1.93),h_tripletmass_mc->FindFixBin(2.01));
        //cout<<"n_mc_peak "<<n_mc_peak<<endl;
        //cout<<"scaled to lumi: "<<n_mc_peak*lumi[i]*xsection_mc*BR/N_MC<<endl;
        Dsyield_MC[i] = n_mc_peak*lumi[i]*xsection_mc*BR/N_MC;
        Dsyield_MC_err[i] = sqrt(n_mc_peak)*(lumi[i]*xsection_mc*BR/N_MC); 
        //Dsyield_MC_err[i] = (f1->IntegralError(1.93, 2.01))*(lumi[i]*xsection_mc*BR/N_MC); 
        //Dsyield_MC[i] = h_tripletmass_mc->GetEntries();
        c3->Update();
        fout->WriteObject(c3,"2mu1trk_invmass_mc");

        fout_yield<<Dsyield_data[i]<<"\t"<<Dsyield_data_err[i]<<"\n";

        cout<<run_lable[i]<<" lumi="<<lumi[i]<<endl;
        cout<<"=================\noverall data/MC scale factor:"<<endl;
        cout<<"data: "<<Dsyield_data[i]<<" +- "<<Dsyield_data_err[i]<<"\n";
        cout<<"MC: "<<Dsyield_MC[i]<<" +- "<<Dsyield_MC_err[i]<<endl;
        cout<<"scale factor: "<<Dsyield_data[i]/Dsyield_MC[i]<<" +- "<<sqrt(pow((Dsyield_data_err[i]/Dsyield_MC[i]),2.0) + pow((Dsyield_data[i]/(Dsyield_MC[i]*Dsyield_MC[i]))*Dsyield_MC_err[i],2.0))<<endl;
    }

    //Drawing control plots using yields from full 2017
    //List of variables
    TString var[] = {
        //"Ptmu1",
       // "Ptmu2",
       // //"Ptmu3",
       // "Etamu1","Etamu2",
       // //"Etamu3",
       // "muIDflag_mu1", "muIDflag_mu2",
       //  "Pmu3","cLP","tKink","segmComp","fv_nC","d0sig",
        "fv_dphi3D",
        "fv_dphi2D",
        "abs(RefVx1 - SVx)",
        "abs(RefVy1 - SVy)",
        "abs(RefVz1 - SVz)",
       //  //"fv_d3D","fv_d3Dsig",
       //  //"bs_sv_d2Dsig",
       //  //"pv_sv_dxy","pv_sv_dxy_sig",
       //    "pv_sv_dxy/pv_sv_dxy_sig",
       //    "mindca_iso","trkRel"
    };
    TString var_names[] = {
        //"Ptmu1",
       // "Ptmu2",
       // //"Ptmu3",
       // "Etamu1","Etamu2",
       // //"Etamu3",
       // "muIDflag_mu1", "muIDflag_mu2",
       //  //"Pmu3","cLP","tKink","segmComp","fv_nC","d0sig","fv_dphi3D",
        "fv_dphi3D",
        "fv_dphi2D",
        "abs(RefVx1 - SVx)",
        "abs(RefVy1 - SVy)",
        "abs(RefVz1 - SVz)",
       //  //"fv_d3D","fv_d3Dsig",
       //  //"bs_sv_d2Dsig",
       //  //"pv_sv_dxy","pv_sv_dxy_sig",
       //    "pv_sv_dxy/pv_sv_dxy_sig",
       //    "mindca_iso","trkRel"
    };
    const int n = sizeof(var)/sizeof(var[0]);

    TH1F *hdata_bkg[n];
    TH1F *hdata_bkg_plus[n];
    TH1F *hdata_bkg_minus[n];
    TH1F *hdata_sgn[n];
    TH1F *hdata_sgn_plus[n];
    TH1F *hdata_sgn_minus[n];
    TH1F *hmc_sgn[n];
    
    TString binning = "";

    for(int k = 0; k<n; k++){

        TString varname = var[k];
        TString varlable = var_names[k];
        cout<<run_lable[nrun-1]<<" "<<varname<<endl;
        TString s = std::to_string(k);

        if(varname=="Ptmu1" || varname=="Ptmu2" || varname=="Ptmu3") binning = "(60,0,30)";
        if(varname=="Etamu1" || varname=="Etamu2" || varname=="Etamu3") binning = "(30,0,2.5)";
        if(varname=="fv_dphi3D") binning = "(42,-0.01,0.12)";
        if(varname=="fv_dphi2D") binning = "(42,-0.01,0.12)";
        if(varname.Contains("Ref")) binning = "(51,-0.01,1.0)";

        tdata_all->Draw(varname+">>hdata_bkg"+s+binning, invmass_SB);
        tdata_all->Draw(varname+">>hdata_bkg_plus"+s+binning, invmass_SB);
        tdata_all->Draw(varname+">>hdata_bkg_minus"+s+binning, invmass_SB);
        tdata_all->Draw(varname+">>hdata_sgn"+s+binning, invmass_peak);
        tdata_all->Draw(varname+">>hdata_sgn_plus"+s+binning, invmass_peak);
        tdata_all->Draw(varname+">>hdata_sgn_minus"+s+binning, invmass_peak);

        tmc->Draw(varname+">>hmc_sgn"+s+binning, invmass_peak);

        hdata_bkg[k] = (TH1F *)gDirectory->Get("hdata_bkg"+s);
        hdata_bkg_plus[k] = (TH1F *)gDirectory->Get("hdata_bkg_plus"+s);
        hdata_bkg_minus[k] = (TH1F *)gDirectory->Get("hdata_bkg_minus"+s);
        hdata_sgn[k] = (TH1F *)gDirectory->Get("hdata_sgn"+s);
        hdata_sgn_plus[k] = (TH1F *)gDirectory->Get("hdata_sgn_plus"+s);
        hdata_sgn_minus[k] = (TH1F *)gDirectory->Get("hdata_sgn_minus"+s);
        hmc_sgn[k] = (TH1F *)gDirectory->Get("hmc_sgn"+s);

        TCanvas *c2 = new TCanvas("c2","c2",150,10,800,800);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        // Upper plot will be in pad1
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->SetGridx();      // Vertical grid
        pad1->Draw();          // Draw the upper pad: pad1
        pad1->cd();            // pad1 becomes the current pad
        gStyle->SetOptTitle(0);
        hmc_sgn[k]->SetTitle(varname);

        Double_t normMC = hmc_sgn[k]->GetEntries();
        //Normalizing Monte Carlo 
        Double_t wNorm = lumi[nrun-1]*xsection_mc*BR/N_MC;
        //Double_t wNorm = lumi[nrun-1]*xsection_mc*BR/N_MC  *  n_mc_peak/hmc_sgn[k]->GetEntries();
        cout<<"wNorm = lumi[nrun-1]*xsection_mc*BR/N_MC = "<<wNorm<<endl;
        hmc_sgn[k]->Scale(wNorm);

        //scaling the SB distribution to number of background events in 1.93,2.01
        Double_t normSB = hdata_bkg[k]->GetEntries();
        hdata_bkg[k]->Scale(Bkgyield_data[nrun-1]/normSB);
        hdata_bkg_plus[k]->Scale( (Bkgyield_data[nrun-1]/normSB) * 1.10);
        hdata_bkg_minus[k]->Scale( (Bkgyield_data[nrun-1]/normSB) * 0.90);

        cout<<"Entries in  hdata_sgn[k] before SB subtraction "<<hdata_sgn[k]->GetEntries()<<endl;
        hdata_sgn[k]->Add(hdata_bkg[k],-1); //subtract h2 from h1 : h1->Add(h2,-1)
        hdata_sgn_plus[k]->Add(hdata_bkg_plus[k],-1); //subtract h2 from h1 : h1->Add(h2,-1)
        hdata_sgn_minus[k]->Add(hdata_bkg_minus[k],-1); //subtract h2 from h1 : h1->Add(h2,-1)

        //Rescaling to same integral
        hmc_sgn[k]->Scale( 1.0 / hmc_sgn[k]->Integral());
        hdata_sgn[k]->Scale( hmc_sgn[k]->Integral()/hdata_sgn[k]->Integral() );
        hdata_sgn_plus[k]->Scale( hmc_sgn[k]->Integral()/hdata_sgn_plus[k]->Integral() );
        hdata_sgn_minus[k]->Scale( hmc_sgn[k]->Integral()/hdata_sgn_minus[k]->Integral() );

        //plot makeup
        double Y_max = std::max(hmc_sgn[k]->GetMaximum(), hdata_sgn[k]->GetMaximum());
        Y_max = Y_max*1.05;
        hmc_sgn[k]->GetYaxis()->SetRangeUser(0, Y_max);

        hmc_sgn[k]->GetYaxis()->SetTitle("a.u.");
        hmc_sgn[k]->GetYaxis()->SetTitleSize(22);
        hmc_sgn[k]->GetYaxis()->SetTitleFont(43);
        hmc_sgn[k]->GetYaxis()->SetTitleOffset(1.5);

        hmc_sgn[k]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        hmc_sgn[k]->GetXaxis()->SetTitleOffset(1.5);

        hmc_sgn[k]->SetLineColor(kBlue);
        hmc_sgn[k]->SetLineWidth(3);
        hmc_sgn[k]->SetFillStyle(3004);
        hmc_sgn[k]->SetFillColor(kBlue);
        hdata_sgn[k]->SetLineColor(kRed);
        hdata_sgn[k]->SetLineWidth(3);
        hdata_sgn[k]->SetFillStyle(3005);
        hdata_sgn[k]->SetFillColor(kRed);

        hdata_sgn_plus[k]->SetLineColor(kBlack);
        hdata_sgn_plus[k]->SetFillStyle(3001);
        hdata_sgn_plus[k]->SetLineWidth(2);
        hdata_sgn_plus[k]->SetLineStyle(2);

        hdata_sgn_minus[k]->SetLineColor(kBlack);
        hdata_sgn_minus[k]->SetFillStyle(3001);
        hdata_sgn_minus[k]->SetLineWidth(2);
        hdata_sgn_minus[k]->SetLineStyle(3);

        hmc_sgn[k]->Draw("hist");
        hdata_sgn[k]->Draw("hist same");
        hdata_sgn_plus[k]->Draw("hist same");
        hdata_sgn_minus[k]->Draw("hist same"); 

        hmc_sgn[k]->SetStats(0);

        Double_t x_low = 0.1, x_high = 0.45, y_low = 0.65, y_high = 0.90; //top left 
        if(varname=="fv_dphi3D") {x_low = 0.55; x_high = 0.90; y_low = 0.65; y_high = 0.90;} //top right 
        TLegend*leg = new TLegend(x_low, y_low, x_high, y_high);
        leg->AddEntry(hmc_sgn[k],"MC DsPhiPi","f");
        leg->AddEntry(hdata_sgn[k],"data "+run_lable[nrun-1]+" (SB subtracted)","f");
        leg->AddEntry(hdata_sgn_plus[k],"data "+run_lable[nrun-1]+" (SB +10\% subtracted)","f");
        leg->AddEntry(hdata_sgn_minus[k],"data "+run_lable[nrun-1]+" (SB -10\% subtracted)","f");
        leg->Draw();

        //K-S consistency test
        Double_t KS = hdata_sgn[k]->KolmogorovTest(hmc_sgn[k]);
        std::stringstream stream_KS;
        stream_KS << std::fixed << std::setprecision(3) << KS;
        std::string strKS = stream_KS.str();
        TString KStstring = "K-S test = "+strKS;
        TLatex* text_KS = new TLatex(x_low+0.02,y_low-0.05, KStstring);
        text_KS->SetTextSize(0.04);
        text_KS->SetNDC(kTRUE);
        //text_KS->Draw("same");

        // lower plot will be in pad2
        c2->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.2);
        pad2->SetGridx(); // vertical grid
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad    
        // Define the ratio plot
        TH1F *h_x_ratio = (TH1F*)hdata_sgn[k]->Clone("h_x_ratio");
        h_x_ratio->Sumw2();
        h_x_ratio->Divide(hmc_sgn[k]);
        h_x_ratio->SetStats(0);
        // Ratio plot settings
        gStyle->SetLineWidth(2);
        h_x_ratio->SetTitle(""); // Remove the ratio title
        h_x_ratio->GetYaxis()->SetTitle("ratio data/MC");
        h_x_ratio->GetYaxis()->SetRangeUser(-0.5,2);
        h_x_ratio->SetLineColor(kBlack);
        h_x_ratio->GetYaxis()->SetTitleSize(22);
        h_x_ratio->GetYaxis()->SetTitleFont(43);
        h_x_ratio->GetYaxis()->SetTitleOffset(2.0);
        h_x_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_x_ratio->GetYaxis()->SetLabelSize(15);

        // X axis ratio plot settings
        h_x_ratio->GetXaxis()->SetTitle(varlable);
        h_x_ratio->GetXaxis()->SetTitleSize(22);
        h_x_ratio->GetXaxis()->SetTitleFont(43);
        h_x_ratio->GetXaxis()->SetTitleOffset(2.5);
        h_x_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_x_ratio->GetXaxis()->SetLabelSize(15);
        //Compute weighted average ratio
        Double_t mean = 0;
        Double_t std_dev = 0;
        for (int c=1; c<=(h_x_ratio->GetNbinsX()); c++)
        {
            if(h_x_ratio->GetBinContent(c) == 0 || h_x_ratio->GetBinError(c) == 0) continue;
            //cout<<c<<" "<<h_x_ratio->GetBinContent(c)<<" +- "<<h_x_ratio->GetBinError(c)<<endl;
            mean = mean + h_x_ratio->GetBinContent(c) / ( h_x_ratio->GetBinError(c) * h_x_ratio->GetBinError(c) );
            std_dev = std_dev + 1/( h_x_ratio->GetBinError(c) * h_x_ratio->GetBinError(c));
        }
        mean = mean/std_dev;
        std_dev = 1/std_dev;
        //Get mean value and error of ratio plot
        cout<<var[k]+run_lable[nrun-1]+" Mean: "<<mean<<endl;
        cout<<var[k]+run_lable[nrun-1]+" StdDev: "<<std_dev<<endl;
        //Draw line corresponding to mean value on ratio plot
        TLine l;
        h_x_ratio->Draw("ep");
        l.DrawLine(h_x_ratio->GetXaxis()->GetXmin(),mean,h_x_ratio->GetXaxis()->GetXmax(),mean);
        h_x_ratio->Draw("same");

        c2->cd();
        c2->Update();
        varname = varname.ReplaceAll(".", "_");
        varname = varname.ReplaceAll(":", "_");
        varname = varname.ReplaceAll(">", "_");
        varname = varname.ReplaceAll("<", "_");
        varname = varname.ReplaceAll("?", "_");
        varname = varname.ReplaceAll("(", "_");
        varname = varname.ReplaceAll(")", "_");
        h_x_ratio->SetName(varname+"_"+TMVA_inputpath_control+"_"+run_lable[nrun-1]);
        h_x_ratio->Write();
        //c2->SaveAs("plots_MVAvalidation/"+varname+"_"+TMVA_inputpath_control+"_"+run_lable[nrun-1]+".png");
        fout->cd();
        fout->WriteObject(c2,varname+run_lable[nrun-1]);

    }
    fout->Write();
    fout->Close();
    fout_yield.close();
 
}

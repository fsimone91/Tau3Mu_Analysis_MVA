#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAddModel.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "TString.h"
#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace RooFit;

std::string to_string_precision(const float a_value, int n)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

void taumasspeak(TString infilename, TString label, TString categ) 
{
    //open root files where to store all plots!!
    TFile *fout = new TFile("output_"+infilename, "RECREATE");
    fout->cd();

    TFile *f = new TFile(infilename,"READ");
    TTree *t = (TTree*)f->Get("FinalTree");
  
    TString invmass_SB = "puFactor*(tripletMass<1.80 && tripletMass>1.70)";
    TString invmass_peak = "puFactor*(tripletMass<2.01 && tripletMass>1.93)";
    //TString invmass_SB = "puFactor*(tripletMassRef<1.80 && tripletMassRef>1.70)";
    //TString invmass_peak = "puFactor*(tripletMassRef<2.01 && tripletMassRef>1.93)";
    TString binning_mass = "(84, 1.60, 2.02)";

    TCut reso_A = "tripletMassReso < 0.007";
    TCut reso_B = "tripletMassReso >= 0.007 && tripletMassReso <= 0.0105";
    TCut reso_C = "tripletMassReso > 0.0105";
    TCut reso_cat = "tripletMassReso < 0"; //always false    

    if(categ.Contains("A")) reso_cat = reso_cat || reso_A;
    if(categ.Contains("B")) reso_cat = reso_cat || reso_B;
    if(categ.Contains("C")) reso_cat = reso_cat || reso_C;

    if(categ=="") reso_cat = "tripletMassReso > 0"; //always true
    TH1F *h_tripletmass;
    t->Draw("tripletMass>>h_tripletmass"+binning_mass, reso_cat);
    //t->Draw("tripletMassRef>>h_tripletmass"+binning_mass, reso_cat);

    h_tripletmass     = (TH1F *)gDirectory->Get("h_tripletmass");

    // Declare observable x
    TCanvas *c5 = new TCanvas("c5","c5",150,10,990,660);
    c5->Update();
    RooRealVar x("x","3mu inv. mass (GeV)",1.60,2.02);

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,Import(*h_tripletmass)); 

    //import data from TTree instead
    //RooRealVar t3m("tripletMass","tripletMass",1.60,2.02);
    //RooDataSet ds("ds", "ds", RooArgSet(t3m), Import(*t));

    // Make plot of binned dataset showing Poisson error bars (RooFit default)
    RooPlot* frame = x.frame(Title(" "));
    dh.plotOn(frame);
    //dh.statOn(frame,Layout(0.55,0.99,0.8)) ;

    //set ranges
    x.setRange("R",1.63,2.0); //full range    
    x.setRange("R1",1.772,1.779); //tau peak
    x.setRange("R2",1.75,1.79); //tau peak

    RooRealVar mass("mass","Central value of Gaussian",1.777,1.772,1.781);
    RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0,0.1);
    RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0,0.1);

    // Fit background
    RooRealVar a("a", "a", -5, -10, 10.);
    RooExponential exp("exp", "exp", x, a);
    exp.fitTo(dh, Range("R"));
    // Fit second gaussian
    RooGaussian gaus2("gaus2", "The signal distribution", x, mass, sigma2); 
    gaus2.fitTo(dh, Range("R2"));
    // Fit first gaussian
    RooGaussian gaus1("gaus1", "The signal distribution", x, mass, sigma1); 
    gaus1.fitTo(dh, Range("R1"));
    
    // Combine the models
    RooRealVar nsig2("nsig2","#signal2 events",1000,0.,50000);
    RooRealVar nsig1("nsig1","#signal1 events",1000,0.,50000);
    RooRealVar nbkg("nbkg","#bkg events",200,0.,5000);
    RooAddPdf model("model","g+a",RooArgList(exp,gaus1,gaus2),RooArgList(nbkg,nsig1,nsig2));
    RooFitResult * r = model.fitTo(dh, Save(true));
    r->Print();
    model.paramOn(frame,Layout(0.6, 0.9, 0.9));

    //plot
    model.plotOn(frame, LineColor(kRed), LineStyle(kDashed));
    model.plotOn(frame, Components(RooArgSet(gaus1)), LineColor(kBlue), LineStyle(kDashed) );
    model.plotOn(frame, Components(RooArgSet(gaus2)), LineColor(kGreen), LineStyle(kDashed) );

//    //trying with CB
//    RooRealVar alpha("alpha","alpha value CB",1,-20,20);
//    RooRealVar n("n", "n", 2, 0, 20);
//    RooRealVar sigma3("sigma3","Width of Gaussian",0.01,0,0.06);
//    RooCBShape CB("cb", "The signal distribution", x, mass, sigma3, alpha, n);
//    CB.fitTo(dh, Range("R"));
//    CB.plotOn(frame);
//    CB.paramOn(frame,Layout(0.7, 0.9, 0.7));

    frame->Draw();
    //Chi2
    cout<<"frame->chiSquare(3) "<<frame->chiSquare(3)<<endl;
    RooChi2Var chi2("chi2","chi2",model,dh);
    cout << "chi2.getVal() " << chi2.getVal() << endl ;

    Double_t ntot = nsig1.getVal() +  nsig2.getVal() +  nbkg.getVal();
    Double_t ntot_err = sqrt(pow(nsig1.getError(), 2.0) + pow(nsig2.getError(), 2.0) + pow(nbkg.getError(), 2.0));
    Double_t f1 = nsig1.getVal()/ntot;
    Double_t f2 = nsig2.getVal()/ntot;
    Double_t f1_err = 1/ntot * sqrt( pow(nsig1.getError(), 2.0) + pow(f1, 2.0) * pow(ntot_err, 2.0) );
    Double_t f2_err = 1/ntot * sqrt( pow(nsig2.getError(), 2.0) + pow(f2, 2.0) * pow(ntot_err, 2.0) );
    Double_t mass_reso = sqrt( f1 * pow(sigma1.getVal(), 2.0) + f2 * pow(sigma2.getVal(), 2.0) );
    Double_t mass_reso_err = 0.5 * 1/mass_reso * sqrt ( pow(sigma1.getVal(),4.0)*pow(f1_err, 2.0) + pow(f1*2*sigma1.getVal()*sigma1.getError(), 2.0) + pow(sigma2.getVal(),4.0)*pow(f2_err, 2.0) + pow(f2*2*sigma2.getVal()*sigma2.getError(), 2.0) );

    cout<<"mass resolution "<<mass_reso<<" #pm "<<mass_reso_err<<endl;
    TString reso_string = to_string_precision(mass_reso*1000,3);
    TString reso_err_string = to_string_precision(mass_reso_err*1000,3);

    if(categ=="") categ = "all";
    TLatex* text = new TLatex(0.1,0.91,label+" m(3mu) resolution - cat. "+categ+" = "+reso_string+" #pm "+reso_err_string+" MeV");
    text->SetNDC(kTRUE);
    text->Draw("same");

    c5->Update();
    c5->SaveAs(infilename+label+categ+"_plot.png");

    fout->Write();
    fout->Close();

}

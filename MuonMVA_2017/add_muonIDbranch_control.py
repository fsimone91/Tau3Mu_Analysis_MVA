import ROOT
from ROOT import TMVA 
from ROOT import gROOT
from array import array
import numpy as np

def get_muon_SF(_h, MVA, pt ):
    SF = [1.0, 0.0]
    if( MVA < _h.GetXaxis().GetXmax() and MVA > _h.GetXaxis().GetXmin() and pt < _h.GetYaxis().GetXmax() and pt > _h.GetYaxis().GetXmin()):
        iMVA = _h.GetXaxis().FindBin(MVA);
        ipt = _h.GetYaxis().FindBin(pt);
        SF[0] = _h.GetBinContent(iMVA,ipt);
        SF[1] = _h.GetBinError(iMVA,ipt);
    return SF

basedir = "/lustrehome/fsimone/MuonID_study/MuonMVA_2017/"

def add_eval_branch(in_fname='infile.root'):
   
   Muon_cLP = array('f', [-1.0]); 
   Muon_cLM = array('f', [-1.0]);
   Muon_staRelChi2 = array('f', [-1.0]);
   Muon_trkRelChi2 = array('f', [-1.0]); 
   Muon_glbdEP = array('f', [-1.0]);
   Muon_trkKink = array('f', [-1.0]);
   Muon_glbKink = array('f', [-1.0]);
   #Muon_glbTrkP = array('f', [-1.0]);
   Muon_nTVH = array('f', [-1.0]);
   Muon_nVPH = array('f', [-1.0]);
   Muon_nMS = array('f', [-1.0]);
   Muon_segComp = array('f', [-1.0]);
   Muon_glbNChi2 = array('f', [-1.0]);
   Muon_inner_nChi2 = array('f', [-1.0]);
   Muon_outer_nChi2 = array('f', [-1.0]);
   Muon_innner_VF = array('f', [-1.0]);

   mu_eta = array('f', [-1.0])
   mu_pt = array('f', [-1.0])
   mu_phi = array('f', [-1.0])
   mu_simPdgId = array('f', [-1.0])
   mu_simMotherPdgId = array('f', [-1.0])
   mu_SoftMVA = array('f', [-1.0])
   
   # Muon Id 1
   reader_MuonId_barrel = TMVA.Reader("!Color:!Silent");
   reader_MuonId_barrel.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon_cLM);
   reader_MuonId_barrel.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon_cLP);
   reader_MuonId_barrel.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon_staRelChi2); 
   reader_MuonId_barrel.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon_trkRelChi2);
   reader_MuonId_barrel.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon_glbdEP);
   reader_MuonId_barrel.AddVariable("log(mu_combinedQuality_trkKink)", Muon_trkKink);
   reader_MuonId_barrel.AddVariable("log(mu_combinedQuality_glbKink)", Muon_glbKink);
   #reader_MuonId_barrel.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon_glbTrkP);
   reader_MuonId_barrel.AddVariable("mu_Numberofvalidpixelhits", Muon_nVPH);
   reader_MuonId_barrel.AddVariable("mu_trackerLayersWithMeasurement", Muon_nTVH);
   reader_MuonId_barrel.AddVariable("mu_numberOfMatchedStations", Muon_nMS);
   reader_MuonId_barrel.AddVariable("mu_segmentCompatibility", Muon_segComp);
   reader_MuonId_barrel.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon_glbNChi2);
   reader_MuonId_barrel.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon_inner_nChi2);
   reader_MuonId_barrel.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon_outer_nChi2);
   reader_MuonId_barrel.AddVariable("mu_innerTrack_validFraction", Muon_innner_VF);
   reader_MuonId_barrel.AddSpectator("mu_eta",mu_eta);
   reader_MuonId_barrel.AddSpectator("mu_pt",mu_pt);
   reader_MuonId_barrel.AddSpectator("mu_phi",mu_phi);
   reader_MuonId_barrel.AddSpectator("mu_simPdgId",mu_simPdgId);
   reader_MuonId_barrel.AddSpectator("mu_simMotherPdgId",mu_simMotherPdgId);
   reader_MuonId_barrel.AddSpectator("mu_SoftMVA",mu_SoftMVA);
   reader_MuonId_barrel.BookMVA( "BDT", basedir+"MuonMVA_2017_june2021_noGLprob_barrel/weights/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles

   reader_MuonId_endcap = TMVA.Reader("!Color:!Silent");
   reader_MuonId_endcap.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon_cLM);
   reader_MuonId_endcap.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon_cLP);
   reader_MuonId_endcap.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon_staRelChi2); 
   reader_MuonId_endcap.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon_trkRelChi2);
   reader_MuonId_endcap.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon_glbdEP);
   reader_MuonId_endcap.AddVariable("log(mu_combinedQuality_trkKink)", Muon_trkKink);
   reader_MuonId_endcap.AddVariable("log(mu_combinedQuality_glbKink)", Muon_glbKink);
   #reader_MuonId_endcap.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon_glbTrkP);
   reader_MuonId_endcap.AddVariable("mu_Numberofvalidpixelhits", Muon_nVPH);
   reader_MuonId_endcap.AddVariable("mu_trackerLayersWithMeasurement", Muon_nTVH);
   reader_MuonId_endcap.AddVariable("mu_numberOfMatchedStations", Muon_nMS);
   reader_MuonId_endcap.AddVariable("mu_segmentCompatibility", Muon_segComp);
   reader_MuonId_endcap.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon_glbNChi2);
   reader_MuonId_endcap.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon_inner_nChi2);
   reader_MuonId_endcap.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon_outer_nChi2);
   reader_MuonId_endcap.AddVariable("mu_innerTrack_validFraction", Muon_innner_VF);
   reader_MuonId_endcap.AddSpectator("mu_eta",mu_eta);
   reader_MuonId_endcap.AddSpectator("mu_pt",mu_pt);
   reader_MuonId_endcap.AddSpectator("mu_phi",mu_phi);
   reader_MuonId_endcap.AddSpectator("mu_simPdgId",mu_simPdgId);
   reader_MuonId_endcap.AddSpectator("mu_simMotherPdgId",mu_simMotherPdgId);
   reader_MuonId_endcap.AddSpectator("mu_SoftMVA",mu_SoftMVA);
   reader_MuonId_endcap.BookMVA('BDT', basedir+"MuonMVA_2017_june2021_noGLprob_endcap/weights/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles

   file_ = ROOT.TFile(in_fname, "READ")
   old_tree = [ file_.Get('TreeMu1'), file_.Get('TreeMu2') ]

   out_fname = in_fname.split(".root")[0] + '_MuonID.root'
   new_file = ROOT.TFile(out_fname,"RECREATE")
   new_tree = [ ROOT.TTree('MuonIDeval_Mu1', 'MuonIDeval_Mu1'), ROOT.TTree('MuonIDeval_Mu2', 'MuonIDeval_Mu2') ]
   #new_tree = old_tree.CloneTree(0)
 
   globalMuonId = [ array('f', [-99.]), array('f', [-99.]) ]
   scalefactor = [ array('f', [-99.]), array('f', [-99.]) ]
   scalefactor_err = [ array('f', [-99.]), array('f', [-99.]) ]
  
   #opening rootfile containing SF for MVA correction
   in_fname_SF = '/lustrehome/fsimone/MVA_2017/MVA_control/ControlPlots/dsphipi_all_output_MVA_control_2017_UL_27june__MVAvalidation_mu2_SF.root'
   file_SF = ROOT.TFile(in_fname_SF, "READ")
   hbarrel = gROOT.FindObject("h2barrel");
   hendcap = gROOT.FindObject("h2endcap");

   for i in xrange(2): 
       new_tree[i].Branch("MuonID", globalMuonId[i], "MuonID/F")
       new_tree[i].Branch("MuonID_SF", scalefactor[i], "MuonID_SF/F")
       new_tree[i].Branch("MuonID_SFerr", scalefactor_err[i], "MuonID_SFerr/F")
   
   for i in xrange(2):
       nentries = old_tree[i].GetEntriesFast()
       print old_tree[i].GetName()
       for j in xrange(nentries):
          if (j%1000==0): print "Processing ",j,"/",nentries," ..."
          old_tree[i].GetEntry(j)
          muon_score = -99.0
          Muon_cLP[0]         = old_tree[i].mu_combinedQuality_chi2LocalPosition;
          Muon_cLM[0]         = old_tree[i].mu_combinedQuality_chi2LocalMomentum;
          Muon_staRelChi2[0]  = old_tree[i].mu_combinedQuality_staRelChi2;
          Muon_trkRelChi2[0]  = old_tree[i].mu_combinedQuality_trkRelChi2;
          Muon_glbdEP[0]      = old_tree[i].mu_combinedQuality_globalDeltaEtaPhi;
          Muon_trkKink[0]     = np.log(0.01+old_tree[i].mu_combinedQuality_trkKink);
          Muon_glbKink[0]     = np.log(0.01+old_tree[i].mu_combinedQuality_glbKink);
          #Muon_glbTrkP[0]     = old_tree[i].mu_combinedQuality_glbTrackProbability;
          Muon_nTVH[0]        = old_tree[i].mu_trackerLayersWithMeasurement;
          Muon_nVPH[0]        = old_tree[i].mu_Numberofvalidpixelhits;
          Muon_nMS[0]         = old_tree[i].mu_numberOfMatchedStations;
          Muon_segComp[0]     = old_tree[i].mu_segmentCompatibility;
          Muon_glbNChi2[0]    = old_tree[i].mu_GLnormChi2;
          Muon_inner_nChi2[0] = old_tree[i].mu_innerTrack_normalizedChi2;
          Muon_outer_nChi2[0] = old_tree[i].mu_outerTrack_normalizedChi2;
          Muon_innner_VF[0]   = old_tree[i].mu_innerTrack_validFraction;
      
          SF = [0, 1] 
          if (abs(old_tree[i].mu_eta)<1.2):
             muon_score = reader_MuonId_barrel.EvaluateMVA("BDT")
             SF = get_muon_SF(hbarrel, muon_score, old_tree[i].mu_pt)
          else:
             muon_score = reader_MuonId_endcap.EvaluateMVA("BDT")
             SF = get_muon_SF(hendcap, muon_score, old_tree[i].mu_pt)
       
          globalMuonId[i][0] = muon_score
          scalefactor[i][0] = SF[0]
          scalefactor_err[i][0] = SF[1]

          new_tree[i].Fill()
       
   new_file.cd()
   new_tree[0].Write()
   new_tree[1].Write()
   new_file.Close()
   file_.Close()

def main():
   #add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1029/AnalysedTree_data_2017B_control_16june.root")
   #add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1030/AnalysedTree_data_2017C_control_16june.root")
   #add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1031/AnalysedTree_data_2017D_control_16june.root")
   #add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1032/AnalysedTree_data_2017E_control_16june.root")
   #add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1033/AnalysedTree_data_2017F_control_16june.root")
   #add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1028/AnalysedTree_MC_2017DsPhiPi_control_16june.root")
   add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1109/AnalysedTree_data_2017B_control_28july.root")
   add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1110/AnalysedTree_data_2017C_control_28july.root")
   add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1111/AnalysedTree_data_2017D_control_28july.root")
   add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1112/AnalysedTree_data_2017E_control_28july.root")
   add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1113/AnalysedTree_data_2017F_control_28july.root")
   add_eval_branch("/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210728_1108/AnalysedTree_MC_2017DsPhiPi_control_28july.root")

if __name__=='__main__':
    main()

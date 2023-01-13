//#include "ROOT/RDataFrame.h"
/*#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
//#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "RooRealVar.h"
#include "RooBernstein.h"
#include "RooAbsPdf.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include <time.h>
#include "TBox.h"
//#include "TASImage.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TSystem.h"
//#include "TString.h"
//#include "TRandom.h"
//#include "TLorentzVector.h"
//#include "TTree.h"
*/


int Efficiencies(){

	TChain* mc = new TChain("Events");
	mc->Add(("/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/private_BtoDMuN_NtoL_nano_2tracks_Mass3/privateHNLnano_2021Jul28/210728_175518/0000/*.root")/*.c_str()*/);
	mc->Add(("/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/private_BtoDMuN_NtoL_nano_2tracks_Mass3/privateHNLnano_2021Jul28/210728_175518/0001/*.root")/*.c_str()*/);

	ROOT::RDataFrame df(*mc);
	const int nDispl=4;
	std::string lxy_cat[nDispl] = {"_sv_lxy<3","_sv_lxy<10","_sv_lxy<20","_sv_lxy>20"};
//	df.Define("TrigMuon_pt","TriggerMuon_pt[0]");
//	df.Define("TrigMuon_Eta","TriggerMuon_eta[0]");
	auto d1 = df.Define("muon_pt","Muon_pt[BToMuMuPi_sel_mu_idx]");
	auto d2 = d1.Define("muon_eta","Muon_eta[BToMuMuPi_sel_mu_idx]");
	std::string trg_mu_sel = " TriggerMuon_pt[0]>0.3 && fabs(TriggerMuon_eta[0])<2.8";
	std::string mu_sel= trg_mu_sel+" && muon_pt>0.3 && abs(muon_eta)<2";
		/* &&  ProbeTracks_pt[evt.BToMuMuPi_pi_idx[j]]>0.7 
		 &&  fabs(evt.ProbeTracks_eta[evt.BToMuMuPi_pi_idx[j]])<2 
		 &&  fabs(evt.ProbeTracks_dz[evt.BToMuMuPi_pi_idx[j]])>0.005 
		 &&  fabs(evt.ProbeTracks_dxy[evt.BToMuMuPi_pi_idx[j]])>0.005 
		 &&  fabs(evt.ProbeTracks_dxyS[evt.BToMuMuPi_pi_idx[j]])>3 
		 &&  fabs(evt.ProbeTracks_dzS[evt.BToMuMuPi_pi_idx[j]])>1.5 
		 &&  evt.BToMuMuPi_sv_prob[j]>0.001
		 &&  evt.BToMuMuPi_hnl_pt[j]>1
		 &&  fabs(evt.BToMuMuPi_hnl_cos2D[j])>0.9
		 &&  evt.BToMuMuPi_sv_lxy[j]/evt.BToMuMuPi_sv_lxye[j]>3 
		 &&  fabs(evt.Muon_dxy[evt.BToMuMuPi_sel_mu_idx[j]])>0.001
		 &&  fabs(evt.Muon_dz[evt.BToMuMuPi_sel_mu_idx[j]])>0.0015
		 &&  fabs(evt.Muon_dxy[evt.BToMuMuPi_sel_mu_idx[j]]/evt.Muon_dxyErr[evt.BToMuMuPi_sel_mu_idx[j]])>1.5
		 &&  fabs(evt.Muon_dz[evt.BToMuMuPi_sel_mu_idx[j]]/evt.Muon_dzErr[evt.BToMuMuPi_sel_mu_idx[j]])>1
		 &&  DeltaPhi(evt.ProbeTracks_phi[evt.BToMuMuPi_pi_idx[j]],evt.BToMuMuPi_fit_pi_phi[j] )<0.03
		 &&  fabs(evt.ProbeTracks_eta[evt.BToMuMuPi_pi_idx[j]]-evt.BToMuMuPi_fit_pi_eta[j] )<0.015
		 &&  deltaR_trgMu<0.5"; */
//	std::string ele_sel =  trg_mu_sel+" &&  ProbeTracks_pt[BToMuEPi_pi_idx]>0.7 ";
/*		 &&  fabs(evt.ProbeTracks_eta[evt.BToMuEPi_pi_idx[j]])<2 
		 &&  fabs(evt.ProbeTracks_dz[evt.BToMuEPi_pi_idx[j]])>0.005 
		 &&  fabs(evt.ProbeTracks_dxy[evt.BToMuEPi_pi_idx[j]])>0.005 
		 &&  fabs(evt.ProbeTracks_dxyS[evt.BToMuEPi_pi_idx[j]])>3 
		 &&  fabs(evt.ProbeTracks_dzS[evt.BToMuEPi_pi_idx[j]])>1.5 
		 &&  fabs(evt.ProbeTracks_DCASig[evt.BToMuEPi_pi_idx[j]])>5 
		 &&  evt.BToMuEPi_sv_prob[j]>0.001
		 &&  fabs(evt.BToMuEPi_hnl_cos2D[j])>0.99 
		 &&  evt.BToMuEPi_sv_lxy[j]/evt.BToMuEPi_sv_lxye[j]>3 
		 &&  fabs(evt.Electron_dxy[evt.BToMuEPi_sel_e_idx[j]])>0.008
		 &&  fabs(evt.Electron_dz[evt.BToMuEPi_sel_e_idx[j]])>0.008
		 &&  fabs(evt.Electron_dxy[evt.BToMuEPi_sel_e_idx[j]]/evt.Electron_dxyErr[evt.BToMuEPi_sel_e_idx[j]])>1.5
		 &&  fabs(evt.Electron_dz[evt.BToMuEPi_sel_e_idx[j]]/evt.Electron_dzErr[evt.BToMuEPi_sel_e_idx[j]])>1
		 &&  deltaR_trgMuPi<2;*/
//	std::string tt_sel =  trg_mu_sel+" &&  ProbeTracks_pt[BToMuEPi_pi_idx]>0.7 ";
	std::string mu_evt = ("All(BToMuMuPi_isMatched && "+mu_sel+")").c_str();
//	std::string ele_evt = ("(!"+mu_evt +"(BToMuEPi_isMatched && "+ele_sel+"))").c_str();

//	std::string tt_evt = ("((!"+mu_evt+" && !"+ele_evt+") && "+tt_sel+")").c_str();
	TH1D* hnl_mass[nDispl];

/*	std::vector<ROOT::RDataFrame*> d_mu;
	std::vector<ROOT::RDataFrame> d_PF;
	std::vector<ROOT::RDataFrame> d_lPt;
	std::vector<ROOT::RDataFrame> d_tt;*/
	for(int i=0; i<nDispl;i++){
	mu_sel += "BToMuMuPi"+lxy_cat[i];
//	ele_sel += "BToMuEPi"+lxy_cat[i];
//	tt_sel += "BToMuEPiHD"+lxy_cat[i];
	hnl_mass[i] = d2.Filter(mu_evt.c_str()).Histo1D("BToMuMuPi_hnl_mass").GetPtr();
/*	d_PF[i] = df.Filter((ele_evt+" && BToMuEPi_isPF").c_str());
	d_lPt[i] = df.Filter((ele_evt+" && BToMuEPi_isLowPt && BToMuEPi_isPFoverlap == 0").c_str());
	d_tt[i] = df.Filter(tt_evt.c_str());*/
	}
	hnl_mass[0]->Draw();
	
	return 1;	
}

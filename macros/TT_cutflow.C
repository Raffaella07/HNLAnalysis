#include "../src/FlatBToMuLPiClass.cc"
#include <iostream>
//#include <utility>
//#include <vector>
//#include <math.h>
//#include <string>
//#include <iostream>
//#include <fstream>
//#include "TH1.h"
//#include "TF1.h"
//#include "TH2.h"
//#include "TMath.h"
//#include "TFile.h"
//#include "TAxis.h"
//#include "RooRealVar.h"
//#include "RooBernstein.h"
//#include "RooAbsPdf.h"
//#include "RooBinning.h"
//#include "RooWorkspace.h"
//#include "RooFFTConvPdf.h"
//#include "RooFitResult.h"
//#include "RooDataHist.h"
//#include "RooHistPdf.h"
//#include "RooDataSet.h"
//#include "RooKeysPdf.h"
//#include "RooPlot.h"
//#include "TLegend.h"
//#include "TCanvas.h"
//#include "TPad.h"
//#include "TLatex.h"
//#include "TLine.h"
//#include <time.h>
//#include "TBox.h"
////#include "TASImage.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TSystem.h"
//#include "TString.h"
//#include "TRandom.h"
//#include "TLorentzVector.h"
//#include "TTree.h"



int TT_cutflow(){

	TChain* c = new TChain("Events");
	c->Add("../python/BToMuLPiSmall_Channel3.root");

	FlatBToMuLPiClass evt;
	evt.Init(c);	

	int num[16]={0};
	int num_single[16]={0};
	int i;
	int den;
	for(i=0;i<evt.fChain->GetEntries();i++){

	evt.fChain->GetEntry(i);
	if (evt.hnl_l_isEle) continue;
	den++;
	bool sel_array[16];
	sel_array[0]= evt.TrgMu_pt>7;
	sel_array[1]= sel_array[0] && fabs(evt.TrgMu_eta)<1.5; 
	sel_array[2] = sel_array[1] && fabs(evt.hnl_pi_eta)<2;
	sel_array[3] = sel_array[2]&& fabs(evt.hnl_pi_dz)>0.005;
	sel_array[4] = sel_array[3]&& fabs(evt.hnl_pi_dxy)>0.005;
	sel_array[5] = sel_array[4]&& fabs(evt.hnl_pi_dxyS)>3;
	sel_array[6] = sel_array[5]&& fabs(evt.hnl_pi_dzS)>1.5;
	sel_array[7] = sel_array[6]&& fabs(evt.hnl_pi_DCAS)>5;
	sel_array[8] = sel_array[7]&& evt.hnl_vtxProb>0.001;
	sel_array[9] = sel_array[8]&& fabs(evt.hnl_cos2D)>0.99;
	sel_array[10] = sel_array[9]&& evt.hnl_lxy_sig>20;
	sel_array[11] = sel_array[10]&& fabs(evt.hnl_l_dxy)>0.008;
	sel_array[12] = sel_array[11]&& fabs(evt.hnl_l_dz)>0.008;
	sel_array[13] = sel_array[12]&& fabs(evt.hnl_l_dxyS)>1.5;
	sel_array[14] = sel_array[13]&& fabs(evt.hnl_l_dzS)>1;
	sel_array[15] = sel_array[14]&& evt.dr_trgmu_hnl<0.6;


	for (int j=0; j<16;j++){
		std::cout << "trg mu pt " << evt.TrgMu_pt << std::endl;	
		if (sel_array[j]==1) num[j]++;
		if (sel_array[j]==1 && evt.CandIdx ==0 ) num_single[j]++;
		
	}

	}


	for (int j=0; j<16;j++){
	
	if(j==0){
	std::cout << "selection n " << j << "efficiency " << (double)(1.0*num[j]/den) <<std::endl;
	std::cout << "single cand " << j << "efficiency " << (double)(1.0*num_single[j]/den) <<std::endl;
	}
	else{ std::cout << "selection n " << j << "efficiency " << (double)(1.0*num[j]/num[j-1]) <<std::endl;
	std::cout << "selection n " << j << "efficiency " << (double)(1.0*num_single[j]/num_single[j-1]) <<std::endl;
	}
	}
	std::cout << "selection toto efficiency " << (double)(1.0*num[0]/den) <<std::endl;
	std::cout << "selection toto efficiency " << (double)(1.0*num_single[0]/den) <<std::endl;
	return 1; 
}

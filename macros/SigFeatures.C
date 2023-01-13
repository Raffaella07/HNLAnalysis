#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
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
#include "../src/BParkTools.cc"
//#include "TASImage.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TSystem.h"
//#include "TString.h"
//#include "TRandom.h"
//#include "TLorentzVector.h"
//#include "TTree.h"




int SigFeatures(std::string filename){


	TFile* f = TFile::Open(filename.c_str());
	TTree* tree = (TTree*) f->Get("Events");
	
	TH1D* h[4][4];
	int nLxyBin = 4;
	int nCategories = 4;
	int i, j;

	for(i=0;i<nCategories;i++){

		for(j=0;j<nLxyBin;j++){


	h[i][j]=new  TH1D((std::to_string(i)+std::to_string(j)).c_str(),(std::to_string(i)+std::to_string(j)).c_str(),30,2.8,3.1);
	TCanvas* c = new TCanvas("c","c",600,800);

	h[i][j]->GetXaxis()->SetTitle("hnl mass (GeV)");
	h[i][j]->GetYaxis()->SetTitle("entries");
	h[i][j]->SetLineWidth(2);
	h[i][j]->SetLineColor(kBlue);
	tree->Draw((std::string("hnl_mass>>")+h[i][j]->GetName()+"(30,2.8,3.1)").c_str(),("Type=="+std::to_string(i)+" &&LxyBin=="+std::to_string(i)).c_str());

	c->SaveAs((std::to_string(i)+"_"+std::to_string(j)+".pdf").c_str());	
	
		}
	}

	return 1;	
		
}
int TTFeatures(std::string filename){


	TFile* f = TFile::Open(filename.c_str());
	TTree* tree = (TTree*) f->Get("Events");
	TH1D* h[2][4];
	int nLxyBin = 4;
	int nFlav = 2;
	int i, j;
	
	setStyle();

	for(i=0;i<nFlav;i++){
		for(j=0;j<nLxyBin;j++){


	h[i][j]=new TH1D(("flav_"+std::to_string(i)+std::to_string(j)).c_str(),("flav_"+std::to_string(i)+std::to_string(j)).c_str(),30,2.8,3.1);
	h[i][j]->GetXaxis()->SetTitle("hnl mass (GeV)");
	h[i][j]->GetYaxis()->SetTitle("entries");
	h[i][j]->SetLineWidth(3);
	h[i][j]->SetLineColor(kBlue);
	TCanvas *c = new TCanvas("c","c",600,800);

	if (i ==0) tree->Draw((std::string("hnl_mass>>")+h[i][j]->GetName()+"(30,2.5,3.2)").c_str(),("Type==3 && hnl_l_isEle == 1 && LxyBin=="+std::to_string(j)).c_str());
	else tree->Draw((std::string("hnl_mass>>")+h[i][j]->GetName()+"(30,2.5,3.2)").c_str(),("Type==3 && hnl_l_isMu == 1 && LxyBin =="+std::to_string(j)).c_str());

	c->SaveAs(("../plots/CategorizedPlots/flav_"+std::to_string(i)+"_"+std::to_string(j)+".pdf").c_str());	
	delete c;		
		}
	}

	TCanvas *big = new TCanvas("bog","big",2000,1400);
	std::cout << "here" << std::endl;
	big->Divide(4,2);
	std::cout << "here" << std::endl;
	gPad->cd(1);	
	h[0][0]->Draw();
	std::cout << "here" << std::endl;
	big->cd(2);	
	std::cout << "here" << std::endl;
//	big->SaveAs("../plots/CategorizedPlots/flavGrid.pdf");
	return 1;
}

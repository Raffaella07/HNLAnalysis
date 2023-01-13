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
//#include "TASImage.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TSystem.h"
//#include "TString.h"
//#include "TRandom.h"
//#include "TLorentzVector.h"
//#include "TTree.h"



int LxyEff(std::string dataset,int nBin, double min, double max){

	
	TFile* in = TFile::Open(filename.c_str());
	TChain* c= new TChain("Events");

	c->Add((filename+"/*.root").c_str());

	std::string lables[6] = {"gen","muon","Pf","LowPt","MuTrack","EleTrack"}


	TH1D* h[6];
	
	
	



}

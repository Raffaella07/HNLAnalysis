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
#include "/cmshome/ratramon/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/src/RooDoubleCBFast.cc"
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
#include "setStyle.C"
//#include "TASImage.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TSystem.h"
//#include "TString.h"
//#include "TRandom.h"
//#include "TLorentzVector.h"
//#include "TTree.h"



int fit(std::string filename, std::string path, double mass,double ctau, int type, std::string wd, double* sigF,std::string sel,std::string presel_file,std::string sel_file){

//	TFile* in = TFile::Open(filename.c_str());

	TChain* c= new TChain("Events");
        bool doOnlyHistos= false;
	c->Add(filename.c_str());
	
	TChain* Data= new TChain("Events");
	std::string datapath[8] = {path+"/Parking1D0_fullPF",path+"/Parking1D1_fullPF",
	path+"/Parking1D2_fullPF",path+"/Parking1D3_fullPF",path+"/Parking1D4_fullPF",
	path+"/Parking1D5_fullPF",path+"/Parking1D6_fullPF",path+"/Parking1D7_fullPF"};
	for(int i=0; i<8;i++){

	Data->Add((datapath[i]+"/*.root").c_str());
	std::cout << Data->GetEntries() << std::endl;
	} 

	std::string lxy[3] = {"hnl_lxy<1","hnl_lxy>1 && hnl_lxy<5","hnl_lxy>5"};

	RooWorkspace wspace("wS");
	gStyle->SetOptFit(0000);
//	gROOT->SetBatch(true);
	gROOT->SetStyle("Plain");
	gStyle->SetGridStyle(3);
	gStyle->SetOptStat(000000);
	gStyle->SetOptTitle(0);
	float SigMu[4]= {336,538,386,173};
	float SigEle[4]= {133,347,323,284};
	std::string channel[4]= {"Muon","PFele","LoePtele","Track"};
	wspace.factory("r_ele[0.359,0.359,0.359]");
	wspace.factory("r_mu[0.641,0.641,0.641]");
	wspace.factory("hnl_mass[3,0,7]");
	wspace.factory("B_mass[5.3,0,1000000]");
	wspace.factory("LepQProd[0,-3,3]");
	wspace.factory("hnl_charge[0,-2,2]");
	wspace.factory("Type[0,0,4]");
	wspace.factory("LxyBin[0,0,4]");
	wspace.factory("dr_trgmu_hnl[0,-100,100]");
	wspace.factory("dilepton_mass[0,0,1000000]");
	wspace.factory("dilepton_pt[0,0,10000000]");
	wspace.factory("hnl_vtxProb[0,0,1]");
	wspace.factory("hnl_pi_dxyS[0,-100000000,10000000]");
	wspace.factory("hnl_l_dxyS[0,-100000000,100000000]");
	wspace.factory("hnl_l_pt[10,0,100000000]");
	wspace.factory("hnl_l_mvaId[0,-2000,2000]");
	wspace.factory("hnl_pi_pt[10,0,10000000]");
	wspace.factory("hnl_lxy_sig[10,-10000000,10000000]");
	wspace.factory("hnl_lxy[10,-100000000,10000000]");
	wspace.factory("likelihood[0,-1,1]");
	wspace.factory("hnl_cos2D[0.995,0,1]");
	wspace.factory("c0[0.1,-1,1]");
	wspace.factory("c1[0.5,-10,10]");
	wspace.factory("c2[0,-10,10]");
	wspace.factory("weight[1,1,1]");
	wspace.factory("Nbkg[100,0,10000]");
	//wspace.factory(("Nele[10,7,"+std::to_string(hist->GetMaximum())+"]").c_str());
	//wspace.factory(("Nmu[10,7,"+std::to_string(hist->GetMaximum())+"]").c_str());
	wspace.var("r_ele");
	wspace.var("r_mu");
	wspace.var("Nsig");
	wspace.var("Nele");
	wspace.var("Nmu");
	wspace.var("N1");
	wspace.var("N2");
	wspace.var("Nbkg");
	wspace.var("width");
	wspace.var("sigma_ele");
	wspace.var("sigma_mu");
	wspace.var("n_cb");
	wspace.var("alpha_cb");
	wspace.var("hnl_mass");
	wspace.var("hnl_cos2D");
	wspace.var("hnl_vtxProb");
	wspace.var("hnl_charge");
	wspace.var("hnl_lxy_sig");
	wspace.var("hnl_lxy");
	wspace.var("hnl_trgmu_hnl");
	wspace.var("dilepton_mass");
	wspace.var("dilepton_pt");
	wspace.var("hnl_l_dxyS");
	wspace.var("hnl_pi_dxyS");
	wspace.var("hnl_l_pt");
	wspace.var("hnl_pi_pt");
	wspace.var("B_mass");
	wspace.var("LepQProd");
	wspace.var("Type");
	wspace.var("LxyBin");
	wspace.var("c0");
	wspace.var("c1");
	wspace.var("c2");
	wspace.var("weigth");
	wspace.var("weight")->setConstant(kTRUE);

	//Define set of variables to include in dataset

	std::string variables = "B_mass,hnl_mass,hnl_charge,hnl_lxy,Type,LxyBin,dr_trgmu_hnl,dilepton_mass,hnl_vtxProb,hnl_cos2D,hnl_l_pt,hnl_l_mvaId,hnl_pi_pt,hnl_l_dxyS,hnl_pi_dxyS,dilepton_pt,likelihood,hnl_lxy_sig,LepQProd";
        wspace.defineSet("treeSet",variables.c_str());
//	MC dataset initialization
	RooDataSet mcSet("mcSet","mcSet",*wspace.set("treeSet"),RooFit::Import(*c),RooFit::Cut(("Type=="+std::to_string(type)).c_str()));

	double SigEff[3],BackgroundEff[3],nbin[3];
	std::string title[3]= {"lxy0to1","lxy1to5","lxyOver5"};
	
	for (int i =0; i<3; i++){

	SigEff[i] = 1.0* mcSet.sumEntries((lxy[i]+sel+" && likelihood>0.9 ").c_str())/mcSet.sumEntries((lxy[i]+sel).c_str());
	nbin[i] = i;
	std::cout << "sig eff " << lxy[i] << " " << SigEff[i]<< std::endl;
	}


	TMultiGraph* m = new TMultiGraph();
	TGraph* g = new TGraph(presel_file.c_str(),"%lg %lg");
	TGraph* g_l = new TGraph(sel_file.c_str(),"%lg %lg"); 
	g->SetLineWidth(2);
	g_l->SetLineWidth(2);
	g->SetMarkerStyle(8);
	g_l->SetMarkerStyle(8);
	g->SetTitle("Preselection");
	g_l->SetTitle("Likelihood>0.9");
	m->Add(g,"LP");
	m->Add(g_l,"LP");

	TGraph* eff = new TGraph(3,nbin,SigEff);
	eff->SetMarkerStyle(8);
	eff->SetLineStyle(2);
	eff->SetMarkerColor(kBlue);
	eff->SetLineColor(kBlue);
	TCanvas * canva = new TCanvas("c","c",800,600);
	setStyle();
	canva->Divide(1,2);
	canva->cd(1);
	gPad->SetBottomMargin(0.01);//0.13);
	gPad->SetLogy();
	m->Draw("ALP");
	m->GetHistogram()->GetXaxis()->Set(3,-0.5,2.5);
	m->GetHistogram()->GetYaxis()->SetTitle("Limit @ 95% C.L.");
	for (int i=1;i<=3;i++) m->GetHistogram()->GetXaxis()->SetBinLabel(i,title[i-1].c_str());
	m->GetHistogram()->LabelsDeflate();
	m->GetYaxis()->SetLabelSize(20);
	gStyle->SetPalette(91);
	m->Draw("A PMC PLC");

	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(43);
	gStyle->SetLegendTextSize(15);

	TLegend* l = gPad->BuildLegend(0.3,0.3,0.3,0.3);
	l->SetFillStyle(0);
	canva->cd(2);
	gPad->SetTopMargin(0);//0.13);
	TH2D* plotter = new TH2D("plotter","plotter",3,-0.5,2.5,10,0,1);
	plotter->GetYaxis()->SetLabelSize(15);
	for (int i=1;i<=3;i++) plotter->GetXaxis()->SetBinLabel(i,title[i-1].c_str());
	plotter->GetYaxis()->SetTitle("N^{gen}_{likelihood>0.9}/N^{gen}_{preselection}");
	plotter->Draw();
	eff->Draw("sameLP");

	canva->SaveAs("LikelihoodEffect.pdf");	
	
	return 0;

}

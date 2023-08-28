#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>
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
#include "LifetimeReweight.C"
#include "/cmshome/ratramon/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/src/RooDoubleCBFast.cc"
#include "/cmshome/ratramon/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/src/RooMultiPdf.cxx"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include <time.h>
#include "TBox.h"
#include "TPython.h"
//#include "TASImage.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TSystem.h"
//#include "TString.h"
//#include "TRandom.h"
//#include "TLorentzVector.h"
//#include "TTree.h"



int genEventCount(std::string filename){

	long int evtCounter=0;
	long int TotCounter =0;
	TChain* mc = new TChain("Runs");
	mc->Add(filename.c_str());
	mc->SetBranchAddress("genEventCount",&evtCounter);


	for(int i=0;i<mc->GetEntries();i++){
		//		std::cout <<"running on " << data->GetEntries() << " tree" << std::endl;
		mc->GetEntry(i);

		//	std::cout <<"running on " << i << " tree" << std::endl;
		//std::cout <<"single counter" << evtCounter  << std::endl;
		TotCounter+=evtCounter;
	}

	return TotCounter;



}

void single_fit(std::string csvIn,std::string lxyLable,std::string sigLable, bool doNNsel,double mass,double ctau,std::string lxy,std::string sig,std::string sel, bool debug, double nn_cut,double* sigF,std::string wd,int nSigma){

	
	std::vector<double> GenEvent_Bu, GenEvent_Bd, GenEvent_Bs,GenEvtCount,Filter_Bu,Filter_Bd,Filter_Bs,Filter,Mass,Ctau,Sigma,Flav_ratio;
	std::vector<std::string> Files, label;
	std::string B_c = "Bc_";
	std::size_t found = csvIn.find(std::string("Bc"));
	bool Bc;
	if (found == std::string::npos){
		 Bc= false;
		 B_c = "";
	}
	else  Bc = true;
	
	//sig_fromcsv("../data/MC_datasets_InclusiveFilterEff.csv",&Files,&label,&GenEvtCount,&Filter,&Mass,&Ctau,&Sigma,&Flav_ratio);
	sig_fromcsv(csvIn.c_str(),&Files,&label,&GenEvent_Bu,&GenEvent_Bd,&GenEvent_Bs,&GenEvtCount,&Filter_Bu,&Filter_Bd,&Filter_Bs,&Filter,&Mass,&Ctau,&Sigma,&Flav_ratio);
	
//	std::cout << Ctau[0] <<std::endl;
/*	std::string files[5][4] = {{"/cmshome/ratramon/Analysis/data/HNLFlatTuples/NewSignature_1p0_ctau10p0/","../data/HNLFlatTuples/NewSignature_1p0_ctau100p0/","../data/HNLFlatTuples/NewSignature_1p0_ctau1000p0/","placeholder"},
		{"../data/HNLFlatTuples/NewSignature_1p5_ctau10p0/","../data/HNLFlatTuples/NewSignature_1p5_ctau100p0/","../data/HNLFlatTuples/NewSignature_1p5_ctau1000p0/","placeholder"},
		{"../data/HNLFlatTuples/NewSignature_2p0_ctau10p0/","../data/HNLFlatTuples/NewSignature_2p0_ctau100p0/","../data/HNLFlatTuples/NewSignature_2p0_ctau1000p0/","placeholder"},
		{"/cmshome/ratramon/Analysis/data/HNLFlatTuples/NewSignature_3p0_ctau1p0/","../data/HNLFlatTuples/NewSignature_3p0_ctau10p0/","../data/HNLFlatTuples/NewSignature_3p0_ctau100p0/","../data/HNLFlatTuples/NewSignature_3p0_ctau1000p0/"},
		{"../data/HNLFlatTuples/NewSignature_4p5_ctau0p1/","../data/HNLFlatTuples/NewSignature_4p5_ctau1p0/","../data/HNLFlatTuples/NewSignature_4p5_ctau10p0/","../data/HNLFlatTuples/NewSignature_4p5_ctau100p0/"},
	};*/
/*	double genEvtCount[5][4]={{9.080620000000000000e+05,4.532249000000000000e+06,1.312593500000000000e+07,-1},
		{7.857100000000000000e+05,1.955714000000000000e+06,7.546074000000000000e+06,-1},
		{5.673050000000000000e+05,1.916707000000000000e+06,8.421721000000000000e+06,-1},
		{5.252660000000000000e+05,6.562680000000000000e+05,1.314002000000000000e+06,3.282584000000000000e+06},
		{5.914900000000000000e+05,5.065250000000000000e+05,5.962870000000000000e+05,1.195035000000000000e+06}};
*/

	//update with filter efficiencies for specific channel
/*	double filter[5][4]={{9.49e-03,	7.90e-03,2.03e-03,-1},
		{8.50e-03,7.36e-03,2.13e-03,-1},
		{7.99e-03,7.12e-03,2.05e-03,-1},
		{1.69e-02,1.69e-02,1.59e-02,5.78e-03},
		{2.58e-02,2.58e-02,2.58e-02,2.51e-02}};
*/
	TChain* c= new TChain("Events");
	int Ntot =0;
	float FilterEff = 0;
	bool doOnlyHistos= false;
/*	double grid[5][4]={     {10,100,1000,-1},
		{10,100,1000,-1},
		{10,100,1000,-1},
		{1,10,100,1000},
		{0.1,1,10,100},
	};*/

	double hnlMass[6]={1,1.5,2,3,4.5,5.5};
	int flavor = 4;
	std::string f_label[5] = {"_Bu","_Bd","_Bs","_Bc",""};
	std::string flav = "";
	bool reweight = 1;//filename.find(".root") != std::string::npos ;
	//std::cout  << filename << "is central sample? " << no_reweight << std::endl;

//		c->Add((filename).c_str());
		//	std::cout << c->GetEntries() << std::endl;
		std::cout << Files.size() << std::endl;
		for(int i=0;i<Files.size();i++){
		//	std::cout << "idx" << i <<" " << Files[i] << " " <<Ctau[i] << " " << ctau <<" " <<mass << " " << Mass[i] << std::endl;
			if (mass==Mass[i] && ctau==Ctau[i]){
				c->Add((Files[i]).c_str());
				reweight = 0;
				if( std::string(sel).find(std::string("is_B")) == std::string::npos  ){
						Ntot=GenEvtCount[i] * Flav_ratio[i];
						if (Bc){
							  Filter[i] = Filter[i]*0.08244;
							  flavor = 3;
				
						}
						FilterEff = Filter[i];
						
					}
				else if (std::string(sel).find(std::string("&& is_Bu")) != std::string::npos ){ 
							Ntot=GenEvent_Bu[i] * Flav_ratio[i];
							FilterEff = Filter_Bu[i];
							flavor = 0;
							}
				else if (std::string(sel).find(std::string("&& is_Bd")) != std::string::npos ){ 
							Ntot=GenEvent_Bd[i] * Flav_ratio[i];
							FilterEff = Filter_Bu[i];
							flavor = 1;
							}
				else if (std::string(sel).find(std::string("&& is_Bs")) != std::string::npos ){ 
							Ntot=GenEvent_Bs[i] * Flav_ratio[i];
							FilterEff = Filter_Bs[i];
							flavor = 2;

							}	
				}	
		}
			
		
	if (reweight){		
		for(int i=0;i<Files.size();i++){

		//	std::cout <<"______________________________________________" <<sel << std::endl;
			if (mass==Mass[i] && Ctau[i]> 0.){

				//		std::cout << i  << "mass " << mass << std::endl;
				//		std::cout << "in if "<<  i  <<  Ctau[i] << " " << ctau << std::endl;
					//	std::cout << i  <<  Filter[i] <<  std::endl;


						c->Add((Files[i]).c_str());
				//		std::cout <<"second out " <<  i  <<" " <<   Ctau[i] << " " << ctau << std::endl;
						if( std::string(sel).find(std::string("is_B")) == std::string::npos  ){

							if (Bc){
								Filter[i] = Filter[i]*0.08244 ;
				    				flavor = 3;
							}

							Ntot += GenEvtCount[i] * Flav_ratio[i] ;
							FilterEff += GenEvtCount[i]*Filter[i] * Flav_ratio[i];
							std::cout<< GenEvtCount[i] << " " << Filter[i] <<std::endl;
							//std::cout<< "in filter " << i  <<" " <<   Ctau[i] << " " << ctau << " " << Flav_ratio[i]<<  std::endl;
						
						}
						else if (std::string(sel).find(std::string("&& is_Bu")) != std::string::npos ){ 
								flavor = 0;
								Ntot += GenEvent_Bu[i] * Flav_ratio[i];
								FilterEff += GenEvent_Bu[i]*Filter_Bu[i]* Flav_ratio[i];
								}
						else if (std::string(sel).find(std::string("&& is_Bd")) != std::string::npos ){ 
								flavor = 1;
								Ntot += GenEvent_Bd[i] * Flav_ratio[i];
								FilterEff += GenEvent_Bd[i]*Filter_Bd[i]* Flav_ratio[i];
								}
						else if (std::string(sel).find(std::string("&& is_Bs")) != std::string::npos ){ 
								flavor = 2;
								Ntot += GenEvent_Bs[i] * Flav_ratio[i];
								FilterEff += GenEvent_Bs[i]*Filter_Bs[i]* Flav_ratio[i];
	
								}	


	//						std::cout << i  <<" " <<  FilterEff << " " << ctau << std::endl;
					}

				}


			
		
		FilterEff = FilterEff/Ntot;
		std::cout << FilterEff <<" "  <<Ntot << std::endl;
	}

	std::cout << "mass here " << flavor <<  std::endl;
	flav = f_label[flavor] ; 
	std::string modelPath;
	std::string category;
	if (doNNsel){

		category = lxyLable+"_"+sigLable;
		//	std::string modelPath = "/cmshome/ratramon/Analysis/python/model/";
		//if (Bc)	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_Bc_26Jan2023_00h25m43s/";
		if (Bc)	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_Bc_28Jul2023_09h10m22s/";
		else	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_26Jul2023_17h22m57s//";
	//	else	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_12Oct2022_15h25m05s//";
//		modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_muon/";
		std::cout << modelPath << std::endl;
		TPython::LoadMacro("../python/ComputeNNScore.py");
		TPython::Exec(("modelSpecs = NNinput('"+modelPath+"')").c_str());
		TPython::Exec(("modelSpecs = NNinput('"+modelPath+"')").c_str());
		if (debug) std::cout << category << std::endl;
		if (debug) std::cout << modelPath << std::endl;

	}else{

		category = "inclusive";
	}
	std::cout << "mass here" << std::endl;
	std::string var_label = category + B_c; 
	RooWorkspace wspace("wS");
	gStyle->SetOptFit(0000);
	gROOT->SetBatch(true);
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
	wspace.factory("B_pt[1,0,10000000]");
	wspace.factory("is_Bu[0.5,-2,1000000]");
	wspace.factory("is_Bd[0.5,-2,1000000]");
	wspace.factory("is_Bs[0.5,-2,1000000]");
	wspace.factory("LepQProd[0,-3,3]");
	wspace.factory("hnl_charge[0,-2,2]");
	wspace.factory("Type[0,0,4]");
	wspace.factory("LxyBin[0,0,4]");
	wspace.factory("dr_trgmu_hnl[0,-100,100]");
	wspace.factory("dr_trgMu_lep[0,-100,100]");
	wspace.factory("dr_Blep_pi[0,-100,100]");
	wspace.factory("dilepton_mass[0,0,1000000]");
	wspace.factory("dilepton_pt[0,0,10000000]");
	wspace.factory("BlepPi_mass[0,0,1000000]");
	wspace.factory("hnl_vtxProb[0,0,1]");
	wspace.factory("hnl_vtxChi2[0,0,1000000]");
	wspace.factory("hnl_pi_PixelLayers[1,0,10000000]");
	wspace.factory("hnl_pi_TrackerLayers[1,0,10000000]");
	wspace.factory("hnl_pi_DCAS[0,-100000000,10000000]");
	wspace.factory("hnl_l_dxyS[0,-100000000,100000000]");
	wspace.factory("hnl_l_pt[10,0,100000000]");
	wspace.factory("l0_pt[10,0,100000000]");
	wspace.factory("l_pt[10,0,100000000]");
	wspace.factory("hnl_l_eta[0,-5,5]");
	wspace.factory("hnl_l_relIso[10,0,100000000]");
	wspace.factory("hnl_l_mvaId[0,-2000,2000]");
	wspace.factory("hnl_pi_pt[10,0,10000000]");
	wspace.factory("TrgMu_pt[10,0,10000000]");
	wspace.factory("TrgMu_relIso[10,0,10000000]");
	wspace.factory("TrgMu_PixelLayers[1,0,10000000]");
	wspace.factory("TrgMu_TrackerLayers[1,0,10000000]");
	wspace.factory("l_relIso[10,0,10000000]");
	wspace.factory("l0_relIso[10,0,10000000]");
	wspace.factory("hnl_lxy_sig[10,-10000000,10000000]");
	wspace.factory("hnl_lxy[10,-100000000,10000000]");
	wspace.factory("hnl_ct[10,0,100000]");
	wspace.factory("likelihood[0,-100,100]");
	wspace.factory("hnl_cos2D[0.995,0,1]");
	wspace.factory("weight[1,0,10000000000000]");
	wspace.factory("PU_weight[1,0,10000000000000]");
	wspace.factory("Trg_weight[1,0,10000000000000]");
	wspace.factory("muId_weight[1,0,10000000000000]");
	wspace.factory("c0[0.1,-1,1]");
	wspace.factory("c2[0,-10,10]");
	wspace.factory("Nbkg[100,0,10000]");
	wspace.var("r_ele");
	wspace.var("r_mu");
	wspace.var("is_Bu");
	wspace.var("is_Bd");
	wspace.var("is_Bs");
	wspace.var("Nsig");
	wspace.var("Nele");
	wspace.var("Nmu");
	wspace.var("N1");
	wspace.var("N2");
	wspace.var("Nbkg");
	wspace.var("width");
	wspace.var(("sigma_ele_"+var_label).c_str());
	wspace.var("sigma_mu");
	wspace.var("n_cb");
	wspace.var("alpha_cb");
	wspace.var("hnl_mass");
	wspace.var("hnl_cos2D");
	wspace.var("hnl_vtxProb");
	wspace.var("hnl_vtxChi2");
	wspace.var("hnl_charge");
	wspace.var("hnl_lxy_sig");
	wspace.var("hnl_lxy");
	wspace.var("hnl_ct");
	wspace.var("dr_Blep_pi");
	wspace.var("dr_trgMu_lep");
	wspace.var("dilepton_mass");
	wspace.var("dilepton_pt");
	wspace.var("BlepPi_mass");
	wspace.var("hnl_l_dxyS");
	wspace.var("hnl_pi_DCAS");
	wspace.var("hnl_l_pt");
	wspace.var("l0_pt");
	wspace.var("l_pt");
	wspace.var("hnl_l_relIso");
	wspace.var("hnl_pi_pt");
	wspace.var("hnl_pi_PixelLayers");
	wspace.var("hnl_pi_TrackerLayers");
	wspace.var("TrgMu_pt");
	wspace.var("TrgMu_relIso");
	wspace.var("TrgMu_PixelLayers");
	wspace.var("TrgMu_TrackerLayers");
	wspace.var("l_relIso");
	wspace.var("l0_relIso");
	wspace.var("B_mass");
	wspace.var("B_pt");
	wspace.var("LepQProd");
	wspace.var("Type");
	wspace.var("LxyBin");
	wspace.var("c0");
	wspace.var("PU_weight");
	wspace.var("weight");
	wspace.var("Trg_weight");
	wspace.var("muId_weight");

	wspace.var("c1");
	wspace.var("c2");
	int Nevents;
	//Define set of variables to include in dataset

	//	std::string variables = "B_mass,hnl_mass,hnl_charge,hnl_lxy,Type,LxyBin,dr_trgmu_hnl,dilepton_mass,hnl_vtxProb,hnl_cos2D,hnl_l_pt,hnl_l_mvaId,hnl_pi_pt,hnl_l_dxyS,hnl_pi_dxyS,dilepton_pt,likelihood,hnl_lxy_sig,LepQProd,hnl_ct,nn_score,nn_score_single";
	std::string variables = "hnl_mass,hnl_charge,likelihood,hnl_ct,hnl_lxy,LepQProd,hnl_cos2D,hnl_vtxProb,hnl_pi_pt,hnl_pi_DCAS,B_mass,hnl_vtxChi2,l_pt,hnl_l_pt,TrgMu_pt,hnl_l_eta,l0_pt,hnl_lxy_sig,B_pt,dr_Blep_pi,dr_trgMu_lep,dilepton_mass,BlepPi_mass,TrgMu_relIso,hnl_l_relIso,l0_relIso,l_relIso,hnl_l_mvaId,PU_weight,Trg_weight,muId_weight,weight,is_Bu,is_Bd,is_Bs,hnl_pi_PixelLayers,hnl_pi_TrackerLayers,TrgMu_PixelLayers,TrgMu_TrackerLayers";
	wspace.defineSet("treeSet",variables.c_str());
	double sigma = 0.00697693*mass +0.00185739 ;

	//	MC dataset initialization
	RooDataSet mcSet_un("mcSet_un","mcSet_un",c,*wspace.set("treeSet"),(lxy+sig+sel).c_str())/*("Type=="+std::to_string(type)+" && "+lxy +sel).c_str()))*/;
	RooArgSet * entry;
	RooRealVar* nn_score= new RooRealVar("nn_score","nn_score",0.5,0,1000000);
	RooRealVar* weight= new RooRealVar("weight","weight",1,0,1000000000);
	RooRealVar* hnl_mass= new RooRealVar("hnl_mass","hnl_mass",1,0,1000000000);
	RooDataSet* temp = new RooDataSet("temp","",RooArgSet(*hnl_mass,*weight),weight->GetName()); 
	//	temp= (RooDataSet*) mcSet_un.emptyClone();//new  RooDataSet(mcSet_un.GetName(),mcSet_un.GetTitle(),*wspace.set("treeSet"),RooFit::WeightVar("weight"));

	std::cout << c->GetEntries() << std::endl;
	//		RooFormulaVar wFunc("Weight","Weight","LifetimeReweight(mass,ctau,x[0])",*wspace.var("hnl_ct")) ;
	//include weights
	mcSet_un.Print();
	TH1D* h_weights=new TH1D("w","w",100,0,300);
	TH1D* h_nn=new TH1D("nn_score","nn_score",100,0.95,1);
	TH1D* h_ct=new TH1D("hnl_ct","hnl_ct",100,0,ctau*10);
	TH1D* hist_pi_pt=new TH1D("hnl_pi_pt","hnl_pi_pt",100,0,10);
	TH1D* hist_cos2D=new TH1D("hnl_cos2D","hnl_cos2D",100,0.995,1);
	TH1D* hist_vtxProb=new TH1D("hnl_vtxProb","hnl_vtxProb",100,0,0.01);
	TH1D* hist_l_pt=new TH1D("hnl_l_pt","hnl_l_pt",100,0,50);
	TH1D* hist_trgmu_pt=new TH1D("hnl_trgmu_pt","hnl_trgmu_pt",100,0,80);
	TH1D* hist_B_mass=new TH1D("B_mass","B_mass",100,0,8);
	TH1D* hist_B_pt=new TH1D("B_pt","B_pt",100,0,100);
	TH1D* hist_dilepton_mass=new TH1D("dilepton_mass","dilepton_mass",100,0,8);
	TH1D* hist_BlepPi_mass=new TH1D("BlepPi_mass","BlepPi_mass",100,0,8);
	Nevents = mcSet_un.sumEntries();
	for (int i =0; i<mcSet_un.sumEntries();i++){
		RooArgSet* entry;
		RooRealVar* target = (RooRealVar* ) mcSet_un.get(i)->find("hnl_ct");
		//	RooRealVar* hnlM = (RooRealVar* ) mcSet_un.get(i)->find("hnl_mass");
		RooRealVar* puWeight = (RooRealVar* ) mcSet_un.get(i)->find("PU_weight");
		RooRealVar* trgWeight = (RooRealVar* ) mcSet_un.get(i)->find("Trg_weight");
		RooRealVar* muIdWeight = (RooRealVar* ) mcSet_un.get(i)->find("muId_weight");
		RooRealVar* hnlcos2D = (RooRealVar* ) mcSet_un.get(i)->find("hnl_cos2D");
		RooRealVar* hnlvtxProb = (RooRealVar* ) mcSet_un.get(i)->find("hnl_vtxProb");
		RooRealVar* hnl_piPt = (RooRealVar* ) mcSet_un.get(i)->find("hnl_pi_pt");
		RooRealVar* hnl_pi_DCAS = (RooRealVar* ) mcSet_un.get(i)->find("hnl_pi_DCAS");
		RooRealVar* bMass = (RooRealVar* ) mcSet_un.get(i)->find("B_mass");
		RooRealVar* hnlchi2 = (RooRealVar* ) mcSet_un.get(i)->find("hnl_vtxChi2");
		RooRealVar* l0_Pt = (RooRealVar* ) mcSet_un.get(i)->find("l0_pt");
		RooRealVar* hnl_lPt = (RooRealVar* ) mcSet_un.get(i)->find("hnl_l_pt");
		RooRealVar* hnltrgmuPt = (RooRealVar* ) mcSet_un.get(i)->find("TrgMu_pt");
		RooRealVar* l_Pt = (RooRealVar* ) mcSet_un.get(i)->find("l_pt");
		RooRealVar* hnl_lxySig = (RooRealVar* ) mcSet_un.get(i)->find("hnl_lxy_sig");
		RooRealVar* B_pt = (RooRealVar* ) mcSet_un.get(i)->find("B_pt");
		RooRealVar* dilep_mass = (RooRealVar* ) mcSet_un.get(i)->find("dilepton_mass");
		RooRealVar* Blep_pi_mass = (RooRealVar* ) mcSet_un.get(i)->find("BlepPi_mass");
		RooRealVar* dr_Blep_pi = (RooRealVar* ) mcSet_un.get(i)->find("dr_Blep_pi");
		RooRealVar* dr_trgMu_lep = (RooRealVar* ) mcSet_un.get(i)->find("dr_trgMu_lep");
		RooRealVar* l0_iso = (RooRealVar* ) mcSet_un.get(i)->find("TrgMu_relIso");
		RooRealVar* l_iso = (RooRealVar* ) mcSet_un.get(i)->find("hnl_l_relIso");
		RooRealVar* pi_pixelLayers = (RooRealVar* ) mcSet_un.get(i)->find("'hnl_pi_PixelLayers");
		RooRealVar* pi_TrkLayers = (RooRealVar* ) mcSet_un.get(i)->find("hnl_pi_TrackerLayers");
		RooRealVar* trgmu_pxLayers = (RooRealVar* ) mcSet_un.get(i)->find("'TrgMu_PixelLayers");
		RooRealVar* trgmu_TrkLayers = (RooRealVar* ) mcSet_un.get(i)->find("TrgMu_TrackerLayers");

		hnl_mass->setVal(mcSet_un.get(i)->getRealValue("hnl_mass"));
		h_ct->Fill(target->getVal(),weight->getVal());
		if (reweight){
			//		std::cout << target->getVal() << " " << LifetimeReweight(mass,ctau,target->getVal())<< std::endl;
			//	std::cout << target->getVal() << " " << LifetimeReweight(mass,ctau,target->getVal())<< std::endl;i
	//	std::cout << "flavoooorrr____________________" << flavor << std::endl;
			double w = LifetimeReweight(csvIn,mass,ctau,target->getVal(),flavor)*puWeight->getVal()*trgWeight->getVal()*muIdWeight->getVal(); 
			if (w<20){
			weight->setVal(w);
			h_weights->Fill(weight->getVal());
			}
			else weight->setVal(0);
			h_ct->Fill(target->getVal(),weight->getVal());
		}
		else{
			weight->setVal(puWeight->getVal()*trgWeight->getVal()*muIdWeight->getVal());
			//	std::cout << weight->getVal() << std::endl;
		}
		//	std::cout << weight->getVal() <<std::endl;
		double var;
		if (doNNsel){
/*		if( Bc )	var = TPython::Eval(("ComputeNNScore("+
						std::to_string(hnl_piPt->getVal())+","+
						std::to_string(hnl_lPt->getVal())+","+
						std::to_string(hnltrgmuPt->getVal())+","+
						std::to_string(bMass->getVal())+","+
						std::to_string(hnlcos2D->getVal())+","+
						std::to_string(hnl_lxySig->getVal())+","+
						std::to_string(hnlvtxProb->getVal())+","+
						std::to_string(hnlchi2->getVal())+","+
						std::to_string(B_pt->getVal())+","+
						std::to_string(dilep_mass->getVal())+","+
						std::to_string(Blep_pi_mass->getVal())+","+
						std::to_string(dr_trgMu_lep->getVal())+","+
						std::to_string(dr_Blep_pi->getVal())+","+
						std::to_string(l0_iso->getVal())+","+
						std::to_string(l_iso->getVal())+","+
						std::to_string(mass)+")").c_str()
					);

		else var = TPython::Eval(("ComputeNNScore_10("+
						std::to_string(hnlcos2D->getVal())+","+
						std::to_string(hnlvtxProb->getVal())+","+
						std::to_string(hnl_piPt->getVal())+","+
						std::to_string(bMass->getVal())+","+
						std::to_string(hnlchi2->getVal())+","+
						std::to_string(hnl_lPt->getVal())+","+
						std::to_string(hnltrgmuPt->getVal())+","+
						std::to_string(hnl_lxySig->getVal())+","+
						std::to_string(B_pt->getVal())+","+
						std::to_string(mass)+")").c_str()
					);*/
		var = TPython::Eval(("ComputeNNScore_hits("+
						std::to_string(hnl_piPt->getVal())+","+
						std::to_string(l0_Pt->getVal())+","+
						std::to_string(l_Pt->getVal())+","+
						std::to_string(bMass->getVal())+","+
						std::to_string(hnlcos2D->getVal())+","+
						std::to_string(hnl_lxySig->getVal())+","+
						std::to_string(hnl_pi_DCAS->getVal())+","+
						std::to_string(hnlvtxProb->getVal())+","+
						std::to_string(B_pt->getVal())+","+
						std::to_string(dilep_mass->getVal())+","+
						std::to_string(Blep_pi_mass->getVal())+","+
						std::to_string(dr_trgMu_lep->getVal())+","+
						std::to_string(dr_Blep_pi->getVal())+","+
						std::to_string(l0_iso->getVal())+","+
						std::to_string(l_iso->getVal())+","+
						std::to_string(mcSet_un.get(i)->getRealValue("hnl_pi_PixelLayers"))+","+
						std::to_string(mcSet_un.get(i)->getRealValue("hnl_pi_TrackerLayers"))+","+
						std::to_string(mcSet_un.get(i)->getRealValue("TrgMu_PixelLayers"))+","+
						std::to_string(mcSet_un.get(i)->getRealValue("TrgMu_TrackerLayers"))+","+
						std::to_string(mass)+")").c_str()
					);

					
//`		std::cout << var << std::endl;
			h_nn->Fill(var);
			if (Float_t(var)>nn_cut){	
				nn_score->setVal(Float_t(var));

				hist_pi_pt->Fill(hnl_piPt->getVal());
				hist_l_pt->Fill(l_Pt->getVal());
				hist_trgmu_pt->Fill(l0_Pt->getVal());
				hist_B_mass->Fill(bMass->getVal());
				hist_B_pt->Fill(B_pt->getVal());
				hist_cos2D->Fill(hnlcos2D->getVal());
				hist_vtxProb->Fill(hnlvtxProb->getVal());
				hist_dilepton_mass->Fill(dilep_mass->getVal());
				hist_BlepPi_mass->Fill(Blep_pi_mass->getVal());
		//		hist_bMass->Fill(bMass->getVal());
					//entry->Print("V");
				//		entry->add(*wspace.var("hnl_mass"));
				//		entry->add(*weight);
				//		entry->add(*nn_score);
				temp->add(RooArgList(*hnl_mass,*weight),weight->getVal());
			}
		}else{

			nn_score->setVal(-1.0);
			//	entry->add(*wspace.var("hnl_mass"));
			//	entry->add(*weight);
			//entry->add(*nn_score);
			temp->add(RooArgList(*hnl_mass,*weight),weight->getVal());

		}
	}		
	h_weights->SaveAs((wd+"/hweights"+category+".root").c_str());
	h_ct->SaveAs((wd+"/h_ct_weighted_"+category+".root").c_str());
	if (doNNsel)h_nn->SaveAs((wd+"/h_nn"+category+".root").c_str());	

	RooDataSet mcSet(*temp,"mcSet");
	mcSet.Print();
/*	TFile f("tree.root", "update");
	hist_pi_pt->Write();
	hist_l_pt->Write();
	hist_trgmu_pt->Write();
	hist_B_mass->Write();
	hist_B_pt->Write();
	hist_cos2D->Write();
	hist_vtxProb->Write();
	hist_dilepton_mass->Write();
	hist_BlepPi_mass->Write();
	f.Close();*/
	
	
	std::cout << "is dataset weighted: " << mcSet.isWeighted() <<std::endl;
	std:: cout << (lxy+sig+sel).c_str() << " " << mcSet.sumEntries()<< " " << mcSet.createHistogram("counter",*wspace.var("hnl_mass"))->GetEntries() <<std::endl;
	if (mcSet.sumEntries()<1e-10){
		sigF[0]	= mcSet.sumEntries();
		sigF[1]	= mcSet.sumEntries();
		sigF[2]	= -1;
		sigF[3]	= -1;
		std::ofstream SigFeatures;
		std:: cout << wd+"/Signal"+flav+"_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt" <<std::endl;
		SigFeatures.open(wd+"/Signal"+flav+"_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt",ios_base::app );
		SigFeatures << mass << " " << ctau << "" << " " << sigF[0] <<  " " << sigF[1] << " " << sigF[2] << " " << mcSet.sumEntries() << " " << Ntot <<  std::endl;
		return 0;


	}
	wspace.factory(("Nsig[10,"+std::to_string(mcSet.sumEntries()*0.8)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());

	wspace.Print();


	wspace.var("hnl_mass")->setMin(mass-0.3);
	wspace.var("hnl_mass")->setMax(mass+0.2);
	wspace.factory(("mean_ele_"+var_label+"["+std::to_string(mcSet.mean(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-2*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+2*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
	wspace.factory(("mean_ele_p_"+var_label+"["+std::to_string(mcSet.mean(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-2*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+2*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());

	wspace.factory(("sigma_ele_"+var_label+"["+std::to_string(sigma*0.3)+","+std::to_string(sigma*0.2)+","+std::to_string(sigma*2.0)+"]").c_str());
	wspace.factory(("sigma_ele_"+var_label).c_str())->Print();
	//wspace.factory(("sigma_ele_"+var_label+"["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.3)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.2)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*2.0)+"]").c_str());
	wspace.factory(("sigma_ele_p_"+var_label+"["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.3)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.2)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*2.0)+"]").c_str());
	wspace.factory("sigma_mu[0,0,0]");
	wspace.factory(("n1_cb_"+var_label+"[0.5,0,10]").c_str());
	wspace.factory(("n1_cb_p_"+var_label+"[0.5,0,10]").c_str());
	wspace.factory(("alpha1_cb_"+var_label+"[0.5,0.1,10]").c_str());
	wspace.factory(("alpha1_cb_p_"+var_label+"[0.5,0.1,10]").c_str());
	wspace.factory(("n2_cb_"+var_label+"[2,1,6]").c_str());
	wspace.factory(("n2_cb_p_"+var_label+"[2,1,6]").c_str());
	wspace.factory(("alpha2_cb_"+var_label+"[1,0,10]").c_str());
	wspace.factory(("alpha2_cb_p_"+var_label+"[1,0,10]").c_str());
	wspace.factory("n_cb[0.5,0,10]");
	wspace.factory("alpha_cb[10,-100,1000]");
	wspace.factory("N1[0.25,0,0.5]");
	wspace.factory("N2[0.25,0,0.5]");
	wspace.factory(("Nele["+std::to_string(mcSet.sumEntries()*1)+","+std::to_string(mcSet.sumEntries()*0.5)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());
	wspace.factory(("Nmu[0,0,0]"));
	wspace.factory(("width[0.1,0,1]"));
	//		wspace.factory("Gaussian::signal(hnl_mass,mean_ele,sigma_ele)");
	//		wspace.factory("Voigtian::signal(hnl_mass,mean_ele,sigma_ele,width)");
//	wspace.factory("CBShape::CB1"+category+"_"+B_c+"(hnl_mass,mean_ele,sigma_ele,alpha1_cb,n1_cb)");//	wspace.factory("CBShape::CB2"+category+"_"+B_c+"(hnl_mass,mean_ele,sigma_ele,alpha2_cb,n2_cb)");
	//	wspace.factory("SUM::signal(Nele*CB1,Nmu*CB2)")

	RooDoubleCBFast signal_nominal("signal_nominal","DoubleSidedCB",*wspace.var("hnl_mass"),*wspace.var(("mean_ele_"+var_label).c_str()),*wspace.var(("sigma_ele_"+var_label).c_str()),*wspace.var(("alpha1_cb_"+var_label).c_str()),*wspace.var(("n1_cb_"+var_label).c_str()),*wspace.var(("alpha2_cb_"+var_label).c_str()),*wspace.var(("n2_cb_"+var_label).c_str()));	
	wspace.var(("mean_ele_p_"+var_label).c_str())->setVal(mass);
	wspace.var(("sigma_ele_p_"+var_label).c_str())->setVal(sigma);
	wspace.var(("alpha1_cb_p_"+var_label).c_str())->setVal(1.2);
	wspace.var(("n1_cb_p_"+var_label).c_str())->setVal(2.5);
	wspace.var(("alpha2_cb_p_"+var_label).c_str())->setVal(1.5);
	wspace.var(("n2_cb_p_"+var_label).c_str())->setVal(4);
	RooDoubleCBFast signal("signal","DoubleSidedCB",*wspace.var("hnl_mass"),*wspace.var(("mean_ele_p_"+var_label).c_str()),*wspace.var(("sigma_ele_p_"+var_label).c_str()),*wspace.var(("alpha1_cb_p_"+var_label).c_str()),*wspace.var(("n1_cb_p_"+var_label).c_str()),*wspace.var(("alpha2_cb_p_"+var_label).c_str()),*wspace.var(("n2_cb_p_"+var_label).c_str()));	
	//parametrized shape for signal
	//      signal fit 	
	signal_nominal.fitTo(mcSet,RooFit::Save(),RooFit::Range(mass-3*mcSet.sigma(*wspace.var("hnl_mass")),mass+3*mcSet.sigma(*wspace.var("hnl_mass"))),RooFit::SumW2Error(kTRUE));
	
  //	wspace.var(("n1_cb_"+category).c_str())->setRange(wspace.var(("n1_cb_"+category).c_str())->getVal()*0.7,wspace.var(("n1_cb_"+category).c_str())->getVal()*1.3);
 //	wspace.var(("n2_cb_"+category).c_str())->setRange(wspace.var(("n2_cb_"+category).c_str())->getVal()*0.7,wspace.var(("n2_cb_"+category).c_str())->getVal()*1.3);
 // 	wspace.var(("alpha1_cb_"+category).c_str())->setRange(wspace.var(("alpha1_cb_"+category).c_str())->getVal()*0.7,wspace.var(("alpha1_cb_"+category).c_str())->getVal()*1.3);
 // 	wspace.var(("alpha2_cb_"+category).c_str())->setRange(wspace.var(("alpha2_cb_"+category).c_str())->getVal()*0.7,wspace.var(("alpha2_cb_"+category).c_str())->getVal()*1.3);

  	//RooFitResult *res = signal.fitTo(mcSet,RooFit::Save(),RooFit::Range(mass-3*mcSet.sigma(*wspace.var("hnl_mass")),mass+3*mcSet.sigma(*wspace.var("hnl_mass"))),RooFit::SumW2Error(kTRUE));

	//resize dataset to background fit window (nSigma), for now is 10 
	wspace.var("hnl_mass")->setMin(wspace.var(("mean_ele_"+var_label).c_str())->getVal()-nSigma*wspace.var(("sigma_ele_"+var_label).c_str())->getVal()) ;
	wspace.var("hnl_mass")->setMax(wspace.var(("mean_ele_"+var_label).c_str())->getVal()+nSigma*wspace.var(("sigma_ele_"+var_label).c_str())->getVal()) ;


	//Signal only plots
	TCanvas* c1 =new TCanvas("MC_Fit", "fit", 800, 800);
	RooPlot* mhnl_frame = wspace.var("hnl_mass")->frame();
	wspace.var("hnl_mass")->setBins(50,"plotter");
	mcSet.plotOn(mhnl_frame, RooFit::Name("template mc"),RooFit::Binning("plotter"));
	signal_nominal.plotOn(mhnl_frame,RooFit::Name("signal"),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
	signal.plotOn(mhnl_frame,RooFit::Name("signal_param"),RooFit::LineColor(3),RooFit::MoveToBack()); // this will show fit overlay on canvas
	mhnl_frame->GetYaxis()->SetTitleOffset(0.9);
	mhnl_frame->GetYaxis()->SetTitleFont(42);
	mhnl_frame->GetYaxis()->SetTitleSize(0.05);
	mhnl_frame->GetYaxis()->SetLabelSize(0.065);
	mhnl_frame->GetYaxis()->SetLabelSize(0.04);
	mhnl_frame->GetYaxis()->SetLabelFont(42);
	mhnl_frame->GetXaxis()->SetTitleOffset(0.9);
	mhnl_frame->GetXaxis()->SetTitleFont(42);
	mhnl_frame->GetXaxis()->SetTitleSize(0.05);
	mhnl_frame->GetXaxis()->SetLabelSize(0.065);
	mhnl_frame->GetXaxis()->SetLabelSize(0.04);
	mhnl_frame->GetXaxis()->SetLabelFont(42);

	mhnl_frame->GetYaxis()->SetTitle("Events");
	mhnl_frame->GetXaxis()->SetTitle("m_{HNL} (GeV)");
	mhnl_frame->SetStats(0);
	mhnl_frame->SetMinimum(0);
	mhnl_frame->Draw();
	Double_t chi2 = mhnl_frame->chiSquare(5);
	std::cout << chi2 << std::endl;
	//	TMathText* latex = new TMathText(0.5,0.8,("#Chi^{2}/n.dof="+std::string(Form("%.2f", chi2))).c_str()) ;
	TLatex* latex= new TLatex();
	latex->SetTextAlign(12);
	latex->SetTextFont(43);
	latex->SetTextSize(30);
	/*	TPaveText* t = new TPaveText(0.1,0.8,0.3,0.62,"NB");
		t->AddText("ciao");
		t->AddText(Form("#chi^{2} = %.2f", chi2));
		t->SetTextAlign(22);
		t->SetTextFont(43);
		t->SetTextSize(30);*/
	//latex->Draw();
	mhnl_frame->addObject(latex);
//	std::cout << mass-0.28 << " " << mhnl_frame->GetMaximum()<<std::endl;
	latex->DrawLatex(mass-9*wspace.var(("sigma_ele_"+var_label).c_str())->getVal(),mhnl_frame->GetMaximum()*0.76,("#chi^{2}/n.dof"+std::string(Form(" %.2f", chi2))).c_str());
	latex->DrawLatex(mass-9*wspace.var(("sigma_ele_"+var_label).c_str())->getVal(),mhnl_frame->GetMaximum()*0.97,("HNL Mass"+std::string(Form(" %.2f", mass))+" GeV c#tau = "+std::string(Form(" %.2f", ctau))+ " mm").c_str());
	latex->DrawLatex(mass-9*wspace.var(("sigma_ele_"+var_label).c_str())->getVal(),mhnl_frame->GetMaximum()*0.90,("mean"+std::string(Form(" %.2f", wspace.var(("mean_ele_"+var_label).c_str())->getVal()))+" GeV " ).c_str());
	latex->DrawLatex(mass-9*wspace.var(("sigma_ele_"+var_label).c_str())->getVal(),mhnl_frame->GetMaximum()*0.83,(" #sigma = "+std::string(Form(" %.3f", wspace.var(("sigma_ele_"+var_label).c_str())->getVal()))+" GeV "  ).c_str());
	c1->SaveAs((wd+"/SigFits/Fit_mass"+flav+"_"+std::to_string(mass)+"_ctau"+std::to_string(ctau)+"_"+var_label+".pdf").c_str());


	RooDataHist* sig_hist = mcSet.binnedClone("sig_hist","sig_hist");

	sigF[0]	= mcSet.sumEntries()*FilterEff/Ntot;
	sigF[1]	= sqrt(mcSet.sumEntries());
	sigF[2]	= wspace.var(("sigma_ele_"+var_label).c_str())->getVal();
	sigF[3]	= -1;



	wspace.var("hnl_mass")->setRange("threeSigma",wspace.var(("mean_ele_"+var_label).c_str())->getVal()-3*sigF[2],wspace.var(("mean_ele_"+var_label).c_str())->getVal()+3*sigF[2]) ;

	RooAbsReal* yield = signal.createIntegral(*wspace.var("hnl_mass"));
	std::cout << mcSet.sumEntries()*FilterEff/Ntot  << " " <<mcSet.sumEntries()<< " "<< FilterEff << " "<<Ntot  <<std::endl;

	//prepare ws for signal shapes
	RooWorkspace w("w");
	//RooRealVar signal_norm("signal_norm","signal_norm", 1.0 , 1.0 , 1.0 );	
	w.import(signal_nominal);
	w.import(signal);
	w.import(mcSet);
	//	w.import(signal_norm);
	w.var(("mean_ele_p_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("sigma_ele_p_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("alpha1_cb_p_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("n1_cb_p_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("alpha2_cb_p_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("n2_cb_p_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("mean_ele_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("sigma_ele_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("alpha1_cb_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("n1_cb_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("alpha2_cb_"+var_label).c_str())->setConstant(kTRUE);
	w.var(("n2_cb_"+var_label).c_str())->setConstant(kTRUE);
	w.var("hnl_mass")->setMin(wspace.var(("mean_ele_"+var_label).c_str())->getVal()-10*wspace.var(("sigma_ele_"+var_label).c_str())->getVal()) ;
	w.var("hnl_mass")->setMax(wspace.var(("mean_ele_"+var_label).c_str())->getVal()+10*wspace.var(("sigma_ele_"+var_label).c_str())->getVal()) ;

	//check and save ws
	w.Print();
	char mass_lable[100];
	char ctau_lable[200];
//	std::cout << "check mass lable for mass " << mass << " " <<int(mass) << <<  std::endl;
	sprintf(mass_lable,"%.2f",mass);
	sprintf(ctau_lable,"%.4f",ctau);
//	std::cout <<  mass_lable <<std::endl;
	std::string ctau_label = std::string(ctau_lable).replace(std::string(ctau_lable).find("\."),1,"p");
	std::string mass_label = std::string(mass_lable).replace(std::string(mass_lable).find("\."),1,"p");
//	std::replace_if((std::string(mass_label)).begin(), (std::string(mass_label)).end(), std::ispunct,'p');
	std::cout << wd+"/SigFits/workspace_signal"+flav+"_bhnl_m_"+std::string(mass_label)+"_ctau_"+std::string(ctau_label)+"_cat_"+category+".root" <<std::endl;
	w.writeToFile((wd+"/SigFits/workspace_signal"+flav+"_bhnl_m_"+std::string(mass_label)+"_ctau_"+std::string(ctau_label)+"_cat_"+category+".root").c_str());

	//save signal yields and sigmas to file	
	std::ofstream SigFeatures;
	std:: cout << wd+"/Signal"+flav+"_Mass"+std::string(mass_label)+"_ctau"+std::string(ctau_lable)+"_"+category+"_features.txt" <<std::endl;
	SigFeatures.open(wd+"/Signal"+flav+"_Mass"+std::string(mass_label)+"_ctau"+std::string(ctau_lable)+"_"+category+"_features.txt",ios_base::app);
	std::cout << SigFeatures.is_open() <<std::endl;
	SigFeatures << mass << " " << ctau << "" << " " << sigF[0] <<  " " << sigF[1] << " " << sigF[2] << " " << mcSet.sumEntries()  << " " << Nevents << std::endl;
	SigFeatures.close();
	


}


//needs a bit of cleaning
int ShapeWSProducer(std::string filename, std::string path,double mass,double ctau, std::string wd, std::string sel,double lumiRatio,bool binned, bool blinded,float nn_cut){

	//this has to be substituted with input from .csv with all signal info

	std::string lxy[3]={"hnl_lxy_sig<50 && ","hnl_lxy_sig>50 && hnl_lxy_sig<150 && ","hnl_lxy_sig>150 && "};
	std::string lxyLables[3]={"lxysig0to50","lxysig50to150","lxysiggt150"};
	std::string sig[2]={"LepQProd<0","LepQProd>0"};
	std::string sig_lables[2] = {"OS","SS"};


	double sigF[3];

	for (int lxy_i =0; lxy_i < 3;lxy_i++){

		for (int s_j =0;s_j<2;s_j++){


			single_fit(filename,lxyLables[lxy_i],sig_lables[s_j], 1,mass, ctau, lxy[lxy_i],sig[s_j],sel,1,nn_cut,sigF,wd,10);

		}
	}
	return 0;
}



/*int AllFits(std::string filename,std::string path, double mass,double ctau,int channel, std::string wd,std::string sel,double lumiRatio, bool binned, bool blinded ){


  int  j;
  double sigFeat[3];
//std::string lxy_cuts[4]= {"hnl_lxy<3", "hnl_lxy>3 && hnl_lxy<10", "hnl_lxy>10 && hnl_lxy<20","hnl_lxy>20"};
std::string lxy_cuts[3]= {"hnl_lxy<1", "hnl_lxy>1 && hnl_lxy<5", "hnl_lxy>5"};
std::string category[3]= {"LxyUnder1", "LxyOver1Under5", "LxyOver5"};

std::ofstream SigFeatures;

SigFeatures.open(wd+"/Signal_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt");
std::cout << wd+"/Signal_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt" << std::endl;


std::cout << "in cycle " << std::endl;
for (j=0;j<3;j++){
const int temp_i = channel;
const int temp_j=j;
fit(filename.c_str(),path.c_str(),mass,ctau,channel,lxy_cuts[j],j,category[j],wd+"/SigFits",sigFeat,sel,lumiRatio,binned,blinded,0.9);
channel=temp_i;
j=temp_j;
SigFeatures << channel << " "  << j << " " << sigFeat[0] <<  " " << sigFeat[1] << " " << sigFeat[2] <<  std::endl;
}

SigFeatures.close();

return 0;
}
*/

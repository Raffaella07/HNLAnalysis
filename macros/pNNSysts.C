#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TH1.h"
#include "setStyle.C"
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

TH1D* testModel(std::string filename,std::string lxyLable,std::string sigLable, bool doNNsel,double mass,double ctau,std::string lxy,std::string sig,std::string sel, bool debug, double nn_cut,std::string wd,bool ele_model, TH1D** h_mass,bool data){

	
	std::vector<double> GenEvtCount,Filter,Mass,Ctau,Sigma,Flav_ratio;
	std::vector<std::string> Files, label;

	
	//sig_fromcsv("../data/MC_datasets_InclusiveFilterEff.csv",&Files,&label,&GenEvtCount,&Filter,&Mass,&Ctau,&Sigma,&Flav_ratio);
	sig_fromcsv("../data/MC_datasets_InclusiveFilterEff.csv",&Files,&label,&GenEvtCount,&Filter,&Mass,&Ctau,&Sigma,&Flav_ratio);
	

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

	bool no_reweight = filename.find(".root") != std::string::npos ;
	std::cout  << filename << "is central sample? " << no_reweight << std::endl;

	if (no_reweight){
		c->Add((filename).c_str());
		//	std::cout << c->GetEntries() << std::endl;
		for(int i=0;i<Files.size();i++){
			if (mass==Mass[i] && ctau ==Ctau[i] ){
						Ntot=GenEvtCount[i] * Flav_ratio[i];
						FilterEff = Filter[i];
					}
				}
			
		
	}else {

		for(int i=0;i<Files.size();i++){

			if (mass==Mass[i] && Ctau[i]>0){

						std::cout << "mass " << mass << std::endl;
						std::cout << Ctau[i] << " " << ctau << std::endl;

						c->Add((Files[i]).c_str());
						Ntot += GenEvtCount[i] * Flav_ratio[i] ;
						FilterEff += GenEvtCount[i]*Filter[i] * Flav_ratio[i];

					}


				}


			
		
		FilterEff = FilterEff/Ntot;
		std::cout << FilterEff <<" "  <<Ntot << std::endl;
	}

	std::string modelPath;
	std::string category;
	if (doNNsel){

		category = lxyLable+"_"+sigLable;
		//	std::string modelPath = "/cmshome/ratramon/Analysis/python/model/";
		if (ele_model)modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_12Oct2022_15h25m05s/";
		else modelPath = "/cmshome/ratramon/Analysis/python/model/mu_model/trainNN_"+category+"_mu/";
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
	wspace.factory("LepQProd[0,-3,3]");
	wspace.factory("hnl_charge[0,-2,2]");
	wspace.factory("Type[0,0,4]");
	wspace.factory("LxyBin[0,0,4]");
	wspace.factory("dr_trgmu_hnl[0,-100,100]");
	wspace.factory("dilepton_mass[0,0,1000000]");
	wspace.factory("dilepton_pt[0,0,10000000]");
	wspace.factory("hnl_vtxProb[0,0,1]");
	wspace.factory("hnl_vtxChi2[0,0,1000000]");
	wspace.factory("hnl_pi_DCAS[0,-100000000,10000000]");
	wspace.factory("hnl_l_dxyS[0,-100000000,100000000]");
	wspace.factory("hnl_l_pt[10,0,100000000]");
	wspace.factory("hnl_l_mvaId[0,-2000,2000]");
	wspace.factory("hnl_pi_pt[10,0,10000000]");
	wspace.factory("TrgMu_pt[10,0,10000000]");
	wspace.factory("hnl_lxy_sig[10,-10000000,10000000]");
	wspace.factory("hnl_to_trgmu[10,-10000000,10000000]");
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
	wspace.var("r_ele");
	wspace.var("r_mu");
	wspace.var("Nsig");
	wspace.var("Nele");
	wspace.var("Nmu");
	wspace.var("N1");
	wspace.var("N2");
	wspace.var("width");
	wspace.var("sigma_ele");
	wspace.var("sigma_mu");
	wspace.var("n_cb");
	wspace.var("alpha_cb");
	wspace.var("alpha");
	wspace.var("alpha1");
	wspace.var("alpha2");
	wspace.var("hnl_mass");
	wspace.var("hnl_cos2D");
	wspace.var("hnl_vtxProb");
	wspace.var("hnl_vtxChi2");
	wspace.var("hnl_charge");
	wspace.var("hnl_lxy_sig");
	wspace.var("hnl_to_trgmu");
	wspace.var("hnl_lxy");
	wspace.var("hnl_ct");
	wspace.var("hnl_Blep_hnl");
	wspace.var("dilepton_mass");
	wspace.var("dilepton_pt");
	wspace.var("hnl_l_dxyS");
	wspace.var("hnl_pi_DCAS");
	wspace.var("hnl_l_pt");
	wspace.var("hnl_pi_pt");
	wspace.var("TrgMu_pt");
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
	wspace.var("mean");
	wspace.var("sigma");
	wspace.var("n1_cb");
	wspace.var("alpha1_cb");
	wspace.var("Nsig");
	wspace.var("Nbkg");
	wspace.var("common");

	wspace.var("c1");
	wspace.var("c2");

	//Define set of variables to include in dataset

	//	std::string variables = "B_mass,hnl_mass,hnl_charge,hnl_lxy,Type,LxyBin,dr_trgmu_hnl,dilepton_mass,hnl_vtxProb,hnl_cos2D,hnl_l_pt,hnl_l_mvaId,hnl_pi_pt,hnl_l_dxyS,hnl_pi_dxyS,dilepton_pt,likelihood,hnl_lxy_sig,LepQProd,hnl_ct,nn_score,nn_score_single";
	std::string variables = "hnl_mass,hnl_charge,likelihood,hnl_ct,hnl_lxy,LepQProd,hnl_cos2D,hnl_vtxProb,hnl_pi_pt,B_mass,hnl_vtxChi2,hnl_l_pt,TrgMu_pt,hnl_lxy_sig,B_pt,dilepton_mass,PU_weight,Trg_weight,muId_weight,weight,hnl_to_trgmu";
	wspace.defineSet("treeSet",variables.c_str());
	//	MC dataset initialization
	RooDataSet mcSet_un("mcSet_un","mcSet_un",c,*wspace.set("treeSet"),(lxy+sig+sel).c_str())/*("Type=="+std::to_string(type)+" && "+lxy +sel).c_str()))*/;
	RooArgSet * entry;
	RooRealVar* nn_score= new RooRealVar("nn_score","nn_score",0.5,0,1000000);
	RooRealVar* weight= new RooRealVar("weight","weight",1,0,1000000000);
	RooRealVar* hnl_mass= new RooRealVar("hnl_mass","hnl_mass",1,0,1000000000);
	RooRealVar* B_mass= new RooRealVar("B_mass","B_mass",5.3,0,1000000);
	RooDataSet* temp = new RooDataSet("temp","",RooArgSet(*B_mass,*weight),weight->GetName()); 
	//	temp= (RooDataSet*) mcSet_un.emptyClone();//new  RooDataSet(mcSet_un.GetName(),mcSet_un.GetTitle(),*wspace.set("treeSet"),RooFit::WeightVar("weight"));

	std::cout << c->GetEntries() << std::endl;
	//		RooFormulaVar wFunc("Weight","Weight","LifetimeReweight(mass,ctau,x[0])",*wspace.var("hnl_ct")) ;
	//include weights
	mcSet_un.Print();
	TH1D* h_weights=new TH1D("w","w",100,0,20);
	TH1D* h_nn=new TH1D("nn_score","nn_score",20,0,1);
	TH1D* h_ct=new TH1D("hnl_ct","hnl_ct",100,0,ctau*100);
//	TH1D* h_hnl =(TH1D*) mcSet_un.createHistogram("hnl_mass",*wspace.var("hnl_mass"));
//	h_hnl->SaveAs((wd+"/h_hnl"+category+".root").c_str());
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
		RooRealVar* bMass = (RooRealVar* ) mcSet_un.get(i)->find("B_mass");
		RooRealVar* hnlchi2 = (RooRealVar* ) mcSet_un.get(i)->find("hnl_vtxChi2");
		RooRealVar* hnl_lPt = (RooRealVar* ) mcSet_un.get(i)->find("hnl_l_pt");
		RooRealVar* hnltrgmuPt = (RooRealVar* ) mcSet_un.get(i)->find("TrgMu_pt");
		RooRealVar* hnl_lxySig = (RooRealVar* ) mcSet_un.get(i)->find("hnl_lxy_sig");
		RooRealVar* hnlpiDCAS = (RooRealVar* ) mcSet_un.get(i)->find("B_pt");
		RooRealVar* hnl_toTrgMu = (RooRealVar* ) mcSet_un.get(i)->find("hnl_to_trgmu");
	
		RooRealVar* lep0_pt; 
		RooRealVar* lep_pt;

		if ( hnl_toTrgMu->getVal()==0){ 
		lep0_pt = (RooRealVar* ) mcSet_un.get(i)->find("TrgMu_pt");
		lep_pt = (RooRealVar* ) mcSet_un.get(i)->find("hnl_l_pt");
		}
		else if ( hnl_toTrgMu->getVal()==1){
		lep_pt = (RooRealVar* ) mcSet_un.get(i)->find("TrgMu_pt");
		lep0_pt = (RooRealVar* ) mcSet_un.get(i)->find("hnl_l_pt");
		}
		

		hnl_mass->setVal(mcSet_un.get(i)->getRealValue("dilepton_mass"));
		if (!no_reweight){
			//		std::cout << target->getVal() << " " << LifetimeReweight(mass,ctau,target->getVal())<< std::endl;
			//	std::cout << target->getVal() << " " << LifetimeReweight(mass,ctau,target->getVal())<< std::endl;
			weight->setVal(LifetimeReweight(mass,ctau,target->getVal())*puWeight->getVal()*trgWeight->getVal()*muIdWeight->getVal());
			h_weights->Fill(weight->getVal());
			h_ct->Fill(target->getVal(),weight->getVal());
		}
		else{
			weight->setVal(puWeight->getVal()/**trgWeight->getVal()*/*muIdWeight->getVal());
			//	std::cout << weight->getVal() << std::endl;
		}
		//	std::cout << weight->getVal() <<std::endl;

		if (doNNsel){
		double var;
		if (ele_model)	var = TPython::Eval(("ComputeNNScore_10("+
						std::to_string(hnlcos2D->getVal())+","+
						std::to_string(hnlvtxProb->getVal())+","+
						std::to_string(hnl_piPt->getVal())+","+
						std::to_string(bMass->getVal())+","+
						std::to_string(hnlchi2->getVal())+","+
						std::to_string(hnl_lPt->getVal())+","+
						std::to_string(hnltrgmuPt->getVal())+","+
						std::to_string(hnl_lxySig->getVal())+","+
						std::to_string(hnlpiDCAS->getVal())+","+
						std::to_string(mass)+")").c_str()
					);
		else	var = TPython::Eval(("ComputeNNScore_10("+
						std::to_string(hnl_piPt->getVal())+","+
						std::to_string(lep_pt->getVal())+","+
						std::to_string(lep0_pt->getVal())+","+
						std::to_string(bMass->getVal())+","+
						std::to_string(hnlcos2D->getVal())+","+
						std::to_string(hnl_lxySig->getVal())+","+
						std::to_string(hnlvtxProb->getVal())+","+
						std::to_string(hnlchi2->getVal())+","+
						std::to_string(hnlpiDCAS->getVal())+","+
						std::to_string(mass)+")").c_str()
					);
			h_nn->Fill(var,weight->getVal());
			h_mass[0]->Fill(hnl_mass->getVal(),weight->getVal());
			h_mass[1]->Fill(bMass->getVal(),weight->getVal());
			nn_score->setVal(Float_t(var));
				//			entry->Print("V");
				//		entry->add(*wspace.var("hnl_mass"));
				//		entry->add(*weight);
				//		entry->add(*nn_score);
			temp->add(RooArgList(*bMass,*target,*weight),weight->getVal());
			
		}else{

			nn_score->setVal(-1.0);
			//	entry->add(*wspace.var("hnl_mass"));
			//	entry->add(*weight);
			//entry->add(*nn_score);
			temp->add(RooArgList(*bMass,*weight),weight->getVal());

		}
			
	}
//	h_weights->SaveAs((wd+"/hweights"+category+".root").c_str());
//	h_ct->SaveAs((wd+"/h_ct_weighted_"+category+".root").c_str());
	if (doNNsel){
		
			if(ele_model){
				h_nn->SaveAs((wd+"/h_nn"+category+"_eleModel.root").c_str());	
				h_mass[1]->SaveAs((wd+"/B_mass"+category+"_eleModel.root").c_str());
				}
			else{
				 h_nn->SaveAs((wd+"/h_nn"+category+"_muModel.root").c_str());	
				 h_mass[1]->SaveAs((wd+"/B_mass"+category+"_muModel.root").c_str());
				}

//	RooDataSet mcSet(*temp,"mcSet");
	if (data){
	temp->Print();
	wspace.factory("mean[5.29,4.9,5.4]");
	wspace.factory("sigma[0.1,0.01,1]");
	wspace.factory("n1_cb[5,0,1000]");
	wspace.factory("alpha1_cb[0.1,0,10000]");
	wspace.factory("alpha[-0.01,-0.5,1]");
	wspace.factory("alpha1[1,0,10]");
	wspace.factory("alpha2[1,0,10]");
	wspace.factory("Nsig[0.25,0,10]");
	wspace.factory("Nbkg[1,0,10]");
	wspace.factory("common[3,0,10]");
	RooLinearVar c1("c1", "c1", *wspace.var("Nsig"), *wspace.var("common"), RooFit::RooConst(0.));
	RooLinearVar c2("c2", "c2", *wspace.var("Nbkg"), *wspace.var("common"), RooFit::RooConst(0.));
	wspace.factory("CBShape::CB1(B_mass,mean,sigma,alpha1_cb,n1_cb)");
	wspace.factory("Bernstein::Bkg(B_mass,{alpha,alpha1,alpha2})");
	wspace.factory("SUM::model(Nsig*CB1,Nbkg*Bkg)");
	RooAbsPdf* mod = wspace.pdf("model");
	mod->fitTo(*temp, RooFit::SumW2Error(true));
	RooStats::SPlot splot("splot", "splot", *temp, mod, RooArgSet(*wspace.var("Nsig"),*wspace.var("Nbkg")));
	std::cout << std::endl
             << "Yield of B is\t" << wspace.var("Nsig")->getVal() << ".  From sWeights it is "
             << splot.GetYieldFromSWeight("Nsig") << std::endl;
	std::cout << std::endl
             << "Yield of B is\t" << wspace.var("Nbkg")->getVal() << ".  From sWeights it is "
             << splot.GetYieldFromSWeight("Nbkg") << std::endl;
	RooDataSet dataw{temp->GetName(), temp->GetTitle(), temp, *temp->get(), nullptr, "Nsig_sw"};
	dataw.Print();
//	auto B_mass_w = ->createHistogram("B_mass",*B_mass);
	//B_mass_w->SaveAs("Bmass_weighted.root");
	}
	}
	return h_nn;
	}


void saveRatio(TH1D* h1,TH1D* h2,std::string label, std::string legendString,bool isEff){

	TCanvas*  c= new TCanvas("c","c",600,800);
	if (!isEff){
	h1->Scale(1.0/h1->Integral());
	h2->Scale(1.0/h2->Integral());
	}
	setStyle();
	

	TRatioPlot * ratio = new TRatioPlot(h1,h2,"pois");	
//	std::cout << muModel->GetName() << std::endl;
	
	ratio->Draw("");

	
	gPad->Modified();
	gPad->Update(); // make sure it’s really (re)drawn
	TPad *p = ratio->GetUpperPad();
//	std::string signalLabel = "mass "+ std::to_string(int(mass))+"."+std::to_string(int(mass-int(mass))*10)+" GeV  c#tau "+ std::to_string(int(ctau))+"."+std::to_string(int(ctau-int(ctau))*10) + " mm";
	gStyle->SetLegendTextSize(0.03);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	TLegend *l = p->BuildLegend(0.1,0.7,0.5,0.9,(legendString).c_str());


//	l->SetHeader(("Signal"+signalLabel).c_str());
//	l->AddEntry(eleModel,"e channel pNN model");
//	l->AddEntry(eleModel,"#mu channel pNN model");
//	l->Draw("same");
	ratio->GetLowerRefYaxis()->SetRangeUser(0.1,3);
	if (isEff) {
	ratio->GetUpperRefYaxis()->SetRangeUser(0.0,1);
	ratio->GetLowerRefXaxis()->SetBinLabel(1,"lxysig0to50");
	ratio->GetLowerRefXaxis()->SetBinLabel(2,"lxysig50to150");
	ratio->GetLowerRefXaxis()->SetBinLabel(3,"lxysiggt150");
	ratio->GetLowerRefYaxis()->SetTitle("Data/MC");
	}
	std::cout << "____________________get number of points" << ratio->GetLowerRefGraph()->GetN() << std::endl;
	ratio->GetLowerRefGraph()->SetMarkerStyle(8);
	ratio->GetLowerRefGraph()->SetMarkerColor(9);
	ratio->GetLowerRefGraph()->SetLineColor(9);
	ratio->GetLowerRefGraph()->SetLineWidth(2);
	ratio->GetLowerRefXaxis()->Draw("same");
	gPad->Modified();
	gPad->Update(); // make sure it’s really (re)drawn
	
	c->SaveAs(("../plots/Systematics/"+label+"_JPsiMuMu_eleModel.pdf").c_str());
}


void compareModels(double* eff,std::string filename,std::string lxyLable,std::string sigLable, bool doNNsel,double mass,double ctau,std::string lxy,std::string sig,std::string sel, bool debug, double nn_cut,std::string wd){

	
	TH1D* mc_mass[2];
	mc_mass[0]=new TH1D("dilep_mass","dilep_mass",60,2.9,3.3);
	mc_mass[1]=new TH1D("B_mass","B_mass",60,4.5,5.4);
	mc_mass[0]->SetTitle("; mass_{#mu#mu}(GeV); normalized entries");
	mc_mass[1]->SetTitle("; mass_{#pi#mu#mu}(GeV); normalized entries");
	mc_mass[0]->SetLineWidth(2);
	mc_mass[0]->SetLineColor(8);
	mc_mass[0]->SetFillColorAlpha(8,0.4);
	mc_mass[1]->SetTitle(";B mass(GeV); normalized entries");
	mc_mass[1]->SetLineWidth(2);
	mc_mass[1]->SetLineColor(8);
	mc_mass[1]->SetFillColorAlpha(8,0.4);

	TH1D* data_mass[2];
	data_mass[0]=new TH1D("dilepton_mass","dilepton_mass",60,2.9,3.3);
	data_mass[1]=new TH1D("B_mass","B_mass",60,4.5,5.4);
	data_mass[0]->SetTitle("; mass_{#mu#mu}(GeV); normalized entries");
	data_mass[1]->SetTitle("; mass_{#pi#mu#mu}(GeV); normalized entries");
	data_mass[0]->SetLineWidth(2);
	data_mass[0]->SetLineColor(9);
	data_mass[0]->SetMarkerColor(9);
	data_mass[0]->SetMarkerStyle(8);
	data_mass[1]->SetLineWidth(2);
	data_mass[1]->SetLineColor(9);
	data_mass[1]->SetMarkerColor(9);
	data_mass[1]->SetMarkerStyle(8);
//	TH1D* eleModel =  testModel(filename,lxyLable,sigLable,1,mass,ctau,lxy,sig,sel,debug, nn_cut, wd,1);
//	TH1D* muModel = testModel(filename,lxyLable,sigLable,1,mass,ctau,lxy,sig,sel,debug, nn_cut, wd,0);

	TH1D* data_score = testModel("../python/test_JPsi_data_1D.root",lxyLable,sigLable,1,mass,ctau,lxy,sig,sel,debug, nn_cut, wd,1,data_mass,1);
	data_score->SetTitle("Data; normalized entries");
	data_score->SetName(" Data");
	data_score->SetLineWidth(2);
	data_score->SetLineColor(9);
	data_score->SetMarkerColor(9);
	data_score->SetMarkerStyle(8);
	TH1D* mc_score = testModel("../python/test_mc_BuToKJpsi_*.root",lxyLable,sigLable,1,mass,ctau,lxy,sig,sel,debug, nn_cut, wd,1,mc_mass,0);
	mc_score->SetTitle("MC ; normalized entries");
	mc_score->SetName(" MC");
	mc_score->SetLineWidth(2);
	mc_score->SetLineColor(8);
	mc_score->SetFillColorAlpha(8,0.4);


        std::string category = lxyLable+"_"+sigLable;
	eff[0] = mc_score->IntegralAndError(mc_score->FindBin(nn_cut),100,eff[1])/mc_score->Integral();
	eff[2] = data_score->IntegralAndError(data_score->FindBin(nn_cut),100,eff[3])/data_score->Integral();
	eff[1] = eff[1]/mc_score->Integral();
	eff[3] = eff[2]/data_score->Integral();

	saveRatio(mc_mass[0],data_mass[0],"hnl_mass_"+category,"#ell#pi mass (GeV) - "+category,0);
	saveRatio(mc_mass[1],data_mass[1],"B_mass_"+category,"#ell#pi mass (GeV) - "+category,0);
	saveRatio(mc_score,data_score,"score_"+category,"pNN score - "+category,0);

//	std:: cout << "____________________________________________________________eff ele " << eff[0] << " eff mu " << eff[1] << std::endl;
		

} 


void wrapper(){

	std::string lxyBins[3] = {"lxysig0to50","lxysig50to150","lxysiggt150"};
	std::string sig[2] = {"OS","SS"};
	std::string lxySel[3] = {"hnl_lxy_sig<50 && ","hnl_lxy_sig>50 && hnl_lxy_sig<150 &&","hnl_lxy_sig>150 && "};
	std::string sigSel[2] = {"LepQProd<0","LepQProd>0"};

	std::vector<double> GenEvtCount,Filter,Mass,Ctau,Sigma,Flav_ratio;
	std::vector<std::string> Files, label;
	sig_fromcsv("../data/MC_datasets_InclusiveFilterEff.csv",&Files,&label,&GenEvtCount,&Filter,&Mass,&Ctau,&Sigma,&Flav_ratio);

	int i,j,k,ic;
	double eff[4];
	double cut[4] = {0.3,0.5,0.8,0.95};
	std::vector<std::vector<double>> eff_ele;
	std::vector<std::vector<double>> eff_mu;

	TH1D* eff_hist[4][2];

	for (ic =0; ic<4;ic++){
		eff_hist[ic][0] = new TH1D("MC",";",3,0,3);
		eff_hist[ic][1] = new TH1D("Data",";",3,0,3);
		eff_hist[ic][0]->SetLineColor(kRed);
		eff_hist[ic][0]->SetLineStyle(7);
		eff_hist[ic][0]->SetLineWidth(2);
		eff_hist[ic][1]->SetLineColor(kRed);
		eff_hist[ic][1]->SetLineWidth(2);
	for (j=0; j<3 ;j++){	
			//eff_hist[ic][0]->GetXaxis()->SetBinLabel(j+1,lxyBins[j].c_str());
			//eff_hist[ic][1]->GetXaxis()->SetBinLabel(j+1,lxyBins[j].c_str());
		for (k=0; k<1 ;k++){	
			std::vector<double> temp_ele;
			std::vector<double> temp_mu;
			for (i=0;i<1;i++){
        			std::string category = lxyBins[j]+" "+sig[k];
			
				compareModels(eff,"../python/test_BuToKJPsi_data_*.root",lxyBins[j],sig[k],1,3.1,1.0,lxySel[j],sigSel[k],"", 1, cut[ic],".");
				temp_ele.push_back(eff[0]);
				temp_mu.push_back(eff[1]);
				std::cout << category << " " << eff[0] << " " << eff[1] << std::endl;
				eff_hist[ic][0]->SetBinContent(j+1,eff[0]);
				eff_hist[ic][0]->SetBinError(j+1,eff[1]);
				eff_hist[ic][1]->SetBinContent(j+1,eff[2]);
				eff_hist[ic][1]->SetBinError(j+1,eff[3]);


			}	
			eff_ele.push_back(temp_ele);
			eff_mu.push_back(temp_mu);
		}
	}
	saveRatio(eff_hist[ic][1],eff_hist[ic][0],"effVsLxyS_"+std::to_string(ic),"#epsilon_{pNN} in LxyS categories, BuToKJPsiMuMu, pNN_{#mu} score > "+std::to_string(cut[ic]),1);

	}
	
/*	for (ic =2; ic<3;ic++){


		
		
	for (j=0; j<3 ;j++){	
			for (k=0; k<1 ;k++){	
				std::cout << category << std::endl;
		//		for (i=0;i<Files.size();i++){
	
				

		//		}
				std::cout << " "  << std::endl;

				}
			}
	
	}*/
}


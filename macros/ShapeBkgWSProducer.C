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


/*int genEventCount(std::string filename){

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
*/


void readFile(const char* filename, std::vector<string>& lines)
{
    lines.clear();
    ifstream file(filename);
    string s;
    while (getline(file, s))
        lines.push_back(s);
}


//needs a bit of cleaning
int ShapeBkg(std::string filename,std::string lxyLable,std::string sigLable,std::string lxy,std::string sig,int i,std::string path,float mass, float sigma,int nSigma, std::string sel, std::string wd,bool binned, bool blinded,float nn_cut){

	//this has to be substituted with input from .csv with all signal info
	std::vector<string> files;
        readFile("../data/input_41p6_hits.txt",files);
/*	std::vector<string> files = {"../data/HNLFlatTuples/Parking_1D_0_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_0_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_0_2_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_1_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_1_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_2_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_2_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_3_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_3_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_4_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_4_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_5_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_5_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_5_2_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_6_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_6_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_7_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_7_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_8_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_8_1_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_9_0_valid/",
		 		"../data/HNLFlatTuples/Parking_1D_9_1_valid/"};*/
//	std::string files[24] = {"../data/HNLFlatTuples/Parking_1D_section0_0_PF","../data/HNLFlatTuples/Parking1D0_newSig_nn_section0_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section2_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section2_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section3_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section4_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section4_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section5_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section6_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section7_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section7_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section8_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section9_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section9_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section10_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section10_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section11_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section11_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section7_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section12_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section12_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section13_0","../data/HNLFlatTuples/Parking1D0_newSig_nn_section13_1","../data/HNLFlatTuples/Parking1D0_newSig_nn_section14_0	TChain* c= new TChain("Events");
	int Ntot =0;
	float FilterEff = 0;
        bool doOnlyHistos= false;
	double hnlMass[6]={1,1.5,2,3,4.5,5.5};
	TChain* c = new TChain("Events");
//	for (int i=0;i<files.size();i++){
	for (int i=0;i<files.size();i++){
	//std::cout << files[i]  << std::endl;
	c->Add((files[i]+"/*.root").c_str());
	
	std::cout << c->GetEntries() << std::endl;
	}	
	std::string B_c = "_Bc";
	std::size_t found = filename.find(B_c);
	bool Bc;
	if (found == std::string::npos){
	Bc= false;
	B_c="";
	}
	else  Bc = true;
	
	//SM resonance edges 
	double low_veto[3] = {1.814,3.047,3.636};
	double high_veto[3] = {1.914,3.147,3.736};


/*	std::string lxy[3]={"hnl_lxy_sig<50","hnl_lxy_sig>50 && hnl_lxy_sig<150","hnl_lxy_sig>150"};
	std::string lxyLables[3]={"lxysig0to50","lxysig50to150","lxysiggt150"};
	std::string sig[2]={"LepQProd<0","LepQProd>0"};
	std::string sig_lables[2] = {"OS","SS"};
*/

	//Define set of variables to include in dataset

//	std::string variables = "B_mass,hnl_mass,hnl_charge,hnl_lxy,Type,LxyBin,dr_trgmu_hnl,dilepton_mass,hnl_vtxProb,hnl_cos2D,hnl_l_pt,hnl_l_mvaId,hnl_pi_pt,hnl_l_dxyS,hnl_pi_dxyS,dilepton_pt,likelihood,hnl_lxy_sig,LepQProd,hnl_ct,nn_score,nn_score_single";
//	MC dataset initialization

	std::string sig_label;
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
	std::string category = lxyLable+"_"+sigLable;

	
	TPython::LoadMacro("../python/ComputeNNScore.py");
	std::string modelPath;
	// = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_12Oct2022_15h25m05s/";
	if (Bc)	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_Bc_28Jul2023_09h10m22s//";
	//if (Bc)	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_Bc_26Jan2023_00h25m43s//";
//	else	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_12Oct2022_15h25m05s//";
	else	modelPath = "/cmshome/ratramon/Analysis/python/model/trainNN_"+category+"_26Jul2023_17h22m57s//";
	std::cout << modelPath <<std::endl;
	TPython::Exec(("modelSpecs = NNinput('"+modelPath+"')").c_str());
	TPython::Exec(("modelSpecs = NNinput('"+modelPath+"')").c_str());
	double processedLumi = 41.6;
	category = lxyLable+"_"+sigLable+"_"+B_c;

	wspace.factory("r_ele[0.359,0.359,0.359]");
	wspace.factory("r_mu[0.641,0.641,0.641]");
	//if mass window upper of lower edge fall under SM resonance, resize the window avoiding the SM resonance
	
	double low_edge = mass-nSigma* sigma;
	double high_edge = mass+nSigma* sigma;
	bool noVeto = true;
	for (int idx=0;idx<3;idx++){
	
	if (mass>0){
		if (low_edge < low_veto[idx] &&  high_edge>low_veto[idx] && high_edge<high_veto[idx]){

	         wspace.factory(("hnl_mass["+std::to_string(mass)+","+std::to_string(mass-nSigma* sigma)+","+std::to_string(low_veto[idx])+"]").c_str());
		noVeto = false;
		}
		else if(low_edge>low_veto[idx] && low_edge<high_veto[idx] and high_edge>high_veto[idx]){

	        wspace.factory(("hnl_mass["+std::to_string(mass)+","+std::to_string(high_veto[idx])+","+std::to_string(mass+nSigma* sigma)+"]").c_str());
		noVeto = false;

		}


	}
	}
//	if (mass> 0 ) wspace.factory(("hnl_mass["+std::to_string(mass)+","+std::to_string(mass-nSigma* sigma)+","+std::to_string(mass+nSigma*sigma)+"]").c_str());
		
	if (noVeto) wspace.factory(("hnl_mass["+std::to_string(mass)+","+std::to_string(mass-nSigma* sigma)+","+std::to_string(mass+nSigma* sigma)+"]").c_str());
	wspace.var("hnl_mass")->Print();
	wspace.factory("B_mass[5.3,0,1000000]");
	wspace.factory("B_pt[1,0,1000000]");
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
	wspace.factory("l0_relIso[10,0,10000000]");
	wspace.factory("l_relIso[10,0,10000000]");
	wspace.factory("hnl_lxy_sig[10,-10000000,10000000]");
	wspace.factory("hnl_lxy[10,-100000000,10000000]");
	wspace.factory("hnl_ct[10,0,100000]");
	wspace.factory("likelihood[0,-100,100]");
	wspace.factory("hnl_cos2D[0.995,0,1]");
	wspace.factory("weight[1,0,10000000000000]");
	wspace.factory("c0[0.1,-1,1]");
	wspace.factory("c2[0,-10,10]");
	wspace.factory(("lumiscale["+std::to_string(41.6/processedLumi)+","+std::to_string(41.6/processedLumi)+","+std::to_string(41.6/processedLumi)+"]").c_str());
	wspace.factory("Nbkg[100,0,10000]");
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
	wspace.var("l_pt");
	wspace.var("l0_pt");
	wspace.var("hnl_l_eta");
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
	wspace.var("weight");
	wspace.var("nn_score");
	wspace.var("c1");
	wspace.var("c2");
	wspace.var("lumiscale");

	//Define set of variables to include in dataset

	std::string variables = "hnl_mass,hnl_charge,likelihood,hnl_ct,hnl_lxy,LepQProd,hnl_cos2D,hnl_vtxProb,hnl_pi_pt,hnl_pi_DCAS,B_mass,hnl_vtxChi2,l_pt,hnl_l_pt,TrgMu_pt,hnl_l_eta,l0_pt,hnl_lxy_sig,B_pt,dr_Blep_pi,dr_trgMu_lep,dilepton_mass,BlepPi_mass,TrgMu_relIso,hnl_l_relIso,l0_relIso,l_relIso,hnl_l_mvaId,weight,hnl_pi_PixelLayers,hnl_pi_TrackerLayers,TrgMu_PixelLayers,TrgMu_TrackerLayers";
        wspace.defineSet("treeSet",variables.c_str());
//	MC dataset initialization
	RooDataSet* data_un =new  RooDataSet("data_un","data_un",c,*wspace.set("treeSet"),(lxy+" && " +sig +sel).c_str());
	RooArgSet * entry;
	RooRealVar* nn_score= new RooRealVar("nn_score","nn_score",0.5,0,1000000);
//	RooRealVar* lumiScale= new RooRealVar("lumi_scale","lumi_scale",41.6/processedLumi,41.6/processedLumi,41.6/processedLumi);
        RooDataSet* temp; 
	temp=new  RooDataSet(data_un->GetName(),data_un->GetTitle(),*wspace.set("treeSet"),"lumiscale");

		TH1D* h_nn=new TH1D("nn_score","nn_score",100,0,1);
		std::cout << c->GetEntries() << std::endl;
//		RooFormulaVar wFunc("Weight","Weight","LifetimeReweight(mass,ctau,x[0])",*wspace.var("hnl_ct")) ;
    		//include weights
		data_un->Print();
		for (int i =0; i<data_un->sumEntries();i++){
//			std::cout << i << std::endl;
			entry = (RooArgSet* ) data_un->get(i);
			if (nn_cut>0){
			entry->add(*nn_score);
//			entry->add(*lumiScale);
			
		RooRealVar* hnlcos2D = (RooRealVar* ) data_un->get(i)->find("hnl_cos2D");
		RooRealVar* hnlvtxProb = (RooRealVar* ) data_un->get(i)->find("hnl_vtxProb");
		RooRealVar* hnl_piPt = (RooRealVar* ) data_un->get(i)->find("hnl_pi_pt");
		RooRealVar* bMass = (RooRealVar* ) data_un->get(i)->find("B_mass");
		RooRealVar* hnl_pi_DCAS = (RooRealVar* ) data_un->get(i)->find("hnl_pi_DCAS");
		RooRealVar* hnlchi2 = (RooRealVar* ) data_un->get(i)->find("hnl_vtxChi2");
		RooRealVar* hnl_lPt = (RooRealVar* ) data_un->get(i)->find("hnl_l_pt");
		RooRealVar* hnltrgmuPt = (RooRealVar* ) data_un->get(i)->find("TrgMu_pt");
		RooRealVar* l_Pt = (RooRealVar* ) data_un->get(i)->find("l_pt");
		RooRealVar* l0_Pt = (RooRealVar* ) data_un->get(i)->find("l0_pt");
		RooRealVar* hnl_lxySig = (RooRealVar* ) data_un->get(i)->find("hnl_lxy_sig");
		RooRealVar* B_pt = (RooRealVar* ) data_un->get(i)->find("B_pt");
		RooRealVar* dilep_mass = (RooRealVar* ) data_un->get(i)->find("dilepton_mass");
		RooRealVar* Blep_pi_mass = (RooRealVar* ) data_un->get(i)->find("BlepPi_mass");
		RooRealVar* dr_Blep_pi = (RooRealVar* ) data_un->get(i)->find("dr_Blep_pi");
		RooRealVar* dr_trgMu_lep = (RooRealVar* ) data_un->get(i)->find("dr_trgMu_lep");
		RooRealVar* l_iso = (RooRealVar* ) data_un->get(i)->find("hnl_l_relIso");
		RooRealVar* l0_iso = (RooRealVar* ) data_un->get(i)->find("TrgMu_relIso");

	//	RooRealVar* hnl_mass->setVal(data_un->get(i)->getRealValue("hnl_mass"));
			
		double var;	
	/*	if (Bc)	var = TPython::Eval(("ComputeNNScore("+
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
						std::to_string(data_un->get(i)->getRealValue("hnl_pi_PixelLayers"))+","+
						std::to_string(data_un->get(i)->getRealValue("hnl_pi_TrackerLayers"))+","+
						std::to_string(data_un->get(i)->getRealValue("TrgMu_PixelLayers"))+","+
						std::to_string(data_un->get(i)->getRealValue("TrgMu_TrackerLayers"))+","+
						std::to_string(mass)+")").c_str()
					);

			h_nn->Fill(var);
			if (Float_t(var)>nn_cut){	
				nn_score->setVal(Float_t(var));
		//		std::cout << Float_t(var) << std::endl;
	//			entry->Print("V");
	//		
				temp->add(*entry);
			}
			
	}else temp->add(*entry);
	}		
		
        RooDataSet data_obs(*temp,"data_obs");
	data_obs.Print();
	std::cout << "is dataset weighted: " << data_obs.isWeighted() <<std::endl;
	h_nn->SaveAs((wd+"/"+category+"h_nn.root").c_str());	
	wspace.factory(("Nsig[10,"+std::to_string(data_obs.sumEntries()*0.8)+","+std::to_string(data_obs.sumEntries()*1.5)+"]").c_str());
        wspace.import(data_obs);	
	wspace.Print();
	wspace.writeToFile((wd+"/Data_"+channel[1]+"_"+std::to_string(i)+"_"+sigLable+B_c+".root").c_str());
	
	
	return 0;
}


/*
int AllFits(std::string filename,std::string path, double mass,double ctau,int channel, std::string wd,std::string sel,double lumiRatio, bool binned, bool blinded ){
 

	int  j;
	double sigFeat[3];
	//std::string lxy_cuts[4]= {"hnl_lxy<3", "hnl_lxy>3 && hnl_lxy<10", "hnl_lxy>10 && hnl_lxy<20","hnl_lxy>20"};
	std::string lxy_cuts[4]= {"hnl_lxy<1", "hnl_lxy>1 && hnl_lxy<5", "hnl_lxy>5"};
	std::ofstream SigFeatures;
	
	SigFeatures.open(wd+"/Signal_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt");
	std::cout << wd+"/Signal_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt" << std::endl;
	

		std::cout << "in cycle " << std::endl;
		for (j=0;j<3;j++){
		const int temp_i = channel;
		const int temp_j=j;
		fit(filename.c_str(),path.c_str(),mass,ctau,channel,lxy_cuts[j],j,wd+"/SigFits",sigFeat,sel,lumiRatio,binned,blinded);
		channel=temp_i;
		j=temp_j;
	 	SigFeatures << channel << " "  << j << " " << sigFeat[0] <<  " " << sigFeat[1] << " " << sigFeat[2] <<  std::endl;
		}

	SigFeatures.close();

	return 0;
	}
*/

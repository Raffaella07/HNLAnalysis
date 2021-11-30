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



int fit(std::string filename, double mass, int type,std::string lxy,int lxyBin, std::string wd, double* sigF,std::string sel){

//	TFile* in = TFile::Open(filename.c_str());
	TChain* c= new TChain("Events");

	c->Add(filename.c_str());



	RooWorkspace wspace("w");
	gStyle->SetOptFit(0000);
	gROOT->SetBatch(true);
	gROOT->SetStyle("Plain");
	gStyle->SetGridStyle(3);
	gStyle->SetOptStat(000000);
	gStyle->SetOptTitle(0);
	float SigMu[4]= {336,538,386,173};
	float SigEle[4]= {133,347,323,284};
	std::string channel[4]= {"Muon","PFele","LoePtele","Track"};
//	std::string lxy[4]= {"lxyUnder3","lxy_Over3_Under10","lxy_over10_Under20","lxy_Over20"};
//	wspace.factory(("mean_ele["+std::to_string(hist->GetMean())+","+std::to_string(hist->GetMean()-hist->GetStdDev())+","+std::to_string(hist->GetMean()-hist->GetStdDev())+"]").c_str());
//	wspace.factory(("mean_mu["+std::to_string(hist->GetMean())+","+std::to_string(hist->GetMean()-hist->GetStdDev())+","+std::to_string(hist->GetMean()-hist->GetStdDev())+"]").c_str());
//	wspace.factory(("sigma_ele[0.1,0,"+std::to_string(hist->GetStdDev()*1.5)+"]").c_str());
//	wspace.factory(("sigma_mu[0.1,0,"+std::to_string(hist->GetStdDev()*1.5)+"]").c_str());
//	wspace.factory("n_cb[0,-10,10]");
//	wspace.factory("alpha_cb[0,-10,10]");
//	wspace.factory(("Nele["+std::to_string(SigEle[type])+","+std::to_string(SigEle[type]-sqrt(SigEle[type]))+","+std::to_string(SigEle[type]+sqrt(SigEle[type]))+"]").c_str());
	wspace.factory("r_ele[0.359,0.359,0.359]");
	wspace.factory("r_mu[0.641,0.641,0.641]");
	if (type ==0 )wspace.factory("hnl_mass[3,0.7,5]");
	else wspace.factory("hnl_mass[3,0.7,5]");
	wspace.factory("B_mass[5.3,2,6]");
	wspace.factory("Type[0,0,4]");
	wspace.factory("LxyBin[0,0,4]");
	wspace.factory("dr_trgmu_hnl[0,-3,3]");
	wspace.factory("dilepton_mass[0,0,7]");
	wspace.factory("dilepton_pt[0,0,100]");
	wspace.factory("hnl_vtxProb[0,0,1]");
	wspace.factory("hnl_pi_dxyS[0,-100,100]");
	wspace.factory("hnl_l_dxyS[0,-100,100]");
	wspace.factory("hnl_l_pt[10,0,100]");
	wspace.factory("hnl_pi_pt[10,0,100]");
	wspace.factory("hnl_lxy[10,0,200]");
	wspace.factory("likelihood[0,-1,1]");
	wspace.factory("LepQProd[0,-1,1]");
	//wspace.factory(("Nele[10,7,"+std::to_string(hist->GetMaximum())+"]").c_str());
	//wspace.factory(("Nmu[10,7,"+std::to_string(hist->GetMaximum())+"]").c_str());
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
	wspace.var("hnl_mass");
	wspace.var("hnl_trgmu_hnl");
	wspace.var("dilepton_mass");
	wspace.var("dilepton_pt");
	wspace.var("hnl_vtxProb");
	wspace.var("hnl_l_dxyS");
	wspace.var("hnl_pi_dxyS");
	wspace.var("B_mass");
	wspace.var("LepQProd");
	wspace.var("Type");
	wspace.var("LxyBin");
	//wspace.var("r_mu")->setVal(1-wspace.var("r_ele")->getVal());
	std::string variables = "B_mass,hnl_mass,Type,LxyBin,dr_trgmu_hnl,dilepton_mass,hnl_vtxProb,hnl_l_pt,hnl_pi_pt,dilepton_pt,likelihood,hnl_lxy,LepQProd";
        wspace.defineSet("treeSet",variables.c_str());
	//RooDataSet mcSet("mcSet","mcSet",*wspace.set("treeSet"),RooFit::Import(*c),RooFit::Cut(("Type=="+std::to_string(type)+" && LxyBin=="+std::to_string(lxyBin)+sel).c_str()));
	RooDataSet mcSet("mcSet","mcSet",*wspace.set("treeSet"),RooFit::Import(*c),RooFit::Cut(("Type=="+std::to_string(type)+" && "+lxy +sel).c_str()));
	std:: cout << ("Type=="+std::to_string(type)+" && "+lxy+sel).c_str() << " " << mcSet.sumEntries()<<std::endl;
	if (mcSet.sumEntries()<10){
	sigF[0]	= mcSet.sumEntries();
	sigF[1]	= mcSet.sumEntries();
	sigF[2]	= -99;
	sigF[3]	= -99;
	return 0;
	

	}else{
	wspace.var("hnl_mass")->setMin(mcSet.mean(*wspace.var("hnl_mass"))-4*mcSet.sigma(*wspace.var("hnl_mass")));
	wspace.var("hnl_mass")->setMax(mcSet.mean(*wspace.var("hnl_mass"))+4*mcSet.sigma(*wspace.var("hnl_mass")));
	}

	wspace.factory(("Nsig[10,"+std::to_string(mcSet.sumEntries()*0.8)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());
	wspace.Print();
	if (type==3){

		wspace.factory(("mean_ele[3.00,"+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory(("mean_mu[3.00,"+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory(("sigma_ele["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.08)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.05)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.5)+"]").c_str());
		wspace.factory(("sigma_mu["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.08)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.05)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.5)+"]").c_str());
		wspace.factory("n1_cb[7,0,10]");
		wspace.factory("alpha1_cb[3,0,10]");
		wspace.factory("n2_cb[5,0,10]");
		wspace.factory("alpha2_cb[-3,-10,0]");
//		wspace.factory(("Nele["+std::to_string(SigEle[lxyBin])+","+std::to_string(SigEle[lxyBin]-sqrt(SigEle[lxyBin]))+","+std::to_string(SigEle[lxyBin]+sqrt(SigEle[lxyBin]))+"]").c_str());
		wspace.factory("Nele[10,1,10000]");
		wspace.factory("N1[10,1,1000000]");
		wspace.factory("N2[10,1,1000000]");
		wspace.factory(("Nmu["+std::to_string(SigMu[lxyBin])+","+std::to_string(SigMu[lxyBin]-sqrt(SigMu[lxyBin]))+","+std::to_string(SigMu[lxyBin]+sqrt(SigMu[lxyBin]))+"]").c_str());
		wspace.factory("CBShape::CB1(hnl_mass,mean_ele,sigma_ele,alpha1_cb,n1_cb)");
		wspace.factory("CBShape::CB2(hnl_mass,mean_ele,sigma_ele,alpha2_cb,n2_cb)");
		wspace.factory("Gaussian::gaus(hnl_mass,mean_ele,sigma_mu)");
	//	wspace.factory("SUM::shape(r_ele*CB1,r_mu*gaus)");
	//	wspace.factory("SUM::shapeGaus(r_mu*gaus)");


	}else
	 if (type==1 || type==2){


		wspace.factory(("mean_ele["+std::to_string(mcSet.mean(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory(("sigma_ele["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.08)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.05)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.5)+"]").c_str());
		wspace.factory(("sigma_mu[0,0,0]"));
		wspace.factory("n1_cb[7,0,10]");
		wspace.factory("alpha1_cb[3,0,10]");
		wspace.factory("n2_cb[5,0,10]");
		wspace.factory("alpha2_cb[-3,-10,0]");
		wspace.factory("n_cb[7,0,10]");
		wspace.factory("alpha_cb[3,0,10]");
		wspace.factory("N1[10,1,1000000]");
		wspace.factory("N2[10,1,1000000]");
 		wspace.factory(("Nele["+std::to_string(mcSet.sumEntries()*1)+","+std::to_string(mcSet.sumEntries()*0.5)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());
		wspace.factory(("Nmu[0,0,0]"));
		wspace.factory("CBShape::CB1(hnl_mass,mean_ele,sigma_ele,alpha1_cb,n1_cb)");
		wspace.factory("CBShape::CB2(hnl_mass,mean_ele,sigma_ele,alpha2_cb,n2_cb)");
	//	wspace.factory("SUM::signal(Nele*CB1,Nmu*CB2)");


	}else if (type==0){


		wspace.factory(("mean_ele["+std::to_string(mcSet.mean(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory(("sigma_ele["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.08)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.05)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.5)+"]").c_str());
	//	wspace.factory(("sigma_ele[0,0,0]"));
		wspace.factory("n1_cb[7,0,10]");
		wspace.factory("alpha1_cb[3,0,10]");
		wspace.factory("n2_cb[5,0,10]");
		wspace.factory("alpha2_cb[-3,-10,0]");
		wspace.factory(("width[0.1,0,1]"));
		wspace.factory("N1[10,1,1000000]");
		wspace.factory("N2[10,1,1000000]");
 		wspace.factory(("Nmu["+std::to_string(mcSet.sumEntries()*1)+","+std::to_string(mcSet.sumEntries()*0.5)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());
		wspace.factory(("Nele[0,0,0]"));
		wspace.factory("Voigtian::voigt(hnl_mass,mean_mu,sigma_mu,width)");
		wspace.factory("Gaussian::gaus(hnl_mass,mean_mu,sigma_mu)");
		wspace.factory("CBShape::CB1(hnl_mass,mean_ele,sigma_ele,alpha1_cb,n1_cb)");
		wspace.factory("CBShape::CB2(hnl_mass,mean_ele,sigma_ele,alpha2_cb,n2_cb)");
	//	wspace.factory("SUM::signal(Nmu*gaus)");


	}

	
	wspace.factory("SUM::signal(N1*CB1,N1*CB2)");
//	RooFFTConvPdf signal("signal","CB1 (X) CB2",*wspace.var("hnl_mass"),*wspace.pdf("CB1"),*wspace.pdf("CB2"));
	
//	if (type ==3 )wspace.factory("SUM::signal(Nele*CB1,Nmu*CB2)");
	RooAbsPdf* signal = wspace.pdf("signal");
	RooFitResult* res = signal->fitTo(mcSet,RooFit::Range(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")),mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass"))),RooFit::Save());
	TCanvas* c1 =new TCanvas("MC_Fit", "fit", 800, 1200);
	RooPlot* mhnl_frame = wspace.var("hnl_mass")->frame();
	mcSet.plotOn(mhnl_frame, RooFit::Name("template mc"));
	signal->plotOn(mhnl_frame,RooFit::Name("signal"),RooFit::Range(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")),mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass"))),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
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
	sigF[2]	= sqrt(pow(wspace.var("sigma_ele")->getVal(),2));//+pow(wspace.var("sigma_mu")->getVal(),2));
	sigF[3]	= sqrt(pow(wspace.var("sigma_ele")->getVal(),2));//+pow(wspace.var("sigma_mu")->getVal(),2));
	
	
	RooAbsReal*  FullSpectrum = wspace.pdf("signal")->createIntegral(*wspace.var("hnl_mass"));
	if (type==0)wspace.var("hnl_mass")->setRange("twoSigma",wspace.var("mean_ele")->getVal()-2*sigF[2],wspace.var("mean_ele")->getVal()+2*sigF[2]) ;
	else wspace.var("hnl_mass")->setRange("twoSigma",wspace.var("mean_ele")->getVal()-3*sigF[2],wspace.var("mean_ele")->getVal()+2*sigF[2]) ;
	RooAbsReal*  yield = wspace.pdf("signal")->createIntegral(*wspace.var("hnl_mass"),RooFit::Range("twoSigma"));
//	RooAbsReal*  yield1 = wspace.pdf("CB1")->createIntegral(*wspace.var("hnl_mass"),RooFit::Range("twoSigma"));
//	RooAbsReal*  yield2 = wspace.pdf("CB2")->createIntegral(*wspace.var("hnl_mass"),RooFit::Range("twoSigma"));
	std:cout <<  yield->getVal()/FullSpectrum->getVal()*mcSet.sumEntries() << std::endl;
	sigF[0]	= yield->getVal()/FullSpectrum->getVal()*mcSet.sumEntries();
	sigF[1]	= sqrt(pow(yield->getPropagatedError(*res)/FullSpectrum->getVal(),2)+pow(FullSpectrum->getPropagatedError(*res)*yield->getVal()/pow(FullSpectrum->getVal(),2),2))*mcSet.sumEntries();
	if (std::isnan(sigF[1])) sigF[1]= -99;
//	sigF[2]	= sqrt(pow(wspace.var("sigma_ele")->getVal(),2)+pow(wspace.var("sigma_mu")->getVal(),2));
//	sigF[3]	= sqrt(pow(wspace.var("sigma_ele")->getVal(),2)+pow(wspace.var("sigma_mu")->getVal(),2));

	c1->SaveAs((wd+"/Fit_"+channel[type]+"_"+std::to_string(lxyBin)+".pdf").c_str());
	std::cout << "type " << type << " lxyBin "<< lxyBin<< std::endl;
	return 0;
}


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


int AllFits(std::string filename, double mass,double ctau, std::string wd,std::string sel ){
 

	int i, j;
	double sigFeat[3];
	std::string lxy_cuts[3]= {"hnl_lxy<1", "hnl_lxy>1 && hnl_lxy<5","hnl_lxy>5"};
	std::ofstream SigFeatures;
	
	SigFeatures.open(wd+"/Signal_mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt");
	std::cout << genEventCount(filename) << std::endl;
	for (i=0;i<4;i++){
	

		std::cout << "in cycle " << std::endl;
		for (j=0;j<3;j++){
		const int temp_i = i;
		const int temp_j=j;
		fit(filename.c_str(),mass,i,lxy_cuts[j],j,wd+"/SigFits",sigFeat,sel);
		i=temp_i;
		j=temp_j;
	 	SigFeatures << i << " "  << j << " " << sigFeat[0] <<  " " << sigFeat[1] << " " << sigFeat[2] <<  std::endl;
		std::cout << "after print  " << i << " " << j << std::endl;
		}

	}
	SigFeatures.close();

	return 0;
	}


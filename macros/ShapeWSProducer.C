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
#include "LifetimeReweight.C"
#include "/cmshome/ratramon/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/src/RooDoubleCBFast.cc"
#include "/cmshome/ratramon/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/src/RooMultiPdf.cxx"
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


//needs a bit of cleaning
int fit(std::string filename, std::string path, double mass,double ctau, int type,std::string lxy,int lxyBin, std::string wd, double* sigF,std::string sel,double lumiRatio,bool binned, bool blinded){

	//this has to be substituted with input from .csv with all signal info
	std::string files[5][4] = {{"../data/HNLFlatTuples/NewSignature_1p0_ctau10p0/","../data/HNLFlatTuples/NewSignature_1p0_ctau100p0/","../data/HNLFlatTuples/NewSignature_1p0_ctau1000p0/","placeholder"},
{"../data/HNLFlatTuples/NewSignature_1p5_ctau10p0/","../data/HNLFlatTuples/NewSignature_1p5_ctau100p0/","../data/HNLFlatTuples/NewSignature_1p5_ctau1000p0/","placeholder"},
{"../data/HNLFlatTuples/NewSignature_2p0_ctau10p0/","../data/HNLFlatTuples/NewSignature_2p0_ctau100p0/","../data/HNLFlatTuples/NewSignature_2p0_ctau1000p0/","placeholder"},
{"../data/HNLFlatTuples/NewSignature_3p0_ctau1p0/","../data/HNLFlatTuples/NewSignature_3p0_ctau10p0/","../data/HNLFlatTuples/NewSignature_3p0_ctau100p0/","../data/HNLFlatTuples/NewSignature_3p0_ctau1000p0/"},
{"../data/HNLFlatTuples/NewSignature_4p5_ctau0p1/","../data/HNLFlatTuples/NewSignature_4p5_ctau1p0/","../data/HNLFlatTuples/NewSignature_4p5_ctau10p0/","../data/HNLFlatTuples/NewSignature_4p5_ctau100p0/"},
};
	double genEvtCount[5][4]={{9.080620000000000000e+05,4.532249000000000000e+06,9.950223000000000000e+06,-1},
				{1.666740000000000000e+05,1.955714000000000000e+06,-1,-1},
                                {5.673050000000000000e+05,1.916707000000000000e+06,-1},
                                {4.411610000000000000e+05,6.562680000000000000e+05,1.314002000000000000e+06,3.282584000000000000e+06},
                                {1.653590000000000000e+05,5.065250000000000000e+05,5.106050000000000000e+05,1.195035000000000000e+06}};

	//update with filter efficiencies for specific channel
	double filter[5][4]={{5.85e-03,	5.17e-03,1.59e-03,-1},
				{3.72e-03,3.35e-03,1.03e-03,-1},
                                {1.78e-03,1.57e-03,8.95e-04,-1},
                                {5.87e-03,5.58e-03,2.40e-03,1.78e-03},
                                {8.62e-04,8.35e-04,8.00e-04,6.10e-04}};

	TChain* c= new TChain("Events");
	int Ntot =0;
	float FilterEff = 0;
        bool doOnlyHistos= false;
	double grid[5][4]={     {10,100,1000,-1},
				{10,100,1000,-1},
				{10,100,-1,-1},
				{1,10,100,1000},
				{0.1,1,10,100},
			};
	double hnlMass[5]={1,1.5,2,3,4.5};
	if (filename!= ""){
	c->Add((filename+"HNLFlat_0.root").c_str());
	for(int i=0;i<5;i++){
		if (mass==hnlMass[i]){
		for(int j=0;j<4;j++){
		    if (ctau == grid[i][j]){
			Ntot=genEvtCount[i][j];
			FilterEff = filter[i][j];
		    }
	}
	}
	}
	}
	else {
		
	for(int i=0;i<5;i++){

	if (mass==hnlMass[i]){
		
		for(int j=0;j<4;j++){
		    if (grid[i][j]>0 && ctau < grid[i][j]){
		    std::cout << "mass " << mass << std::endl;
		    std::cout << grid[i][j] << " " << ctau << std::endl;
			
		    c->Add((files[i][j]+"HNLFlat_0.root").c_str());
		    Ntot += genEvtCount[i][j];
		    FilterEff += genEvtCount[i][j]*filter[i][j];

							}


				}


	  		}		
		}
	FilterEff = FilterEff/Ntot;
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
	wspace.factory("hnl_ct[10,0,100000]");
	wspace.factory("likelihood[0,-100,100]");
	wspace.factory("hnl_cos2D[0.995,0,1]");
	wspace.factory("c0[0.1,-1,1]");
	wspace.factory("c1[0.5,-10,10]");
	wspace.factory("c2[0,-10,10]");
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
	wspace.var("hnl_charge");
	wspace.var("hnl_lxy_sig");
	wspace.var("hnl_lxy");
	wspace.var("hnl_ct");
	wspace.var("hnl_Blep_hnl");
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

	//Define set of variables to include in dataset

	std::string variables = "B_mass,hnl_mass,hnl_charge,hnl_lxy,Type,LxyBin,dr_trgmu_hnl,dilepton_mass,hnl_vtxProb,hnl_cos2D,hnl_l_pt,hnl_l_mvaId,hnl_pi_pt,hnl_l_dxyS,hnl_pi_dxyS,dilepton_pt,likelihood,hnl_lxy_sig,LepQProd,hnl_ct";
        wspace.defineSet("treeSet",variables.c_str());
//	MC dataset initialization
        RooDataSet* temp; 
	if (filename==""){
		std::cout << c->GetEntries() << std::endl;
//		RooFormulaVar wFunc("Weight","Weight","LifetimeReweight(mass,ctau,x[0])",*wspace.var("hnl_ct")) ;
		RooDataSet mcSet_un("mcSet_un","mcSet_un",c,*wspace.set("treeSet"),("Type=="+std::to_string(type)+" && "+lxy +sel).c_str())/*("Type=="+std::to_string(type)+" && "+lxy +sel).c_str()))*/;
    		//include weights
	    	RooDataSet* Addweights; 
		mcSet_un.Print();
		RooRealVar* weight= new RooRealVar("weight","weight",1,0,1000000);
		RooArgSet entry;
		entry.add(*wspace.var("hnl_mass"));
		entry.add(*weight);
                	
		temp=new  RooDataSet(mcSet_un.GetName(),mcSet_un.GetTitle(),entry,RooFit::WeightVar("weight"));
		TH1D* h_weights=new TH1D("w","w",100,0,20);
		for (int i =0; i<mcSet_un.sumEntries();i++){
		 	RooRealVar* target = (RooRealVar* ) mcSet_un.get(i)->find("hnl_ct");
			RooRealVar* hnlM = (RooRealVar* ) mcSet_un.get(i)->find("hnl_mass");
			std::cout << target->getVal() << " " << LifetimeReweight(mass,ctau,target->getVal())<< std::endl;
			wspace.var("hnl_mass")->setVal(hnlM->getVal());
			entry.Print("V");
			h_weights->Fill(LifetimeReweight(mass,ctau,target->getVal()));
			temp->add(entry,LifetimeReweight(mass,ctau,target->getVal()));
		}
		TH1D* h_weightsCopy = (TH1D*)temp->createHistogram("weight_copy",*weight,RooFit::Binning(100,0,20));
		h_weightsCopy->SaveAs("hweights_copy.root");
        	h_weights->SaveAs("hweights.root");	
		
	}else{	
		temp =new  RooDataSet("temp","temp",c,*wspace.set("treeSet"),("Type=="+std::to_string(type)+" && "+lxy +sel).c_str());
	//	temp->Print();
	}
        RooDataSet mcSet(*temp,"mcSet");
	std:: cout << ("Type=="+std::to_string(type)+" && "+lxy+sel).c_str() << " " << mcSet.sumEntries() <<std::endl;
	wspace.factory(("Nsig[10,"+std::to_string(mcSet.sumEntries()*0.8)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());
	
	wspace.Print();

	//Define signal fit function, might shrink the section but for now lets be redundant
	//might be useful for signal shape syst studies 
	if (type==3){

		wspace.var("hnl_mass")->setMin(mass-0.3);
		wspace.var("hnl_mass")->setMax(mass+0.2);
		wspace.factory(("mean_ele[3,"+std::to_string(mass-0.5*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+0.5*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory(("mean_mu[3.00,"+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory(("sigma_ele["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.07)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.03)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.5)+"]").c_str());
		wspace.factory(("sigma_mu["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.08)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.05)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.5)+"]").c_str());
		wspace.factory("n1_cb[3,0,10]");
		wspace.factory("alpha1_cb[3,0,100]");
		wspace.factory("n2_cb[7,0,10]");
		wspace.factory("alpha2_cb[-3,-100,0]");
//		wspace.factory(("Nele["+std::to_string(SigEle[lxyBin])+","+std::to_string(SigEle[lxyBin]-sqrt(SigEle[lxyBin]))+","+std::to_string(SigEle[lxyBin]+sqrt(SigEle[lxyBin]))+"]").c_str());
		wspace.factory("Nele[10,1,10000]");
		wspace.factory("N1[100,0,1000000]");
		wspace.factory("N2[100,0,1000000]");
		wspace.factory(("Nmu["+std::to_string(SigMu[lxyBin])+","+std::to_string(SigMu[lxyBin]-sqrt(SigMu[lxyBin]))+","+std::to_string(SigMu[lxyBin]+sqrt(SigMu[lxyBin]))+"]").c_str());
		wspace.factory("CBShape::CB1(hnl_mass,mean_ele,sigma_ele,alpha1_cb,n1_cb)");
		wspace.factory("CBShape::CB2(hnl_mass,mean_ele,sigma_ele,alpha2_cb,n2_cb)");
		wspace.factory("Gaussian::gaus(hnl_mass,mean_ele,sigma_mu)");
	//	wspace.factory("SUM::shape(r_ele*CB1,r_mu*gaus)");
	//	wspace.factory("SUM::shapeGaus(r_mu*gaus)");


	}else
	 if (type==1 || type==2|| type ==0){


		wspace.var("hnl_mass")->setMin(mass-0.3);
		wspace.var("hnl_mass")->setMax(mass+0.2);
		wspace.factory(("mean_ele["+std::to_string(mcSet.mean(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-2*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+2*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory(("sigma_ele["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.17)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.13)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.5)+"]").c_str());
		wspace.factory(("sigma_mu[0,0,0]"));
		wspace.factory("n1_cb[10,0,1000]");
		wspace.factory("alpha1_cb[0.1,0,10000]");
		wspace.factory("n2_cb[10,0,10000]");
		wspace.factory("alpha2_cb[0.1,0,10000]");
		wspace.factory("n_cb[0.5,0,1]");
		wspace.factory("alpha_cb[10,0,1000]");
		wspace.factory("N1[0.25,0,0.5]");
		wspace.factory("N2[0.25,0,0.5]");
 		wspace.factory(("Nele["+std::to_string(mcSet.sumEntries()*1)+","+std::to_string(mcSet.sumEntries()*0.5)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());
		wspace.factory(("Nmu[0,0,0]"));
		wspace.factory(("width[0.1,0,1]"));
//		wspace.factory("Gaussian::signal(hnl_mass,mean_ele,sigma_ele)");
//		wspace.factory("Voigtian::signal(hnl_mass,mean_ele,sigma_ele,width)");
		wspace.factory("CBShape::CB1(hnl_mass,mean_ele,sigma_ele,alpha1_cb,n1_cb)");
		wspace.factory("CBShape::CB2(hnl_mass,mean_ele,sigma_ele,alpha2_cb,n2_cb)");
	//	wspace.factory("SUM::signal(Nele*CB1,Nmu*CB2)")


	}else if (type==4){

		wspace.var("hnl_mass")->setMin(mcSet.mean(*wspace.var("hnl_mass"))-0.3);
		wspace.var("hnl_mass")->setMax(mcSet.mean(*wspace.var("hnl_mass"))+1*mcSet.sigma(*wspace.var("hnl_mass")));

		wspace.factory(("mean_ele["+std::to_string(mcSet.mean(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))-3*mcSet.sigma(*wspace.var("hnl_mass")))+","+std::to_string(mcSet.mean(*wspace.var("hnl_mass"))+3*mcSet.sigma(*wspace.var("hnl_mass")))+"]").c_str());
		wspace.factory("mean_conv[0,-3,3]");
		wspace.factory(("sigma_ele["+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.03)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*0.02)+","+std::to_string(mcSet.sigma(*wspace.var("hnl_mass"))*1)+"]").c_str());
		wspace.factory("sigma_conv[0.1,0.0001,1]");
	//	wspace.factory(("sigma_ele[0,0,0]"));
		wspace.factory("n1_cb[1,0.5,1]");
		wspace.factory("alpha1_cb[7,0,10]");
		wspace.factory("n2_cb[1,0.5,1.5]");
		wspace.factory("alpha2_cb[-7,-10,0]");
		wspace.factory(("width[0.1,0,1]"));
		wspace.factory("N1[10,1,1000000]");
		wspace.factory("N2[10,1,1000000]");
 		wspace.factory(("Nmu["+std::to_string(mcSet.sumEntries()*1)+","+std::to_string(mcSet.sumEntries()*0.5)+","+std::to_string(mcSet.sumEntries()*1.5)+"]").c_str());
		wspace.factory(("Nele[0,0,0]"));
		wspace.factory("Voigtian::voigt(hnl_mass,mean_mu,sigma_mu,width)");
		wspace.factory("Gaussian::signal(hnl_mass,mean_mu,sigma_mu)");
		wspace.factory("CBShape::CB1(hnl_mass,mean_ele,sigma_ele,alpha1_cb,n1_cb)");
		wspace.factory("CBShape::CB2(hnl_mass,mean_ele,sigma_ele,alpha2_cb,n2_cb)");


	}
//      actual signal function used:	
	RooDoubleCBFast signal("signal","DoubleSidedCB",*wspace.var("hnl_mass"),*wspace.var("mean_ele"),*wspace.var("sigma_ele"),*wspace.var("alpha1_cb"),*wspace.var("n1_cb"),*wspace.var("alpha2_cb"),*wspace.var("n2_cb"));	
//      signal fit 	
	RooFitResult* res = signal.fitTo(mcSet,RooFit::Range(mcSet.mean(*wspace.var("hnl_mass"))-2*mcSet.sigma(*wspace.var("hnl_mass")),mcSet.mean(*wspace.var("hnl_mass"))+2*mcSet.sigma(*wspace.var("hnl_mass"))),RooFit::Save(),RooFit::SumW2Error(kTRUE));
	
	//resize dataset to background fit window (nSigma), for now is 10 

	wspace.var("hnl_mass")->setMin(wspace.var("mean_ele")->getVal()-10*wspace.var("sigma_ele")->getVal()) ;
	wspace.var("hnl_mass")->setMax(wspace.var("mean_ele")->getVal()+10*wspace.var("sigma_ele")->getVal()) ;
	
	//Signal only plots


	TCanvas* c1 =new TCanvas("MC_Fit", "fit", 1200, 800);
	RooPlot* mhnl_frame = wspace.var("hnl_mass")->frame();
	mcSet.plotOn(mhnl_frame, RooFit::Name("template mc"),RooFit::Binning(40,wspace.var("mean_ele")->getVal()-6*wspace.var("sigma_ele")->getVal(),wspace.var("mean_ele")->getVal()+6*wspace.var("sigma_ele")->getVal()));
	signal.plotOn(mhnl_frame,RooFit::Name("signal"),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
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
	Double_t chi2 = mhnl_frame->chiSquare(5);
	std::cout << chi2 << std::endl;
	mhnl_frame->Draw();
	TLatex* latex= new TLatex();
	latex->SetTextAlign(12);
   	latex->SetTextFont(43);
   	latex->SetTextSize(30);
	mhnl_frame->addObject(latex);
  	latex->DrawLatex(wspace.var("mean_ele")->getVal()-5*wspace.var("sigma_ele")->getVal(),mhnl_frame->GetMaximum()*0.92,("#chi^{2}/n.dof"+std::string(Form(" %.2f", chi2))).c_str());
  	latex->DrawLatex(wspace.var("mean_ele")->getVal()-5*wspace.var("sigma_ele")->getVal(),mhnl_frame->GetMaximum()*0.97,("lxy bin "+std::to_string(lxyBin)+", HNL Mass"+std::string(Form(" %.2f", mass))+" GeV c#tau = "+std::string(Form(" %.2f", ctau))+ " mm").c_str());

	//export sigma and its error
	sigF[2]	= sqrt(pow(wspace.var("sigma_ele")->getVal(),2));//+pow(wspace.var("sigma_mu")->getVal(),2));
	sigF[3]	= sqrt(pow(wspace.var("sigma_ele")->getVal(),2));//+pow(wspace.var("sigma_mu")->getVal(),2));

	c1->SaveAs((wd+"/Fit_"+channel[type]+"_"+std::to_string(lxyBin)+".pdf").c_str());
	std::cout <<"entries (weighted) "<< mcSet.sumEntries()<< " NgenTot "<<  Ntot << "filter eff " << FilterEff << "eff * acc " << mcSet.sumEntries()*FilterEff/Ntot <<std::endl;
	//save also in output  vector Nsig + its error
	sigF[0]	= yield->getVal()/FullSpectrum->getVal()*mcSet.sumEntries()*FilterEff/Ntot;
	sigF[1]	= sqrt(pow(yield->getPropagatedError(*res)/FullSpectrum->getVal(),2)+pow(FullSpectrum->getPropagatedError(*res)*yield->getVal()/pow(FullSpectrum->getVal(),2),2))*mcSet.sumEntries()*FilterEff/Ntot;
	if (std::isnan(sigF[1])) sigF[1]= -99;
	RooAbsReal*  FullSpectrum = signal.createIntegral(*wspace.var("hnl_mass"));
	wspace.var("hnl_mass")->setRange("threeSigma",wspace.var("mean_ele")->getVal()-3*sigF[2],wspace.var("mean_ele")->getVal()+3*sigF[2]) ;
	RooAbsReal*  yield = signal.createIntegral(*wspace.var("hnl_mass"),RooFit::Range("sixSigma"));

	//prepare ws for signal shapes
	RooWorkspace w("w");
	RooRealVar signal_norm("signal_norm","signal_norm", 1.0 , 1.0 , 1.0 );	
	w.import(signal);
	w.import(signal_norm);
	w.var("mean_ele")->setConstant(kTRUE);
	w.var("sigma_ele")->setConstant(kTRUE);
	w.var("alpha1_cb")->setConstant(kTRUE);
	w.var("n1_cb")->setConstant(kTRUE);
	w.var("alpha2_cb")->setConstant(kTRUE);
	w.var("n2_cb")->setConstant(kTRUE);
	w.var("hnl_mass")->setMin(wspace.var("mean_ele")->getVal()-10*wspace.var("sigma_ele")->getVal()) ;
	w.var("hnl_mass")->setMax(wspace.var("mean_ele")->getVal()+10*wspace.var("sigma_ele")->getVal()) ;

	//check and save ws
	w.Print();
	w.writeToFile((wd+"/shapes_"+channel[type]+"_"+std::to_string(lxyBin)+".root").c_str());

	//save signal yields and sigmas to file	
	SigFeatures.open(wd+"/Signal_mass"+std::to_string((int)(mass))+/*"_ctau"+std::to_string((int)(ctau))+*/"_features.txt",ios_base::app );
	SigFeatures << mass << " " << ctau << " " << sigF[0] <<  " " << sigF[1] << " " << sigF[2] <<  std::endl;
	return 0;
}



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


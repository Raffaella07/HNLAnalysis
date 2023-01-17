//#include "/interface/GBRForestTools.h"
#include "../interface/BParkTools.h"
#include "../interface/Utils.h"
#include "../interface/SmallBToMuLPiClass.h"
#include "TTree.h"
//#include <ROOT/RDataFrame.hxx>
//#include <ROOT/RVec.hxx>
//#include "../../../XGBoost-FastForest/include/fastforest.h"
#include <stdint.h>
//#include "TStopwatch.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <math.h>
#include <dirent.h>
#include <string>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TTree.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBDT.h"


void FillHists(SmallBToMuLPiClass* ,std::vector<TH1D*>* ,int,double );
void genEventCount(std::string* ,int ,int* );

int main(int argc, char **argv){

	//getting all the bkg samples in input 
	int nDataset=9;
	TChain* qcd_in[nDataset];
	TChain* data_in=new TChain("Events");

		data_in->Add(("/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/SkimTuples/crab_HNL_BParking1A_000*10.root")/*.c_str()*/);




	SmallBToMuLPiClass qcd[nDataset];
	SmallBToMuLPiClass data;

	SmallBToMuLPiClass* qcd_ptr[nDataset];
	SmallBToMuLPiClass* data_ptr;
//iQCD_samples/crab_HNLMC_QCD_Pt80-120_ext_Chunk2.root
	data.Init(data_in);
	data_ptr = &data;	

	std::string lables[nDataset]={/*"QCD_#muEnriched_#hat{p}15To20",*/"QCD_#muEnriched_#hat{p}20To30","QCD_#muEnriched_#hat{p}30To50","QCD_#muEnriched_#hat{p}50To80","QCD_#muEnriched_#hat{p}80To120","QCD_#muEnriched_#hat{p}80To120_ext","QCD_#muEnriched_#hat{p}120To170","QCD_#muEnriched_#hat{p}120To170_ext","QCD_#muEnriched_#hat{p}170To300"};
	std::string paths[nDataset]={
 		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt15to20_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt20to30_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt30to50_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt50to80_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt80to120_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt80to120_ext_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt120to170_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt120to170_ext_Chunk*.root",
		"/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/Skim_QCDsamples/crab_HNLMC_QCD_Pt170to300_Chunk*.root"
	};
	double XCrossWeight[nDataset];
	double XCross[nDataset]={2.806e+06,2.514e+06,1.368e+06,3.779e+05,8.838e+04,8.838e+04,2.125e+04,2.125e+04,7.011e+03}; // croos section in fb
//	double XCross[nDataset]={8.490e+08,3.978e+08,1.070e+08,1.573e+07,2.343e+06,2.343e+06,4.080e+05,4.080e+05,1.037e+05}; // croos section in fb
//	double XCross[nDataset]={1273190000,558528000,139803000,19222500,2758420,2758420,469797,469797,117989};
	double FilterEff[nDataset]={0.00300,0.00530,0.01182,0.02276,0.03844,0.03844,0.05362,0.05362,0.073358};
	int MiniEvts[nDataset]={0};
	int nVar=11;

	std::vector<std::vector<TH1D*>> vars;
		
	vars.resize(nDataset+1);	
	std::string variables[nVar] = {"mass_{HNL+#mu}(GeV)","mass_{HNL}(GeV)","p_{T}^{HNL+#mu}","p_{T}^{HNL}","#hat{p}_{T}(GeV)","nElectron","PV_npvs","PFMvaId","charge_{HNL}","HNL_Lxy_SR","HNL_Lxy_CR"};
	int nBin[nVar]= {50, 50, 50,50,50 ,50,50 , 50, 6,10,10};
	double min[nVar]= { 0, 0.2, 0, 0, 0 , 0 , 0 ,-10,0,0,0};
	double max[nVar]= {8, 7 ,20,70,50,100,100, 8 ,3,30,30};
	double min_array[5] = {0,3,10,20,30};
	int i, j;
	double dataLumi =1000*0.776;
	std::cout <<"before generator stuff"<< vars.size()<< std::endl;
	vars.at(0).resize(nVar);	
	for(i=0;i<nVar;i++){
		if (i<(nVar-2))vars.at(0).at(i) = new TH1D(("data_"+variables[i]).c_str(),"data-BParkingD5Part1",nBin[i],min[i],max[i]);
		else vars.at(0).at(i) = new TH1D(("data_"+variables[i]).c_str(),"data-BParkingD5Part1",4,min_array);
		}
	
/*
	TH1D* Bmass[nDataset+1]; 
	TH1D* HNLmass[nDataset+1]; 
	TH1D* HNLPt[nDataset+1]; 
	TH1D* CandPt[nDataset+1]; 
*/
	genEventCount(paths,nDataset,MiniEvts);
/*	Bmass[0] = new TH1D("data_B","data-BParkingD5Part1",50,0,8);
	HNLmass[0] = new TH1D("data_hnl","data-BParkingD5Part1",50,0.2,7);
	HNLPt[0] = new TH1D("data_B","data-BParkingD5Part1",50,0,20);
	CandPt[0] = new TH1D("data_hnl","data-BParkingD5Part1",50,0,70);*/

	FillHists(data_ptr,&vars.at(0),nVar,1);
	std::ofstream QCDWeightsFile;
	QCDWeightsFile.open("QCDWeights.txt");
	std::cout << "data entries"<< vars.at(0).at(0)->GetBinContent(25) << std::endl;  

	double SumWeightsMC =0;
	//std::cout << "before dataset cycle"<< std::endl;
	  
	for (i=0;i<nDataset;i++){
	vars.at(i+1).resize(nVar);	
//		std::cout << "initialising dataset"<< std::endl;  
		qcd_in[i] = new TChain("Events");
		qcd_in[i]->Add(paths[i].c_str());
		std::cout << " dataset added to chain"<< paths[i]  << " "<<  qcd_in[i]->GetEntries() << std::endl;  
		qcd[i].Init(qcd_in[i]);
		qcd_ptr[i] = &qcd[i];
		if(i==3 || i==5 )XCrossWeight[i] =dataLumi*XCross[i]/(MiniEvts[i]+MiniEvts[i+1]);
		else XCrossWeight[i] = dataLumi* XCross[i]/MiniEvts[i];
		std::cout << XCross[i] <<" " <<  FilterEff[i] <<" " << MiniEvts[i]<< std::endl; 
		QCDWeightsFile << XCrossWeight[i] <<std::endl; 
//		if(i==4 || i==6)SumWeightsMC += (MiniEvts[i]+MiniEvts[i+1])/(XCross[i]*1000);
//		else SumWeightsMC += MiniEvts[i]/(XCross[i]*1000);
		std::cout << SumWeightsMC <<std::endl; 
		for (j=0;j<nVar;j++){
		if (j<(nVar-2))vars.at(i+1).at(j) = new TH1D((variables[j]+lables[i]).c_str(),lables[i].c_str(),nBin[j],min[j],max[j]);
		else vars.at(i+1).at(j) = new TH1D((variables[j]+lables[i]).c_str(),(variables[j]+lables[i]).c_str(),4,min_array);
		//std::cout << vars.at(i+1).at(j)->GetName() <<std::endl; 
	/*	Bmass[i+1] = new TH1D(("Bmass"+lables[i]).c_str(),lables[i].c_str(),50,0,8);
		HNLmass[i+1] = new TH1D(("HNLmass"+lables[i]).c_str(),lables[i].c_str(),50,0.2,7);
		HNLPt[i+1] = new TH1D(("HNLPt"+lables[i]).c_str(),lables[i].c_str(),50,0,20);
		CandPt[i+1] = new TH1D(("CandPt"+lables[i]).c_str(),lables[i].c_str(),50,0,70);
	*/
//	if (i==5 || i==7)FillHists(qcd_ptr[i],Bmass[i],HNLmass[i],HNLPt[i],CandPt[i],XCrossWeight[i]);
	}
	if (i==4 || i==6)FillHists(qcd_ptr[i],&vars.at(i),nVar,XCrossWeight[i-1]);
	else FillHists(qcd_ptr[i],&vars.at(i+1),nVar,XCrossWeight[i]);
		
	}
	QCDWeightsFile.close();
	std::cout <<"sum of weights" << SumWeightsMC << std::endl;
	for(i=0;i<nVar;i++){
	TH1D* tempVar[nDataset];
	double sumInt =  vars.at(1).at(i)->Integral()+vars.at(2).at(i)->Integral()+vars.at(3).at(i)->Integral()+vars.at(5).at(i)->Integral()+vars.at(7).at(i)->Integral()+vars.at(4).at(i)->Integral()+vars.at(8).at(i)->Integral(); 
	for (j=0;j<nDataset+1;j++){
		tempVar[j]=(TH1D*)vars.at(j).at(i)->Clone();
		if (j==0)tempVar[j]->Scale(1/tempVar[j]->Integral());
		else tempVar[j]->Scale(1/sumInt);
		
	}
//		if (j!=0)tempVar[j]->Scale(1/SumWeightsMC);
	
	RatioPlot(variables[i].c_str(),tempVar,nDataset+1,(variables[i]+"_DataVsQCD").c_str(),false,lables);
	}
	
/*	RatioPlot("m_{HNL} (GeV)",HNLmass,10,"HNLmass_DataVsQCD",false,lables);	
	RatioPlot("p_{T}^{HNL} (GeV)",HNLPt,10,"HNLpt_DataVsQCD",false,lables);	
	RatioPlot("p_{T}^{HNL+#mu} (GeV)",CandPt,10,"CandPt_DataVsQCD",false,lables);	
*/

	return 0;
}


void FillHists(SmallBToMuLPiClass* evt,std::vector<TH1D*>* hists,int nVars,double weight ){

	int i;
	TFile* f_puW = TFile::Open("PUweights.root");
	TH1D* PUw = (TH1D*) f_puW->Get("h");
	for(i=0;i<evt->fChain->GetEntries();i++){

		evt->fChain->GetEntry(i);
		double PUweight;
		if (weight != 1 ){ 
		for (int j =1; j<PUw->GetNbinsX();j++){
			
			if(evt->PV_npvs >= PUw->GetXaxis()->GetBinLowEdge(j) && evt->PV_npvs < PUw->GetXaxis()->GetBinUpEdge(j) )PUweight = PUw->GetBinContent(j);
	//		std::cout << PUweight << std::endl;

		}
		weight  = weight * PUweight;
		}
		if(i%100000==0)	std::cout << "on entry " << i <<" event " << evt->event << std::endl;

		if (evt->nBToMuEPi ==0) continue;
		//	if (evt.nBToKEE ==0) continue;
		if(evt->nElectron==0) continue;
		if(evt->HLT_Mu9_IP6==0 && evt->HLT_Mu12_IP6==0) continue;
		hists->at(4)->Fill(evt->TriggerMuon_pt[0],weight);
		int nHNL = evt->nBToMuEPi;
		//		if(nHNL>0) den++;
		//for (int j =0; j<evt.nBToKEE; j++){}
		//		if((int)(evt.event) ==  1132)
		//	std::cout << nHNL << " "  <<evt->nElectron << std::endl;
		//hists->at(9)->Fill(evt->Ge,weight);
		int PFEleCount=0;
		std::vector<double> PreSel_HNLpt;
		for (int j =0; j<nHNL; j++){

		double deltaR_trgMu,deltaR_trgMuPi;
		TLorentzVector HNL, TrgMu,Pi;
		HNL.SetPtEtaPhiM( evt->BToMuEPi_pt[j], evt->BToMuEPi_eta[j], evt->BToMuEPi_phi[j], evt->BToMuEPi_hnl_mass[j]);
		//	std::cout << "check 1  " <<std::endl;
		Pi.SetPtEtaPhiM( evt->BToMuEPi_pi_pt[j], evt->BToMuEPi_pi_eta[j], evt->BToMuEPi_pi_phi[j], 0.139526);
		//	std::cout << "check 2  " <<std::endl;
		TrgMu.SetPtEtaPhiM( evt->TriggerMuon_pt[0], evt->TriggerMuon_eta[0], evt->TriggerMuon_phi[0], evt->TriggerMuon_mass[0]);
		//	std::cout << "check 3  " <<std::endl;
			
		deltaR_trgMu = HNL.DeltaR(TrgMu);
		deltaR_trgMuPi = Pi.DeltaR(TrgMu);
			
			//std::cout << "here " << evt->BToMuEPi_l_isPF[j] << " " << fabs(evt->BToMuEPi_hnl_cos2D[j]) <<" " << evt->BToMuEPi_sv_prob[j] << std::endl; 
			bool preSel = evt->TriggerMuon_pt[0]>7
		 &&  fabs(evt->TriggerMuon_eta[0])<1.5
		 &&  evt->BToMuEPi_pi_pt[j]>0.7 
		 &&  fabs(evt->BToMuEPi_pi_eta[j])<2 
		 &&  fabs(evt->BToMuEPi_pi_dz[j])>0.005 
		 &&  fabs(evt->BToMuEPi_pi_dxy[j])>0.005 
		 &&  fabs(evt->BToMuEPi_pi_dxySig[j])>3 
		 &&  fabs(evt->BToMuEPi_pi_dzSig[j])>1.5 
		 &&  evt->BToMuEPi_sv_prob[j]>0.001
		 &&  fabs(evt->BToMuEPi_hnl_cos2D[j])>0.99 
		 &&  evt->BToMuEPi_sv_lxy[j]/evt->BToMuEPi_sv_lxye[j]>3
		 &&  fabs(evt->BToMuEPi_l_dxy[j])>0.008
		 &&  fabs(evt->BToMuEPi_l_dz[j])>0.008
		 &&  fabs(evt->BToMuEPi_l_dxy[j]/evt->BToMuEPi_l_dxyErr[j])>3
		 &&  fabs(evt->BToMuEPi_l_dz[j]/evt->BToMuEPi_l_dzErr[j])>1.5
		 &&  deltaR_trgMu<0.6;



		
	//		std::cout << "here 1" << std::endl; 
			if (!preSel) continue; // HNL candidates are here ranked by pt, if first passes selection, ok, else, discard event
			
		//	std::cout << "pushing back entry " << std::endl; 
				PreSel_HNLpt.push_back(evt->BToMuEPi_hnl_pt[j]);
			}				
	//		Int_t ele_idx = evt->Electron_genPartIdx[evt->BToMuEPi_sel_e_idx[j]];
	//		int trk_idx = evt->ProbeTracks_genPartIdx[evt->BToMuEPi_pi_idx[j]];
	//		bool EleMatch =  ele_idx != -1 && ele_idx <Ngen;
	//		bool TrackMatch =  ele_idx != -1 && trk_idx <Ngen;
			int theBIdx = -1;
			if(PreSel_HNLpt.size()>0){
			std::sort(PreSel_HNLpt.begin(), PreSel_HNLpt.end(), std::greater <>());
				for (int k =0; k<nHNL; k++){
					
				if(evt->BToMuEPi_hnl_pt[k]==PreSel_HNLpt.at(0)){
				
			//std::cout << "in index if" << std::endl; 
				theBIdx = k;

					
				}
				}
			}
			
			hists->at(8)->Fill(evt->BToMuEPi_hnl_charge[theBIdx],weight);
			if(evt->BToMuEPi_l_isPF[theBIdx]) PFEleCount++;
			bool ControlRegion,SignalRegion; 
	//		if (weight!=1 ) ControlRegion = evt->BToMuEPi_hnl_charge[theBIdx]!=0;
			//else 
			ControlRegion = evt->BToMuEPi_mass[theBIdx]>5.7  && evt->BToMuEPi_l_isPF[theBIdx] && evt->TriggerMuon_pt[theBIdx]>7;
			SignalRegion = evt->BToMuEPi_mass[theBIdx]>4.5 &&  evt->BToMuEPi_mass[theBIdx]<5.7  && evt->BToMuEPi_l_isPF[theBIdx] && evt->TriggerMuon_pt[theBIdx]>7;
	//			std::cout << "in for " << evt->BToMuEPi_hnl_charge[theBIdx]<< " "<<  evt->Electron_isPF[evt->BToMuEPi_sel_e_idx[theBIdx]] <<std::endl;
			if (SignalRegion) hists->at(9)->Fill(evt->BToMuEPi_sv_lxy[theBIdx],weight);
	
			if (!ControlRegion) continue;

			hists->at(0)->Fill(evt->BToMuEPi_mass[theBIdx],weight);
			hists->at(1)->Fill(evt->BToMuEPi_hnl_mass[theBIdx],weight);
	//		std::cout <<"before generator stuff"<< hists->at(1)->GetName()<< std::endl;
			hists->at(2)->Fill(evt->BToMuEPi_hnl_pt[theBIdx],weight);
			hists->at(3)->Fill(evt->BToMuEPi_pt[theBIdx],weight);
			hists->at(7)->Fill(evt->BToMuEPi_l_pfmvaId[theBIdx],weight);
			hists->at(10)->Fill(evt->BToMuEPi_sv_lxy[theBIdx],weight);
	//		std::cout <<"before generator stuff"<< hists->at(0)->GetName()<< std::endl;
	//		std::cout <<"after generator stuff"<<std::endl;
			hists->at(6)->Fill(evt->PV_npvs,weight);
			}
			//hists->at(5)->Fill(PFEleCount,weight);
		
		
	}
	//return 0;


void genEventCount(std::string* paths,int nDatasets, int* TotCounter){

   long int evtCounter;
 //  int TotCounter[nDatasets]={0};
   int j;

   for(j=0;j<nDatasets;j++){

	   TChain* data=new TChain("Runs");

	   data->Add(paths[j].c_str()); 

   	//   TBranch *EventCount  = data->GetBranch("genEventCount");
	   data->SetBranchAddress("genEventCount",&evtCounter);
   	   int i;

  	   for(i=0;i<data->GetEntries();i++){
//		std::cout <<"running on " << data->GetEntries() << " tree" << std::endl;
		data->GetEntry(i);
		
	//	std::cout <<"running on " << i << " tree" << std::endl;
	//	std::cout <<"single counter" << evtCounter  << std::endl;
		TotCounter[j]+=evtCounter;
   }
	std::cout << TotCounter[j] << std::endl;
	delete data;
//	delete EventCount;
	}

   

}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

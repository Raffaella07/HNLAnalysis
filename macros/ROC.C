#include "setStyle.C"

void ROC(std::string signal, std::string background, std::string RocBDT, std::string RocBDTAll,std::string testBDTsig, std::string testBDTbkg){

	int i;
	setStyle();
	TFile* fSig = TFile::Open(signal.c_str());
	TFile* fBdt = TFile::Open(RocBDT.c_str());
	TFile* fBdtAll = TFile::Open(RocBDTAll.c_str());
	TFile* fBDTSig = TFile::Open(testBDTsig.c_str());
	TFile* fBDTBkg = TFile::Open(testBDTbkg.c_str());
	
	TH1::AddDirectory(kFALSE);
	//ROC BDT scelto 
	TH1D * bdt = (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_effBvsS")->Clone(); // test sample (Control region Q !=0 data)
	TH1D * bdt_train = (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_trainingEffBvsS")->Clone(); // training sample (Control region Q!=0) 
	//ROC BDT all vars
	TH1D * bdtAll = (TH1D*)fBdtAll->Get("dataset/Method_BDT/BDT/MVA_BDT_effBvsS")->Clone();
	TH1D * sigBDT= (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_S")->Clone();
	TH1D * bkgBDT= (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_B")->Clone();
	TH1D * sig = (TH1D*)fSig->Get("sig")->Clone();
	TFile* fBkg = TFile::Open(background.c_str());
	TH1D * bkg = (TH1D*)fBkg->Get("sig")->Clone();
	TH1D * BDTbkg = (TH1D*)fBDTBkg->Get("sig")->Clone();
	TH1D * BDTsig = (TH1D*)fBDTSig->Get("sig")->Clone();

	std::vector <double> SigEff;
	std::vector <double> BkgReg;
	std::vector <double> SigEffBDT;
	std::vector <double> BkgRegBDT;
	int Nbin = sig->GetXaxis()->GetNbins();
	for(i=1;i<Nbin;i++){
	
//		if(i%4==0){
			
			SigEff.push_back(sig->Integral(i,Nbin));	
			BkgReg.push_back(bkg->Integral(i,Nbin));	
			std::cout << i << " " << sig->Integral(i,Nbin) << " " << bkg->Integral(i,Nbin)<< " " <<sig->GetXaxis()->GetBinCenter(i) << std::endl;

//		}


	}	
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	for(i=5;i<BDTsig->GetXaxis()->GetNbins()-5;i++){
	
//		if(i%4==0){
			
			SigEffBDT.push_back(BDTsig->Integral(i,Nbin));	
			BkgRegBDT.push_back(BDTbkg->Integral(i,Nbin));	
			std::cout << i << " " << BDTsig->Integral(i,Nbin)/BDTsig->Integral()<< " " << BDTbkg->Integral(i,Nbin)/BDTbkg->Integral() << " " <<BDTsig->GetXaxis()->GetBinCenter(i) << std::endl;

//		}


	}	
	TGraph* roc = new TGraph((int)(Nbin-1),SigEff.data(),BkgReg.data());
	TGraph* roc_testBDT = new TGraph((int)(BDTsig->GetXaxis()->GetNbins()-11),SigEffBDT.data(),BkgRegBDT.data());
	
	TCanvas* rcanva = new TCanvas("rcanva","rcanva",600,800);
	TH2D* plotter = new TH2D("plotter","plotter",10,0.5,1.00,10,0.00000001,1);
	gStyle->SetLegendFont(42);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendTextSize(0.03);
	TLegend* l = new TLegend(0.2,0.15,0.4,0.3);
	l->SetBorderSize(0);
	rcanva->SetLogy();
	gStyle->SetOptStat(0);
	roc->SetLineColor(kAzure+7);
	roc->SetLineWidth(3);
	plotter->GetXaxis()->SetTitle("signal efficiency");
	//plotter->GetXaxis()->SetNdivisions(-10,kTRUE);
	plotter->GetYaxis()->SetTitle("bkg efficiency");
	bdt->SetLineColor(kRed);
	bdt->SetLineWidth(3);
	bdtAll->SetLineColor(kBlue);
	bdtAll->SetLineWidth(3);
	bdt_train->SetLineColor(kOrange+7);
	bdt_train->SetLineWidth(3);
	roc_testBDT->SetLineColor(kMagenta);
	roc_testBDT->SetLineWidth(3);
	l->AddEntry(roc,"likelihood","l");
	l->AddEntry(bdt,"BDT (12 variables set, test sample Q!=0)","l");
	l->AddEntry(bdt_train,"BDT (12 variables set, training sample Q!=0)","l");
	l->AddEntry(bdtAll,"BDT (all variables, test sample Q!=0)","l");
	l->AddEntry(roc_testBDT,"BDT (test sample Q ==0)","l");
	
	plotter->Draw();
	bdtAll->Draw("sameL");
	bdt->Draw("sameL");
	bdt_train->Draw("sameL");
	roc->Draw("same");
	roc_testBDT->Draw("sameL");
	l->Draw();
	
//	rcanva->BuildLegend();


	rcanva->SaveAs("ROC.pdf");


}

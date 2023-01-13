#include "setStyle.C"

//void ROC(std::string signal, std::string background, std::string RocBDT, std::string RocBDTAll,std::string testBDTsig, std::string testBDTbkg){}
void ROC(int nBuiltROCs, int nHistoROCs,std::string BuiltROCs, std::string HistoROCs, bool debug, std::string nametag){

	int i;
	setStyle();
	std::string ROCFiles[nBuiltROCs];
	std::string HistoFiles[nHistoROCs][2];
	std::string HistoROCDoc[nHistoROCs];
	std::string BuiltROCDoc[nBuiltROCs];
	
	//retrieving filenames for already built ROC curves
	if (nBuiltROCs!=0){
	 ifstream file (BuiltROCs.c_str());
         if (file.is_open())
         { int i =0;
	   std::string line;
           while ( getline (file,line) )
          {
            std::stringstream ss(line);
            std::string r, doc;
            if (ss >> r >> doc )
             { 
          ROCFiles[i] = r;    
	  BuiltROCDoc[i] = doc;
	  i++;
          }
	}
          file.close();
        }
	}
	//retrieving filenames for histos for ROC curves
	if (nHistoROCs!=0){
	ifstream myfile (HistoROCs.c_str());
        if (myfile.is_open())
        {
	  int i =0;
	  std::string line;
          while ( getline (myfile,line) )
         {
	    
	std::cout<< line << std::endl;
            std::stringstream ss(line);
            std::string s, b, doc;
            if (ss >> s >> b >> doc )
             { 
	       HistoFiles[i][0] = s;	
	       HistoFiles[i][1] = b;	
	       HistoROCDoc[i] = doc;

	       i++;
	
                 }
	}
         myfile.close();
        }
	}
	//initialising containers
	TH1::AddDirectory(kFALSE);
	TH1D* builtROCs[nBuiltROCs];
	TH1D* histoROCs[nHistoROCs][2];
	std::vector<double> SigEff[nHistoROCs];
	std::vector<double> BkgEff[nHistoROCs];
	TGraph* roc[nBuiltROCs]; 	

	//plotting 
	TCanvas* rcanva = new TCanvas("rcanva","rcanva",600,800);
	TH2D* plotter = new TH2D("plotter","plotter",10,0.5,1.00,10,0.00000001,1);
	gStyle->SetLegendFont(42);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendTextSize(0.03);
	TLegend* l = new TLegend(0.2,0.17,0.4,0.3);
	l->SetBorderSize(0);
	rcanva->SetLogy();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(91);
	plotter->GetXaxis()->SetTitle("signal efficiency");
	plotter->GetYaxis()->SetTitle("bkg efficiency");
	plotter->Draw();


	for (int i =0; i < nHistoROCs; i++){
				
  		TFile* fsig = TFile::Open(HistoFiles[i][0].c_str());
		histoROCs[i][0] = (TH1D*)fsig->Get("sig")->Clone();
  		TFile* fbkg = TFile::Open(HistoFiles[i][1].c_str());
		histoROCs[i][1] = (TH1D*)fbkg->Get("sig")->Clone();
		int Nbin = histoROCs[i][0]->GetXaxis()->GetNbins();
		for(int j=1;j<Nbin;j++){
	
//		if(i%4==0){
			
			SigEff[i].push_back(histoROCs[i][0]->Integral(j,Nbin));	
			BkgEff[i].push_back(histoROCs[i][1]->Integral(j,Nbin));	
			if(debug) std::cout << j << " " <<histoROCs[i][0]->Integral(j,Nbin) << " " << histoROCs[i][1]->Integral(j,Nbin)<< " " <<histoROCs[i][1]->GetXaxis()->GetBinCenter(j) << std::endl;
//			if(debug) std::cout << j << " " <<SigEff[i].back() << " " << BkgEff[i].back() << " " <<histoROCs[i][1]->GetXaxis()->GetBinCenter(j) << std::endl;

//		}


					}
		roc[i] = new TGraph((int)(Nbin-1),SigEff[i].data(),BkgEff[i].data());
		roc[i]->SetLineWidth(2);	
		roc[i]->Draw("SAMEL PLC");
		l->AddEntry(roc[i],HistoROCDoc[i].c_str(),"l");

	}		
	for (int i =0; i < nBuiltROCs; i++){
				
  		TFile* f = TFile::Open(ROCFiles[i].c_str());
		builtROCs[i] = (TH1D*)f->Get("dataset/Method_BDT/BDT/MVA_BDT_effBvsS")->Clone(); // test sample (Control region Q !=0 data)
		builtROCs[i]->SetLineWidth(2);	
		builtROCs[i]->Draw("SAMEL PLC");
		l->AddEntry(builtROCs[i],BuiltROCDoc[i].c_str());
	}
	l->Draw("same");
	rcanva->SaveAs(("ROC"+nametag+".pdf").c_str());
	delete rcanva;
	return 0;
/*	TFile* fSig = TFile::Open(signal.c_str());
	TFile* fBdt = TFile::Open(RocBDT.c_str());
	TFile* fBdtAll = TFile::Open(RocBDTAll.c_str());
	TFile* fBDTSig = TFile::Open(testBDTsig.c_str());
	TFile* fBDTBkg = TFile::Open(testBDTbkg.c_str());
	
	//ROC BDT scelto 
	TH1D * bdt = (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_effBvsS")->Clone(); // test sample (Control region Q !=0 data)
	TH1D * bdt_train = (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_trainingEffBvsS")->Clone(); // training sample (Control region Q!=0) 
	TH1D * sigBDT= (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_S")->Clone();
	TH1D * bkgBDT= (TH1D*)fBdt->Get("dataset/Method_BDT/BDT/MVA_BDT_B")->Clone();
	//ROC BDT all vars
	TH1D * bdtAll = (TH1D*)fBdtAll->Get("dataset/Method_BDT/BDT/MVA_BDT_effBvsS")->Clone();
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
			std::cout << i << " " << BDTsig->GetBinContent(i)/BDTsig->Integral()<< " " << BDTbkg->GetBinContent(i)/BDTbkg->Integral()<< " " << BDTsig->Integral(i,Nbin)/BDTsig->Integral()<< " " << BDTbkg->Integral(i,Nbin)/BDTbkg->Integral() << " " <<  BDTsig->GetXaxis()->GetBinCenter(i) << std::endl;

//		}


	}	
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	for(i=1;i<sigBDT->GetXaxis()->GetNbins();i++){
	
//		if(i%4==0){
			
			std::cout << i << " " <<sigBDT->GetXaxis()->GetBinCenter(i) << " " << sigBDT->GetBinContent(i)/sigBDT->Integral()<< " " << bkgBDT->GetBinContent(i)/bkgBDT->Integral() << " " << sigBDT->Integral(i,Nbin)/sigBDT->Integral()<< " " << bkgBDT->Integral(i,Nbin)/bkgBDT->Integral()  << std::endl;

//		}


	}	
	TGraph* roc = new TGraph((int)(Nbin-1),SigEff.data(),BkgReg.data());
	TGraph* roc_testBDT = new TGraph((int)(BDTsig->GetXaxis()->GetNbins()-11),SigEffBDT.data(),BkgRegBDT.data());
	
	TCanvas* rcanva = new TCanvas("rcanva","rcanva",600,800);
	TH2D* plotter = new TH2D("plotter","plotter",10,0.5,1.00,10,0.00000001,1);
	gStyle->SetLegendFont(42);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendTextSize(0.03);
	TLegend* l = new TLegend(0.2,0.17,0.4,0.3);
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
	l->AddEntry(roc,"likelihood - 6 uncorrelated vars set","l");
	l->AddEntry(bdtAll,"BDT (12 variables set, test sample Q!=0)","l");
	l->AddEntry(bdt_train,"BDT (12 variables set, training sample Q!=0)","l");
//	l->AddEntry(bdt,"BDT (all variables, test sample Q!=0)","l");
	l->AddEntry(roc_testBDT,"BDT (12 variables set,test sample Q ==0)","l");
	
	plotter->Draw();
//	bdtAll->Draw("sameL");
	bdt->Draw("sameL");
	bdt_train->Draw("sameL");
	roc->Draw("same");
	roc_testBDT->Draw("sameL");
	l->Draw();
	
//	rcanva->BuildLegend();


	rcanva->SaveAs("ROC.pdf");
*/

}

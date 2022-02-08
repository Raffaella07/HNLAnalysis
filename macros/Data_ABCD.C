
void setStyle() {


	// set the TStyle
	TStyle* style = new TStyle("DrawBaseStyle", "");
	style->SetCanvasColor(0);
	style->SetPadColor(0);
	style->SetFrameFillColor(0);
	style->SetStatColor(0);
	style->SetOptStat(0);
	style->SetOptFit(0);
	style->SetTitleFillColor(0);
	style->SetCanvasBorderMode(0);
	style->SetPadBorderMode(0);
	style->SetFrameBorderMode(0);

	style->SetPadLeftMargin(0.12);
	style->cd();
	// For the canvas:
	style->SetCanvasBorderMode(0);
	style->SetCanvasColor(kWhite);
	style->SetCanvasDefH(600); //Height of canvas
	style->SetCanvasDefW(600); //Width of canvas
	style->SetCanvasDefX(0); //POsition on screen
	style->SetCanvasDefY(0);
	// For the Pad:
	style->SetPadBorderMode(0);
	style->SetPadColor(kWhite);
	style->SetPadGridX(false);
	style->SetPadGridY(false);
	style->SetGridColor(0);
	style->SetGridStyle(3);
	style->SetGridWidth(1);
	// For the frame:
	style->SetFrameBorderMode(0);
	style->SetFrameBorderSize(1);
	style->SetFrameFillColor(0);
	style->SetFrameFillStyle(0);
	style->SetFrameLineColor(1);
	style->SetFrameLineStyle(1);
	style->SetFrameLineWidth(1);
	// Margins:
	style->SetPadTopMargin(0.10);
	style->SetPadBottomMargin(0.14);//0.13);
	style->SetPadLeftMargin(0.16);//0.16);
	style->SetPadRightMargin(0.2);//0.02);
	// For the Global title:
	style->SetOptTitle(0);
	style->SetTitleFont(42);
	style->SetTitleColor(1);
	style->SetTitleTextColor(1);
	style->SetTitleFillColor(10);
	style->SetTitleFontSize(0.05);
	// For the axis titles:
	style->SetTitleColor(1, "XYZ");
	style->SetTitleFont(42, "XYZ");
	style->SetTitleSize(0.05, "XYZ");
	style->SetTitleXOffset(1.15);//0.9);
	style->SetTitleYOffset(1.5); // => 1.15 if exponents
	// For the axis labels:
	style->SetLabelColor(1, "XYZ");
	style->SetLabelFont(42, "XYZ");
	style->SetLabelOffset(0.007, "XYZ");
	style->SetLabelSize(0.045, "XYZ");
	// For the axis:
	style->SetAxisColor(1, "XYZ");
	style->SetStripDecimals(kTRUE);
	style->SetTickLength(0.03, "XYZ");
	style->SetNdivisions(510, "XYZ");
	style->SetPadTickX(1); // To get tick marks on the opposite side of the frame
	style->SetPadTickY(1);
	// for histograms:
	style->SetHistLineColor(1);
	// for the pallete
	Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red  [5] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue [5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(5, stops, red, green, blue, 100);
	style->SetNumberContours(100);

	style->cd();

}
void RatioPlot(std::string titlestring,TH1D* hist, TH1D* histMock, const char* filename,bool log,double plotScale){

	std::cout << "in function "  << std::endl;
	setStyle();
//	TH1D* h[nHist];
	TCanvas* canvas = new TCanvas("canva","canva",650,700);
	TPad *p_graph = new TPad("p_graph","p_graph",0,0.2,1,1); //creazione ratio plot
	TH1D* empty = new TH1D("empty","empty",10,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
	std::string PNGPATH = "/eos/home-r/ratramon/www/HNL/";
	TLatex l;
	THStack s;
	gStyle->SetPalette(kPastel);
	Color_t sel_colors[5]= {kAzure+7,kRed-7,kOrange+7,kGreen-3,kMagenta};
	double x_lable,y_lable;
	y_lable= -1;
	x_lable= 0.60;

	//double integral = 1.0*sum->Integral(sum->GetXaxis()->FindBin(4.5),sum->GetXaxis()->FindBin(5.7))/sum->Integral(sum->GetXaxis()->FindBin(5.7),50);
//	std::cout << "mass integrals " << integral <<  std::endl;
	
//	p6_graph->SetTopMargin(0.3);
	TLegend* leg= new TLegend(0.4,0.85-0.25,0.4+0.3,0.85);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
//	sum->SetNameTitle("h", "my normalized stack");
//	sum->Scale(1./sum->Integral());
//	sum->SetLineWidth(3);
//	sum->SetLineColor(kGray+3);
//	sum->SetFillStyle(0);
//	std::cout << "compare entries " << h1->GetEntries() << "reco "  << h2->GetEntries() << std::endl;
	empty->GetXaxis()->SetMaxDigits(2);
	if (log){ p_graph->SetLogy();
	//	empty->GetYaxis()->SetRangeUser(1,std::max(h1->GetMaximum()*2.5,h2->GetMaximum()*2.5));
		empty->GetYaxis()->SetRangeUser(0.0001,y_lable);
	}else {
		empty->GetYaxis()->SetRangeUser(0,std::max(y_lable,hist->GetMaximum())*plotScale);
		
	}	

	l.SetTextSize(0.04);
	l.SetTextAlign(13);
//	if(lables != NULL){
//	l.SetTextColor(kBlack);
//	l.DrawLatex(h[0]->GetXaxis()->GetXmin()+0.05*(h[0]->GetXaxis()->GetXmax()-h[0]->GetXaxis()->GetXmin()),y_lable,lable.c_str());
//	}
	canvas->cd();
	p_graph->SetBottomMargin(0.03);
	empty->GetXaxis()->SetTickLength(0.06);
	p_graph->SetTicks();
	//p_graph->SetTickSize();
//	p_graph->SetGridy();
	p_graph->Draw();
	p_graph->cd(); 
//	empty->GetXaxis()->SetTitle(titlestring.c_str());
	empty->GetYaxis()->SetTitle("entries");
	empty->GetYaxis()->SetTitleSize(0.2);
	empty->GetXaxis()->SetLabelSize(0);
	empty->GetYaxis()->SetLabelSize(0.03);
	empty->Draw("");
	hist->SetLineWidth(2);
	hist->SetMarkerStyle(8);
	hist->SetLineColor(sel_colors[0]);
	hist->SetMarkerColor(sel_colors[0]);
	hist->SetFillStyle(3003);
	histMock->SetLineWidth(2);
	histMock->SetMarkerStyle(8);
	histMock->SetLineColor(sel_colors[1]);
	histMock->SetMarkerColor(sel_colors[1]);
	histMock->SetFillStyle(0);
	//s.SetLneStyle(0);
//	s.Scale(1/s.Integral())
	hist->Draw("EPsame");
	histMock->Draw("EPsame");
//	if(nHist>2)sum->Draw("HISTsame");
// I suppose you have a vector of entities you want to draw here
// for (int i=0;i<cols.GetSize() || i >= myObjects.size() ;i++) {
//   myObjects[i].SetMarkerColor(cols.At(i));
//   }
	leg->AddEntry(hist,"Signal region QCD MC","lep");
 	leg->AddEntry(histMock,"BParking 1A, 0.76 fb^{-1}","lep");
	
//	l.DrawLatex(h[i]->GetXaxis()->GetXmin()+x_lable*(h[i]->GetXaxis()->GetXmax()-h[i]->GetXaxis()->GetXmin()),.95*y_lable-i*y_lable/20,(std::string(h[i]->GetTitle())+" "+std::to_string((int)h[i]->GetEntries())).c_str());
	//leg->AddEntry(sum,("Total Bkg "+std::to_string((int)hist->GetEntries())).c_str(),"f");
	leg->Draw();
//	s.Draw("HISTsame");
/*	l.SetTextColor(kRed-6);
	l.DrawLatex(h1->GetXaxis()->GetXmin()+0.55*(h1->GetXaxis()->GetXmax()-h1->GetXaxis()->GetXmin()),.98*h1->GetMaximum(),(std::string("Signal ")+std::to_string((int)h1->GetEntries())).c_str());
	l.SetTextColor(kAzure+7);
	l.DrawLatex(h1->GetXaxis()->GetXmin()+0.55*(h1->GetXaxis()->GetXmax()-h1->GetXaxis()->GetXmin()),.88*h1->GetMaximum(),"Combinatorial ");
	l.DrawLatex(h1->GetXaxis()->GetXmin()+0.55*(h1->GetXaxis()->GetXmax()-h1->GetXaxis()->GetXmin()),.82*h1->GetMaximum(),(std::string("bkg ")+std::to_string((int)h2->GetEntries())).c_str());
//	if (lable ) lables1D(canvas,h1);*/
	
	std::cout << "clones for ratioplots "  << std::endl;
	TH1D* sigCopy =  dynamic_cast<TH1D *>(hist->Clone()); 
	TH1D* TotBkgCopy =  dynamic_cast<TH1D *>(histMock->Clone());
	
	std::cout << "after clones for ratioplots "  << std::endl;
	sigCopy->Divide(TotBkgCopy);
	sigCopy->SetName("hist");
	TPad *p6_graph = new TPad("p6_graph","p6_graph",0,0,1,0.2); //creazione ratio plot

	p6_graph->SetTicks();
	sigCopy->GetYaxis()->SetRangeUser(0.3,2.5);
	sigCopy->GetXaxis()->SetTitle(titlestring.c_str());

	sigCopy->GetYaxis()->SetTitle("Data / Bkg");

	sigCopy->GetYaxis()->SetLabelSize(0.1);
	sigCopy->GetYaxis()->SetNdivisions(5);
	sigCopy->GetYaxis()->SetTitleSize(0.15);
	sigCopy->GetYaxis()->SetTitleOffset(0.6);
	sigCopy->GetYaxis()->SetLabelOffset(0.03);

	sigCopy->GetXaxis()->SetLabelSize(0.1);
	sigCopy->GetXaxis()->SetTitleSize(0.15);

	sigCopy->SetMarkerStyle(8);
	sigCopy->SetMarkerSize(.8);
	sigCopy->SetMarkerColor(sel_colors[0]);

	p6_graph->SetTopMargin(0);
	p6_graph->SetBottomMargin(0.5);
	p6_graph->SetTickx();
	p6_graph->SetTicky();
	p6_graph->SetGridy();
	canvas->cd();
	//setStyle();
	p6_graph->Draw("same");
	p6_graph->cd(); 
	sigCopy->GetXaxis()->SetTickLength(0.06);
	sigCopy->Draw("EP");
	sigCopy->GetYaxis()->Draw("same");


	std::cout << " ratioplots  draw"  << std::endl;
	canvas->SaveAs((std::string(filename)+".pdf").c_str());
//	if(nHist<5)sigCopy->SaveAs((std::string(filename)+".root").c_str());
//	canvas->SaveAs((PNGPATH+std::string(filename)+".png").c_str());
	//canvas->Clear();

	std::cout << "after integrals "  << std::endl;
	delete empty;
	delete canvas;
}

void SavePlot (std::string titlestring, TH1D * histo, const char * filename, bool log, TF1* fit, bool lable){

	TCanvas* canvas = new TCanvas("canva","canva",600,550);
	//TFile* histos = new TFile("gausslog.root","UPDATE","",0);
//	std::string PNGPATH = "/eos/home-r/ratramon/www/";
//	std::string PDFPATH = "../plots/";
	//const char * temptitle = titlestring.c_str();
	setStyle();
	if (log) canvas->SetLogy();
	histo->GetXaxis()->SetMaxDigits(2);
	histo->SetLineWidth(1);
	histo->SetMarkerStyle(8);
	histo->SetMarkerSize(1);
	histo->SetLineColor(kAzure+7);
	histo->SetMarkerColor(kAzure+7);
	histo->GetXaxis()->SetTitle(titlestring.c_str());
	histo->GetYaxis()->SetTitle("entries");
	std::cout << "axis title" << titlestring.c_str() << std::endl;
	histo->Draw("PE1");
	if(fit != NULL) fit->Draw("samePE1");
//	if(lable)lables1D(canvas,histo)
	canvas->SaveAs((std::string(filename)+".pdf").c_str());
//	canvas->SaveAs((PNGPATH+std::string(filename)+".png").c_str());
	canvas->Clear();
//	histos->Close();


	delete canvas;
}

void AddQCDToHist(int SetIdx,double QCDweight,TChain* tree, std::string filter, std::string leaf,TH1D * hist, std::string name, float lumi, float mass,double sigma, bool window, float NsigmaUp, float NsigmaDown){
//void AddQCDToHist(int SetIdx,double QCDweight,ROOT::RDataFrame d, std::string filter, std::string leaf,ROOT::RDF::TH1DModel h_model,TH1D * hist, std::string name)
	//ROOT::EnableImplicitMT(4);
	if(SetIdx ==0){
			tree->Draw((leaf+">>"+std::string(hist->GetName())).c_str(),filter.c_str());
//			d.Filter(filter.c_str()).Histo1D(h_model,leaf.c_str()).GetValue().Copy(*hist);
			hist->Scale(QCDweight*lumi);		
			hist->SetTitle(name.c_str());		
			hist->SetName(name.c_str());		
	}else{
			TH1D * temp;
			//std::cout << "mass " << mass << " sigma " << sigma << std::endl;
			if(window) temp = new TH1D("temp","temp",1,mass-NsigmaDown*sigma,mass+NsigmaUp*sigma);
			else temp = new TH1D("temp","temp",25,2.5,3.5);
			tree->Draw((leaf+">>"+std::string(temp->GetName())).c_str(),filter.c_str());
		//	d.Filter(filter.c_str()).Histo1D(h_model,leaf.c_str()).GetValue().Copy(*temp);
			temp->Scale(QCDweight*lumi);		
			temp->SetTitle(name.c_str());		
			temp->SetName(name.c_str());		
			hist->Add(temp);



		}
}
int Data_ABCD(std::string path,int SameSign, bool SuperQCD,bool MassWindow,double mass,double ctau,int channel,std::string wd, std::string selection, bool unblinded){

	gROOT->SetBatch(kTRUE);	
	std::cout << wd+"/BkgYields_.txt"<< std::endl;
	const int nDataset = 9;
	const int LxyBin = 4;
	const int Channels = 4;
	const int ABCD = 6;
	int i;
	double BMassUpLimit = 5.7;
	double BMassLowLimit = 4.5;
	std::string ChannelLable[Channels]= {"Muon","PF","LowPt","Track"};
	std::string cutX = "(hnl_vtxProb >0.05)";
	std::string cutY = "(hnl_cos2D >0.996)";
	std::string RegCut[ABCD];
	std::string RegLable[ABCD] = {"AMock","B","C","D","TF"};
	
	RegCut[0] = cutX + " && " + cutY; //A
	RegCut[1] = cutX + " && !" + cutY; //B
	RegCut[2] = "!"+cutX + " && !" + cutY; //C
	RegCut[3] = "!"+cutX + " && " + cutY; //D
	RegCut[4] = "hnl_charge==0"; //Global TF on A 
	RegCut[5] = " "; // recomputed A 
	//ROOT::RDF::TH1DModel *  h_models[2];
//	const int nRegion,nChannel, nLxyBin;
	const int nRegion = 5; //B,C,D with A blinded and mocked
	const int nChannel = 4;
	const int nLxyBin = 3;
//	std::string lxy_cuts[4]= {"hnl_lxy<3", "hnl_lxy>3 && hnl_lxy<10", " hnl_lxy<20 && hnl_lxy>10","hnl_lxy>20"};

	std::string lxy_cuts[3]= {"hnl_lxy<1", "hnl_lxy>1 && hnl_lxy<5", "hnl_lxy>5"};
	double sigma[nChannel][nLxyBin];
	std::ifstream sigFeatures;
	sigFeatures.open((wd+"/Signal_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt").c_str());
	std::cout << wd+"/Signal_Mass"+std::to_string((int)(mass))+"_ctau"+std::to_string((int)(ctau))+"_features.txt" << std::endl;
	if (sigFeatures.is_open())
  	{
		int ch,lxy;
		double sig,sigErr,Sigma;
		while(sigFeatures>>ch>>lxy>>sig>>sigErr>>Sigma){
		if (Sigma != -99)sigma[ch][lxy]= Sigma;
		else sigma[ch][lxy] = 0.07; 
		std::cout << "sigma" << ch << " " << lxy << " " << sigma[ch][lxy] << std::endl;

		}
	}
	TChain* c = new TChain("Events");
	std::string datapath[10] = {path+"/Parking1D0_LowPt_Ldata",path+"/Parking1D1_LowPt_Ldata",path+"/Parking1D2_LowPt_Ldata",path+"/Parking1D3_LowPt_Ldata",path+"/Parking1D4_LowPt_Ldata",path+"/Parking1D5_LowPt_Ldata",path+"/Parking1D6_LowPt_Ldata",path+"/Parking1D7_LowPt_Ldata",path+"/Parking1D8_LowPt_Ldata",path+"/Parking1D9_LowPt"};///,"~/Analysis/data/HNLFlatTuples/Parking1A5","~/Analysis/data/HNLFlatTuples/Parking1A6";
//	c->Add((datapath[0]+"/HNLFlat_*7.root").c_str());
	c->Add((datapath[0]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[1]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[2]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[3]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[4]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[5]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[6]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[7]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[8]+"/*.root").c_str());
	std::cout << c->GetEntries() << std::endl;
	c->Add((datapath[9]+"/*.root").c_str());
/*	std::cout << c->GetEntries() << std::endl;

	std::cout << c->GetEntries() << std::endl;
*/	
	std::string INPATH = "/cmshome/ratramon/Analysis/data/HNLFlatTuples/";
	//std::string InDataset[nDataset] = {"QCD_Pt15_20","QCD_Pt20_30","QCD_Pt30_50","QCD_Pt50_80","QCD_Pt80_120","QCD_Pt80_120_ext","QCD_Pt120_170","QCD_Pt120_170_ext","QCD_Pt170_300"};
	std::string InDataset[nDataset] = {"Pt-15to20_noTrk","Pt-20to30_noTrk","Pt-30to50_noTrk","Pt-50to80_noTrk","Pt-80to120_noTrk","Pt-80to120_ext_noTrk","Pt-120to170_noTrk","Pt-120to170_ext_noTrk","Pt-170to300_noTrk"};
	double QCDweights[nDataset];
	std::ifstream infile("/cmshome/ratramon/Analysis/plugins/QCDWeights.txt");
	if (infile.is_open()){
		i =0;
		while ( !infile.eof()){
            
		infile>>QCDweights[i];
		i++;
		}
	}

	std::vector<TChain*> tree;
	for (i=0;i<nDataset;i++){
	        std::cout << INPATH+InDataset[i]+"/HNLFlat*.root" << std::endl;
		TChain* ctemp = new TChain("Events");
		ctemp->Add((INPATH+InDataset[i]+"/HNLFlat*.root").c_str());
		//ROOT::RDataFrame dtemp("Events", (INPATH+InDataset[i]+"/HNLFlat*.root").c_str());	
		std::cout << QCDweights[i] << std::endl;
		//d.push_back(dtemp);
		tree.push_back(ctemp);
		
	}
	TH1D* hnlMass[nRegion][nChannel][nLxyBin];
	TH1D* hnlMass_MC[nChannel][nLxyBin];
	std::string sign;
	std::string signSel;
	int j =  channel;
	
	if (SameSign==0){
	sign = "";
	signSel = "&&  hnl_charge==0 && B_mass<7";
	}
	else if (SameSign==1){
	sign = "SS";
	signSel = "&& LepQProd>0 && hnl_charge==0 && B_mass<7 ";
	}
	else if (SameSign==2){
	sign = "OS";
	signSel= " &&  LepQProd<0 && hnl_charge==0 && B_mass<7 ";
	}
	double NsigmaUp, NsigmaDown;
	for(int i = 0; i< nRegion; i++){
			if (j==0){

			NsigmaUp = 2;
			NsigmaDown = 2;
	
			}else{


			NsigmaUp = 2;
			NsigmaDown = 2;

			}
			for(int k = 0; k< nLxyBin; k++){
			std::string Masstitle = ChannelLable[j]+"HNLMass_LxyBin"+std::to_string(k)+"_"+"Region_"+RegLable[i]+"_"+sign;
			std::string hName = RegLable[i]+"_"+ChannelLable[j]+"_LxyBin"+std::to_string(k)+sign;
			if (MassWindow)hnlMass[i][j][k]=new TH1D(("window_"+Masstitle).c_str(),("window_"+Masstitle).c_str(),1,mass-NsigmaDown*sigma[j][k],mass+NsigmaUp*sigma[j][k]);
			else hnlMass[i][j][k]=new TH1D(("window_"+Masstitle).c_str(),("window_"+Masstitle).c_str(),25,2.5,3.5);
			std::cout << "mass " << mass << " sigma " << sigma[j][k] << " interval " <<  mass-2*sigma[j][k]<< std::endl;
			std::string Massfilter;
			//if(j==0) Massfilter =" Type =="+std::to_string(j)+" && LxyBin =="+std::to_string(lxy_cuts[k])+" && "+RegCut[i]+signSel+" && fabs(hnl_l_mvaId)==1"+selection;
			if( i==0 || i ==3) 
			Massfilter =" Type =="+std::to_string(j)+" && "+/*std::to_string*/(lxy_cuts[k])+" && "+RegCut[i]+signSel+selection;
			else if(i==1 || i ==2 ) Massfilter =" Type =="+std::to_string(j)+" && "+/*std::to_string*/(lxy_cuts[k])+" && "+RegCut[i]+signSel;
			
			if(j==0) Massfilter +=" && fabs(hnl_l_mvaId)==1";
			std::cout << Massfilter<< std::endl;
	//		if(i!=0){
			c->Draw(("hnl_mass>>"+("window_"+Masstitle)).c_str(),Massfilter.c_str());
		//	std::cout << "channel" << j  << "lxybin"<<k  << " " << mass-3*sigma[j][k] << " " <<mass+3*sigma[j][k] <<" "  <<hnlMass[i][j][k]->Integral()<<std::endl;
	//		}
			if (SuperQCD && i ==0){
			for (int id=0;id<nDataset;id++){
		//	std::cout << "MC"<< std::endl;
			
			if (id==0){
			if(MassWindow)hnlMass_MC[j][k]=new TH1D(Masstitle.c_str(),Masstitle.c_str(),1,mass-NsigmaDown*sigma[j][k],mass+NsigmaUp*sigma[j][k]);
			else hnlMass_MC[j][k]=new TH1D(Masstitle.c_str(),Masstitle.c_str(),25,2.5,3.5);
			}
		//	std::cout <<"ch" << j << "lxybin" << k <<  "mass " << mass << " sigma " << sigma[j][k] << " interval " <<  mass-2*sigma[j][k]<< std::endl;
			AddQCDToHist(id,QCDweights[id],tree.at(id),Massfilter,"hnl_mass",hnlMass_MC[j][k],Masstitle, 0.76, mass, sigma[j][k],MassWindow,NsigmaUp, NsigmaDown);
			
			}
			}
			}
		}
	

		if (unblinded){
		for(int nBin = 1; nBin< hnlMass[0][j][0]->GetXaxis()->GetNbins()+1; nBin++){
			for(int k = 0; k< nLxyBin; k++){
		std::cout << hnlMass[0][j][k]->GetBinContent(nBin) << " "  << hnlMass[2][j][k]->GetBinContent(nBin)<< std::endl;
		if(hnlMass[2][j][k]->GetBinContent(nBin)!=0 )hnlMass[0][j][k]->SetBinContent(nBin,1.0*hnlMass[3][j][k]->GetBinContent(nBin)*hnlMass[1][j][k]->GetBinContent(nBin)/(hnlMass[2][j][k]->GetBinContent(nBin)));	
		std::cout << hnlMass[0][j][k]->GetBinContent(nBin)<< std::endl;
		if(hnlMass[2][j][k]->GetBinContent(nBin)!=0 )hnlMass[4][j][k]->SetBinContent(nBin,1.0*hnlMass[1][j][k]->GetBinContent(nBin)/(hnlMass[2][j][k]->GetBinContent(nBin)));	
		if(hnlMass[2][j][k]->GetBinContent(nBin)!=0 )hnlMass[0][j][k]->SetBinError(nBin,sqrt(pow(hnlMass[1][j][k]->GetBinError(nBin)*hnlMass[3][j][k]->GetBinContent(nBin)/hnlMass[2][j][k]->GetBinContent(nBin),2)+pow(hnlMass[3][j][k]->GetBinError(nBin)*hnlMass[1][j][k]->GetBinContent(nBin)/hnlMass[2][j][k]->GetBinContent(nBin),2)+pow(hnlMass[2][j][k]->GetBinError(nBin)*hnlMass[1][j][k]->GetBinContent(nBin)*hnlMass[3][j][k]->GetBinError(nBin)/pow(hnlMass[2][j][k]->GetBinContent(nBin),2),2)));	
		if(hnlMass[2][j][k]->GetBinContent(nBin)!=0 )hnlMass[4][j][k]->SetBinError(nBin,sqrt(pow(hnlMass[1][j][k]->GetBinError(nBin)/hnlMass[2][j][k]->GetBinContent(nBin),2)+pow(hnlMass[2][j][k]->GetBinError(nBin)*hnlMass[1][j][k]->GetBinContent(nBin)/pow(hnlMass[2][j][k]->GetBinContent(nBin),2),2)));	
			

			}	
	
	}	
	}
	std::ofstream bkg;
	bkg.open((wd+"/BkgYields_.txt").c_str());
	for(int i = 0; i< nRegion; i++){
		for(int k = 0; k< nLxyBin; k++){
//	hnlMass[i][j][k]->Scale(41.6/1.79);
	SavePlot ("HNL Mass(GeV)", hnlMass[i][j][k], (wd+"/ABCDplots/"+std::string(hnlMass[i][j][k]->GetName())).c_str(), false, NULL, false);

	if (i==4 || i==0 ) hnlMass[i][j][k]->SaveAs((wd+"/ABCDplots/"+std::string(hnlMass[i][j][k]->GetName())+"TF_freezed.root").c_str());
		
	if (i==0){
	double integ;
	double error;
	integ = hnlMass[i][j][k]->IntegralAndError(1,1,error);
	bkg << j  << " " << k  << " " << integ << " " << error <<std::endl;
	hnlMass_MC[j][k]->Scale(1.0*hnlMass[i][j][k]->Integral()/hnlMass_MC[j][k]->Integral());
	RatioPlot ("HNL Mass(GeV)", hnlMass_MC[j][k],hnlMass[i][j][k],(wd+"/ABCDplots/DataVsQCD_"+std::string(hnlMass[i][j][k]->GetName())).c_str(), false, 1.1);	
	}
	}
	
	}

	return 1;
}



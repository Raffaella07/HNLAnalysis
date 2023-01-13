//#include "/interface/GBRForestTools.h"
//#include "../interface/BParkTools.h"
//#include "../interface/Utils.h"
//#include "../interface/SmallBToMuLPiClass.h"
//#include "TTree.h"
////#include <ROOT/RDataFrame.hxx>
////#include <ROOT/RVec.hxx>
////#include "../../../XGBoost-FastForest/include/fastforest.h"
//#include <stdint.h>
////#include "TStopwatch.h"
//#include <iostream>
//#include <fstream>
//#include <utility>
//#include <vector>
//#include <math.h>
//#include <dirent.h>
//#include "TFile.h"
//#include "TAxis.h"
//#include "TTree.h"
//#include "TLatex.h"
//#include "TCanvas.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TString.h"
//#include "TLorentzVector.h"
//#include "TMVA/Tools.h"
//#include "TMVA/Reader.h"
//#include "TMVA/MethodBDT.h"
//#include "ROOT/RDataFrame.hxx"


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
	hist->SetFillStyle(0);
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
 	leg->AddEntry(histMock,"ABCD QCD MC","lep");
	
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
void SuperPlot (std::string titlestring, TH1D * hist1,TH1D * hist2,  const char * filename, bool log, TF1* fit, bool lable){

	TCanvas* canvas = new TCanvas("canva","canva",600,550);
	//TFile* histos = new TFile("gausslog.root","UPDATE","",0);
//	std::string PNGPATH = "/eos/home-r/ratramon/www/";
//	std::string PDFPATH = "../plots/";
	//const char * temptitle = titlestring.c_str();
	setStyle();
	if (log) canvas->SetLogy();
	hist1->GetXaxis()->SetMaxDigits(2);
	hist1->SetLineWidth(1);
	hist1->SetMarkerStyle(8);
	hist1->SetMarkerSize(1);
	hist1->SetLineColor(kAzure+7);
	hist1->SetMarkerColor(kAzure+7);
	hist2->GetXaxis()->SetMaxDigits(2);
	hist2->SetLineWidth(1);
	hist2->SetMarkerStyle(8);
	hist2->SetMarkerSize(1);
	hist2->SetLineColor(kRed-7);
	hist2->SetMarkerColor(kRed-7);
	hist1->GetXaxis()->SetTitle(titlestring.c_str());
	hist1->GetYaxis()->SetTitle("entries");
	std::cout << "axis title" << titlestring.c_str() << std::endl;
	hist1->Draw("PE");
	hist2->Draw("samePE");
	if(fit != NULL) fit->Draw("samePE1");
//	if(lable)lables1D(canvas,histo)
	canvas->SaveAs((std::string(filename)+".pdf").c_str());
//	canvas->SaveAs((PNGPATH+std::string(filename)+".png").c_str());
	canvas->Clear();
//	histos->Close();

/*	TCanvas* canvaNorm = new TCanvas("canva","canva",600,550);

	hist1->DrawNormalized("PE1");
	hist2->DrawNormalized("samePE1");
	
	
	canvaNorm->SaveAs((std::string(filename)+"_norm.pdf").c_str());*/
	delete canvas;
//	delete canvaNorm;
}


void AddQCDToHist(int SetIdx,double QCDweight,TChain* tree, std::string filter, std::string leaf,ROOT::RDF::TH1DModel h_model,TH1D * hist, std::string name){
//void AddQCDToHist(int SetIdx,double QCDweight,ROOT::RDataFrame d, std::string filter, std::string leaf,ROOT::RDF::TH1DModel h_model,TH1D * hist, std::string name)
	//ROOT::EnableImplicitMT(4);
	if(SetIdx ==0){
			tree->Draw((leaf+">>"+std::string(hist->GetName())).c_str(),filter.c_str());
//			d.Filter(filter.c_str()).Histo1D(h_model,leaf.c_str()).GetValue().Copy(*hist);
			hist->Scale(QCDweight);		
			hist->SetTitle(name.c_str());		
			hist->SetName(name.c_str());		
	}else{ 
			TH1D * temp = new TH1D("temp","temp",40,1,7);
			tree->Draw((leaf+">>"+std::string(temp->GetName())).c_str(),filter.c_str());
		//	d.Filter(filter.c_str()).Histo1D(h_model,leaf.c_str()).GetValue().Copy(*temp);
			temp->Scale(QCDweight);		
			temp->SetTitle(name.c_str());		
			temp->SetName(name.c_str());		
			hist->Add(temp);



		}
}
int QCD_histos(){

	gROOT->SetBatch(kTRUE);	
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
	std::string RegLable[ABCD] = {"A_SR","B","C","D","TF","MockA"};
	std::string SignLable [2] = {"OS","SS"};
	RegCut[0] = cutX + " && " + cutY; //A
	RegCut[1] = cutX + " && !" + cutY; //B
	RegCut[2] = "!"+cutX + " && !" + cutY; //C
	RegCut[3] = "!"+cutX + " && " + cutY; //D
	RegCut[4] = " "; //Global TF on A 
	RegCut[5] = " "; // recomputed A 
	ROOT::RDF::TH1DModel *  h_models[2];
	h_models[0] = new ROOT::RDF::TH1DModel("h0","h1",40, 1, 7);
	h_models[1] = new ROOT::RDF::TH1DModel("h1","h1",4,0,4);
	std::string INPATH = "/cmshome/ratramon/Analysis/data/HNLFlatTuples/";
	//std::string InDataset[nDataset] = {"QCD_Pt15_20","QCD_Pt20_30","QCD_Pt30_50","QCD_Pt50_80","QCD_Pt80_120","QCD_Pt80_120_ext","QCD_Pt120_170","QCD_Pt120_170_ext","QCD_Pt170_300"};
	std::string InDataset[nDataset] = {"Pt-15to20","Pt-20to30","Pt-30to50","Pt-50to80","Pt-80to120","Pt-80to120_ext","Pt-120to170","Pt-120to170_ext","Pt-170to300"};

	//qcd sample weights 
	double QCDweights[nDataset];
	std::ifstream infile("/cmshome/ratramon/Analysis/plugins/QCDWeights.txt");
	if (infile.is_open()){
		i =0;
		while ( !infile.eof()){
            
		infile>>QCDweights[i];
		i++;
		}
	}

	
	TH1D* hnlMass[Channels][LxyBin][ABCD][2];
//	TH1D* hnlLxy[Channels][ABCD][2];
//	std::vector<ROOT::RDataFrame> d;
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
	
	
	for(int sidx=0; sidx<2;sidx++){
	for(int ch=0; ch<Channels;ch++){
	for(int reg=0; reg<4;reg++){

  //              hnlLxy[ch][reg][sidx]= new TH1D();
			for (i=0;i<nDataset;i++){
		/*	std::string Lxyfilter;
			if (sidx >0)Lxyfilter = "LepQProd>0 &&  Type =="+std::to_string(ch)+" && "+RegCut[reg];
			else Lxyfilter = "LepQProd<0 &&  Type =="+std::to_string(ch)+" && "+RegCut[reg];
			std::string Lxytitle = ChannelLable[ch]+"HNL_LxyBin_"+"Region"+RegLable[reg]+SignLable[sidx]+RegLable[reg];
			AddQCDToHist(i,QCDweights[i],d.at(i),Lxyfilter,"hnl_lxy",*h_models[1],hnlLxy[ch][reg][sidx],Lxytitle);
		*/
			for(int l=0; l<LxyBin;l++){
			std::string Massfilter;
			if (sidx>0) Massfilter ="LepQProd>0 && Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+RegCut[reg];
			else  Massfilter ="LepQProd<0 && Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+RegCut[reg];
			std::string Masstitle = ChannelLable[ch]+"HNLMass_LxyBin"+std::to_string(l)+"_"+"Region"+RegLable[reg]+SignLable[sidx];
		       	if (i==0)	hnlMass[ch][l][reg][sidx]= new TH1D(Masstitle.c_str(),Masstitle.c_str(),40,1,7);
		//	if (i ==0 ) std::cout << Masstitle << std::endl;
		//	AddQCDToHist(i,QCDweights[i],d.at(i),Massfilter,"hnl_mass",*h_models[0],hnlMass[ch][l][reg][sidx],Masstitle);
			AddQCDToHist(i,QCDweights[i],tree.at(i),Massfilter,"hnl_mass",*h_models[0],hnlMass[ch][l][reg][sidx],Masstitle);
		/*	if (i==0){	
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+RegCut[reg]).c_str()).Histo1D(*h_models[0],"hnl_mass").GetValue().Copy(*hnlMass[ch][l][reg]);
			hnlMass[ch][l][reg]->Scale(QCDweights[i]);		
			std::cout << CR_hnlMass[ch][l]->GetEntries()<<std::endl;		
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+CR_sel).c_str()).Histo1D({("CR_hnlLxy"+label).c_str(),("CR_hnlLxy"+label).c_str(),4u,0,4},"LxyBin").GetValue().Copy(*CR_hnlLxy[ch][l]);
			CR_hnlLxy[ch][l]->Scale(QCDweights[i]);		
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && " +SR_sel).c_str()).Histo1D({("SR_hnlMass"+label).c_str(),("SR_hnlMass"+label).c_str(),100u,0,7},"hnl_mass").GetValue().Copy(*SR_hnlMass[ch][l]);
			SR_hnlMass[ch][l]->Scale(QCDweights[i]);		
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && " +SR_sel).c_str()).Histo1D({("SR_hnlLxy"+label).c_str(),("SR_hnlLxy"+label).c_str(),4u,0,4},"LxyBin").GetValue().Copy(*SR_hnlLxy[ch][l]);
			SR_hnlLxy[ch][l]->Scale(QCDweights[i]);		

		}else if (i>0){
			TH1D * temp = new TH1D();
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+CR_sel).c_str()).Histo1D({("CR_hnlMass"+label).c_str(),("CR_hnlMass"+label).c_str(),100u,0,7},"hnl_mass").GetValue().Copy(*temp);
			temp->Scale(QCDweights[i]);		
			CR_hnlMass[ch][l]->Add(temp);
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+ CR_sel).c_str()).Histo1D({("CR_hnlLxy"+label).c_str(),("CR_hnlLxy"+label).c_str(),4u,0,4},"LxyBin").GetValue().Copy(*temp);
			temp->Scale(QCDweights[i]);		
			CR_hnlLxy[ch][l]->Add(temp);
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+ SR_sel).c_str()).Histo1D({("SR_hnlMass"+label).c_str(),("SR_hnlMass"+label).c_str(),100u,0,7},"hnl_mass").GetValue().Copy(*temp);
			temp->Scale(QCDweights[i]);		
			SR_hnlMass[ch][l]->Add(temp);
			d.at(i).Filter(("Type =="+std::to_string(ch)+" && LxyBin =="+std::to_string(l)+" && "+ SR_sel).c_str()).Histo1D({("SR_hnlLxy"+label).c_str(),("SR_hnlLxy"+label).c_str(),4u,0,4},"LxyBin").GetValue().Copy(*temp);
			temp->Scale(QCDweights[i]);		
			SR_hnlLxy[ch][l]->Add(temp);

		


		
		}*/
		}
		
		}

	
	for(int l=0; l<LxyBin;l++){
	SavePlot ("HNL mass(GeV)", hnlMass[ch][l][reg][sidx], ("ABCDplots/"+std::string(hnlMass[ch][l][reg][sidx]->GetName())).c_str(), false, NULL, false);	
	}
//	SavePlot ("HNL Lxy(GeV)", hnlLxy[ch][reg][sidx], ("ABCDplots/"+std::string(hnlLxy[ch][reg][sidx]->GetName())).c_str(), false, NULL, false);	
//	std::cout<< "______________filename "<< hnlLxy[ch][reg][sidx]->GetName()<< std::endl;	
	}
	
	for(int l=0; l<LxyBin;l++){
		hnlMass[ch][l][4][sidx]= new TH1D("","",40,1,7);
		hnlMass[ch][l][4][sidx]->SetName((ChannelLable[ch]+"HNLMass_LxyBin"+std::to_string(l)+"_"+"Region"+RegLable[4]+"_"+SignLable[sidx]).c_str());
		hnlMass[ch][l][4][sidx]->SetTitle((ChannelLable[ch]+"HNLMass_LxyBin"+std::to_string(l)+"_"+"Region"+RegLable[4]+"_"+SignLable[sidx]).c_str());
		for (int nBin=1; nBin<hnlMass[ch][l][0][sidx]->GetXaxis()->GetNbins();nBin++){
		if(hnlMass[ch][l][2][sidx]->GetBinContent(nBin)!=0 )hnlMass[ch][l][4][sidx]->SetBinContent(nBin,1.0*hnlMass[ch][l][1][sidx]->GetBinContent(nBin)*hnlMass[ch][l][3][sidx]->GetBinContent(nBin)/hnlMass[ch][l][2][sidx]->GetBinContent(nBin));	
		if(hnlMass[ch][l][2][sidx]->GetBinContent(nBin)!=0 )hnlMass[ch][l][4][sidx]->SetBinError(nBin,sqrt(pow(hnlMass[ch][l][1][sidx]->GetBinError(nBin)*hnlMass[ch][l][3][sidx]->GetBinContent(nBin)/hnlMass[ch][l][2][sidx]->GetBinContent(nBin),2)+pow(hnlMass[ch][l][3][sidx]->GetBinError(nBin)*hnlMass[ch][l][1][sidx]->GetBinContent(nBin)/hnlMass[ch][l][2][sidx]->GetBinContent(nBin),2)+pow(hnlMass[ch][l][1][sidx]->GetBinError(nBin)*hnlMass[ch][l][3][sidx]->GetBinContent(nBin)*hnlMass[ch][l][2][sidx]->GetBinError(nBin)/pow(hnlMass[ch][l][2][sidx]->GetBinContent(nBin),2),2)));	
	}
	std::cout << "________________________bin content " << hnlMass[ch][l][4][sidx]->GetBinContent(1) <<" numerator " << hnlMass[ch][l][1][sidx]->GetBinContent(1)*hnlMass[ch][l][3][sidx]->GetBinContent(1) << " denominator " << hnlMass[ch][l][2][sidx]->GetBinContent(1) <<  " " << (float)1.0*hnlMass[ch][l][1][sidx]->GetBinContent(1)*hnlMass[ch][l][3][sidx]->GetBinContent(1)/(float)hnlMass[ch][l][2][sidx]->GetBinContent(1) <<  std::endl;
	std::cout << hnlMass[ch][l][4][sidx]->GetEntries() <<std::endl;
	SavePlot ("HNL mass(GeV)", hnlMass[ch][l][4][sidx], ("ABCDplots/"+std::string(hnlMass[ch][l][4][sidx]->GetName())).c_str(), false, NULL, false);	
	}

/*	hnlLxy[ch][4][sidx]= new TH1D("","",4,0,4);
	for (int nBin=1; nBin<hnlLxy[0][4][sidx]->GetXaxis()->GetNbins();nBin++){
		hnlLxy[ch][4][sidx]->SetBinContent(nBin,1.0*hnlLxy[ch][1][sidx]->GetBinContent(nBin)*hnlLxy[ch][3][sidx]->GetBinContent(nBin)/hnlLxy[ch][2][sidx]->GetBinContent(nBin));	
		if(hnlLxy[ch][2][sidx]->GetBinContent(nBin)!=0 )hnlLxy[ch][4][sidx]->SetBinError(nBin,sqrt(pow(hnlLxy[ch][1][sidx]->GetBinError(nBin)*hnlLxy[ch][3][sidx]->GetBinContent(nBin)/hnlLxy[ch][2][sidx]->GetBinContent(nBin),2)+pow(hnlLxy[ch][3][sidx]->GetBinError(nBin)*hnlLxy[ch][1][sidx]->GetBinContent(nBin)/hnlLxy[ch][2][sidx]->GetBinContent(nBin),2)+pow(hnlLxy[ch][1][sidx]->GetBinError(nBin)*hnlLxy[ch][3][sidx]->GetBinContent(nBin)*hnlLxy[ch][2][sidx]->GetBinError(nBin)/pow(hnlLxy[ch][2][sidx]->GetBinContent(nBin),2),2)));	
	}
	hnlLxy[ch][4][sidx]->SetName((ChannelLable[ch]+"HNL_LxyBin_"+"Region"+RegLable[4]+"_"+SignLable[sidx]).c_str());
	hnlLxy[ch][4][sidx]->SetTitle((ChannelLable[ch]+"HNL_LxyBin_"+"Region"+RegLable[4]+"_"+SignLable[sidx][sidx]).c_str());
	SavePlot ("HNL Lxy(GeV)", hnlLxy[ch][4][sidx], ("ABCDplots/"+std::string(hnlLxy[ch][4][sidx]->GetName())).c_str(), false, NULL, false);	
	
	RatioPlot ("HNL Lxy(GeV)",hnlLxy[ch][0][sidx],hnlLxy[ch][4][sidx], ("ABCDplots/super_"+std::string(hnlLxy[ch][4][sidx]->GetName())).c_str(), false,1.1);	
	*/
	double scale;
	for(int l=0; l<LxyBin;l++){
	std::cout << hnlMass[ch][l][0][sidx]->GetEntries()<< std::endl;
	if (l==0 || l==3) scale = 1.2;
	else scale = 3;
	RatioPlot ("HNL Mass(GeV)",hnlMass[ch][l][0][sidx],hnlMass[ch][l][4][sidx], ("ABCDplots/super_"+std::string(hnlMass[ch][l][4][sidx]->GetName())).c_str(), false, scale);	
	}


}
}
	//hnlMass[0][0][0]->SaveAs("histTest.root");
//	hnlMass[0][0][4][sidx]->SaveAs("histTest_Mock.root");
return 1;
}




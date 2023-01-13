#include "setStyle.C"



void EffPlotter(){

	TFile * mc_mu = TFile::Open("../plugins/HNL_mu.root");
	TFile * mc_e = TFile::Open("../plugins/HNL_e.root");
	TFile * input = TFile::Open("../python/data_noTrgMuClean_loose.root");

	TH1D* h_mu = (TH1D*) mc_mu->Get("HNL_mu");
	TH1D* h_e = (TH1D*) mc_e->Get("HNL_e");
	TTree* looseTree  = (TTree*) input->Get("Events");

	TH1D * h_PF = new TH1D("h_PF","PF electrons",50,0,50);
	TH1D * h_LowPt = new TH1D("h_LowPt","LowPt electrons",50,0,50);
	TH1D * h_Mu = new TH1D("h_mu","muons",50,0,50);
	TH1D * h_trke = new TH1D("h_trke","trk electrons",50,0,50);
	TH1D * h_trkmu = new TH1D("h_trkmu","trk muons",50,0,50);
	
	looseTree->Draw("hnl_lxy>>h_mu","Type==0");
	looseTree->Draw("hnl_lxy>>h_PF","Type==1");
	looseTree->Draw("hnl_lxy>>h_LowPt","Type==2");
	looseTree->Draw("hnl_lxy>>h_trke","Type==3 && toMu==1");
	looseTree->Draw("hnl_lxy>>h_trkmu","Type==3 && toEle==1");


	TEfficiency* eff_PF = new TEfficiency(*h_PF,*h_e);
	h_LowPt->Add(h_PF);
	TEfficiency* eff_LowPt = new TEfficiency(*h_LowPt,*h_e);
	//eff_LowPt->Add(*eff_PF);
	h_trke->Add(h_LowPt);
	TEfficiency* eff_trke = new TEfficiency(*h_trke,*h_e);
	TEfficiency* eff_Mu = new TEfficiency(*h_Mu,*h_mu);
	h_trkmu->Add(h_Mu);
	TEfficiency* eff_trkmu = new TEfficiency(*h_trkmu,*h_mu);

	setStyle();
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(42);

	TCanvas* canva_e = new TCanvas("canvaE","canvaE",800,600);
	TLegend* l_e = new TLegend(0.5,0.7,0.87,0.87);
	l_e->SetFillColor(0);
	TH2D* plotter = new TH2D("plotter","plotter",10,0,50,10,0,1);	

	plotter->GetXaxis()->SetTitle("HNL_lxy (cm)");
	plotter->GetYaxis()->SetTitle("#varepsilon_{reconstruction}");
	
	eff_PF->SetMarkerStyle(8);
	eff_LowPt->SetMarkerStyle(8);
	eff_trke->SetMarkerStyle(8);
	eff_Mu->SetMarkerStyle(8);
	eff_trkmu->SetMarkerStyle(8);

	l_e->SetHeader("#mu HNL #rightarrow e#pi signatures");
	l_e->AddEntry(eff_PF,"PF electrons","lp");
	l_e->AddEntry(eff_LowPt,"PF+LowPt electrons","lp");
	l_e->AddEntry(eff_trke,"PF+LowPt+trk electrons","lp");

	eff_PF->SetLineColor(kAzure+7);
	eff_PF->SetMarkerColor(kAzure+7);
	eff_LowPt->SetLineColor(kBlue-7);
	eff_LowPt->SetMarkerColor(kBlue-7);
	eff_trke->SetLineColor(kBlue-3);
	eff_trke->SetMarkerColor(kBlue-3);
	eff_Mu->SetLineColor(kRed-7);
	eff_Mu->SetMarkerColor(kRed-7);
	eff_trkmu->SetLineColor(kRed-4);
	eff_trkmu->SetMarkerColor(kRed-4);

	plotter->Draw();
	eff_PF->Draw("sameP");
	eff_LowPt->Draw("sameP");
	eff_trke->Draw("sameP");

	l_e->Draw("same");
	canva_e->SaveAs("RecoEff_ele.pdf");

	TCanvas* canva_mu = new TCanvas("canvaMu","canvaMu",800,600);
	
	TLegend* l_mu = new TLegend(0.5,0.7,0.87,0.87);
	l_mu->SetFillColor(0);
	l_mu->SetHeader("#mu HNL #rightarrow #mu#pi signatures");
	l_mu->AddEntry(eff_Mu,"Muons","lp");
	l_mu->AddEntry(eff_trkmu,"Muons+MuonTrk","lp");
	plotter->Draw();
	eff_Mu->Draw("sameP");
	eff_trkmu->Draw("sameP");

	l_mu->Draw("same");
	canva_mu->SaveAs("RecoEff_mu.pdf");
	//canva->BuildLegend();	


}

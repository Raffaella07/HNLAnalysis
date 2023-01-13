#include "setStyle.C"

int plotBkgYields(std::string file1, std::string file2, std::string file3 ){

	setStyle();

	TCanvas* c = new TCanvas("c","c",1000,600);
	TH2D* plotter = new TH2D("plotter","plotter",6,0,6,10,0.1,1000000);
	plotter->GetXaxis()->SetBinLabel(1,"lxy0to1_OS");
	plotter->GetXaxis()->SetBinLabel(2,"lxy1to5_OS");
	plotter->GetXaxis()->SetBinLabel(2,"lxygtr5_OS");
	plotter->GetXaxis()->SetBinLabel(3,"lxy0to1_OS");
	plotter->GetXaxis()->SetBinLabel(3,"lxy0to1_SS");
	plotter->GetXaxis()->SetBinLabel(2,"lxy1to5_OS");
	plotter->GetXaxis()->SetBinLabel(3,"lxygtr5_OS");
	plotter->GetXaxis()->SetBinLabel(4,"lxy0to1_SS");
	plotter->GetXaxis()->SetBinLabel(4,"lxy1to5_SS");
	plotter->GetXaxis()->SetBinLabel(5,"lxy1to5_SS");
	plotter->GetXaxis()->SetBinLabel(6,"lxygtr5_SS");
	plotter->Draw();
	c->SetLogy();
//	gStyle->SetPalette(91);
	TGraphErrors * g[3];
	gPad->SetGrid();
	plotter->GetYaxis()->SetTitle("background yields");
	g[0]= new TGraphErrors("YieldsFrom5p3ToFullLumi_Mass1.txt","%lg %lg %lg %lg");
	g[1]= new TGraphErrors("YieldsFrom5p3ToFullLumi_Mass3_cut.txt","%lg %lg %lg %lg");
	g[2]= new TGraphErrors("YieldsFrom5p3ToFullLumi_Mass4p5.txt","%lg %lg %lg %lg");
	g[0]->SetTitle("2#sigma window around 1 GeV");
	g[1]->SetTitle("2#sigma window around 3 GeV");
	g[2]->SetTitle("2#sigma window around 4.5 GeV");
	g[0]->SetMarkerStyle(8);
	g[0]->SetMarkerColor(kOrange);
	g[0]->SetLineColor(kOrange);
	g[0]->SetMarkerSize(1.5);
	g[0]->SetLineWidth(1);
	g[1]->SetMarkerStyle(8);
	g[1]->SetLineWidth(1);
	g[1]->SetMarkerColor(kRed+1);
	g[1]->SetLineColor(kRed+1);
	g[1]->SetMarkerSize(1.5);
	g[2]->SetMarkerStyle(8);
	g[2]->SetMarkerColor(kRed+4);
	g[2]->SetLineColor(kRed+4);
	g[2]->SetLineWidth(1);
	g[2]->SetMarkerSize(1.5);
//	g[0]->Draw("PEsame");
	g[1]->Draw("PEsame");
//	g[2]->Draw("PEsame");
	c->BuildLegend();
	c->SaveAs("BkgCutCat.pdf");

	return 1;
}

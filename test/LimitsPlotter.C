#include "../macros/setStyle.C"





void limPlotter(std::string limFile,double mass ,double xlow, double xup, double ylow, double yup, std::string tag, std::string channel){

	setStyle();
	TH2D * plotter = new TH2D("plotter","plotter",10,xlow,xup,10,ylow,yup);
	
	TGraph* line = new TGraph(limFile.c_str(),"%lg %lg lg lg lg lg  ");
	TGraphAsymmErrors* g = new TGraphAsymmErrors(limFile.c_str(),"%lg %lg %lg %lg lg lg ");
	TGraphAsymmErrors* o = new TGraphAsymmErrors(limFile.c_str(),"%lg %lg %lg %lg %lg %lg ");

	TCanvas* c = new TCanvas("c","c",600,800);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendFont(42);
	TLegend* l = new TLegend();
	c->SetLogy();
	c->SetLogx();
	plotter->GetXaxis()->SetTitle("|V|^{2}");
	plotter->GetYaxis()->SetTitle(" limit @ 95% C.L ");
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	g->SetFillColor(kGreen-3);
	o->SetFillColor(kOrange);
	l->AddEntry(line,"expected limit","l");
	l->AddEntry(g,"#pm 1 #sigma","f");
	l->AddEntry(o,"#pm 2 #sigma","f");
	plotter->Draw();
	o->Draw("same3");
	g->Draw("same3");
	line->Draw("sameL");
	l->Draw("same");


	c->SaveAs(("CombinedLimits/Mass"+std::to_string(mass)+"_"+tag+"_"+channel+".pdf").c_str());
	



}

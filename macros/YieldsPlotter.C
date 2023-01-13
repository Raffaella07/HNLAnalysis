#include "setStyle.C"

void YieldsPlotter(std::string wd, std::string sign ){

	TMultiGraph* m = new TMultiGraph();
	
	
	std::string lables[3] = {"Yields_LxyUnder1","Yields_LxyOver1_Under5","Yields_LxyOver5"};

	TGraph* signal[3];
	TGraph* background[3];
	int colors[3] = {kOrange+7,kMagenta+2,kAzure-3};
	std::string lxy_string[3] = { " lxy [0,1] cm","lxy [1,5] cm ","lxy over 5 cm" };
	
	int i; 
	setStyle();

	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendFont(42);
	TLegend* l = new TLegend(0.51,0.20,0.85,0.40);
//	l->SetFillStyle(0);
	l->SetHeader(("#mu e-#pi "+sign+" channel expected yields").c_str());

	for (i=0; i<3;i++){


		signal[i] = new TGraph((wd+"/"+lables[i]+"_"+sign+".txt").c_str(),"%lg %lg %*lg");
		background[i] = new TGraph((wd+"/"+lables[i]+"_"+sign+".txt").c_str(),"%lg %*lg %lg");
		signal[i]->SetLineWidth(2);
		background[i]->SetLineWidth(2);
		background[i]->SetLineStyle(7);
		signal[i]->SetLineColor(colors[i]);
		background[i]->SetLineColor(colors[i]);
		m->Add(signal[i],"l");	
//		m->Add(background[i],"l");	
		l->AddEntry(signal[i],("signal, m = 3 GeV, "+lxy_string[i]).c_str(),"l");
//		l->AddEntry(background[i],("background, "+lxy_string[i]).c_str(),"l");
		

	}

	TCanvas * c = new TCanvas("c","c",800,600);
	TH2D * plotter = new TH2D("plotter","plotter",10,0.000001,0.005,10,0.00001,10000000);
	plotter->GetXaxis()->SetTitle("#frac{|V_{#mu} #upoint V*_{e}|^{2}}{|V|^{2}}");
	plotter->GetYaxis()->SetTitle("expected yields @ 41.6 fb^{-1}");

	plotter->Draw();
	c->SetLogy();
	c->SetLogx();
	c->SetGrid();
	m->Draw("l");
	l->Draw("same");

	c->SaveAs((wd+"/ExpectedYields_"+sign+".pdf").c_str());
		
	
	





}

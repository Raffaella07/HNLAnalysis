#include "setStyle.C"
#include <vector>

void LikelihoodScan( std::string inputList, bool debug, std::string nametag){

	int i;
	setStyle();
	std::vector<std::string> LimFile;	
	std::vector<std::string> LimTitle;	
	std::vector<int> LimType;	
	gStyle->SetPalette(91);


	//retrieving files with limits 
	 ifstream file (inputList.c_str());
         if (file.is_open())
         { 
	   std::string line;
           while ( getline (file,line) )
          {
          {
            std::stringstream ss(line);
            std::string filename, title, Ltype;
            if (ss >> filename >> title >> Ltype )
             { 
	  LimFile.push_back(filename);
	  LimTitle.push_back(title);
	  LimType.push_back(std::atoi(Ltype.c_str()));
	  
	  }
          }
        }
	}
        file.close();

	const int nFiles = LimFile.size();
	TGraph* g[nFiles];
	TMultiGraph *mg = new TMultiGraph();
	TH2D* plotter = new TH2D("plotter","plotter",10,0.5,5,10,0.004,100);


	for (i=0; i<nFiles;i++){

		g[i] = new TGraph(LimFile.at(i).c_str(),"%lg %lg %*lg %*lg %*lg %*lg");
		if (LimType.at(i)==0 ) g[i]->SetMarkerStyle(8);
		else g[i]->SetMarkerStyle(4);
		g[i]->SetLineWidth(2);
		g[i]->SetTitle(LimTitle.at(i).c_str());
		mg->Add(g[i], "PL");
		
	}

	TCanvas* c = new TCanvas ("c","c",600,800);
	c->SetLogy();
	mg->Draw("A plc pmc");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(42);
	gStyle->SetLegendFillColor(0);
	mg->GetXaxis()->SetTitle("HNL Mass (GeV)");
	mg->GetYaxis()->SetTitle("95% C.L limit median");
	//mg->GetXaxis()->SetLimits(1.5,3.5);
//	mg->GetYaxis()->SetLimits(0.003,100);
	c->BuildLegend(0.65,0.20,0.85,0.45);
	c->SaveAs("Lim_LikelihoodScan.pdf");
	return 0;
}

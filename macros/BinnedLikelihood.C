#include <sstream>
//script for likelihoood computation, independent from n and types of likelihood variables 

//Input arguments:
//ch: analysis channel: Mu, PF..
//mass: mass of MC signal used in sig distributions
//ctau: ctau of MC signal used in sig distributions
//point: vector containing the ordinated values on which the likelihood is to be computed
//order has to be the same as variables order in the pdf input lists
//
//ref for likelihood computation: slide 4 in https://indico.cern.ch/event/1120716/contributions/4705640/attachments/2379972/4066298/PFHigheEff_Likelihood%20studies.pdf
double BinnedLikelihood(int ch, int mass, int ctau,double* point){ 

	bool debug = false;
	std::ifstream inHistos;
	std::vector<TH1D*> sigHistos;
	std::vector<TH1D*> bkgHistos;
	
	//input list example in /macros/L_histos_ch1_Mass3_ctau10.txt
	inHistos.open("/cmshome/ratramon/Analysis/macros/L_histos_ch"+std::to_string(ch)+"_Mass"+std::to_string(mass)+"_ctau"+std::to_string(ctau)+".txt");
		
		
	if (inHistos.is_open())
  	{
		bool signal;
		std::string filename, histName;
		while(inHistos>>filename>>signal>> histName){
		TFile* f = TFile::Open(filename.c_str());
		TH1::AddDirectory(kFALSE);
		TH1D * h = new TH1D();
		h = (TH1D*)f->Get(histName.c_str())->Clone();
		if(signal)  sigHistos.push_back(h);
		else  bkgHistos.push_back(h);
		delete f;
		}
	}
	inHistos.close();
	double PSig = 1;
	double PBkg = 1;
	if (debug)std::cout << sigHistos.size() << std::endl;
	for (int nH=0; nH<sigHistos.size();nH++){
	if (debug) std::cout << point[nH] << std::endl;
	 
	for (int i=0; i<sigHistos[nH]->GetXaxis()->GetNbins()+1;i++){
		if ((sigHistos[nH]->GetXaxis()->FindBin(point[nH])==i) || (i==sigHistos[nH]->GetXaxis()->GetNbins() && point[nH]>sigHistos[nH]->GetXaxis()->GetBinLowEdge(i)) ){		

  		if(debug)std::cout<< point[nH] << " " << i << " bin content "<< sigHistos[nH]->GetBinContent(i) << std::endl;
        	PSig  *= sigHistos[nH]->GetBinContent(i);
  		if(debug)std::cout <<"Psig " << PSig<< std::endl;
        	}	
        	if (bkgHistos[nH]->GetXaxis()->FindBin(point[nH])==i|| (i==bkgHistos[nH]->GetXaxis()->GetNbins() && point[nH]>bkgHistos[nH]->GetXaxis()->GetBinLowEdge(i))){
  		if(debug)std::cout << point[nH] << " " << i << " bin content "<< bkgHistos[nH]->GetBinContent(i) << std::endl;
        	PBkg *= bkgHistos[nH]->GetBinContent(i);
  		if(debug)std::cout << "Pbkg " <<  PBkg << std::endl;

        	}	

         
        		


        	}
	
	if (debug) std::cout << " " << std::endl;
        }

  	if(debug)std::cout << PSig<< std::endl;
  	if(debug)std::cout << PBkg<< std::endl;
        double L = PSig/(PSig+PBkg);
  	if(debug){
	std::cout << L<< std::endl;
  	if(L==1)std::cout << "______________________________________________________________________________________________________________ "<< std::endl;
  	else std::cout << " "<< std::endl;
	}
	return L;	

	}




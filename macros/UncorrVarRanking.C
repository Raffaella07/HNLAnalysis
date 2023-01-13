
void UncorrVarRanker(std::string infile, std::string corrFile){

  ifstream myfile ("Ranking.txt");
  int i;
  TFile* fCorr = TFile::Open(corrFile.c_str());
  TH2D * sig = (TH2D*)fCorr->Get("dataset/CorrelationMatrixS")->Clone();
  TH2D * bkg = (TH2D*)fCorr->Get("dataset/CorrelationMatrixB")->Clone();
   
  if (myfile.is_open())
  { 
    int counter=0;
    while ( getline (myfile,line) )
    {
      cout << line << '\n';
      counter++;
      for (i =1; i<sig->GetXaxis()->GetNbins();i++){

	if (file == sig->GetXaxis()->GetBinLabel()):

      		for (int j =1; j<sig->GetXaxis()->GetNbins();j++){

		}
	
		
	}


	
    }
    myfile.close();
  }






}

#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"

void sig_fromcsv(std::string filename, std::vector<std::string>* flatIn,std::vector<std::string>* label, std::vector<double>* genEvent_Bu,std::vector<double>* genEvent_Bd,std::vector<double>* genEvent_Bs,std::vector<double>* genEvent,std::vector<double>* filterEff_Bu,std::vector<double>* filterEff_Bd,std::vector<double>* filterEff_Bs,std::vector<double>* filterEff,std::vector<double>* mass,std::vector<double>* ctau, std::vector<double>* sigma,std::vector<double>* Flav_ratio){

	fstream fin;

	// Open an existing file
	fin.open(filename.c_str(), ios::in);


	// Read the Data from the file
	// as String Vector
	vector<string> row;
	string line, word, temp;

	while (fin >> temp) {

		row.clear();

		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);

		// used for breaking words
		stringstream s(line);

		// read every column data of a row and
		// store it in a string variable, 'word'
		int counter = 0;
		if (line.find('+')== std::string::npos ){
		 continue;
		}
		while (getline(s, word,' ')) {

			// add all the column data
			// of a row to a vector
			counter++;
			row.push_back(word);
		}
		flatIn->push_back(row[15].c_str());
		label->push_back(row[1]);
//		std:: cout << flatIn->at(0)<< std::endl;
		genEvent_Bu->push_back(atof(row[5].c_str()));
		genEvent_Bd->push_back(atof(row[6].c_str()));
		genEvent_Bs->push_back(atof(row[7].c_str()));
		genEvent->push_back(atof(row[8].c_str()));
		filterEff_Bu->push_back(atof(row[9].c_str()));
		filterEff_Bd->push_back(atof(row[10].c_str()));
		filterEff_Bs->push_back(atof(row[11].c_str()));
		filterEff->push_back(atof(row[12].c_str()));
		//std:: cout << row[12]<< std::endl;
		mass->push_back(atof(row[2].c_str()));
		ctau->push_back(atof(row[3].c_str()));
		sigma->push_back(atof(row[13].c_str()));
		Flav_ratio->push_back(atof(row[14].c_str()));


		// Compare the roll number
/*		if (roll2 == rollnum) {

			// Print the found data
			count = 1;
			cout << "Details of Roll " << row[0] << " : \n";
			cout << "Name: " << row[1] << "\n";
			cout << "Maths: " << row[2] << "\n";
			cout << "Physics: " << row[3] << "\n";
			cout << "Chemistry: " << row[4] << "\n";
			cout << "Biology: " << row[5] << "\n";
			break;
		}*/
	}
//	if (count == 0)
//		cout << "Record not found\n";*/
}



double LifetimeReweight(std::string csvIn,float mass, float ctau,float ct,int flavor){


	std::vector<double> GenEvent_Bu, GenEvent_Bd, GenEvent_Bs,GenEvtCount,Filter_Bu,Filter_Bd,Filter_Bs,Filter,Mass,Ctau,Sigma,Flav_ratio;
	std::vector<std::string> Files, label;
	
	//sig_fromcsv("../data/MC_datasets_InclusiveFilterEff.csv",&Files,&label,&GenEvtCount,&Filter,&Mass,&Ctau,&Sigma,&Flav_ratio);
	sig_fromcsv(csvIn.c_str(),&Files,&label,&GenEvent_Bu,&GenEvent_Bd,&GenEvent_Bs,&GenEvtCount,&Filter_Bu,&Filter_Bd,&Filter_Bs,&Filter,&Mass,&Ctau,&Sigma,&Flav_ratio);

//dummy initialization for miniAOD genevt count numbers 
/*	double genEvtCount[5][4]={{9.080620000000000000e+05,4.532249000000000000e+06,1.312593500000000000e+07,-1},
				{7.857100000000000000e+05,1.955714000000000000e+06,7.546074000000000000e+06,-1},
                                {5.673050000000000000e+05,1.916707000000000000e+06,8.421721000000000000e+06,-1},
                                {5.252660000000000000e+05,6.562680000000000000e+05,1.314002000000000000e+06,3.282584000000000000e+06},
                                {5.914900000000000000e+05,5.065250000000000000e+05,5.962870000000000000e+05,1.195035000000000000e+06}};


//dummy initialization for the signal grid

	double grid[5][4]={     {10,100,1000,-1},
				{10,100,1000,-1},
				{10,100,1000,-1},
				{1,10,100,1000},
				{0.1,1,10,100},
			};
	double hnlMass[5]={1,1.5,2,3,4.5};
	
	double filter[5][4]={{9.49e-03,	7.90e-03,2.03e-03,-1},
				{8.50e-03,7.36e-03,2.13e-03,-1},
                                {7.99e-03,7.12e-03,2.05e-03,-1},
                                {1.69e-02,1.69e-02,1.59e-02,5.78e-03},
                                {2.58e-02,2.58e-02,2.58e-02,2.51e-02}};
*/

	
	long double den=0.;
	long double GenCountTot=0;
	double avgFilter=0;

	int f;
		for(int i=0;i<Files.size();i++){

//			std::cout << "mass " << mass<< " " << Mass[i]-mass <<  " ctau" << Ctau[i] <<  std::endl;
			if (fabs(mass-Mass[i])/mass<0.000001 && Ctau[i]>0. ){
				if (flavor ==4 || flavor ==3){	

				den += GenEvtCount[i]/(Ctau[i]*Filter[i]) * exp(-ct/(Ctau[i]));
				GenCountTot+=GenEvtCount[i];
				avgFilter+=GenEvtCount[i]*Filter[i];

				}else if (flavor == 0){

				den += GenEvent_Bu[i]/(Ctau[i]*Filter_Bu[i]) * exp(-ct/(Ctau[i]));
				GenCountTot+=GenEvent_Bu[i];
				avgFilter+=GenEvent_Bu[i]*Filter_Bu[i];
		
				}else if (flavor == 1){

				den += GenEvent_Bd[i]/(Ctau[i]*Filter_Bd[i]) * exp(-ct/(Ctau[i]));
				GenCountTot+=GenEvent_Bd[i];
				avgFilter+=GenEvent_Bd[i]*Filter_Bd[i];

				}else if (flavor == 2){

				den += GenEvent_Bs[i]/(Ctau[i]*Filter_Bs[i]) * exp(-ct/(Ctau[i]));
				GenCountTot+=GenEvent_Bs[i];
				avgFilter+=GenEvent_Bs[i]*Filter_Bs[i];

				}
	

	  		}
	}
		
	avgFilter = avgFilter/GenCountTot;
	long double num = GenCountTot*1.0/(ctau*avgFilter) * exp(-ct/(ctau));
/*	if (num !=0){
  	std::cout << "ctau " << ctau<<" ct " << ct  << std::endl;
  	std::cout << "avgFilter" << avgFilter<< std::endl;
  	std::cout << "GenEvtCount" << GenCountTot << std::endl;
  	std::cout << "suppression" << TMath::Exp(-ct/(ctau))<< std::endl;
	}*/
/*	std::cout << "num " << num << std::endl;
 	std::cout << "den " << den << std::endl;*/
//  	std::cout << "weight " << num*1.0/den << std::endl;
	return num/den;	

}

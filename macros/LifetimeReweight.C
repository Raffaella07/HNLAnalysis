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

double LifetimeReweight(float mass, float ctau,float ct){


//dummy initialization for miniAOD genevt count numbers 
	double genEvtCount[5][4]={{9.080620000000000000e+05,4.532249000000000000e+06,9.950223000000000000e+06,-1},
				{1.666740000000000000e+05,1.955714000000000000e+06,-1,-1},
                                {5.673050000000000000e+05,1.916707000000000000e+06,-1},
                                {4.411610000000000000e+05,6.562680000000000000e+05,1.314002000000000000e+06,3.282584000000000000e+06},
                                {1.653590000000000000e+05,5.065250000000000000e+05,5.106050000000000000e+05,1.195035000000000000e+06}};


//dummy initialization for the signal grid

	double grid[5][4]={     {10,100,1000,-1},
				{10,100,-1,-1},
				{10,100,-1,-1},
				{1,10,100,1000},
				{0.1,1,10,100},
			};
	double hnlMass[5]={1,1.5,2,3,4.5};
	


	
	double den=0.;
	int GenCountTot=0;


	for(int i=0;i<5;i++){

	if (mass==hnlMass[i]){
		
		for(int j=0;j<4;j++){

		    if (grid[i][j]>0 && ctau < grid[i][j]){
				den += genEvtCount[i][j]/grid[i][j] * exp(-ct/grid[i][j]);
				GenCountTot+=genEvtCount[i][j];

							}


				}


	  		}		
		}
	double num = GenCountTot*1.0/ctau * exp(-ct/ctau);
	std::cout << "num " << num << std::endl;
	std::cout << "den " << den << std::endl;
	return num/den;	

}

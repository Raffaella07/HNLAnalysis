#Datacards handler

#can be used to modify luminosity @ which limits are computed
#taking care of scaling for multiflavor coupling

import sys
from ctypes import *
import ROOT
import numpy as np 
import glob
import math
from ROOT import TH1D
from operator import itemgetter
#import multiprocessing
from HNLUtils import BR_HNLmupion,getVV,BToMuX


def scaleYields(datacard,scale,couplingScaling,tag,start):

	
	with open(datacard) as d:
	
		lines = d.readlines()
	linesArray = [l.split(" ") for l in lines]
	for i in range(0,len(linesArray)):
		if linesArray[i][start]=='rate':
		#	print linesArray[i]
			counter =0
			for j in range(start+1, len(linesArray[i])):
				if (linesArray[i][j] != '') and (linesArray[i][j] != '\n'):
					print(linesArray[i][j])
					print(counter)
					
					if couplingScaling and counter%2==0:
						linesArray[i][j] = np.str(np.float(linesArray[i][j])*scale)
					
					print(linesArray[i][j])
					counter+=1
					if not couplingScaling:
							linesArray[i][j] = np.str(np.float(linesArray[i][j])*scale)
		#	print linesArray[i]
			#print(i)	
		#	print linesArray[i-1]
		
	########for j in range(0, len(linesArray[i])):
	########	if (linesArray[i][j] == ''): 
	########		linesArray[i][j]==' '
	########print linesArray[i]
	fout = open(("../"+datacard.strip(".txt")+tag+".txt"), "w")
	
	#np.savetxt(fout,linesArray,fmt="%s")
	for j in range(0,len(linesArray)):
	    #print(linesArray[j])
	    np.savetxt(fout,linesArray[j],fmt="%s",newline=' ')
	
	fout.close()


def MoveSysts(datacard,Oldvalue,newValue):


	with open(datacard, 'r') as file :
	  filedata = file.read()
	
	# Replace the target string
	filedata = filedata.replace(np.str(Oldvalue)+" ", np.str(newValue)+" ")
	
	# Write the file out again
	with open(datacard, 'w') as file:
	  file.write(filedata)	


#def CouplingScaling(idatacard,f):

MuonCards = ["../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_1_v2_5p4em02_cat_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_1_v2_5p4em03_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_1p5_v2_7p1em03_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_1p5_v2_7p1em04_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_1p5_v2_7p1em05_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_2_v2_1p7em03_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_2_v2_1p7em04_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_2_v2_1p7em05_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_3_v2_2p2em03_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_3_v2_2p2em04_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_3_v2_2p2em05_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_3_v2_2p2em06_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_4p5_v2_2p9em03_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_4p5_v2_2p9em04_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_4p5_v2_2p9em05_combined.txt","../output/counting_experiment_realData_lumiD1/datacard_bhnl_m_4p5_v2_2p9em06_combined.txt"]
ElectronCards = ["../output/PFCards/HNL_M1_ctau1000_PF_OSxSS_combined.txt","../output/PFCards/HNL_M1_ctau100_PF_OSxSS_combined.txt","../output/PFCards/HNL_M1_ctau10_PF_OSxSS_combined.txt","../output/PFCards/HNL_M2_ctau1000_PF_OSxSS_combined.txt","../output/PFCards/HNL_M2_ctau100_PF_OSxSS_combined.txt","../output/PFCards/HNL_M2_ctau10_PF_OSxSS_combined.txt","../output/PFCards/HNL_M3_ctau1000_PF_OSxSS_combined.txt","../output/PFCards/HNL_M3_ctau100_PF_OSxSS_combined.txt","../output/PFCards/HNL_M3_ctau10_PF_OSxSS_combined.txt","../output/PFCards/HNL_M3_ctau1_PF_OSxSS_combined.txt","../output/PFCards/HNL_M4_ctau0_PF_OSxSS_combined.txt","../output/PFCards/HNL_M4_ctau100_PF_OSxSS_combined.txt","../output/PFCards/HNL_M4_ctau10_PF_OSxSS_combined.txt","../output/PFCards/HNL_M4_ctau1_PF_OSxSS_combined.txt"]
#for card in MuonCards:
#	MoveSysts(card,1.3,1.1)	
#	scaleYields(card,(41.6/5.12),False,"ScaledFullLumi",0)
#	scaleYields("../"+card.strip(".txt")+"ScaledFullLumi.txt",0.5*0.5,True,"_fmu0p5",1)
for card in ElectronCards:
	scaleYields(card,0.5*0.5,True,"_fmu0p5",0)

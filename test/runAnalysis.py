import os
import re 
import sys
import ROOT
from ROOT import gROOT,TChain 
import numpy as np
import glob

outdir = sys.argv[1]
selection = sys.argv[2]
part = "1A"
Adir = '/cmshome/ratramon/Analysis/'
#Bkgpath = "/cmshome/ratramon/Analysis/datasets/Data/1A"
Sigpath = "/cmshome/ratramon/Analysis/datasets/Signal/"
#prepare signals
with open("Signals.txt") as fSig:
	hnl_signal = fSig.readlines()

hnl_signal = [l.strip("\n") for l in hnl_signal]
hnl_sigFeatures = [l.strip(".root") for l in hnl_signal]
print(hnl_sigFeatures)
#lable for output files
sig_FileLables = ["HNL_%s"%(l) for l in hnl_sigFeatures]
hnl_sigFeatures = [l.split("_") for l in hnl_sigFeatures]
#masses
Mass = [hnl_sigFeatures[i][0].strip("Mass") for i in range(0,len(hnl_sigFeatures))] # strip string structure
print(Mass)
Mass = [np.float(Mass[i].replace("p",".")) for i in range(0,len(hnl_sigFeatures)) ] #correct format +  cast to float
#ctaus
ctau = [hnl_sigFeatures[i][1].strip("ctau") for i in range(0,len(hnl_sigFeatures)) ] # strip string structure
print(ctau)
ctau = [np.float(ctau[i].replace("p",".")) for i in range(0,len(hnl_sigFeatures)) ]#correct format + cast to float

#formatted lables
sig_OutLables = ["HNL Mass = %.1f GeV, c$tau = %.1f mm"%(Mass[i],ctau[i]) for i in range(0,len(Mass)) ] 
sig = np.zeros((4,4,5))
LimComb = np.zeros((4,5))
#prepare data
fullLumi = 41.6
if part == "1A":

	lumi = 0.76
	background = Adir+"/datasets/Data/1A_LikelihoodV1" 

# cp file to working area
for i in range(0, len(hnl_signal)):
	#Create output directory
	workDir = Adir+"/output/"+outdir
	os.system("mkdir %s"%(workDir))
	os.system("cd %s"%(workDir))
	#directory for signal fits
	os.system("mkdir %s/SigFits"%workDir)
	#directory for ABCD plots 
	os.system("mkdir %s/ABCDplots"%workDir)
	#Directory for datacards
	os.system("mkdir %s/Datacards"%workDir)
	
	gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/AllFits.C')
	
	gROOT.SetBatch(True)
	#fits to signal: produces a Signal_mass"+np.int(mass)+"_"+np.int(ctau)+"_features.txt file
	#file is positioned in working dir
	#file format is 'ChannelIdx LxyBinIdx SignalYield SignalYield_error Signal_sigma' 
	ROOT.AllFits(Sigpath+hnl_signal[i],Mass[i],ctau[i],workDir,selection)
	gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/Data_ABCD.C')
	#ABCD method for background estimates: performs ABCD in bins of channel and displacement and produces an output .txt file
	ROOT.Data_ABCD(background,0,True, True, Mass[i],ctau[i],workDir,selection)
	sigOutFile ="/Signal_mass%d_ctau%d_features.txt"%(np.int(Mass[i]),np.int(ctau[i]))	
	os.system("cd "+Adir+"/python")
	os.system("pwd")
	os.system("python "+Adir+"/python/Card_producer.py %s %s %f %f %s"%(sigOutFile,"/BkgYields_.txt",Mass[i],ctau[i],workDir))
#	execfile(Adir+"/python Card_producer.py %s %s %f %f %s"%(sigOutFile,"/Bkg_Yields.txt",Mass[i],ctau[i],workDir))
#	os.system("cp %s",%sigFileLable)
	os.system("cd ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/")
	Cards = glob.glob(workDir+"/Datacards/*.txt")	
	Cards_lable = [ l.strip(".txt") for l in Cards]
	Cards_lable = [ l.split("/") for l in Cards_lable]
	print(Cards_lable)
	Cards_lable = [Cards_lable[l][8] for l in range(0,len(Cards))]
	print(Cards_lable)
	lxy = []
	combinedString = []
	combinedString.append(" ")
	combinedString.append(" ")
	combinedString.append(" ")
	combinedString.append(" ")
	
	chLables = ["Mu","PF","LowPt","Track"]
	#lxyLables = ["LxyUnder3","LxyOver3_Under10","LxyOver10_Under20","LxyOver20"]
	lxyLables = ["LxyUnder1","LxyOver1_Under5","LxyOver5"]
	for k in range(0,len(Cards)):
		os.system(" combine -M AsymptoticLimits --run blind %s  "%(Cards[k]))
		os.system(" mv higgsCombineTest.AsymptoticLimits.mH120.root %s/Lim_%s.root"%(workDir,Cards_lable[k])) 
		c= TChain("limit")
		c.Add("%s/Lim_%s.root"%(workDir,Cards_lable[k]))
		for s in range(0,len(chLables)):
			match =Cards_lable[k].find(chLables[s])
	#		print(Cards_lable[k])
	#		print(match)
			if match!=-1:
				chIdx = s
		#		print(chIdx)
				combinedString[s] = combinedString[s] + " " +Cards[k]
				
		for s in range(0,len(lxyLables)):
			print(Cards_lable[k])
			print(lxyLables[s])
			match = Cards_lable[k].find(lxyLables[s])
			if match!=-1:
				lxyIdx = s	
		temp = []
		for idx in range(0,c.GetEntries()):
			c.GetEntry(idx)
			sig[chIdx][lxyIdx][idx]=c.GetLeaf("limit").GetValue(0)
			#print(temp[idx])
	#	print(sig[k/4][k%4])
	print(combinedString)						
	for k in range(0, len(combinedString)):
		os.system("mkdir "+ workDir+"/Datacards/CombinedCards")
		combinedCard = workDir+"/Datacards/CombinedCards/HNL_M"+np.str(np.int(Mass[i])) +"_ctau"+np.str(np.int(ctau[i]))+"_"+chLables[k]+"_combined.txt"
		print(combinedString[k])
		os.system("python ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "+ combinedString[k]+" > " +combinedCard)
		os.system(" combine -M AsymptoticLimits --run blind  %s "%(combinedCard))
		os.system(" mv higgsCombineTest.AsymptoticLimits.mH120.root %s.root"%(combinedCard.strip(".txt"))) 
		c= TChain("limit")
		c.Add("%s.root"%(combinedCard.strip(".txt")))
		for idx in range(0,c.GetEntries()):
			c.GetEntry(idx)
			LimComb[k][idx]=c.GetLeaf("limit").GetValue(0)
	with open("%s/Limits.txt"%workDir, "w") as limFile:
		for i in range(0,len(sig)):
			
			for j in range(0,len(sig[0])):
			
    				limFile.write("%d %d %f \n"%(i,j,sig[i][j][2]))
				
	
	with open("%s/CombinedLimits.txt"%workDir, "w") as limFile:
		for i in range(0,len(sig)):
			
    			limFile.write("%d  %f \n"%(i,LimComb[i][2]))

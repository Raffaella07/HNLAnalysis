import os
import re 
import sys
import ROOT
from ROOT import gROOT,TChain,gSystem 
import numpy as np
import glob
import time 
import HNLUtils
print(time)

def stepMonitoring(step):
	
	
	print("\n")
	print("\n")
	print("----------------------------------------------------------------------------------------------------------------------------------------")
	print("____________________________________________________process Time: %f____________________________________________________________________"%(time.time()))
	print("----------------------------------------------------------------------------------------------------------------------------------------")
	print("________________________________________________________%s__________________________________________"%step)
	print("________________________________________________________________________________________________________________________________________")
	print("\n")
	print("\n")




def LxyBinnedLimits(Sigpath,background,workDir,hnl_signal,sel,channel,lumi):

  chLables = ["Mu","PF","LowPt","Track"]
  chargeLables = ["OS","SS"]
  charge = ["&& LepQProd<0"," && LepQProd>0"]
  print(channel)
  chIdx = chLables.index(str(channel))
  print(chIdx)
  hnl_sigFeatures = hnl_signal.strip("_LV*.root") #prepare tag for signals
  sub= "HNL_%s_%s"%(hnl_sigFeatures,channel) 
  print(hnl_sigFeatures)
  hnl_sigFeatures = hnl_sigFeatures.split("_") 
  print(hnl_sigFeatures)
  #masses
  Mass = hnl_sigFeatures[0].strip("Mass")  # strip string structure
  print(Mass)
  Mass = np.float(Mass.replace("p","."))  #correct format +  cast to float
  #ctaus
  ctau = hnl_sigFeatures[1].strip("ctau")  # strip string structure
  print(ctau)
  ctau = np.float(ctau.replace("p",".")) #correct format + cast to float
 
  for idx in range(0,2):

     selection = sel + charge[idx]
     subdir = sub+"_"+chargeLables[idx]
 
     stepMonitoring("Preparing output directory")
     os.system("mkdir %s/%s"%(workDir,subdir))
     #directory for signal fits
     os.system("mkdir %s/%s/SigFits"%(workDir,subdir))
     #directory for ABCD plots 
     os.system("mkdir %s/%s/ABCDplots"%(workDir,subdir))
     #Directory for datacards
     os.system("mkdir %s/%s/Datacards"%(workDir,subdir))
     outSubDir = workDir +"/"+subdir
     sig = np.zeros((3,5)) #
     LimComb = np.zeros((5))
     
     #fits to signal: produces a Signal_mass"+np.int(mass)+"_"+np.int(ctau)+"_features.txt file
     #file is positioned in working dir
     #file format is 'ChannelIdx LxyBinIdx SignalYield SignalYield_error Signal_sigma' 
     stepMonitoring("Fit to signal")
     gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/AllFits.C')
     ROOT.AllFits(Sigpath+hnl_signal,Mass,ctau,chIdx,outSubDir,selection)
     stepMonitoring("bakcground estimates: ABCD")
     gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/Data_ABCD.C')
    # ABCD method for background estimates: performs ABCD in bins of channel and displacement and produces an output .txt file
     ROOT.Data_ABCD(background,0,True, True, Mass,ctau,int(chIdx),outSubDir,selection,False)
     stepMonitoring("producing Datacards...")
     sigOutFile ="Signal_"+hnl_sigFeatures[0]+"_"+hnl_sigFeatures[1]+"_features.txt"
     print(sigOutFile)
     os.system("python "+Adir+"/python/Card_producer.py %s %s %f %f %s %s %s %s"%(sigOutFile,"/BkgYields_.txt",Mass,ctau,outSubDir,channel,lumi,chargeLables[idx])) 
     print(chargeLables[idx])
   #  execfile(Adir+"/python Card_producer.py %s %s %f %f %s"%(sigOutFile,"/Bkg_Yields.txt",Mass,ctau,workDir))
     Cards = glob.glob(outSubDir+"/Datacards/*.txt")	
     Cards_lable = [ l.strip(".txt") for l in Cards]
     Cards_lable = [ l.split("/") for l in Cards_lable]
     print(Cards_lable)
     Cards_lable = [Cards_lable[l][9] for l in range(0,len(Cards))]
     stepMonitoring("Datacards summary")
     print(Cards_lable)
     print(chargeLables[idx])
     for k in range(0,len(Cards)):
     	os.system(" combine -M AsymptoticLimits --run blind %s  "%(Cards[k]))
     	os.system(" mv higgsCombineTest.AsymptoticLimits.mH120.root %s/Lim_%s.root"%(outSubDir,Cards_lable[k])) 
     	c= TChain("limit")
     #	print(chargeLables[idx])
     	c.Add("%s/Lim_%s.root"%(outSubDir,Cards_lable[k]))
   	print("%s/Lim_%s.root"%(outSubDir,Cards_lable[k]))
     	for s in range(0,len(lxyLables)):
     	#	print(Cards_lable[k])
     	#	print(lxyLables[s])
     		match = Cards_lable[k].find(lxyLables[s])
     		if match!=-1:
     			lxyIdx = s	
    			for LimIdx in range(0,c.GetEntries()):
   
     				 c.GetEntry(LimIdx)
   				 print("Idx of displ bin %d "%lxyIdx)
     				 sig[lxyIdx][LimIdx]=c.GetLeaf("limit").GetValue(0)
     print(chargeLables[idx])
     with open("%s/Limits.txt"%outSubDir, "a") as limFile:
           for i in range(0,len(sig)):

           	limFile.write("%d %f %f %f %f %f\n"%(i,sig[i][2],sig[i][2]-sig[i][1],sig[i][3]-sig[i][2],sig[i][2]-sig[i][0],sig[i][4]-sig[i][2]))
     
     combinedString= ""
     print(chargeLables[idx])
     for k in range(0,len(Cards_lable)): # prepare string containing datacards for different lxy bins to be combined
           	 # prepare string containing datacards for different lxy bins to be combined
           	print(Cards_lable[k])
           	print(channel)
           	match =Cards_lable[k].find(channel)
           #	print(Cards_lable[k])
           	print(match)
           	if match!=-1:
           		combinedString =  combinedString+Cards[k]+" " 
     print(combinedString)
     print(idx)
     #combine datacards and extract combined limit for a singl mass and |V|^2 point 
     os.system("mkdir "+ workDir+"/CombinedCards")
     combinedCard = workDir+"/CombinedCards/HNL_M"+np.str(np.int(Mass)) +"_ctau"+np.str(np.int(ctau))+"_"+channel+"_"+chargeLables[idx]+"_combined.txt"
     print(chargeLables[idx])
     print("python ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "+ combinedString+" >> " +combinedCard)
     os.system(" python ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "+ combinedString+" >> " +combinedCard)
     os.system(" combine -M AsymptoticLimits --run blind  %s "%(combinedCard))
     os.system(" mv higgsCombineTest.AsymptoticLimits.mH120.root %s.root"%(combinedCard.strip(".txt"))) 
     c= TChain("limit")
     c.Add("%s.root"%(combinedCard.strip(".txt")))
     for limIdx in range(0,c.GetEntries()):
     	c.GetEntry(LimIdx)
     	LimComb[LimIdx]=c.GetLeaf("limit").GetValue(0)
     with open("%s/LimitsCombined.txt"%workDir, "a") as limFile:
     
     	limFile.write("%f %f %f %f %f %f \n"%(HNLUtils.getVV(Mass,ctau,True),LimComb[2],LimComb[2]-LimComb[1],LimComb[3]-LimComb[2],(LimComb[2]-LimComb[0]),LimComb[4]-LimComb[2]))


chLables = ["Mu","PF","LowPt","Track"]
#	lxyLables = ["LxyUnder3","LxyOver3_Under10","LxyOver10_Under20","LxyOver20"]
lxyLables = ["LxyUnder1","LxyOver1_Under5","LxyOver5"]
outdir = sys.argv[1] #name of the directory containing all the info for the analysis for a single mass point
sel = sys.argv[2] #to be revised 
channel = sys.argv[3] #to be revised 
part = "1D"
Adir = '/cmshome/ratramon/Analysis/' #parent directory for the analysis
#Bkgpath = "/cmshome/ratramon/Analysis/datasets/Data/1A"
Sigpath = "/cmshome/ratramon/Analysis/datasets/Signal/" #Signal samples directory 
#prepare signals
with open("Signals.txt") as fSig:
	hnl_signal = fSig.readlines() #

hnl_signal = [l.strip("\n") for l in hnl_signal]
hnl_sigFeatures = [l.strip(".root") for l in hnl_signal] #prepare tag for signals
print(hnl_sigFeatures)
#lable for output files
sig_FileLables = ["HNL_%s"%(l) for l in hnl_sigFeatures]
hnl_sigFeatures = [l.split("_") for l in hnl_sigFeatures]
#masses

Mass = [hnl_sigFeatures[i][0].strip("Mass") for i in range(0,len(hnl_sigFeatures)) ] # strip string structure
print(Mass)
Mass = [np.float(Mass[i].replace("p",".")) for i in range(0,len(hnl_sigFeatures)) ] #correct format +  cast to float
#ctaus
ctau = [hnl_sigFeatures[i][1].strip("ctau") for i in range(0,len(hnl_sigFeatures)) ] # strip string structure
print(ctau)
ctau = [np.float(ctau[i].replace("p",".")) for i in range(0,len(hnl_sigFeatures)) ]#correct format + cast to float

#formatted lables
sig_OutLables = ["HNL Mass = %.1f GeV, c$tau = %.1f mm"%(Mass[i],ctau[i]) for i in range(0,len(Mass)) ] 
sig = np.zeros((3,5)) #
LimComb = np.zeros((5))
#prepare data
fullLumi = 41.6
if part == "1A":

	lumi = 0.76
	background = Adir+"/datasets/Data/1A_LikelihoodV6" 
if part == "1D":

	lumi = 1.793843944 
	background = Adir+"/datasets/Data/1D" 

# cp file to working area
#Create output directory
workDir = Adir+"/output/"+outdir
os.system("mkdir %s"%(workDir))
os.system("cd %s"%(workDir))
gROOT.SetBatch(True)
chLables = ["Mu","PF","LowPt","Track"]
#	lxyLables = ["LxyUnder3","LxyOver3_Under10","LxyOver10_Under20","LxyOver20"]
lxyLables = ["LxyUnder1","LxyOver1_Under5","LxyOver5"]
#subdirectories for different |V|^2 for signal fits
#for subdir in hnl_signal:
#	LxyBinnedLimits(Sigpath,background,workDir,subdir,sel,channel,lumi) #computes datacards+limits for the n-displacement bins (no specified coupling scenario hypothesis) 
#gSystem.Unload('/cmshome/ratramon/Analysis/macros/Data_ABCD.C')								   
#gSystem.Unload('/cmshome/ratramon/Analysis/macros/AllFits.C')								   
gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/YieldsPlotter.C')
ROOT.YieldsPlotter(workDir,"OS")
ROOT.YieldsPlotter(workDir,"SS")



  #combine datacards to obtain combined limit of the n-displacement categories 
#	for k in range(0,len(Cards)):
#		os.system(" combine -M AsymptoticLimits --run blind %s  "%(Cards[k]))
#		os.system(" mv higgsCombineTest.AsymptoticLimits.mH120.root %s/Lim_%s.root"%(outSubDir,Cards_lable[k])) 
#		c= TChain("limit")
#		c.Add("%s/Lim_%s.root"%(outSubDir,Cards_lable[k]))
#		for s in range(0,1):#len(chLables)):
#			match =Cards_lable[k].find(chLables[s])
#	#		print(Cards_lable[k])
#	#		print(match)
#			if match!=-1:
#				chIdx = s
#		#		print(chIdx)
#				combinedString[s] = combinedString[s] + " " +Cards[k]
#				
#		for s in range(0,len(lxyLables)):
#			print(Cards_lable[k])
#			print(lxyLables[s])
#			match = Cards_lable[k].find(lxyLables[s])
#			if match!=-1:
#				lxyIdx = s	
#		temp = []
#		for idx in range(0,c.GetEntries()):
#			c.GetEntry(idx)
#			sig[chIdx][lxyIdx][idx]=c.GetLeaf("limit").GetValue(0)
#			#print(temp[idx])
#	#	print(sig[k/4][k%4])
#	print(combinedString)						
#	for k in range(0, len(combinedString)):
#		os.system("mkdir "+ workDir+"/Datacards/CombinedCards")
#		combinedCard = workDir+"/Datacards/CombinedCards/HNL_M"+np.str(np.int(Mass[i])) +"_ctau"+np.str(np.int(ctau[i]))+"_"+chLables[k]+"_combined.txt"
#		print(combinedString[k])
#		os.system("python ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "+ combinedString[k]+" > " +combinedCard)
#		os.system(" combine -M AsymptoticLimits --run blind  %s "%(combinedCard))
#		os.system(" mv higgsCombineTest.AsymptoticLimits.mH120.root %s.root"%(combinedCard.strip(".txt"))) 
#		c= TChain("limit")
#		c.Add("%s.root"%(combinedCard.strip(".txt")))
#		for idx in range(0,c.GetEntries()):
#			c.GetEntry(idx)
#			LimComb[k][idx]=c.GetLeaf("limit").GetValue(0)
#	with open("%s/Limits.txt"%workDir, "w") as limFile:
#		for i in range(0,len(sig)):
#			
#			for j in range(0,len(sig[0])):
#			
#    				limFile.write("%d %d %f \n"%(i,j,sig[i][j][2]))
#				
#	
#	with open("%s/CombinedLimits.txt"%workDir, "w") as limFile:
#		for i in range(0,len(sig)):
#			
#    			limFile.write("%d  %f \n"%(i,LimComb[i][2]))

#from optparse import OptionParser
#...
#parser = OptionParser()
#parser.add_option("-C", "--CouplingDependent",
#                  action="store", type="bool", dest="explicitCoupling"
#		  help="project limits to specific coupling fractions scenarios")
#parser.add_option("-c", "--combine",
#                  action="store", type="bool", dest="doCombination"
#		  help="combine OS ad SS card and produce limits plot")
#parser.add_option("-y", "--yields",
#                  action="store", type="bool", dest="doYieldsPlots"
#		  help="produce plots for signal and bkg yields superimposed")
#parser.add_option("-l", "--lumi",
#                  action="store", type="float", dest="lumi"
#		  help="processed luminosity")
#parser.add_option("-d", "--directory",
#                  action="store", type="string", dest="outdir"
#		  help="name of output directory (all output dirs are under Analysis/output)")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")
#
#(options, args) = parser.parse_args()
#

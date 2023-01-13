import os
import re 
import sys
import ROOT
from ROOT import gROOT,TChain,gSystem 
import numpy as np
import glob
import time 
import HNLUtils
import csv
from multiprocessing import Process,Value,Array, Pool
print(time)

def stepMonitoring(step):
	
	
	print("\n")
	print("\n")
	print("***********************************************************")
	print("Process Time: %f              "%(time.time()))
	print("***********************************************************")
	print(" %s          "%step)
	print("***********************************************************")
	print("\n")
	print("\n")

#split csv input by mass independently of mass values 
def splitByMass(signals):
	
	massSplitted = [] 
	 
	singleMass=[]
	comment = 0 
	for i,feat in enumerate(signals):
#		print feat[0]
	 	if i==0:
	 		singleMass.append(feat)
	 		continue
#		back_check =
		back_check = signals[i-1]
#		print feat 
		if (feat[2] == back_check[2]):
			 singleMass.append(feat)
			
		if (feat[2] != back_check[2]):
			massSplitted.append(singleMass)
			singleMass = []	
	 		singleMass.append(feat)
		if (i==len(signals)-1):
			massSplitted.append(singleMass)

#	for mass in massSplitted:
#		print mass
#		print "\n"
	return massSplitted 

def LxyBinnedLimits(Sigpath,background,workDir,signal,sel,channel,lumi,counting,NNcut):

# chLables = ["Mu","PF","LowPt","Track"]
# chargeLables = ["OS","SS"]
# charge = ["&& LepQProd<0"," && LepQProd>0"]
## print "_____________________________________",Sigpath+signal
# chIdx = chLables.index(str(channel))
## print(chIdx)
  hnl_sigFeatures = ""
  Bc = False
  if 'Bc' in signal:
        Bc = True 
 # if ".root" in signal:
  print "signal",signal
  hnl_sigFeatures = signal.split("/") #prepare tag for signal
  print hnl_sigFeatures
  for string in hnl_sigFeatures:
	if "NewSignature_" in string:
	 	hnl_sigFeatures = string.strip("NewSignature_") #prepare tag for signal

  subdir= "HNL_%s_%s"%(hnl_sigFeatures,channel) 
#  print(hnl_sigFeatures)
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

 
 #  selection = sel + charge[idx]
 # subdir = sub+"_"+chargeLables[idx]
 
  stepMonitoring("Preparing output directory")
  os.system("mkdir %s/%s"%(workDir,subdir))

  #directory for signal fits
  os.system("mkdir %s/%s/SigFits"%(workDir,subdir))

  #directory for ABCD plots 
  os.system("mkdir %s/%s/ABCDplots"%(workDir,subdir))

  #Directory for datacards
  os.system("mkdir %s/%s/Datacards"%(workDir,subdir))
  outSubDir = workDir +"/"+subdir
  sig = np.zeros((6,5)) #
  LimComb = np.zeros((5))
 # print outSubDir
  
# #fits to signal: produces a Signal_mass"+np.int(mass)+"_"+np.int(ctau)+"_features.txt file
# #file is positioned in working dir
# #file format is 'ChannelIdx LxyBinIdx SignalYield SignalYield_error Signal_sigma' 


  # do counting experiment - deprecated 
# if counting == 1:   
# 	 stepMonitoring("Fit to signal- signal yield extraction")
# 	 ROOT.AllFits(Sigpath+signal,Mass,ctau,chIdx,outSubDir,selection)
# 	 stepMonitoring("Background yield extraction  ")
#        ROOT.Data_ABCD(background,0,True, True, Mass,ctau,int(chIdx),outSubDir,selection,False)


# #do shape analysis
# else:
#      stepMonitoring("Signal and background shape + yields extraction")
#      if ".root" in signal:
#       	sigIn = Sigpath+signal
#      else:
#       	sigIn = ""
#      ROOT.ShapeWSProducer(sigIn,".", Mass,ctau,outSubDir,sel,1,1, 0,NNcut)
#

  #datacards
  stepMonitoring("producing Datacards...")
  #sigOutFile ="Signal_"+hnl_sigFeatures[0]+"_"+hnl_sigFeatures[1]+"_features.txt"
  #sigOutFile ="Signal_Mass"+np.str(np.int(Mass))+"_ctau"+np.str(np.int(ctau))+"_features.txt"
  sigOutFile ="Signal_Mass"+('{:.1f}'.format(Mass)).replace('.','p')+"_ctau"+'{}p{}'.format(int(ctau),int((ctau-int(ctau))*1000))+"_features.txt"
 # print(sigOutFile)
  os.system("python "+Adir+"/python/Card_producer.py %s %s %f %f %s %s %s %d"%(sigOutFile,"/BkgYields_.txt",Mass,ctau,outSubDir,channel,lumi,counting)) 

  print outSubDir+"/Datacards/*.txt"
  Cards = glob.glob(outSubDir+"/Datacards/*.txt")	
  Cards_lable = [ l.strip(".txt") for l in Cards]
  Cards_lable = [ l.split("/") for l in Cards_lable]
  Cards_lable = [ Cards_lable[l][9] for l in range(0,len(Cards_lable))]
  print Cards_lable 

  stepMonitoring("Datacards summary")
  print(Cards_lable)
  lxyIdx = -1 
  for k in range(0,len(Cards)):
  	os.system("combine -M AsymptoticLimits --run blind %s  "%(Cards[k]))
  #	os.system("echo mv higgsCombineTest.AsymptoticLimits.mH120.root %s/Lim_%s.root"%(outSubDir,Cards_lable[k])) 
  	os.system("mv higgsCombineTest.AsymptoticLimits.mH120.root %s/Lim_%s.root"%(outSubDir,Cards_lable[k])) 
  for k in range(0,len(Cards)):
  	c= TChain("limit")
  #	print(chargeLables[idx])
  	c.Add("%s/Lim_%s.root"%(outSubDir,Cards_lable[k]))
 #       print("%s/Lim_%s.root"%(outSubDir,Cards_lable[k]))
  	for s in range(0,len(categories_lower)):
  	#	print(Cards_lable[k])
  	#	print(categories[s])
  		match = Cards_lable[k].find(categories_lower[s])
  		if match!=-1:
  			lxyIdx = s
        		print categories_lower[s],Cards_lable[k],lxyIdx	
      	for LimIdx in range(0,c.GetEntries()):
  
  			 c.GetEntry(LimIdx)
  #    			         print("Idx of displ bin %d "%lxyIdx)
        		 sig[lxyIdx][LimIdx]=c.GetLeaf("limit").GetValue(0)
        		 print lxyIdx,LimIdx,sig[lxyIdx][LimIdx]
        			 
  print sig
  print('opening limits file for writing ___________________________')
  with open("%s/Limits.txt"%outSubDir, "a") as limFile:
        for i in range(0,len(sig)):
        	print("%s %f %f %f %f %f\n"%(categories_lower[i],sig[i][2],sig[i][2]-sig[i][1],sig[i][3]-sig[i][2],sig[i][2]-sig[i][0],sig[i][4]-sig[i][2])) 
        	limFile.write("%s %f %f %f %f %f\n"%(categories_lower[i],sig[i][2],sig[i][2]-sig[i][1],sig[i][3]-sig[i][2],sig[i][2]-sig[i][0],sig[i][4]-sig[i][2]))
 

  #categories combinations per mass point 
  combinedString= ""
  for k in range(0,len(Cards_lable)): # prepare string containing datacards for different lxy bins to be combined
        	 # prepare string containing datacards for different lxy bins to be combined
   #     	print(Cards_lable[k])
    #    	print(channel)
        	match =Cards_lable[k].find('HNL_')
        #	print(Cards_lable[k])
     #   	print(match)
        	if match!=-1:
        		combinedString =  combinedString+Cards[k]+" " 
  print(combinedString)

  #combine datacards and extract combined limit for a singl mass and |V|^2 point 
# os.system("mkdir "+ workDir+"/CombinedCards")
# combinedCard = workDir+"/CombinedCards/HNL_m_"+"{:.2f}".format(Mass).replace('.','p')+"_ctau"+"{:.3f}".format(ctau).replace('.','p')+"_"+channel+"_combined.txt"
# LxyCombinedCards = combinedCard

# if combinedString!= "":
# #print("python ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "+ combinedString+" >> " +combinedCard)
#         print(" python ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "+ combinedString+" >> " +combinedCard)
#         os.system(" python ~/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "+ combinedString+" >> " +combinedCard)
#         os.system(" combine -M AsymptoticLimits --run blind  %s "%(combinedCard))
#         os.system(" mv higgsCombineTest.AsymptoticLimits.mH120.root %s.root"%(combinedCard.strip(".txt"))) 
#         c= TChain("limit")
#         c.Add("%s.root"%(combinedCard.strip(".txt")))


#         for LimIdx in range(0,c.GetEntries()):
# 		c.GetEntry(LimIdx)
# 		LimComb[LimIdx]=c.GetLeaf("limit").GetValue(0)

#         #results as a function of the couplings
#	  with open("%s/LimitsCombined.txt"%workDir, "a") as limFile:
# 
# 		limFile.write("%f %f %f %f %f %f \n"%(HNLUtils.getVV(Mass,ctau,True),LimComb[2],LimComb[2]-LimComb[1],LimComb[3]-LimComb[2],(LimComb[2]-LimComb[0]),LimComb[4]-LimComb[2]))

def singleMassPoint(signal):

        #lable for output files
        for ct_w in ct_reweights:
        	ct_formatted = str(ct_w).replace('.','p')
        	print signal[0]
        	temp_sig = signal[0][:]
        	temp_sig[1] = temp_sig[1].split('_')[0] + "_ctau"+ct_formatted
        	temp_sig[3] = str(ct_w)
        	temp_sig[4] = ' '
        	temp_sig[5] = '-99'
        	temp_sig[6] = '-99'
        	temp_sig[8] = temp_sig[8].replace(signal[0][1],temp_sig[1].split('_')[0] + "_ctau"+ct_formatted) 
        	temp_sig[9] = signal[0][8]
        	temp_sig[8] = temp_sig[9].replace('.root','')		
 		print temp_sig,signal[0][1]
 		print signal[0]
        	signal.append(temp_sig)
        	
        sig_FileLables = ["HNL_%s"%(l[1]) for l in signal]
        hnl_sigFeatures = [l[1].split("_") for l in signal]
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
        
        #Create output directory
        
 	print 'counting', counting
        workDir = Adir+"/output/Mass"+hnl_sigFeatures[0][0]+"_"+outdir
        os.system("mkdir %s"%(workDir))
        os.system("cd %s"%(workDir))
 #       print Mass	
        #background + discrete profiling - once per mass point
        os.system("mkdir "+workDir+"/Background")
        tag = "DiscreteProfiling"
# 	ROOT.ShapeBkg("", ".",Mass[0],float(signal[0][7]),10,"&& hnl_charge==0",workDir+"/Background",1, 0,NNcut)
 # 	os.system("python /cmshome/ratramon/CMSSW_10_2_13/src/flashggFinalFit/Background/submitFtest.py "+str(Mass[0])+" '"+workDir+"/Background/' "+tag)
        for subdir in signal:
 		LxyBinnedLimits("",background,workDir,subdir[9],sel,channel,lumi,counting,NNcut) #computes datacards+limits for the n-displacement bins (no specified coupling scenario hypothesis) 
	


 
lxyLables = ["LxySUnder50","LxySOver50Under150","LxySOver150"]
categories = ["LxySUnder50_OS","LxySOver50Under150_OS","LxySOver150_OS","LxySUnder50_SS","LxySOver50Under150_SS","LxySOver150_SS"]
categories_lower = ["lxysig0to50_OS","lxysig0to50_SS","lxysig50to150_OS","lxysig50to150_SS","lxysiggt150_OS","lxysiggt150_SS"]
outdir = sys.argv[1] #name of the directory containing all the info for the analysis for a single mass point
sel = sys.argv[2] #to be revised 
channel = sys.argv[3] #to be revised 
counting = np.int(sys.argv[4]) #to be revised 
NNcut = np.float(sys.argv[5]) #to be revised
doBc = np.int(sys.argv[6])
#sigma = np.float(sys.argv[6]) 

part = "1D"
Adir = '/cmshome/ratramon/Analysis/' #parent directory for the analysis
#Bkgpath = "/cmshome/ratramon/Analysis/datasets/Data/1A"
Sigpath = "/cmshome/ratramon/Analysis/datasets/Signal/" #Signal samples directory 
ct_reweights = [
# 0.001,
# 0.0015,
# 0.002,
# 0.003,
# 0.005,
# 0.007,
#  0.01,
# 0.015,
# 0.02,
# 0.03,
# 0.04,
# 0.05,
# 0.07,
# 0.1,
# 0.15,
# 0.2,
# 0.3,
# 0.4,
# 0.5,
# 0.7,
# 1., 
# 1.5,
# 2., 
# 3., 
# 4., 
# 5., 
# 7., 
# 10., 
# 15., 
# 20., 
# 30., 
# 40., 
# 50., 
# 70., 
# 100., 
# 150., 
# 200., 
# 300., 
# 400., 
# 500., 
# 700., 
# 1000., 
# 5000.,
# 10000.,
# 50000.,
 ]

#prepare signals
if doBc:
	csvIn = '../data/MC_Bc_datasets_InclusiveFilterEff.csv'
else: 
	csvIn = '../data/MC_datasets_InclusiveFilterEff.csv'

with open(csvIn) as f:
    reader = csv.reader(f,delimiter=' ')
    next(reader)
    next(reader)
    signals = list(reader)
#rint(signals)
#split signal by mass hypothesis
splitMasses = splitByMass(signals)
#print splitMasses
#with open("Signals.txt") as fSig:
#	hnl_signal = fSig.readlines() #
#print hnl_signal
#hnl_signal = [l.strip("\n") for l in hnl_signal]
#hnl_sigFeatures = [l.strip(".root") for l in hnl_signal] #prepare tag for signals
#print(hnl_sigFeatures)

#prepare data
fullLumi = 41.6
if part == "1A":

	lumi = 0.75
	background = Adir+"/datasets/Data/1A_nn_5vars" 
if part == "1D":

	lumi = 4.91 
	background = Adir+"/datasets/Data/1D" 

# cp file to working area
gROOT.SetBatch(True)
if counting == 1:
    gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/Data_ABCD.C')
    gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/AllFits.C')
    print("Will perform counting experiment limits")
else: 
    gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/ShapeWSProducer.C')
    gROOT.LoadMacro('/cmshome/ratramon/Analysis/macros/ShapeBkgWSProducer.C')
    print("Will perform parametric binned shape analysis ")

#    for sigs in splitMasses:
 # 	 singleMassPoint(sigs)

#    Tpool = Pool(6)
    Tpool = Pool(6)
    Tpool.map(singleMassPoint,splitMasses) 	
    Tpool.close()
    Tpool.join()


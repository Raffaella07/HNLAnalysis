#Datacards producer
import os
import sys
from ctypes import *
import ROOT
import numpy as np 
import glob
import math
from ROOT import TH1D
from operator import itemgetter
#import multiprocessing
from HNLUtils import BR_HNLmupion,BR_HNLelepion,getVV,BToMuX,BToEX,BcToMuX,BcToEX

def prefit_fromCard(Nsig,mass,ctau,sig,bkg,sigma,category,Bc):

	c = ROOT.TCanvas()
#	w = ROOT.RooWorkSpace()
	
	print sig
	w_sig = sig.Get("w")
	w_bkg= bkg.Get("multipdf")
		
	x = w_bkg.var("hnl_mass")
	x.setBins(50,'plotter')
	x.setBins(50,'plot_sig')
	
	sig_hist = w_sig.data("mcSet")
	bkg_hist = w_bkg.data("data_obs")
	sig_fit = w_sig.pdf("signal") 
	bkg_fit = w_bkg.pdf("env_pdf_0_2018_13TeV_bern1") 

	frame = x.frame()
	x.setRange("left", float(mass) - 10*float(sigma), float(mass)-2*float(sigma));
 	x.setRange("right", float(mass)+ 2*float(sigma), float(mass)+10*float(sigma));
# 	x.setRange("sig", float(mass)- 3*float(sigma), float(mass)+ 3*float(sigma));
	print "n_entries",sig_hist.sumEntries()
	bkg_hist.plotOn(frame, ROOT.RooFit.Name("data_left"),ROOT.RooFit.Binning("plotter"), ROOT.RooFit.CutRange("left"))
	bkg_hist.plotOn(frame, ROOT.RooFit.Name("data_right"),ROOT.RooFit.Binning("plotter"), ROOT.RooFit.CutRange("right"))
#	bkg_hist.plotOn(frame, ROOT.RooFit.Name("data"),ROOT.RooFit.CutRange("left,right"),ROOT.RooFit.Binning("plotter"))
	sig_hist.plotOn(frame, ROOT.RooFit.Name("sig_hist"),ROOT.RooFit.Binning("plot_sig"),ROOT.RooFit.CutRange("sig"),ROOT.RooFit.FillColor(ROOT.kOrange+7),ROOT.RooFit.FillStyle(3013),ROOT.RooFit.DrawOption("HIST"))

	sig_fit.plotOn(frame, ROOT.RooFit.Name("signal_fit"),ROOT.RooFit.LineColor(ROOT.kOrange+7),ROOT.RooFit.Normalization(Nsig, ROOT.RooAbsReal.NumEvent))
	bkg_fit.plotOn(frame, ROOT.RooFit.Name("data_fit"),ROOT.RooFit.Binning("plotter"), ROOT.RooFit.Range("left"), ROOT.RooFit.LineColor(ROOT.kBlue))
	bkg_fit.plotOn(frame, ROOT.RooFit.Name("data_fit"),ROOT.RooFit.Binning("plotter"), ROOT.RooFit.Range("right"), ROOT.RooFit.LineColor(ROOT.kBlue))
	
	
	frame.GetYaxis().SetTitleOffset(0.9)
	frame.GetYaxis().SetTitleFont(42)
	frame.GetYaxis().SetTitleSize(0.05)
	frame.GetYaxis().SetLabelSize(0.065)
	frame.GetYaxis().SetLabelSize(0.04)
	frame.GetYaxis().SetLabelFont(42)
	frame.GetXaxis().SetTitleOffset(0.9)
	frame.GetXaxis().SetTitleFont(42)
	frame.GetXaxis().SetTitleSize(0.05)
	frame.GetXaxis().SetLabelSize(0.065)
	frame.GetXaxis().SetLabelSize(0.04)
	frame.GetXaxis().SetLabelFont(42)

	frame.GetYaxis().SetTitle("Events")
	frame.GetXaxis().SetTitle("m_{HNL} (GeV)")
	frame.SetStats(0)
	frame.SetMinimum(0)
	frame.Draw()
	label = '{}/prefit_m_{}_ctau_{}_{}.pdf'.format(wd,'{:.2f}'.format(mass).replace('.','p'),'{:.3f}'.format(ctau).replace('.','p'),category)

	print label
	c.SaveAs(label)

HNLMass = np.float(sys.argv[2])
HNLctau = np.float(sys.argv[3])
sigmas = ["0.010","0.013","0.017","0.025","0.035","0.035"]
masses = [1,1.5,2,2.2,3,4.5,5]
enhanceFactors = [1,1,1,1,1,1]
sigma = 0.0014 +0.0074 * HNLMass
channel = sys.argv[5]
lumi =np.float( sys.argv[6])
lxy_in = sys.argv[8]
sig_in = sys.argv[9]
qsign = ["OS","SS"]
wd = sys.argv[4]
n = len(sys.argv)
fname = sys.argv[11:n]
fbkg =sys.argv[1]
counting =np.int(sys.argv[7]) 
#doDirac =np.int(sys.argv[10])
 
print(counting)
Bc = False
blind = False
genEvts= "/cmshome/ratramon/Analysis/python/genEvts.txt" 
#print(sigFeatures)
# array contains 4 columns:  signal channel, displacement category,MC signal yield,MC yield error,MC mean, MC sigma
if counting==1:
 bkgYields = np.loadtxt(wd+fbkg)
else:
 bkgYields = np.zeros((3,4)) #

channels = ["Mu","PF","LowPt","Track"]
lxy = ["LxySUnder50","LxySOver50Under150","LxySOver150"]
otherlxy = ["lxysig0to50","lxysig50to150","lxysiggt150"]


with open(genEvts) as fGen:
	hnl_gen = fGen.readlines()
bkgYields = np.zeros((3,4)) #


channels = ["Mu","PF","LowPt","Track"]
lxy = ["LxySUnder50","LxySOver50Under150","LxySOver150"]
otherlxy = ["lxysig0to50","lxysig50to150","lxysiggt150"]
corr_pNN = [1.,1,1]
syst_pNN = [1.15,1.05,1.10]

with open(genEvts) as fGen:
	hnl_gen = fGen.readlines()
			
		
Xsec_b = 572.0 * pow(10,9) #fb
f_u = 0.4
f_c = 0.00261 

v2 = getVV(HNLMass,HNLctau,True)
br_NtoEPi =BR_HNLmupion(HNLMass) 
br_NtoMuPi =BR_HNLelepion(HNLMass) 

if 'Bc' in wd:
	Bc = True 
#d = ROOT.RDataFrame("Events","../data/Bsidebands/*.root")
#d_sig = ROOT.RDataFrame("Events","../data/SigRegion/*.root")
#TF = 1 #compute transfer factor from CR to SR
FullLumi = 41.6 #fb -1
ProcessedLumi = lumi #fb -1
lhe_eff = 0.08244
oldNorm = False
if Bc:
 B_c = 'Bc_'
 Bs = ['Bc_']
 ff = f_u 
 br_MuNuInclusive = BcToMuX(HNLMass,1)
 br_ENuInclusive = BcToEX(HNLMass,1)
else:
 B_c = ''
 Bs = ['Bu_','Bd_','Bs_']
 ff = f_u
if oldNorm:
	Bs = ['']
print "do Bc samples",Bc
f = open("plot.txt", "a")
if HNLMass == 3.0:
	start =0 #skip low displacement category mass 3 GeV - no jpsi 
else:
	start = 0

for j in range(start,len(lxy)):
  for k in range(0,len(qsign)):
		#selects each single channel and lxy combination and extracts the background yields in the signal region (mu +- 2 sigma )
	#ABCD 
	
	if otherlxy[j] != lxy_in:
		continue
	if qsign[k] != sig_in:
		continue
#	br_MuNuInclusive = BToMuX(3,1)
#	br_ENuInclusive = BToEX(3,1)
#	br_NtoEPi =BR_HNLmupion(3) 
#	br_NtoMuPi =BR_HNLelepion(3) 
#	v2 = getVV(HNLMass,sigFeatures[j][1],True)	
	category =lxy_in+"_"+sig_in 
   #  	filename =(wd+"/Datacards/HNL_m_"+str(HNLMass)+"_ctau_"+str(HNLctau)+"_"+str(channel)+"_"+str(lxy[j])+"_"+qsign[k]+".txt")
     	filename ='{}/Datacards/HNL_m_{}_ctau_{}_cat_{}_{}.txt'.format(wd,'{:.2f}'.format(HNLMass).replace('.','p'),('{:.4f}'.format(HNLctau)).replace('.','p'),otherlxy[j],qsign[k])
	print filename
#	bkgfile = wd+"/../Background/DiscreteProfiling/window_m"+str(int(HNLMass))+".0_s"+str(sigma)+"_ns10/ws/CMS-BHNL_multipdf_"+otherlxy[j]+"_"+qsign[k]+".root"
	if Bc:
		bkgfile = '{}/../Background/DiscreteProfiling/window_m{:.2f}_s{:.3f}_ns10/ws/workspace_multipdf_bhnl_m_{}_cat_{}_{}{}.root'.format(wd,float(HNLMass),float(sigma),'{:.2f}'.format(HNLMass).replace('.','p'),otherlxy[j],qsign[k],"Bc_")
	else:
		bkgfile = '{}/../Background/DiscreteProfiling/window_m{:.2f}_s{:.3f}_ns10/ws/workspace_multipdf_bhnl_m_{}_cat_{}_{}{}.root'.format(wd,float(HNLMass),float(sigma),'{:.2f}'.format(HNLMass).replace('.','p'),otherlxy[j],qsign[k],"")
	bkg = ROOT.TFile.Open(bkgfile)	
	flag  = str(bkg.Get("multipdf"))
	if "multipdf" in flag:
		print  "multipdf found!"
	else:
		print "multipdf not found - discrete profiling failed, will not build datacard for this bin"
		continue
	multi = bkg.Get('multipdf')
	veto = bool(multi.pdf('gaus_veto'))
	print "is there SM resonance? ", veto
	print category	
	if counting==1:
		Yield = bkgYields[j][2]
	else:
		Yield = 1
#	print(Yield)
	Nsig = 0

	Nevents =0  
	for idx,B in enumerate(Bs):
		br_MuNu = BToMuX(HNLMass,1,B)
		br_ENu = BToEX(HNLMass,1,B)
		print br_MuNu,br_ENu
		print BcToMuX(HNLMass,1),BcToEX(HNLMass,1)
		if Bc:
			sigfile = '{}/SigFits/workspace_signal_{}bhnl_m_{}_ctau_{}_cat_{}_{}.root'.format(wd,'Bc_','{:.2f}'.format(HNLMass).replace('.','p'),('{:.4f}'.format(HNLctau)).replace('.','p'),otherlxy[j],qsign[k]) 
		else:
			sigfile = '{}/SigFits/workspace_signal_{}bhnl_m_{}_ctau_{}_cat_{}_{}.root'.format(wd,'Bu_','{:.2f}'.format(HNLMass).replace('.','p'),('{:.4f}'.format(HNLctau)).replace('.','p'),otherlxy[j],qsign[k]) 
		sig = ROOT.TFile.Open(sigfile)

		if "signal" in str(sig):
			print  "sig shape found!"
		else:
			print "signal not found - low stats bin  will not build datacard for this bin"
			continue
 
		featPath =  wd+"/"+str(fname[idx].strip("[,]"))
		print fname,featPath
#		else: 
#			print fname
#			featPath =  wd+"/"+str(fname.strip("[,]"))
#			print fname,featPath
		if os.path.isfile(featPath):
			sigFeatures = np.loadtxt(featPath) 
			print featPath
		else:
			continue
		if sigFeatures[3]==-1 or sigFeatures[2]==0: 
			continue
	#signal yield is computed here WITHOUT CHOOSING A SPECIFIC COUPLING SCENARIO--for the coupling scenario choice see..

		sigEntries =1.2* sigFeatures[2]*corr_pNN[j] #1.2 for MC matching corrections + pNN correction applied flat for different B flavs 
		
		print sigEntries			 
		Nsig = Nsig + sigEntries * FullLumi * Xsec_b/ff* (br_MuNu *   v2   *  br_NtoEPi+br_ENu *   v2   *  br_NtoMuPi)
		Nevents = Nevents+sigFeatures[6]	
	print "Nsig ", Nsig 
 	rate = Nsig/Nevents
	#prefit_fromCard(Nsig,HNLMass,HNLctau,sig,bkg,sigma,otherlxy[j]+'_'+qsign[k],B_c) 
	rsh = open(filename, "w")
 	if counting==1:
 	  rsh.write("#Counting experiment,Signal M %d GeV, ctau %d cm, datacard for %s channel, %s \n"%(HNLMass,HNLctau,channel,lxy[j]))
 	else:
 	  rsh.write("#Shape analysis,Signal M %d GeV, ctau %d cm, datacard for %s channel, %s \n"%(HNLMass,HNLctau,channel,lxy[j]))
 	rsh.write("imax 1  number of channels\n")
 	rsh.write("jmax * number of backgrounds\n")
 	rsh.write("kmax * number of nuisance parameters\n")
 	if  counting==0 :
 	  rsh.write("----------------------------------------------------------------------------------------------------------------------------------------------------\n")
 	  rsh.write("shapes signal m%s_ctau%s_%s_%s_%s_%s %s w:signal\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,sigfile))
 	  rsh.write("shapes background m%s_ctau%s_%s_%s_%s_%s %s multipdf:CMS_hgg_%s_%s%s_2018_13TeV_bkgshape\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,bkgfile,otherlxy[j],qsign[k],B_c))
	  if veto:
 	 	 rsh.write("shapes veto m%s_ctau%s_%s_%s_%s_%s %s multipdf:gaus_veto\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,bkgfile))

 	  rsh.write("shapes data_obs m%s_ctau%s_%s_%s_%s_%s %s multipdf:data_obs\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,bkgfile))
 	rsh.write("----------------------------------------------------------------------------------------------------------------------------------------------------\n")
 	rsh.write("bin			m%s_ctau%s_%s_%s_%s_%s\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
 	if blind:
		rsh.write("observation		-1\n") #blinded
	else:
	 	rsh.write("observation		%f\n"%multi.data('data_obs').sumEntries()) #blinded
 	rsh.write("----------------------------------------------------------------------------------------------------------------------------------------------------\n")
	if veto:
	 	rsh.write("bin  						m%s_ctau%s_%s_%s_%s_%s			m%s_ctau%s_%s_%s_%s_%s		m%s_ctau%s_%s_%s_%s_%s\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
 		rsh.write("process						signal				 		background		veto\n")
 		rsh.write("process						-1			 			1			 2\n")
 		rsh.write("rate   						%f			 		%f		%f\n"%(Nsig,Yield,Yield))
 		rsh.write("-----------------------------------------------------------------------------------------------------------------------------------------------------\n")
   #	rsh.write("lumi				  		lnN		1.025	 			 		-\n")
  		rsh.write("syst_sig_trigger_sf_m%s_ctau%s_%s_%s_%s_%s				lnN		1.05	 			 		-				-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_muid_sf_m%s_ctau%s_%s_%s_%s_%s				lnN		1.01	 			 		-				-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_eleid_sf_m%s_ctau%s_%s_%s_%s_%s				lnN		1.03	 			 		-				-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_track_eff_m%s_ctau%s_%s_%s_%s_%s				lnN		1.05	 			 		-				-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_shape_m%s_ctau%s_%s_%s_%s_%s				lnN		1.15	 	-		 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_stat_m%s_ctau%s_%s_%s_%s_%s		gmN %d		%f	 			 		-				-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,Nevents,rate))
  		rsh.write("syst_sig_genmatch_m%s_ctau%s_%s_%s_%s_%s				lnN		1.05 			 		-		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_sel_m%s_ctau%s_%s_%s_%s_%s		lnN		1.15	 		-	-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_norm		lnN		1.15	 			 		-				1.15\n")
  		rsh.write("syst_sig_fc   		lnN		%f	 			 		-				-\n"%(1.0 if not Bc else 1.24))
 # 	rsh.write("syst_bkg_m%s_ctau%s_%s_%s_%s		lnN		-	 			 		1.1\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k]))
  	else:
	 	rsh.write("bin  						m%s_ctau%s_%s_%s_%s_%s			m%s_ctau%s_%s_%s_%s_%s\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
 		rsh.write("process						signal				 		background\n")
 		rsh.write("process						-1			 			1\n")
 		rsh.write("rate   						%f			 		%f\n"%(Nsig,Yield))
 		rsh.write("-----------------------------------------------------------------------------------------------------------------------------------------------------\n")
   #	rsh.write("lumi				  		lnN		1.025	 			 		-\n")
  		rsh.write("syst_sig_trigger_sf_m%s_ctau%s_%s_%s_%s_%s				lnN		1.05	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_muid_sf_m%s_ctau%s_%s_%s_%s_%s				lnN		1.01	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_eleid_sf_m%s_ctau%s_%s_%s_%s_%s				lnN		1.03	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_track_eff_m%s_ctau%s_%s_%s_%s_%s				lnN		1.05	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_shape_m%s_ctau%s_%s_%s_%s_%s				lnN		1.15	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_stat_m%s_ctau%s_%s_%s_%s_%s					gmN %d		%f	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c,Nevents,rate))
  		rsh.write("syst_sig_genmatch_m%s_ctau%s_%s_%s_%s_%s				lnN		1.05	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_sel_m%s_ctau%s_%s_%s_%s_%s					lnN		1.15	 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],B_c))
  		rsh.write("syst_sig_norm							lnN		1.15	 		-\n")
  		rsh.write("syst_sig_fc						   		lnN		%f	 		-\n"%(1.0 if not Bc else 1.24))
	rsh.write("-----------------------------------------------------------------------------------------------------------------------------------------------------\n")
 	if counting==1:
  	  rsh.write("m%s_ctau%s_%s_%s_%s autoMCStats 0 0 1\n"%(HNLMass,HNLctau,channel,lxy[j], qsign[k]))
 	else:
  #	  rsh.write("lumiScale_m%s_ctau%s_%s_%s_%s rateParam * * %f [%f,%f]\n"%(HNLMass,HNLctau,channel,otherlxy[j], qsign[k],FullLumi/ProcessedLumi,FullLumi/ProcessedLumi,FullLumi/ProcessedLumi))
  	  rsh.write("pdfindex_%s_%s%s_2018_13TeV discrete \n"%(otherlxy[j],qsign[k],B_c))
 	
 	rsh.close()	

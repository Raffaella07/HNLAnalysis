#Datacards producer

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
import CMSlumi

def prefit_fromCard(Nsig,mass,ctau,sig,bkg,sigma,category):

	ROOT.gROOT.SetBatch()	
	c = ROOT.TCanvas()
#	w = ROOT.RooWorkSpace()
	print sig
	w_sig = sig.Get("w")
	w_bkg= bkg.Get("multipdf")
		
	x = w_bkg.var("hnl_mass")
	x.setBins(100,'plotter')
	x.setBins(50,'plot_sig')
	
	sig_hist = w_sig.data("mcSet")
	bkg_hist = w_bkg.data("data")
	sig_fit = w_sig.pdf("signal") 
	bkg_fitlist = w_bkg.allVars()
	bkg_fitlist.Print()
#	print type(bkg_fitlist)
#	for fit in bkg_fitlist:
	#	if "2018_13TeV_bern1" in fit.GetTitle():
	bkg_fit = w_bkg.pdf('CMS_hgg_'+category+'_2018_13TeV_bkgshape')
	x.setRange("left", float(mass) - 10*float(sigma), float(mass)-3*float(sigma));
 	x.setRange("right", float(mass)+ 3*float(sigma), float(mass)+10*float(sigma));
 	x.setRange("sig", float(mass)- 3*float(sigma), float(mass)+ 3*float(sigma));

	frame = x.frame()
	print "n_entries",sig_hist.sumEntries(),Nsig
#	bkg_hist.plotOn(frame, ROOT.RooFit.Name("data_left"),ROOT.RooFit.Binning("plotter"), ROOT.RooFit.Range("left"))
	bkg_hist.plotOn(frame, ROOT.RooFit.Name("data"),ROOT.RooFit.CutRange("left,right"),ROOT.RooFit.Binning("plotter"))
	sig_hist.plotOn(frame, ROOT.RooFit.Name("sig_hist"),ROOT.RooFit.Binning("plot_sig"),ROOT.RooFit.Rescale(Nsig*1.0/sig_hist.sumEntries()),ROOT.RooFit.CutRange("sig"),ROOT.RooFit.FillColor(ROOT.kOrange+7),ROOT.RooFit.LineColor(ROOT.kOrange+7),ROOT.RooFit.FillStyle(3003),ROOT.RooFit.DrawOption("F"))

	sig_fit.plotOn(frame, ROOT.RooFit.Name("signal_fit"),ROOT.RooFit.LineColor(ROOT.kOrange+7),ROOT.RooFit.Normalization(Nsig,ROOT.RooAbsReal.NumEvent))
	bkg_fit.plotOn(frame, ROOT.RooFit.Name("data_fit"),ROOT.RooFit.Range("left,right"),ROOT.RooFit.Normalization(w_bkg.var('CMS_hgg_'+category+'_2018_13TeV_bkgshape_norm').getVal(),ROOT.RooAbsReal.NumEvent))
#	bkg_fit.plotOn(frame, ROOT.RooFit.Name("data_fit"),ROOT.RooFit.Binning("plotter"), ROOT.RooFit.Range("right"), ROOT.RooFit.LineColor(ROOT.kBlue))
	
	frame.SetTitle("")	
	frame.GetYaxis().SetTitleOffset(0.9)
	frame.GetYaxis().SetTitleFont(42)
	frame.GetYaxis().SetTitleSize(0.0)
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
	frame.SetMaximum(max(w_bkg.var('CMS_hgg_'+category+'_2018_13TeV_bkgshape_norm').getVal()/25,Nsig/8.))
	frame.Draw()
	CMSlumi.CMSlumi(c,iPeriod=5,iPosX=11,lumiText='',extraText='Preliminary '+category)
	label = '{}/prefit_m_{}_ctau_{}_{}.pdf'.format(wd,'{:.1f}'.format(mass).replace('.','p'),'{:.1f}'.format(ctau).replace('.','p'),category)

	print label
	c.SaveAs(label)

HNLMass = np.float(sys.argv[3])
HNLctau = np.float(sys.argv[4])
sigmas = ["0.010","0.013","0.017","0.025","0.035","0.035"]
masses = [1,1.5,2,3,4.5,5]
enhanceFactors = [0.01,1,1,1,1,1]
print(sigmas)
for i in range(0,len(sigmas)):
	print(i)
	print(sigmas[i])
	if masses[i]==HNLMass:
		sigma = sigmas[i] 
		eF = enhanceFactors[i]
channel = sys.argv[6]
lumi =np.float( sys.argv[7])
qsign = ["OS","SS"]
wd = sys.argv[5]
fname =sys.argv[1]
fbkg =sys.argv[2]
counting =np.int(sys.argv[8]) 
print(counting)
Bc = False 
genEvts= "/cmshome/ratramon/Analysis/python/genEvts.txt" 
sigFeatures = np.loadtxt(wd+"/"+fname) 
print(sigFeatures)
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
NgenArray = [l.split(" ") for l in hnl_gen]
print(NgenArray)
Ngen = -1
for i in range(0,len(NgenArray)):
	if (np.float(NgenArray[i][0]) == HNLMass):
       	  if (np.float(NgenArray[i][1]) == HNLctau):
		Ngen = np.float(NgenArray[i][2])
			
		
Xsec_b = 472.8 * pow(10,9) #fb
f_u = 0.4
f_c = 0.00261 
br_MuNuInclusive = BToMuX(HNLMass,1)
br_ENuInclusive = BToEX(HNLMass,1)
print(br_MuNuInclusive)
v2 = getVV(HNLMass,HNLctau,True)
print(v2)
br_NtoEPi =BR_HNLmupion(HNLMass) 
br_NtoMuPi =BR_HNLelepion(HNLMass) 
if 'Bc' in wd:
	Bc = True 
#d = ROOT.RDataFrame("Events","../data/Bsidebands/*.root")
#d_sig = ROOT.RDataFrame("Events","../data/SigRegion/*.root")
#TF = 1 #compute transfer factor from CR to SR
FullLumi = 4.91 #fb -1
ProcessedLumi = lumi #fb -1
lhe_eff = 0.08244
if Bc:
 ff = f_u/(f_c*lhe_eff) 
 br_MuNuInclusive = BcToMuX(HNLMass,1)
 br_ENuInclusive = BcToEX(HNLMass,1)
else:
 ff = f_u
print "do Bc samples",Bc
f = open("plot.txt", "a")
if HNLMass == 3.0:
	start =2 #skip low displacement category mass 3 GeV - no jpsi 
else:
	start = 0
for j in range(start,len(lxy)):
  for k in range(0,len(qsign)):
		#selects each single channel and lxy combination and extracts the background yields in the signal region (mu +- 2 sigma )
	#ABCD 
#	br_MuNuInclusive = BToMuX(3,1)
#	br_ENuInclusive = BToEX(3,1)
#	br_NtoEPi =BR_HNLmupion(3) 
#	br_NtoMuPi =BR_HNLelepion(3) 
#	v2 = getVV(HNLMass,sigFeatures[j][1],True)
	category =lxy[j]+"_"+qsign[k] 
   #  	filename =(wd+"/Datacards/HNL_m_"+str(HNLMass)+"_ctau_"+str(HNLctau)+"_"+str(channel)+"_"+str(lxy[j])+"_"+qsign[k]+".txt")
     	filename ='{}/Datacards/HNL_m_{}_ctau_{}p{}_cat_{}_{}.txt'.format(wd,'{:.1f}'.format(HNLMass).replace('.','p'),int(HNLctau),int((HNLctau-int(HNLctau))*1000),otherlxy[j],qsign[k])
	print filename
#	bkgfile = wd+"/../Background/DiscreteProfiling/window_m"+str(int(HNLMass))+".0_s"+str(sigma)+"_ns10/ws/CMS-BHNL_multipdf_"+otherlxy[j]+"_"+qsign[k]+".root"
	bkgfile = '{}/../Background/DiscreteProfiling/window_m{:.2f}_s{:.3f}_ns10/ws/workspace_multipdf_bhnl_m_{}_cat_{}_{}.root'.format(wd,float(HNLMass),float(sigma),'{:.1f}'.format(HNLMass).replace('.','p'),otherlxy[j],qsign[k])
	sigfile = '{}/SigFits/workspace_signal_bhnl_m_{}_ctau_{}p{}_cat_{}_{}.root'.format(wd,'{:.1f}'.format(HNLMass).replace('.','p'),int(HNLctau),int((HNLctau-int(HNLctau))*1000),otherlxy[j],qsign[k]) 
	bkg = ROOT.TFile.Open(bkgfile)	
	sig = ROOT.TFile.Open(sigfile)
	
	if "signal" in str(sig):
		print  "sig shape found!"
	else:
		print "signal not found - low stats bin  will not build datacard for this bin"
		continue 
	flag  = str(bkg.Get("multipdf"))
	if "multipdf" in flag:
		print  "multipdf found!"
	else:
		print "multipdf not found - discrete profiling failed, will not build datacard for this bin"
		continue 
	print category	
	if counting==1:
		Yield = bkgYields[j][2]
	else:
		Yield = 1
#	print(Yield)
	if sigFeatures[2*j+k][3]==-1 or sigFeatures[j+k][2]==0: 
		continue
#	#signal yield is computed here WITHOUT CHOOSING A SPECIFIC COUPLING SCENARIO--for the coupling scenario choice see..
	if channel == "Tracks":
		#        channel reco efficiency    *   lumi   *   sigmaB  *  BR(B->MuNuX)   * |V|eff *  BR(N->epi) * filter_eff
		Nsig = sigFeatures[2*j+k][2]/(Ngen)     * ProcessedLumi * Xsec_b/f_u*br_MuNuInclusive *   v2   * br_NtoEPi   * filter_eff
	else:
		#        channel reco efficiency    *   lumi   *   sigmaB  *  BR(B->MuNuX)   * |V|eff *  BR(N->epi) * filter_eff
		Nsig = eF * sigFeatures[2*j+k][2] * FullLumi * Xsec_b/ff* (br_MuNuInclusive *   v2   *  br_NtoEPi+br_ENuInclusive *   v2   *  br_NtoMuPi)

		print "Nsig ", Nsig 
#	print ("%f %f"%(v2,Nsig))
#	f.write("%f %f\n"%(v2,Nsig))
# 	print (wd+"/Yields_"+lxy[j]+".txt")
# 	with open(wd+"/../Yields_"+lxy[j]+"_"+qsign+".txt","a") as file:
 #		l = np.str(v2) + " " + np.str(Nsig) + " "+ np.str(Yield) +"\n"
 #		file.writelines(l)
 	
	prefit_fromCard(Nsig,HNLMass,HNLctau,sig,bkg,sigma,otherlxy[j]+'_'+qsign[k])
 	rsh = open(filename, "w")
 	if counting==1:
 	  rsh.write("#Counting experiment,Signal M %d GeV, ctau %d cm, datacard for %s channel, %s \n"%(HNLMass,HNLctau,channel,lxy[j]))
 	else:
 	  rsh.write("#Shape analysis,Signal M %d GeV, ctau %d cm, datacard for %s channel, %s \n"%(HNLMass,HNLctau,channel,lxy[j]))
 	rsh.write("imax 1  number of channels\n")
 	rsh.write("jmax 1 number of backgrounds\n")
 	rsh.write("kmax * number of nuisance parameters\n")
 	if  counting==0 :
 	  rsh.write("----------------------------------------------------------------------------------------------------------------------------------------------------\n")
 	  rsh.write("shapes signal m%s_ctau%s_%s_%s_%s %s w:signal\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],sigfile))
          
 	  #rsh.write("shapes background m%s_ctau%s_%s_%s_%s %s/../Background/DiscreteProfiling/window_m%.2f_s%s_ns10/ws/CMS-BHNL_multipdf_%s_%s.root multipdf:CMS_hgg_%s_%s_2018_13TeV_bkgshape\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],wd,HNLMass,sigma,otherlxy[j],qsign[k],otherlxy[j],qsign[k]))
 	  rsh.write("shapes background m%s_ctau%s_%s_%s_%s %s multipdf:CMS_hgg_%s_%s_2018_13TeV_bkgshape\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],bkgfile,otherlxy[j],qsign[k]))
 	  rsh.write("shapes data_obs m%s_ctau%s_%s_%s_%s %s multipdf:data\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],bkgfile))
 	rsh.write("----------------------------------------------------------------------------------------------------------------------------------------------------\n")
 	rsh.write("bin			m%s_ctau%s_%s_%s_%s\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k]))
 	rsh.write("observation		-1\n") #blinded
 	rsh.write("----------------------------------------------------------------------------------------------------------------------------------------------------\n")
 	rsh.write("bin  						m%s_ctau%s_%s_%s_%s			m%s_ctau%s_%s_%s_%s\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k],HNLMass,HNLctau,channel,lxy[j],qsign[k]))
 	rsh.write("process						signal				 		background\n")
 	rsh.write("process						-1			 			1\n")
 	rsh.write("rate   						%f			 		%f\n"%(Nsig,Yield))
 	rsh.write("-----------------------------------------------------------------------------------------------------------------------------------------------------\n")
#   	rsh.write("lumi				  		lnN		1.025	 			 		-\n")
  	rsh.write("syst_sig_trigger_sf				lnN		1.05	 			 		-\n")
  	rsh.write("syst_sig_muid_sf				lnN		1.01	 			 		-\n")
  	rsh.write("syst_sig_eleid_sf				lnN		1.03	 			 		-\n")
  	rsh.write("syst_sig_track_eff				lnN		1.05	 			 		-\n")
  	rsh.write("syst_sig_pheno				lnN		1.15	 			 		-\n")
  	rsh.write("syst_sig_norm_m%s_ctau%s_%s_%s_%s		lnN		1.15	 			 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k]))
  	rsh.write("syst_sig_sel_m%s_ctau%s_%s_%s_%s		lnN		1.2	 			 		-\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k]))
 # 	rsh.write("syst_bkg_m%s_ctau%s_%s_%s_%s		lnN		-	 			 		1.1\n"%(HNLMass,HNLctau,channel,lxy[j],qsign[k]))
  	rsh.write("-----------------------------------------------------------------------------------------------------------------------------------------------------\n")
 	if counting==1:
  	  rsh.write("m%s_ctau%s_%s_%s_%s autoMCStats 0 0 1\n"%(HNLMass,HNLctau,channel,lxy[j], qsign[k]))
 	else:
  #	  rsh.write("lumiScale_m%s_ctau%s_%s_%s_%s rateParam * * %f [%f,%f]\n"%(HNLMass,HNLctau,channel,otherlxy[j], qsign[k],FullLumi/ProcessedLumi,FullLumi/ProcessedLumi,FullLumi/ProcessedLumi))
  	  rsh.write("pdfindex_%s_%s_2018_13TeV discrete \n"%(otherlxy[j],qsign[k]))
 	
 	rsh.close()	

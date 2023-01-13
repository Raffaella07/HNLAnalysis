import os 
import numpy as np
import ROOT
import glob
from multiprocessing import Process,Value,Array, Pool
from multiprocessing.pool import ThreadPool as threadPool
from HNLUtils import BR_HNLmupion,BR_HNLelepion,getVV,BToMuX,BToEX

def SignificanceMax(sig):
	
	#define likelihood cut scan
	Lcut = np.arange(start=0.0,stop=1,step=0.02,dtype="float")
	significance = np.zeros(len(Lcut),dtype="float") 
	#sig and bkg rdataframes
	pathToFlat ="/cmshome/ratramon/Analysis/data/HNLFlatTuples/NewSignature_" 
        signal = ROOT.RDataFrame("Events",pathToFlat+str(sig[0])+"/HNLFlat_0.root")
	data_files =[ "/cmshome/ratramon/Analysis/data/HNLFlatTuples/Parking1A0_newSig_nn_section0/*.root","/cmshome/ratramon/Analysis/data/HNLFlatTuples/Parking1A0_newSig_nn_section1/*.root","/cmshome/ratramon/Analysis/data/HNLFlatTuples/Parking1A0_newSig_nn_section2/*.root","/cmshome/ratramon/Analysis/data/HNLFlatTuples/Parking1A0_newSig_nn_section3/*.root"]
	dataVec = ROOT.std.vector(ROOT.std.string)()
	for file_ in data_files:
		dataVec.push_back(file_)
		
	data = ROOT.RDataFrame("Events",dataVec)
	#filter data around given signal hypothesis
	print("%s %s"%(sig[0],"hnl_mass>"+str(sig[1]-10*sig[6])+ " && "+"hnl_mass>"+str(sig[1]+10*sig[6])))
	data = data.Filter("hnl_mass>"+str(sig[1]-10*sig[6])+ " && "+"hnl_mass>"+str(sig[1]+10*sig[6]))
	#utilities to compute Nsig and Nbkg
	Xsec_b = 472.8 * pow(10,9) #fb
	f_u = 0.4
	br_MuNuInclusive = BToMuX(sig[1],1)
	br_ENuInclusive = BToEX(sig[1],1)
	v2 = getVV(sig[1],sig[2],True)
	br_NtoEPi =BR_HNLmupion(sig[1]) 
	br_NtoMuPi =BR_HNLelepion(sig[1]) 
	ProcessedLumi = 41.6
	lumiscale = 41.6/0.75 #processed in used datasets 

	#cut based selection implementation
	mu_lowMass_OS = "( hnl_mass<2 && LepQProd<0\
							 && ( (fabs(hnl_lxy)<1 &&  hnl_pi_pt>0.7 &&  fabs(hnl_lxy_sig)>30 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>0 &&  fabs(1-hnl_cos2D)<2e-03  )\
							 ||  (fabs(hnl_lxy)>1  && fabs(hnl_lxy)<5 && hnl_pi_pt>0.7 &&  fabs(hnl_lxy_sig)>150 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>10 &&  fabs(1-hnl_cos2D)<1e-04  )\
							 || (fabs(hnl_lxy)>5 && hnl_pi_pt>2 &&  fabs(hnl_lxy_sig)>150 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>40 &&  fabs(1-hnl_cos2D)<3e-05  )\
							))" 

	mu_lowMass_SS = " (hnl_mass<2 && LepQProd>0\
							 && ( (fabs(hnl_lxy)<1 && hnl_pi_pt>0.7 &&  fabs(hnl_lxy_sig)>50 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>0 &&  fabs(1-hnl_cos2D)<2e-03  )\
							 ||  (fabs(hnl_lxy)>1  && fabs(hnl_lxy)<5 && hnl_pi_pt>0.7 &&  fabs(hnl_lxy_sig)>150 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>10 &&  fabs(1-hnl_cos2D)<3e-04  )\
							 || (fabs(hnl_lxy)>5 && hnl_pi_pt>2 &&  fabs(hnl_lxy_sig)>150 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>40 &&  fabs(1-hnl_cos2D)<3e-05  )\
							))" 
	mu_mediumMass_OS = " (hnl_mass>2 && hnl_mass<4.5 && LepQProd<0\
						  	 && ( (fabs(hnl_lxy)<1 && hnl_pi_pt>0.8 &&  fabs(hnl_lxy_sig)>90 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>5 &&  fabs(1-hnl_cos2D)<2e-03  )\
							 ||  (fabs(hnl_lxy)>1  && fabs(hnl_lxy)<5 && hnl_pi_pt>1.5 &&  fabs(hnl_lxy_sig)>300 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>40 &&  fabs(1-hnl_cos2D)<2e-04  )\
							 || (fabs(hnl_lxy)>5 && hnl_pi_pt>3.5 &&  fabs(hnl_lxy_sig)>300 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>40 &&  fabs(1-hnl_cos2D)<2e-05  )\
							))" 

	mu_mediumMass_SS = "  (hnl_mass>2 && hnl_mass<4.5 && LepQProd>0\
							 &&( (fabs(hnl_lxy)<1 && hnl_pi_pt>0.8 &&  fabs(hnl_lxy_sig)>100 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>5 &&  fabs(1-hnl_cos2D)<2e-03  )\
							 ||  (fabs(hnl_lxy)>1  && fabs(hnl_lxy)<5 && hnl_pi_pt>1.5 &&  fabs(hnl_lxy_sig)>300 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>60 &&  fabs(1-hnl_cos2D)<2e-04  )\
							 || (fabs(hnl_lxy)>5 && hnl_pi_pt>3.5 &&  fabs(hnl_lxy_sig)>300 && min(fabs(hnl_pi_dxyS),fabs(hnl_l_dxyS))>40 &&  fabs(1-hnl_cos2D)<3e-05  )\
							))" 
	mu_highMass_OS = " (hnl_mass>4.5 && LepQProd<0\
							 &&( (fabs(hnl_lxy)<1 && hnl_pi_pt>2 &&  fabs(hnl_lxy_sig)>100 )\
							 ||  (fabs(hnl_lxy)>1  && fabs(hnl_lxy)<5 && hnl_pi_pt>4 &&  fabs(hnl_lxy_sig)>100  )\
							 || (fabs(hnl_lxy)>5 && hnl_pi_pt>5 &&  fabs(hnl_lxy_sig)>200  )\
							))" 

	mu_highMass_SS = " (hnl_mass>4.5 && LepQProd>0\
							 && ( (fabs(hnl_lxy)<1 &&  hnl_pi_pt>2 &&  fabs(hnl_lxy_sig)>70 )\
							 ||  (fabs(hnl_lxy)>1  && fabs(hnl_lxy)<5 && hnl_pi_pt>4 &&  fabs(hnl_lxy_sig)>100  )\
							 || (fabs(hnl_lxy)>5 && hnl_pi_pt>5 &&  fabs(hnl_lxy_sig)>200  )\
							))"
	perCat_sel= selection +  " &&  hnl_charge==0 && ("+mu_lowMass_OS +"||"+mu_lowMass_SS +"||"+mu_mediumMass_OS +"||"+mu_mediumMass_SS +"||"+mu_highMass_OS +"||"+mu_highMass_SS +")" 
	N_mc = signal.Filter(perCat_sel).Count()
	Nsig = N_mc.GetValue() * 1.0/sig[4] * sig[5] * 2 * ProcessedLumi * Xsec_b/f_u*(br_MuNuInclusive *   v2   *  br_NtoEPi+br_ENuInclusive *   v2   *  br_NtoMuPi) 
	Ndata = data.Filter(perCat_sel).Count()
	Nbkg = Ndata.GetValue()* lumiscale
	print(np.sqrt(abs(2*(((Nsig+Nbkg)* np.log(1+Nsig/Nbkg))-Nsig)))) 
	sig_cutBased = np.sqrt(abs(2*(((Nsig+Nbkg)* np.log(1+Nsig/Nbkg))-Nsig)))
	
	for i,cut in enumerate(Lcut):
		N_mc = signal.Filter(selection+" &&  hnl_charge==0 && nn_score>"+str(cut)).Count()
		Nsig = N_mc.GetValue() * 1.0/sig[4] * sig[5] * 2 * ProcessedLumi * Xsec_b/f_u*(br_MuNuInclusive *   v2   *  br_NtoEPi+br_ENuInclusive *   v2   *  br_NtoMuPi) 
		Ndata = data.Filter(selection+" && hnl_charge==0 && nn_score>"+str(cut)).Count()
		Nbkg = Ndata.GetValue()* lumiscale
#		print(np.sqrt(abs(2*(((Nsig.GetValue()+Nbkg.GetValue())* np.log(1+Nsig.GetValue()/Nbkg.GetValue()))-Nsig.GetValue())))) 
		significance[i] = np.sqrt(abs(2*(((Nsig+Nbkg)* np.log(1+Nsig/Nbkg))-Nsig)))
	sig_cutBased = (sig_cutBased-significance[0])/significance[0] *100
	significance = (significance-significance[0])/significance[0] *100
	print(" %s cutB %f"%(sig[0],sig_cutBased))
	print(significance)
	print(Lcut)
	graph = ROOT.TGraph(len(Lcut-2),Lcut,significance)
	graph_point = ROOT.TGraph(1,np.ones(1),np.array(sig_cutBased))
	graph.SetTitle(sig[0])
	graph_point.SetTitle(sig[0]+"_cutBased")
	graph.SetMarkerStyle(8)
	graph_point.SetMarkerStyle(29)
	graph_point.SetMarkerSize(2)
	graph.SetLineWidth(2)
	return graph,graph_point
	


def ratioplot(branch):

	#qcd_samples=glob.glob('../data/HNLFlatTuples/Pt-*_newSig_genAnalysis/*1.root')
	data = glob.glob('../data/HNLFlatTuples/Parking1D0_newSig_nn_section*/*.root')

#	chain_qcd = ROOT.TChain('Events')
	chain_data = ROOT.TChain('Events')
#	for sample in qcd_samples:
#		chain_qcd.Add(sample)
	for sample in data:
		chain_data.Add(sample)

#	qcd = ROOT.RDataFrame(chain_qcd)
	data = ROOT.RDataFrame(chain_data);


	display = qcd.Display()
	display.Print()

#	qcd = qcd.Filter('hnl_charge!=0')
	data = data.Filter('hnl_charge!=0')

#	h_qcd = qcd.Histo1D(branch[0],branch[1],"QCDweight")
	h_data = data.Histo1D(branch[0],branch[1])
#	h_qcdNorm = ROOT.TH1D( h_qcd.DrawNormalized(""))
	h_dataNorm = h_data.DrawNormalized("")
#	h_qcd.SetTitle('qcd;'+branch_lable[]+';normalized to unity')
#	h_data.SetTitle('qcd;'+branch_lable+';normalized to unity')

	#ROOT.gROOT.LoadMacro("../macros/setStyle.C")
	#ROOT.setStyle()
	h_dataNorm.SetMarkerStyle(8)
	h_qcdNorm.SetFillStyle(1001)
	h_dataNorm.SetFillStyle(1001)
	h_dataNorm.SetLineWidth(2)
	h_qcdNorm.SetFillColor(ROOT.kBlue-7)
	h_qcdNorm.SetLineColor(ROOT.kBlue-7)
	
	
	rp = ROOT.TRatioPlot(h_qcdNorm, h_dataNorm);	
	c = ROOT.TCanvas('c','c',800,600)	
	rp.Draw()
	rp.SetLeftMargin(0.15)
#	ROOT.gStyle.SetTitleYOffset(1.5)# // => 1.15 if exponents
	rp.GetLowerRefGraph().SetMarkerStyle(8)	
	rp.GetLowerRefGraph().GetYaxis().SetRangeUser(0.2,2)	
	c.SaveAs("../plots/QCDvsData_"+branch[1]+".pdf")
	



def cutflow(sig):

	pathToNano="/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/"
	signal_lables = ["Mass 1.0 GeV c#tau 1000 mm","Mass 3.0 GeV c#tau 100 mm","Mass 4.5 GeV c#tau 1p0 mm","data"]
	
	presel1= ["trgmu_pt>7","abs(trgmu_eta)<1.5","sel_lep_pt>1.5","abs(sel_lep_eta)<2",\
			channel+"_pi_pt>0.7","abs("+channel+"_pi_eta)<2","abs("+channel+"_pi_dz)>0.005","abs("+channel+"_pi_dxy)>0.005",\
			"abs("+channel+"_pi_dzS)>1.5","abs("+channel+"_pi_dxyS)>3","abs("+channel+"_pi_DCASig)>5",\
			]
	presel_displaced_lep=["_dz>0.0015","_dxy>0.001","_dzS>1","_dxyS>1.5"]
	presel2 = [channel+"_sv_prob>0.001",channel+"_hnl_cos2D>0.995",channel+"_sv_lxy_sig>15",channel+"_mass<8"]
	baseline = [channel+"_hnl_charge==0","((trgmu_soft && "+channel+"_hnl_to_trgmu==0) ||(trgmu_loose && "+channel+"_hnl_to_trgmu==1))","sel_lep_mvaId>-3",channel+ "_pi_highPurityFlag",channel+"_sv_lxy_sig>0"]
	loss = np.zeros((30,2),dtype='float')
#	loss_cumulative = np.zeros((30,2),dtype='float')
	print(loss)

        df = ROOT.RDataFrame("Events",pathToNano+sig+"*11.root")
        df = df.Define("trgmu_pt","Take(Muon_pt,"+channel+"_trg_mu_idx)")
        df = df.Define("trgmu_eta","Take(Muon_eta,"+channel+"_trg_mu_idx)")
        df = df.Define("trgmu_dxy","Take(abs(Muon_dxy),"+channel+"_trg_mu_idx)")
        df = df.Define("trgmu_dxyS","Take(abs(Muon_dxyS),"+channel+"_trg_mu_idx)")
        df = df.Define("trgmu_dz","Take(abs(Muon_dz),"+channel+"_trg_mu_idx)")
        df = df.Define("trgmu_dzS","Take(abs(Muon_dzS),"+channel+"_trg_mu_idx)")
        df = df.Define("trgmu_soft","Take(abs(Muon_softId),"+channel+"_trg_mu_idx)")
        df = df.Define("trgmu_loose","Take(abs(Muon_looseId),"+channel+"_trg_mu_idx)")
	
	if (channel=="BToMuMuPi"):
		lep = "mu"
		collection = "Muon"
	elif (channel=="BToMuEPi"):
		lep = "ele"
		collection = "Electron"
		df = df.Define("sel_lep_mvaId",   "Take("+collection+"_pfmvaId,"  +channel+"_sel_"+lep+"_idx)")
	df = df.Define("sel_lep_pt",   "Take("+collection+"_pt,"  +channel+"_sel_"+lep+"_idx)")
	df = df.Define("sel_lep_eta",  "Take("+collection+"_eta," +channel+"_sel_"+lep+"_idx)")
	df = df.Define("sel_lep_dxy",  "Take(abs("+collection+"_dxy),"  +channel+"_sel_"+lep+"_idx)")
	df = df.Define("sel_lep_dxyS", "Take(abs("+collection+"_dxy/"+collection+"_dxyErr)," +channel+"_sel_"+lep+"_idx)")
	df = df.Define("sel_lep_dz",   "Take(abs("+collection+"_dz),"  +channel+"_sel_"+lep+"_idx)")
	df = df.Define("sel_lep_dzS",  "Take(abs("+collection+"_dz/"+collection+"_dzErr)," +channel+"_sel_"+lep+"_idx)")
#       df = df.Define("sel_"+b[1],channel+b[1]+"["+presel_loose+mc_match+"]")
	channel_sel = ["((sel_lep"+i+" && "+channel+"_hnl_to_trgmu==0) || (trgmu"+i+" && "+channel+"_hnl_to_trgmu==1))"for i in presel_displaced_lep]
	preselection = d = np.array([ presel1 + channel_sel+presel2+baseline]).T.flatten()
	print len(preselection), len(lables)
	#actual cutflow

	if ("BToHNL" in sig):
		match = channel+"_isMatched==1"
	else:
		match = channel+"_hnl_mass"
	for i in range (0,len(lables)):
		    if (i==0):
			 sel = preselection[i]
		   	 den = df.Define("den_size",channel+"_hnl_mass["+match+"].size()").Sum("den_size")
			 den0 = den
		   	 num = df.Define("num_size",channel+"_hnl_mass["+match+" && "+preselection[i]+"].size()").Sum("num_size")
	            else:

		   	 den = df.Define("den_size",channel+"_hnl_mass["+match+" && " +sel+"].size()").Sum("den_size")
			 sel += " && "+preselection[i]
		   	 num = df.Define("num_size"," "+channel+"_hnl_mass["+match+" && "+sel+"].size()").Sum("num_size")
		    loss[i][0] = (num.GetValue()*1.0/den.GetValue())*100
		    loss[i][1] =  (num.GetValue()*1.0/den0.GetValue())*100
#	print("total effieciency loss %f"%(num.GetValue()*1.0/den0.GetValue()-1)*100)
	return loss
	

def plotter(b):
	pathToNano="/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/"
	signals=["ParkingBPH1/crab_data_Run2018A_part1_section1/220809_102720/0000/",
		 "BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass1p0_ctau1000p0/220810_150846/0000/",
		 "BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass3p0_ctau100p0/220810_150913/0000/",
		 "BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass4p5_ctau1p0/220810_151253/0000/",
		]
	signal_lables = ["background (unblinded data)","Mass 1.0 GeV c#tau 1000 mm","Mass 3.0 GeV c#tau 100 mm","Mass 4.5 GeV c#tau 1 mm"]
	# build matrix of type 'histo model','var' for easy plotting
	
	presel_loose = channel+"_pi_pt>0.7 \
			 && abs("+channel+"_pi_eta)<2.5 \
		  	&& "+channel+"_trgmu_pt>7 \
			 && abs("+channel+"_trgmu_eta)<2 \
			 && "+channel+"_sel_lep_pt> 1.5 \
			 && abs("+channel+"_sel_lep_eta)<2.5"
	
	ROOT.gStyle.SetPalette(91)
        c = ROOT.TCanvas("c","c")
        l = ROOT.TLegend(0.4,0.6,0.7,0.9)
	l.SetHeader("Signal region, inclusive")
	Max = -99
        for i,sig in enumerate(signals):
        	
		if i!=0:
			mc_match=" && "+ channel+"_isMatched==1 "
		else:
			
			mc_match="  "
        	df = ROOT.RDataFrame("Events",pathToNano+sig+"*.root")
                df = df.Define(channel+"_trgmu_pt","Take(Muon_pt,"+channel+"_trg_mu_idx)")
        	df = df.Define(channel+"_trgmu_eta","Take(Muon_eta,"+channel+"_trg_mu_idx)")
                df = df.Define(channel+"_trgmu_dxy","Take(Muon_dxy,"+channel+"_trg_mu_idx)")
        	df = df.Define(channel+"_trgmu_dxyS","Take(Muon_dxyS,"+channel+"_trg_mu_idx)")
                df = df.Define(channel+"_trgmu_dz","Take(Muon_dz,"+channel+"_trg_mu_idx)")
        	df = df.Define(channel+"_trgmu_dzS","Take(Muon_dzS,"+channel+"_trg_mu_idx)")
		if (channel=="BToMuMuPi"):
			lep = "mu"
			collection = "Muon"
		elif (channel=="BToMuEPi"):
			lep = "ele"
			collection = "Electron"
		df = df.Define(channel+"_sel_lep_pt",   "Take("+collection+"_pt,"  +channel+"_sel_"+lep+"_idx)")
		df = df.Define(channel+"_sel_lep_eta",  "Take("+collection+"_eta," +channel+"_sel_"+lep+"_idx)")
		df = df.Define(channel+"_sel_lep_dxy",  "Take("+collection+"_dxy,"  +channel+"_sel_"+lep+"_idx)")
		df = df.Define(channel+"_sel_lep_pfmvaId",  "Take("+collection+"_pfmvaId,"  +channel+"_sel_"+lep+"_idx)")
		df = df.Define(channel+"_sel_lep_dxyS", "Take("+collection+"_dxy/"+collection+"_dxyErr," +channel+"_sel_"+lep+"_idx)")
		df = df.Define(channel+"_sel_lep_dz",   "Take("+collection+"_dz,"  +channel+"_sel_"+lep+"_idx)")
		df = df.Define(channel+"_sel_lep_dzS",  "Take("+collection+"_dz/"+collection+"_dzErr," +channel+"_sel_"+lep+"_idx)")
        	print(b[0])
        	print(presel_loose+mc_match)
        	df = df.Define(channel+"_sel_"+b[1],channel+b[1]+"["+presel_loose+mc_match+"]")
		n_endcaps = df.Filter("! "+channel+"_hnl_mass[abs("+channel+"_sel_lep_eta)>1.5"+mc_match+"].empty()").Count()
		n_tot = df.Filter("! "+channel+"_hnl_mass[abs("+channel+"_sel_lep_eta)>-1"+mc_match+"].empty()").Count()
		print("%d %s %f"%(i,signal_lables[i],n_endcaps.GetValue()*1.0/n_tot.GetValue()))
        	h = df.Histo1D(b[0],channel+"_sel_"+b[1])
		h.SetBinContent(h.GetNbinsX(),h.GetBinContent(h.GetNbinsX())+h.GetBinContent(h.GetNbinsX()+1))
        	if (i==0):	
			
        		h_0 = h.DrawNormalized("HISTF")
			h_0.SetLineWidth(2)
			h_0.SetFillStyle(3005)
			h_0.SetFillColor(ROOT.kGray+2)
			h_0.SetLineColor(ROOT.kGray+2)
			if (Max<h_0.GetMaximum()):
				Max = h_0.GetMaximum()
       			l.AddEntry(h_0,signal_lables[i],"fl")
  #      		h_0 = h.DrawNormalized("HISTPLC")
#			h_0.SetLineWidth(2)
#			Max =h_0.GetMaximum()
 #       		l.AddEntry(h_0,signal_lables[i],"l")
        	else:
        		h_copy = h.DrawNormalized("HISTsamePLC")
			h_copy.SetLineWidth(2)
			if (Max<h_copy.GetMaximum()):
				Max = h_copy.GetMaximum()
       			l.AddEntry(h_copy,signal_lables[i],"l")
        l.Draw()
	l.SetBorderSize(0)		
	l.SetFillStyle(0)
        l.SetTextSize(0.03)	
	h_0.SetMaximum(Max*1.5)
	h_0.GetXaxis().SetLabelSize(18)
	#c.SetLogy()
	ROOT.gStyle.SetLegendBorderSize(0)
	ROOT.gStyle.SetLegendFont(41)
	cms = ROOT.TLatex()
	cms.SetTextSize(0.04)
	cms.SetTextAlign(12)
	cms.DrawLatex(h_0.GetXaxis().GetBinLowEdge(1),Max*1.55,"CMS #bf{Preliminary}")
        #print("../plots/NewCentralValidation/loose"+str(b[1])+".pdf")
        c.SaveAs("../plots/NewCentralValidation/Zoom_"+channel+str(b[1])+"_log.pdf")

def FillRecoEffBins(channel,df_gen,df,h_num,pt_bin,lxy_bin,i,j,den_sel,num_sel,binValue,binError):

	ptBin_sel = channel+"_hnl_pt>"+str(pt_bin[j])+ " && "+channel+"_hnl_pt<"+str(pt_bin[j+1])
	lxyBin_sel = channel+"_sv_lxy>"+str(lxy_bin[i])+ " && "+channel+"_sv_lxy<"+str(lxy_bin[i+1])
	ptBin_gen_sel = "hnl_pt>"+str(pt_bin[j])+ " && hnl_pt<"+str(pt_bin[j+1])
	lxyBin_gen_sel = "Lxy>"+str(lxy_bin[i]*10)+ " && Lxy<"+str(lxy_bin[i+1]*10)
	denBin = df_gen.Filter(den_sel+ " && "+ptBin_gen_sel+" && "+ lxyBin_gen_sel).Count()
	numBin = df.Filter("! "+channel+"_hnl_mass["+num_sel+" && "+ptBin_sel+" && "+lxyBin_sel+"].empty()").Count()
	print(numBin.GetValue())
	print(denBin.GetValue())
	binValue[i] = numBin.GetValue()*1.0/denBin.GetValue()
	binError[i] = np.sqrt(pow(np.sqrt(numBin.GetValue())*1.0/denBin.GetValue(),2)+pow(numBin.GetValue()*np.sqrt(denBin.GetValue())/pow(denBin.GetValue(),2),2))


def RecoEfficiencies(path,genFile,channel):

	#reco efficiency binning
	
	lxy_bin= np.array([0.,1.,3.,5.,10.,15.,30.,50,100.],dtype=float)
	pt_bin= np.array([1,7,10,15,30,100],dtype=float)
	#channel depending selections 
	#define dataset
		   #&& (mu_fromHNL_pt>7 || mu_fromB_pt>7 )\
	den_sel = "pi_fromHNL_pt>0.7\
		   && abs(pi_fromHNL_eta)<2.5\
		   && mu_fromHNL_pt>1.5\
		   && mu_fromB_pt>1.5\
    		 && abs(mu_fromB_eta)<2.5"
	if (channel=="BToMuMuPi"):
		lepton_sel=" && abs(mu_fromB_pdgid)==13 && abs(mu_fromHNL_pdgid)==13"
	elif (channel=="BToMuEPi"):
		lepton_sel= "&& (( abs(mu_fromB_pdgid)==11 && abs(mu_fromHNL_pdgid)==13) || (abs(mu_fromB_pdgid)==13 && abs(mu_fromHNL_pdgid)==11))"
	den_sel = den_sel + lepton_sel
	num_sel = channel+"_pi_pt>0.7 \
		 && abs("+channel+"_pi_eta)<2.5 \
		 && "+channel+"_trgmu_pt>7 \
		 && abs("+channel+"_trgmu_eta)<2.5 \
		 && "+channel+"_sel_lep_pt> 1.5 \
		 && abs("+channel+"_sel_lep_eta)<2.5\
		 && "+channel+"_isMatched==1"

	ROOT.ROOT.EnableImplicitMT()
	ROOT.gStyle.SetOptStat(0000)
	ROOT.gROOT.LoadMacro("../macros/setStyle.C")
	ROOT.setStyle()
	df = ROOT.RDataFrame("Events",path+"*.root")
	df_gen = ROOT.RDataFrame("tree","/cmshome/ratramon/CMSSW_10_2_15_HNLsGen/src/HNLsGen/genLevelAnalysis/"+genFile)
	df = df.Define(""+channel+"_trgmu_pt","Take(Muon_pt,"+channel+"_trg_mu_idx)")
	df = df.Define(""+channel+"_trgmu_eta","Take(Muon_eta,"+channel+"_trg_mu_idx)")
	if (channel=="BToMuMuPi"):
		df = df.Define(""+channel+"_sel_lep_pt","Take(Muon_pt,"+channel+"_sel_mu_idx)")
		df = df.Define(""+channel+"_sel_lep_eta","Take(Muon_eta,"+channel+"_sel_mu_idx)")
	elif (channel=="BToMuEPi"):
		df = df.Define(""+channel+"_sel_lep_pt","Take(Electron_pt,"+channel+"_sel_ele_idx)")
		df = df.Define(""+channel+"_sel_lep_eta","Take(Electron_eta,"+channel+"_sel_ele_idx)")
	h_num = ROOT.TH2D("num","model;Lxy_{lep#pi}(cm);lep#pi p_{T}(GeV)",8,0,7,5,0,4)
	h_den = ROOT.TH2D("den","model;Lxy_{lep#pi}(cm);lep#pi p_{T}(GeV)",8,0,7,5,0,4)
	for j in range(0,len(pt_bin)-1):
		binValue = Array('f', len(lxy_bin)-1)
		binError = Array('f', len(lxy_bin)-1)
		pool = ROOT.TProcessExecutor(15);
		pxes = (Process(target=FillRecoEffBins,args=(channel,df_gen,df,h_num,pt_bin,lxy_bin,i,j,den_sel,num_sel,binValue,binError) ) for i in range(0,len(lxy_bin)-1))
		for p in pxes:
			p.start()
			p.join() 
		for i in range(0,len(lxy_bin)-1):

			h_num.GetXaxis().SetBinLabel(i+1,"("+str(int(lxy_bin[i]))+","+str(int(lxy_bin[i+1]))+")")
			h_num.GetYaxis().SetBinLabel(j+1,"("+str(int(pt_bin[j]))+","+str(int(pt_bin[j+1]))+")")
			h_num.SetBinContent(i+1,j+1,binValue[i])
			h_num.SetBinError(i+1,j+1,binError[i])
#	#histo model for efficiency numerator & denominator binning
	c = ROOT.TCanvas("c","c",800,600)
#	#extract histo with selections
#	h_num.Divide(h_den)
	#plotting
	h_num.Draw("COLZTEXTE")
	h_num.GetZaxis().SetRangeUser(0,1)
	ROOT.gStyle.SetPaintTextFormat(".2f")
	ROOT.gStyle.SetPalette(57)
	c.SaveAs("../plots/Efficiencies/RecoEff_"+channel+"lxyVsPt.pdf")


if __name__ == '__main__':
	channel = "BToMuEPi"
	selection = "hnl_lxy<1 && LepQProd<0" 
	ROOT.gROOT.LoadMacro("../macros/setStyle.C")
	ROOT.setStyle()

#	RecoEfficiencies("/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass3p0_ctau1000p0/220715_154019/0000/","outputfiles/V00_5filters_genStudies_Central/mass3.0_ctau1000_miniGenTree.root","BToMuEPi")

	samples=["BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass1p0_ctau1000p0//221123_165448/0000/",
		# "BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass2p0_ctau100p0/220715_153807/0000/",
		 "BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass3p0_ctau100p0/221123_152213//0000/",
		 "BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass4p5_ctau1p0/221123_152250/0000/",
		 ""
		]
#	InfoSamples= [  ["1p0_ctau10p0",1,10,0.01,1.043740000000000000e+05],["1p0_ctau100p0",1,100,0.01],["1p0_ctau1000p0",1,1000,0.01],
		 	#["2p0_ctau10p0",2,10,0.02],["2p0_ctau100p0",2,100,0.02],["2p0_ctau1000p0",2,1000,0.02],
		 	#["3p0_ctau1p0",3,1,0.025],["3p0_ctau10p0",3,10,0.025],["3p0_ctau100p0",3,100,0.025],["3p0_ctau1000p0",3,1000,0.025],
		 	#["4p5_ctau0p1",4.5,0.1,0.035],["4p5_ctau1p0",4.5,1,0.035],["4p5_ctau10p0",4.5,10,0.035],["4p5_ctau100p0",4.5,100,0.035]
#	]
	import csv
#	InfoSamples = np.genfromtxt('/cmshome/ratramon/Analysis/data/MC_datasets_InclusiveFilterEff.csv',delimiter = ' ',dtype=["U20","float64","float64","U200","float64","float64","float64"])
#	print(InfoSamples)
	cut = np.zeros(30)
	cutflows = [cut for i in range(0,len(samples)) ]
	presel1= ["trgmu_pt>7","abs(trgmu_eta)<1.5","sel_lep_pt>1.5","abs(sel_lep_eta)<2",\
			channel+"_pi_pt>0.7","abs("+channel+"_pi_eta)<2","abs("+channel+"_pi_dz)>0.005","abs("+channel+"_pi_dxy)>0.005",\
			"abs("+channel+"_pi_dzS)>1.5","abs("+channel+"_pi_dxyS)>3","abs("+channel+"_pi_DCASig)>5",\
			]
	presel_displaced_lep=["_dz>0.0015","_dxy>0.001","_dzS>1","_dxyS>1.5"]
	presel2 = [channel+"_sv_prob>0.001",channel+"_hnl_cos2D>0.995",channel+"_sv_lxy_sig>15",channel+"_mass<8"]
	lables = ["trig muon pt >7 GeV",\
		  "trig muon |#eta| < 1.5",\
		  "electron pt >0.7 GeV",\
		  "electron |\eta|<2",\
		  "pion pt >0.7 GeV",\
		  "pion |#eta|<2",\
		  "pion dz >0.005 cm",\
		  "pion dxy>0.005 cm",\
		  "pion dz significance>1.5",\
		  "pion dxy significance>3",\
		  "pion DCA significance >5",\
		  "displaced lepton dz > 0.0015cm",\
		  "displaced lepton dxy >0.001 cm",\
		  "displaced lepton dz significance> 1 ",\
		  "displaced lepton dxy significance >1.5",\
		  "vertex probability>0.001",\
		  "cos2D>0.995",\
		  "hnl lxy significance >20 ",\
		  "#mul#pi mass <8 GeV",
		  "hnl charge = 0",
		  "muonId",
		  "pfmvaId>-3",
		  "pi purity flag",
		  "LxyS>0",
		  ]


#	br =[ 
 	 	 #	[("hnl_mass","hnl_mass;m_{l\pi} (GeV) ;normalized to unity",100,0,5.5),"hnl_mass"],
 	 	 #	[("B_mass","B_mass;m_{\mu l\pi} (GeV) ;normalized to unity",100,0,6.5),"B_mass"],
 #	 	 	[("dilepton_mass","dilepton_mass;m_{\mu l} (GeV) ;normalized to unity",100,0,5.5),"dilepton_mass"],
 #	 	 	[("Blep_pi_mass","Blep_pi_mass;m_{l_{B}\pi} (GeV) ;normalized to unity",100,0,6.5),"Ble_pi_mass"],
		#]
#	ratioplot(br[0])
#	ratioplot(br[1])
#	scanGraph = [ [ROOT.TGraph(),ROOT.TGraph] for i in range(0, len(InfoSamples))]
	#print(len(scanGraph))
	#print(type(scanGraph[0]))
	#print(len(samples))
	branches= [#[("hnl_mass","hnl_mass;mass_{l#pi}(GeV);entries",50,0,7),"_hnl_mass"],
 		 #	[("hnl_pt","hnl_pt;hnl_pt;entries",100,0,100),"_hnl_pt"],
 		 	[("mass","mass;mass_{l#mu#pi}(GeV);entries",100,0,10),"_mass"],
 		 # 	[("hnl_lxy","hnl_lxy;hnl_lxy_{l#pi}(cm);entries",40,10,50),"_sv_lxy"],
 	 	 # 	[("hnl_lxy_sig","hnl_lxy_sig;hnl_lxy_{sigl#pi}(cm);entries",100,0,100),"_sv_lxy_sig"],
 	 	 #	[("hnl_cos2D","hnl_cos2D;hnl_cos2D;entries",100,0.99,1),"_hnl_cos2D"],
 	         #	[("hnl_ct","hnl_ct;hnl_c#tau;entries",50,0,10),"_hnl_ct"],
 	 	 	[("hnl_probVtx","hnl_probVtx;vertex probability;entries",100,0,1),"_sv_prob"],
 	 	 #	[("hnl_trgmu_pt","hnl_trgmu_pt;trg#mu p_{T}(GeV);entries",50,0,30),"_trgmu_pt"],
 	 	 #	[("hnl_trgmu_eta","hnl_trgmu_eta;trg#mu #eta;entries",50,-3,3),"_trgmu_eta"],
 	 	 #	[("hnl_trgmu_dxy","hnl_trgmu_dxy;trg#mu dxy;entries",50,0,100),"_trgmu_dxy"],
 	 	 #       [("hnl_trgmu_dxyS","hnl_trgmu_dxyS;trg#mu dxyS;entries",30,0,20),"_trgmu_dxyS"],
 	 	 #	[("hnl_trgmu_dz","hnl_trgmu_dz;trg#mu dz;entries",50,0,100),"_trgmu_dz"],
 	 	 #	[("hnl_trgmu_dzS","hnl_trgmu_dzS;trg#mu dzS;entries",50,0,100),"_trgmu_dzS"],
 	 	 #	[("hnl_pi_pt","hnl_pi_pt;#pi p_{T}(GeV);entries",50,0,10),"_pi_pt"],
 	 	 #	[("hnl_pi_eta","hnl_pi_eta;#pi #eta;entries",50,-3,3),"_pi_eta"],
 	 	 #	[("hnl_pi_dxyS","hnl_pi_dxyS;#pi dxyS;entries",30,0,20),"_pi_dxyS"],
 	 	 #	[("hnl_sel_lep_pt","hnl_sel_lep_pt;selected lepton p_{T}(GeV);entries",50,0,30),"_sel_lep_pt"],
 	 	 #	[("hnl_sel_lep_eta","hnl_sel_lep_eta;selected lepton #eta;entries",50,-3,3),"_sel_lep_eta"],
 	 	 #	[("hnl_sel_lep_pfmvaId","hnl_sel_lep_pfmvaId;selected lepton PF MVA ID;entries",50,-12,8),"_sel_lep_pfmvaId"],
 	 	 #	[("hnl_sel_lep_dxy","hnl_sel_lep_dxy;selected lepton dxy;entries",50,0,100),"_sel_lep_dxy"],
 	 	 #	[("hnl_sel_lep_dxyS","hnl_sel_lep_dxyS;selected lepton dxyS;entries",30,0,20),"_sel_lep_dxyS"],
 	 	 #	[("hnl_sel_lep_dz","hnl_sel_lep_dz;selected lepton dz;entries",50,0,100),"_sel_lep_dz"],
 	 	 #	[("hnl_sel_lep_dzS","hnl_sel_lep_dzS;selected lepton dzS;entries",50,0,100),"_sel_lep_dzS"],
			]
 #	pool=Pool(processes=14)
  #      pool.map(plotter,branches)
   #     pool.close()
    #    pool.join()
#	plotter(branches[0])
 	Tpool = Pool(processes = 4)
 	cutflows = Tpool.map(cutflow,samples)
#	cutflows[0] = cutflow(samples[0])
#	print cutflows
#	cutflows[1] = cutflow(samples[1])
#	print cutflows
#	cutflows[2] = cutflow(samples[2])
#	print cutflows
#	cutflows[3] = cutflow(samples[3])
# 	Tpool.close()
# 	Tpool.join()
#	print cutflows
	cutflows = np.asarray(cutflows,dtype=np.float32)
 	for i,lable in enumerate(lables):
 	 	for j in range(0,len(cutflows)):
 	 	 	lable = lable+' & {:.2f}%'.format(cutflows[j][i][0])+' & {:.2f}%'.format(cutflows[j][i][1])
 		print lable
 	print(cutflows)
	
#	Tpool = threadPool(20)
# 	scanGraph = Tpool.map(SignificanceMax,InfoSamples)
# 	Tpool.close()
# 	Tpool.join()
# 	c = ROOT.TCanvas("c","c",1000,800)
#	multi = ROOT.TMultiGraph()
#	multi_point = ROOT.TMultiGraph()
#	multi.SetTitle(";likelihood cut;significance")
#	ROOT.gStyle.SetPalette(80)
#	for g in scanGraph:
#
#		multi.Add(g[0])
#		g[1].SetMarkerColor(g[0].GetLineColor())
#	multi.Draw("apl PLC PMC")
#	multi.GetXaxis().SetRangeUser(0,1.2)
#	for g in scanGraph:
#		print(g[0].GetLineColor())
#		g[1].SetMarkerColor(g[0].GetLineColor())
#		multi_point.Add(g[1])
#	multi_point.GetXaxis().SetRangeUser(0,1.2)
#	multi.SetMinimum(min(multi.GetYaxis().GetXmin(),multi_point.GetYaxis().GetXmin()))
#	multi.SetMaximum(1.5*max(multi.GetYaxis().GetXmax(),multi_point.GetYaxis().GetXmax()))
#	multi_point.Draw("p")
#	multi.GetXaxis().SetLimits(-0.1,1.1)
#	multi_point.GetXaxis().SetLimits(-0.1,1.1)
#       l = c.BuildLegend(0.2,0.7,0.5,0.9,selection)	
#	l.SetBorderSize(0)
#	l.SetFillColor(0)
#	selection = selection.replace(">","Over")
#	selection = selection.replace("<","Under")
#	selection = selection.replace("&&","_")
#	selection = selection.replace(" ","")
#	c.SaveAs("/cmshome/ratramon/Analysis/plots/LSigScans_"+str(InfoSamples[0][1])+"_"+selection+".pdf")
#	cutflow()

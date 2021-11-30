from root_numpy import root2array, tree2array, fill_hist,array2tree
from root_numpy import testdata
import sys
sys.path.append("/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/python/2.7.14-omkpbe4/lib/python2.7/site-packages/")
from rootpy.tree import Tree, TreeModel
from rootpy.tree import IntCol, FloatCol, FloatArrayCol, CharCol, CharArrayCol
from rootpy.io import root_open
from ctypes import *
import ROOT
import numpy as np
import glob
import math
from ROOT import TH1D,TLorentzVector,gROOT
from operator import itemgetter
import multiprocessing
from sklearn.externals.joblib import Parallel, delayed
from tqdm import tqdm


class Event(TreeModel):
    evtIdx = IntCol()
    nCand= IntCol()
    Type = IntCol() # idx for channel category : pi-mu :0 , pi-PF :1, pi-lowpT: 2, track track 3 
    LxyBin = IntCol()	# idx for displacement bin : lxy< 3 : 0, lxy <10 : 1 , lxy < 20:2, lxy >20:3
    ABCDReg = IntCol()
    LepQProd = IntCol()
    B_mass= FloatCol()
    B_pt= FloatCol()
    B_eta= FloatCol()
    B_phi= FloatCol()
    TrgMu_pt= FloatCol()
    TrgMu_eta= FloatCol()
    TrgMu_phi= FloatCol()
    dr_trgMu_lep= FloatCol()
    dr_trgMu_pi= FloatCol()
    hnl_mass= FloatCol()
    hnl_pt= FloatCol()
    hnl_eta= FloatCol()
    hnl_phi= FloatCol()
    hnl_charge= FloatCol()
    hnl_lxy= FloatCol()
    hnl_cos2D= FloatCol()
    hnl_lxy_sig= FloatCol()
    hnl_drLepPi = FloatCol()
    dr_trgmu_hnl = FloatCol()
    hnl_vtxProb = FloatCol()
    hnl_vtxChi2 =  FloatCol()
    hnl_vx =  FloatCol()
    hnl_ex =  FloatCol()
    hnl_vy =  FloatCol()
    hnl_ey =  FloatCol()
    hnl_vz =  FloatCol()
    hnl_ez =  FloatCol()
    hnl_l_pt =  FloatCol()
    hnl_l_eta =  FloatCol()
    hnl_l_phi =  FloatCol()
    hnl_l_dz =  FloatCol()
    hnl_l_dzS =  FloatCol()
    hnl_l_dxy =  FloatCol()
    hnl_l_dxyS =  FloatCol()
    hnl_l_DCAS =  FloatCol()
    hnl_l_mvaId =  FloatCol()
    hnl_l_nHits = FloatCol()
    hnl_pi_pt =  FloatCol()
    hnl_pi_eta =  FloatCol()
    hnl_pi_phi =  FloatCol()
    hnl_pi_dz =  FloatCol()
    hnl_pi_dzS =  FloatCol()
    hnl_pi_dxy =  FloatCol()
    hnl_pi_dxyS =  FloatCol()
    hnl_pi_DCAS =  FloatCol()
    hnl_pi_nHits = FloatCol()
    dphi_pi_fitpi = FloatCol()
    dilepton_mass = FloatCol()
    dilepton_pt = FloatCol()
    likelihood = FloatCol()
#    cosTheta_star = FloatCol()
#def createHistos(full, chString):
#  listHist = []
#  #print len(l[0][0])
#  for j in range(0,len(lxy_lables)):
#	temp = []
#   	for i in range(0,len(in_branches)):
#	
##	temp.resize(4)
##	print i
#		temp.append(ROOT.TH1D(chString+str(i)+str(j),chString+in_branches[i]+lxy_lables[j], Nbin[i], binMin[i], binMax[i]))
##	print len(temp)
#	listHist.append(temp)
#  #print len(l[0])	
##  print len(listHist)
#  return listHist   
def dPhi(a, b):

	delta = a-b
	while (delta >= np.pi):
		 delta -= 2*np.pi
  	while (delta < -np.pi):
		 delta += 2*np.pi
        return delta

def dR(a_phi, b_phi,a_eta, b_eta):

	deltaPhi = dPhi(a_phi,b_phi)
	deltaEta = a_eta-b_eta
	return math.sqrt(pow(deltaPhi,2)+pow(deltaEta,2))

def sortBranches(l,full, chString):
	J = 1.0*j
	if (len(l[0])!=0):
 	   
         test = []
#  	 print "BEFORE SORT"
 #    	 print l
 	 test = l
 #	if j % len(rta_charge[0]) == 1000:
 	#	print('{:.3f}'.format(J/len(rta_charge[0])))
	 if not full:
        	 test = zip(*test) 
      	         test = sorted(test,key= itemgetter(8) ,reverse=True)
                 test = zip(*test)
         l = test

def TwoBodyMass(sig,pt1,eta1,phi1,mass1,pt2,eta2,phi2,mass2):

	if sig<3:
	
		p1 = TLorentzVector()
		p2 = TLorentzVector()
		TLorentzVector.SetPtEtaPhiM(p1,pt1,eta1,phi1,mass1)
		TLorentzVector.SetPtEtaPhiM(p2,pt2,eta2,phi2,mass2)
#		print( (p1+p2).M())
		
		return (p1+p2).M()

	else:

		p1_l = TLorentzVector()
		p2_t = TLorentzVector()
		p1_t = TLorentzVector()
		p2_l = TLorentzVector()
		TLorentzVector.SetPtEtaPhiM(p1_l,pt1,eta1,phi1,mass1)
		TLorentzVector.SetPtEtaPhiM(p2_t,pt2,eta2,phi2,mass2)
		TLorentzVector.SetPtEtaPhiM(p1_t,pt1,eta1,phi1,mass2)
		TLorentzVector.SetPtEtaPhiM(p2_l,pt2,eta2,phi2,mass1)

#		print((p1_l+p2_t).M())
	#	print((p1_t+p2_l).M())
		

		return (p1_l+p2_t).M()

def TwoBodyPt(sig,pt1,eta1,phi1,mass1,pt2,eta2,phi2,mass2):

	if sig<3:
	
		p1 = TLorentzVector()
		p2 = TLorentzVector()
		TLorentzVector.SetPtEtaPhiM(p1,pt1,eta1,phi1,mass1)
		TLorentzVector.SetPtEtaPhiM(p2,pt2,eta2,phi2,mass2)
#		print( (p1+p2).M())
		
		return (p1+p2).Pt()

	else:

		p1_l = TLorentzVector()
		p2_t = TLorentzVector()
		p1_t = TLorentzVector()
		p2_l = TLorentzVector()
		TLorentzVector.SetPtEtaPhiM(p1_l,pt1,eta1,phi1,mass1)
		TLorentzVector.SetPtEtaPhiM(p2_t,pt2,eta2,phi2,mass2)
		TLorentzVector.SetPtEtaPhiM(p1_t,pt1,eta1,phi1,mass2)
		TLorentzVector.SetPtEtaPhiM(p2_l,pt2,eta2,phi2,mass1)

		#print((p1_l+p2_t).Pt())
	#	print((p1_t+p2_l).M())
		

		return (p1_l+p2_t).Pt()







K_MASS = 0.493677;
PI_MASS = 0.139571;
LEP_SIGMA = 0.0000001;
K_SIGMA = 0.000016;
PI_SIGMA = 0.000016;
MUON_MASS = 0.10565837;
ELECTRON_MASS = 0.000511;
mass = [MUON_MASS,ELECTRON_MASS,ELECTRON_MASS,ELECTRON_MASS]
n_branches = 8
binMin = [0,0,0,0,0,0,0] 
binMax = [7,20,100,1,10,50,50] 
Nbin = [100,100,100,100,100,100,100]
sig = ["BToMuMuPi","BToMuPFPi","BToMuLowPtPi","BToMutt"]
lxyMu_bins = ["BToMuMuPi_sv_lxy<3","(BToMuMuPi_sv_lxy>3 && BToMuMuPi_sv_lxy<10)","(BToMuMuPi_sv_lxy>10 && BToMuMuPi_sv_lxy<20)","BToMuMuPi_sv_lxy>20"]
lxyEle_bins = ["BToMuEPi_sv_lxy<3","BToMuEPi_sv_lxy>3 && BToMuEPi_sv_lxy<10","BToMuEPi_sv_lxy>10 && BToMuEPi_sv_lxy<20","BToMuEPi_sv_lxy>20"]
lxytt_bins = ["BToMuEPiHD_sv_lxy<3","BToMuEPiHD_sv_lxy>3 && BToMuEPiHD_sv_lxy<10","BToMuEPiHD_sv_lxy>10 && BToMuEPiHD_sv_lxy<20","BToMuEPiHD_sv_lxy>20","ProbeTracks_dz[BToMuLPiHD_pi_idx]"]
lxy_lables = ["lxyUnder3","lxyUnder10","lxyUnder20","lxyOver20"]
leaf_lables = ["Hnl_mass","Hnl_pt","lxy","cos2D","Bmass","Bpt","lxy_sig"]
mu_branches = ["BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_eta","BToMuMuPi_phi",
               "Muon_pt[BToMuMuPi_trg_mu_idx]","Muon_eta[BToMuMuPi_trg_mu_idx]","Muon_phi[BToMuMuPi_trg_mu_idx]",
               "BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_sv_lxy_sig","BToMuMuPi_dr_mu_pi",
               "BToMuMuPi_dr_trgmu_hnl","BToMuMuPi_sv_prob","BToMuMuPi_sv_chi2","BToMuMuPi_sv_x","BToMuMuPi_sv_xe","BToMuMuPi_sv_y","BToMuMuPi_sv_ye","BToMuMuPi_sv_z","BToMuMuPi_sv_ze",
               "Muon_pt[BToMuMuPi_sel_mu_idx]","Muon_eta[BToMuMuPi_sel_mu_idx]","Muon_phi[BToMuMuPi_sel_mu_idx]","Muon_dz[BToMuMuPi_sel_mu_idx]","Muon_dzS[BToMuMuPi_sel_mu_idx]",
               "Muon_dxy[BToMuMuPi_sel_mu_idx]","Muon_dxyS[BToMuMuPi_sel_mu_idx]","Muon_softId[BToMuMuPi_trg_mu_idx]","Muon_looseId[BToMuMuPi_sel_mu_idx]",
               "ProbeTracks_pt[BToMuMuPi_pi_idx]","ProbeTracks_eta[BToMuMuPi_pi_idx]","ProbeTracks_phi[BToMuMuPi_pi_idx]","ProbeTracks_dz[BToMuMuPi_pi_idx]","ProbeTracks_dzS[BToMuMuPi_pi_idx]",
               "ProbeTracks_dxy[BToMuMuPi_pi_idx]","ProbeTracks_dxyS[BToMuMuPi_pi_idx]","ProbeTracks_DCASig[BToMuMuPi_pi_idx]","BToMuMuPi_fit_pi_phi","BToMuMuPi_hnl_charge","Muon_charge[BToMuMuPi_sel_mu_idx]","Muon_charge[BToMuMuPi_trg_mu_idx]","BToMuMuPi_hnl_eta","BToMuMuPi_hnl_phi",
		"BToMuMuPi_dilepton_mass","BToMuMuPi_dilepton_pt","Muon_mediumId[BToMuMuPi_trg_mu_idx]","ProbeTracks_ndof[BToMuMuPi_pi_idx]"
#, "BToMuMuPi_hnl_cos2D_star"
		]

ele_branches = ["BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_eta","BToMuEPi_phi",
               "Muon_pt[BToMuEPi_trg_mu_idx]","Muon_eta[BToMuEPi_trg_mu_idx]","Muon_phi[BToMuEPi_trg_mu_idx]",
               "BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_sv_lxy_sig","BToMuEPi_dr_ele_pi",
               "BToMuEPi_dr_trgmu_hnl","BToMuEPi_sv_prob","BToMuEPi_sv_chi2","BToMuEPi_sv_x","BToMuEPi_sv_xe","BToMuEPi_sv_y","BToMuEPi_sv_ye","BToMuEPi_sv_z","BToMuEPi_sv_ze",
               "Electron_pt[BToMuEPi_sel_ele_idx]","Electron_eta[BToMuEPi_sel_ele_idx]","Electron_phi[BToMuEPi_sel_ele_idx]","Electron_dz[BToMuEPi_sel_ele_idx]","Electron_dzErr[BToMuEPi_sel_ele_idx]",
               "Electron_dxy[BToMuEPi_sel_ele_idx]","Electron_dxyErr[BToMuEPi_sel_ele_idx]","Electron_pfmvaId[BToMuEPi_sel_ele_idx]","Electron_mvaId[BToMuEPi_sel_ele_idx]",
               "ProbeTracks_pt[BToMuEPi_pi_idx]","ProbeTracks_eta[BToMuEPi_pi_idx]","ProbeTracks_phi[BToMuEPi_pi_idx]","ProbeTracks_dz[BToMuEPi_pi_idx]","ProbeTracks_dzS[BToMuEPi_pi_idx]",
               "ProbeTracks_dxy[BToMuEPi_pi_idx]","ProbeTracks_dxyS[BToMuEPi_pi_idx]","ProbeTracks_DCASig[BToMuEPi_pi_idx]","BToMuEPi_fit_pi_phi","BToMuEPi_hnl_charge","Electron_charge[BToMuEPi_sel_ele_idx]","Muon_charge[BToMuEPi_trg_mu_idx]","BToMuEPi_hnl_eta","BToMuEPi_hnl_phi"
,"BToMuEPi_dilepton_mass","BToMuEPi_dilepton_pt","Muon_softId[BToMuEPi_trg_mu_idx]","ProbeTracks_ndof[BToMuEPi_pi_idx]"
#,"BToMuEPi_hnl_cos2D_star"
	]

tt_branches = ["BToMuEPiHD_mass","BToMuEPiHD_pt","BToMuEPiHD_eta","BToMuEPiHD_phi",
               "Muon_pt[BToMuEPiHD_trg_mu_idx]","Muon_eta[BToMuEPiHD_trg_mu_idx]","Muon_phi[BToMuEPiHD_trg_mu_idx]",
               "BToMuEPiHD_hnl_mass","BToMuEPiHD_hnl_pt","BToMuEPiHD_sv_lxy","BToMuEPiHD_hnl_cos2D","BToMuEPiHD_sv_lxy_sig","BToMuEPiHD_dr_mu_pi",
               "BToMuEPiHD_dr_trgmu_hnl","BToMuEPiHD_sv_prob","BToMuEPiHD_sv_chi2","BToMuEPiHD_sv_x","BToMuEPiHD_sv_xe","BToMuEPiHD_sv_y","BToMuEPiHD_sv_ye","BToMuEPiHD_sv_z","BToMuEPiHD_sv_ze",
               "ProbeTracks_pt[BToMuEPiHD_sel_e_idx]","ProbeTracks_eta[BToMuEPiHD_sel_e_idx]","ProbeTracks_phi[BToMuEPiHD_sel_e_idx]","ProbeTracks_dz[BToMuEPiHD_sel_e_idx]","ProbeTracks_dzS[BToMuEPiHD_sel_e_idx]",
               "ProbeTracks_dxy[BToMuEPiHD_sel_e_idx]","ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx]","ProbeTracks_DCASig[BToMuEPiHD_sel_e_idx]","ProbeTracks_ndof[BToMuEPiHD_sel_e_idx]",
               "ProbeTracks_pt[BToMuEPiHD_pi_idx]","ProbeTracks_eta[BToMuEPiHD_pi_idx]","ProbeTracks_phi[BToMuEPiHD_pi_idx]","ProbeTracks_dz[BToMuEPiHD_pi_idx]","ProbeTracks_dzS[BToMuEPiHD_pi_idx]",
               "ProbeTracks_dxy[BToMuEPiHD_pi_idx]","ProbeTracks_dxyS[BToMuEPiHD_pi_idx]","ProbeTracks_DCASig[BToMuEPiHD_pi_idx]","BToMuEPiHD_fit_pi_phi","BToMuEPiHD_hnl_charge", "ProbeTracks_charge[BToMuEPiHD_sel_e_idx]","Muon_charge[BToMuEPiHD_trg_mu_idx]","BToMuEPiHD_hnl_eta","BToMuEPiHD_hnl_phi"
,"BToMuEPiHD_dilepton_mass","BToMuEPiHD_dilepton_pt","Muon_softId[BToMuEPiHD_trg_mu_idx]","ProbeTracks_ndof[BToMuEPiHD_pi_idx]"
#,"BToMuEPiHD_hnl_cos2D_star"
		]
#tt_branches = ["BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_sv_lxy_sig"]
in_branches = ["BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_sv_lxy_sig"]
print str(sys.argv)
isQCD = int(sys.argv[1])
isMC = 1
datatype = isMC
print "initialize list of input"
if datatype == 1:
	n = len(sys.argv)
	filenames = sys.argv[3:n]
	mu_sig= "BToMuMuPi_isMatched ==1 &&"# &&  BToMuMuPi_hnl_charge ==0 "
	pf_sig= "BToMuEPi_isMatched ==1 &&  BToMuEPi_hnl_charge ==0 && Electron_isPF[BToMuEPi_sel_ele_idx] && "
	lowpT_sig= "BToMuEPi_isMatched &&  BToMuEPi_hnl_charge ==0 && Electron_isLowPt[BToMuEPi_sel_ele_idx] && Electron_isPFoverlap[BToMuEPi_sel_ele_idx]==0 && "
	tt_sig= "BToMuEPiHD_hnl_charge ==0 && "
#	BToMuEPiHD_IsMatcheD 
elif datatype == 0:

	mu_sig= "  " 
	pf_sig="  Electron_isPF[BToMuEPi_sel_ele_idx] && " 
	lowpT_sig="  Electron_isLowPt[BToMuEPi_sel_ele_idx] && " 
	tt_sig= "  "
	n = len(sys.argv)
	filenames = sys.argv[3:n]
#	sel= "(BToMuEPi_mass<2.2 || BToMuEPi_mass>7) && BToMuEPi_hnl_charge ==0"
#	filenames = glob.glob('/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/ParkingBPH1/crab_data_Run2018A_part1/210729_073310/0000*/*121.root')#/bparknano_56.root'
print "list of input completed"
#intree = filenames.Get('Events')

# and convert the TTree into an array
rta_mu =[] 
rta_muSel = []

rta_pf =[] 
rta_pfSel = []
rta_lowPt =[] 
rta_lowPtSel = []
rta_tt =[] 
rta_ttSel = []
trg_BPark_mu = "(Muon_fired_HLT_Mu7_IP4[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP3[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP5[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP6[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8p5_IP3p5[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP4[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP5[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP6[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu10p5_IP3p5[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu12_IP6[BToMuMuPi_trg_mu_idx]==1 ) && "

trg_BPark_ele = "(Muon_fired_HLT_Mu7_IP4[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP3[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP6[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8p5_IP3p5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP4[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP6[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu10p5_IP3p5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu12_IP6[BToMuEPi_trg_mu_idx]==1 ) && "
trg_BPark_trk = "(Muon_fired_HLT_Mu7_IP4[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP3[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP6[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8p5_IP3p5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP4[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP6[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu10p5_IP3p5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu12_IP6[BToMuEPiHD_trg_mu_idx]==1 ) && "
mu_sel = " Muon_pt[BToMuMuPi_trg_mu_idx]>7  &&  Muon_softId[BToMuMuPi_trg_mu_idx]==1  && Muon_looseId[BToMuMuPi_sel_mu_idx]==1  && fabs(Muon_eta[BToMuMuPi_trg_mu_idx])<1.5 &&  Muon_pt[BToMuMuPi_sel_mu_idx]>1.5 &&  fabs(Muon_eta[BToMuMuPi_sel_mu_idx])<2 &&  ProbeTracks_pt[BToMuMuPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuMuPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuMuPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuMuPi_pi_idx])>0.005  && fabs(ProbeTracks_dzS[BToMuMuPi_pi_idx])>1.5 && fabs(ProbeTracks_dxyS[BToMuMuPi_pi_idx])>3 &&  fabs(ProbeTracks_DCASig[BToMuMuPi_pi_idx])>5   &&  BToMuMuPi_sv_prob>0.001  &&  BToMuMuPi_hnl_cos2D>0.99 &&  BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye>20  &&  fabs(Muon_dxy[BToMuMuPi_sel_mu_idx])>0.001 &&  fabs(Muon_dz[BToMuMuPi_sel_mu_idx])>0.0015 &&  fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx])>1.5 &&  fabs(Muon_dzS[BToMuMuPi_sel_mu_idx])>1  &&  fabs(ProbeTracks_eta[BToMuMuPi_pi_idx]-BToMuMuPi_fit_pi_eta )<0.015";
pf_sel = "  Muon_pt[BToMuEPi_trg_mu_idx]>7  && Muon_softId[BToMuEPi_trg_mu_idx]==1 && fabs(Muon_eta[BToMuEPi_trg_mu_idx])<1.5 && Electron_pfmvaId[BToMuEPi_sel_ele_idx]> -3  && ProbeTracks_pt[BToMuEPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuEPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPi_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPi_pi_idx])>5  &&  BToMuEPi_sv_prob>0.001 &&  BToMuEPi_hnl_cos2D>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>3  &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx]/Electron_dxyErr[BToMuEPi_sel_ele_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx]/Electron_dzErr[BToMuEPi_sel_ele_idx])>1 "# &&  deltaR_trgMuPi<2;
#lowpt missing cut on lowptmvaId
lowPt_sel = "   Muon_pt[BToMuEPi_trg_mu_idx]>7&&  Muon_softId[BToMuEPi_trg_mu_idx]==1   &&  fabs(Muon_eta[BToMuEPi_trg_mu_idx])<1.5 && Electron_mvaId[BToMuEPi_sel_ele_idx]> -3  && ProbeTracks_pt[BToMuEPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuEPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPi_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPi_pi_idx])>5  &&  BToMuEPi_sv_prob>0.001 &&  BToMuEPi_hnl_cos2D>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>3  &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx]/Electron_dxyErr[BToMuEPi_sel_ele_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx]/Electron_dzErr[BToMuEPi_sel_ele_idx])>1 "# &&  deltaR_trgMuPi<2;

tt_sel = "  Muon_pt[BToMuEPiHD_trg_mu_idx]>7 &&  Muon_softId[BToMuEPiHD_trg_mu_idx]==1 &&  fabs(Muon_eta[BToMuEPiHD_trg_mu_idx])<1.5 &&  fabs(ProbeTracks_eta[BToMuEPiHD_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPiHD_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPiHD_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPiHD_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPiHD_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPiHD_pi_idx])>5  &&  BToMuEPiHD_sv_prob>0.001 && BToMuEPiHD_hnl_cos2D>0.99  &&  BToMuEPiHD_sv_lxy/BToMuEPiHD_sv_lxye>20  &&  fabs(ProbeTracks_dxy[BToMuEPiHD_sel_e_idx])>0.008 &&  fabs(ProbeTracks_dz[BToMuEPiHD_sel_e_idx])>0.008 &&  fabs(ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx])>1.5 &&  fabs(ProbeTracks_dzS[BToMuEPiHD_sel_e_idx])>1"# &&  deltaR_trgMuPi<2;
inputFiles = tqdm(filenames)
#ROOT.EnableImplicitMT(); # enable multi-threading
#	rta_mu.append(root2array(inputFiles,"Events",branches= mu_branches,object_selection={mu_sig+" && "+lxyMu_bins[j]:['BToMuMuPi_hnl_mass','BToMuMuPi_hnl_pt','BToMuMuPi_sv_lxy','BToMuMuPi_hnl_cos2D','BToMuMuPi_mass','BToMuMuPi_pt','BToMuMuPi_sv_lxy_sig']}))

rta = []
rta.append(root2array(inputFiles,"Events",branches= mu_branches,object_selection={mu_sig+trg_BPark_mu+mu_sel:["BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_eta","BToMuMuPi_phi",
        "Muon_pt[BToMuMuPi_trg_mu_idx]","Muon_eta[BToMuMuPi_trg_mu_idx]","Muon_phi[BToMuMuPi_trg_mu_idx]",   
       "BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_sv_lxy_sig","BToMuMuPi_dr_mu_pi",
       "BToMuMuPi_dr_trgmu_hnl","BToMuMuPi_sv_prob","BToMuMuPi_sv_chi2","BToMuMuPi_sv_x","BToMuMuPi_sv_xe","BToMuMuPi_sv_y","BToMuMuPi_sv_ye","BToMuMuPi_sv_z","BToMuMuPi_sv_ze",
       "Muon_pt[BToMuMuPi_sel_mu_idx]","Muon_eta[BToMuMuPi_sel_mu_idx]","Muon_phi[BToMuMuPi_sel_mu_idx]","Muon_dz[BToMuMuPi_sel_mu_idx]","Muon_dzS[BToMuMuPi_sel_mu_idx]",
       "Muon_dxy[BToMuMuPi_sel_mu_idx]","Muon_dxyS[BToMuMuPi_sel_mu_idx]","Muon_softId[BToMuMuPi_trg_mu_idx]","Muon_looseId[BToMuMuPi_sel_mu_idx]",
       "ProbeTracks_pt[BToMuMuPi_pi_idx]","ProbeTracks_eta[BToMuMuPi_pi_idx]","ProbeTracks_phi[BToMuMuPi_pi_idx]","ProbeTracks_dz[BToMuMuPi_pi_idx]","ProbeTracks_dzS[BToMuMuPi_pi_idx]",
       "ProbeTracks_dxy[BToMuMuPi_pi_idx]","ProbeTracks_dxyS[BToMuMuPi_pi_idx]","ProbeTracks_DCASig[BToMuMuPi_pi_idx]","BToMuMuPi_fit_pi_phi","BToMuMuPi_hnl_charge", "Muon_charge[BToMuMuPi_sel_mu_idx]","BToMuMuPi_hnl_eta","BToMuMuPi_hnl_phi",
	"BToMuMuPi_dilepton_mass","BToMuMuPi_dilepton_pt","ProbeTracks_ndof[BToMuMuPi_pi_idx]"
#,"BToMuMuPi_hnl_cos2D_star"
	]
}))
rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel:["BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_eta","BToMuEPi_phi",
       "BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_sv_lxy_sig","BToMuEPi_dr_ele_pi",
       "BToMuEPi_dr_trgmu_hnl","BToMuEPi_sv_prob","BToMuEPi_sv_chi2","BToMuEPi_sv_x","BToMuEPi_sv_xe","BToMuEPi_sv_y","BToMuEPi_sv_ye","BToMuEPi_sv_z","BToMuEPi_sv_ze",
       "Electron_pt[BToMuEPi_sel_ele_idx]","Electron_eta[BToMuEPi_sel_ele_idx]","Electron_phi[BToMuEPi_sel_ele_idx]","Electron_dz[BToMuEPi_sel_ele_idx]","Electron_dzErr[BToMuEPi_sel_ele_idx]",
       "Electron_dxy[BToMuEPi_sel_ele_idx]","Electron_dxyErr[BToMuEPi_sel_ele_idx]","Electron_pfmvaId[BToMuEPi_sel_ele_idx]","Electron_mvaId[BToMuEPi_sel_ele_idx]",
       "ProbeTracks_pt[BToMuEPi_pi_idx]","ProbeTracks_eta[BToMuEPi_pi_idx]","ProbeTracks_phi[BToMuEPi_pi_idx]","ProbeTracks_dz[BToMuEPi_pi_idx]","ProbeTracks_dzS[BToMuEPi_pi_idx]",
       "ProbeTracks_dxy[BToMuEPi_pi_idx]","ProbeTracks_dxyS[BToMuEPi_pi_idx]","ProbeTracks_DCASig[BToMuEPi_pi_idx]","BToMuEPi_fit_pi_phi","BToMuEPi_hnl_charge","Electron_charge[BToMuEPi_sel_ele_idx]","BToMuEPi_hnl_eta","BToMuEPi_hnl_phi","BToMuEPi_dilepton_mass","BToMuEPi_dilepton_pt","Muon_softId[BToMuEPi_trg_mu_idx]","ProbeTracks_ndof[BToMuEPi_pi_idx]"#,"BToMuEPi_hnl_cos2D_star"
	]
}))
 #	rta_pf.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+" && "+lxyEle_bins[j]:["BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_sv_lxy_sig"]}))
#rta_lowPtSel.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={lowpT_sig+" && "+lxyEle_bins[j]:["BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_sv_lxy_sig"]}))
rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={lowpT_sig+trg_BPark_ele+lowPt_sel:["BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_eta","BToMuEPi_phi",
       "BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_sv_lxy_sig","BToMuEPi_dr_ele_pi",
       "BToMuEPi_dr_trgmu_hnl","BToMuEPi_sv_prob","BToMuEPi_sv_chi2","BToMuEPi_sv_x","BToMuEPi_sv_xe","BToMuEPi_sv_y","BToMuEPi_sv_ye","BToMuEPi_sv_z","BToMuEPi_sv_ze",
       "Electron_pt[BToMuEPi_sel_ele_idx]","Electron_eta[BToMuEPi_sel_ele_idx]","Electron_phi[BToMuEPi_sel_ele_idx]","Electron_dz[BToMuEPi_sel_ele_idx]","Electron_dzErr[BToMuEPi_sel_ele_idx]",
       "Electron_dxy[BToMuEPi_sel_ele_idx]","Electron_dxyErr[BToMuEPi_sel_ele_idx]","Electron_pfmvaId[BToMuEPi_sel_ele_idx]","Electron_mvaId[BToMuEPi_sel_ele_idx]",
       "ProbeTracks_pt[BToMuEPi_pi_idx]","ProbeTracks_eta[BToMuEPi_pi_idx]","ProbeTracks_phi[BToMuEPi_pi_idx]","ProbeTracks_dz[BToMuEPi_pi_idx]","ProbeTracks_dzS[BToMuEPi_pi_idx]",
       "ProbeTracks_dxy[BToMuEPi_pi_idx]","ProbeTracks_dxyS[BToMuEPi_pi_idx]","ProbeTracks_DCASig[BToMuEPi_pi_idx]","BToMuEPi_fit_pi_phi","BToMuEPi_hnl_charge","Electron_charge[BToMuEPi_sel_ele_idx]","BToMuEPi_hnl_eta","BToMuEPi_hnl_phi","BToMuEPi_dilepton_mass","BToMuEPi_dilepton_pt","Muon_softId[BToMuEPi_trg_mu_idx]","ProbeTracks_ndof[BToMuEPi_pi_idx]"#,"BToMuEPi_hnl_cos2D_star"
	]
}))

#rta_tt.append(root2array(inputFiles,"Events",branches= tt_branches,object_selection={tt_sig+" && "+lxytt_bins[j]:["BToMuEPiHD_hnl_mass","BToMuEPiHD_hnl_pt","BToMuEPiHD_sv_lxy","BToMuEPiHD_hnl_cos2D","BToMuEPiHD_mass","BToMuEPiHD_pt","BToMuEPiHD_sv_lxy_sig"]}))
rta.append(root2array(inputFiles,"Events",branches= tt_branches,object_selection={tt_sig+trg_BPark_trk+tt_sel:["BToMuEPiHD_mass","BToMuEPiHD_pt","BToMuEPiHD_eta","BToMuEPiHD_phi",
       "BToMuEPiHD_hnl_mass","BToMuEPiHD_hnl_pt","BToMuEPiHD_sv_lxy","BToMuEPiHD_hnl_cos2D","BToMuEPiHD_sv_lxy_sig","BToMuEPiHD_dr_mu_pi",
       "BToMuEPiHD_dr_trgmu_hnl","BToMuEPiHD_sv_prob","BToMuEPiHD_sv_chi2","BToMuEPiHD_sv_x","BToMuEPiHD_sv_xe","BToMuEPiHD_sv_y","BToMuEPiHD_sv_ye","BToMuEPiHD_sv_z","BToMuEPiHD_sv_ze",
       "ProbeTracks_pt[BToMuEPiHD_sel_e_idx]","ProbeTracks_eta[BToMuEPiHD_sel_e_idx]","ProbeTracks_phi[BToMuEPiHD_sel_e_idx]","ProbeTracks_dz[BToMuEPiHD_sel_e_idx]","ProbeTracks_dzS[BToMuEPiHD_sel_e_idx]",
       "ProbeTracks_dxy[BToMuEPiHD_sel_e_idx]","ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx]","ProbeTracks_DCASig[BToMuEPiHD_sel_e_idx]","ProbeTracks_ndof[BToMuEPiHD_sel_e_idx]",
       "ProbeTracks_pt[BToMuEPiHD_pi_idx]","ProbeTracks_eta[BToMuEPiHD_pi_idx]","ProbeTracks_phi[BToMuEPiHD_pi_idx]","ProbeTracks_dz[BToMuEPiHD_pi_idx]","ProbeTracks_dzS[BToMuEPiHD_pi_idx]",
       "ProbeTracks_dxy[BToMuEPiHD_pi_idx]","ProbeTracks_dxyS[BToMuEPiHD_pi_idx]","ProbeTracks_DCASig[BToMuEPiHD_pi_idx]","BToMuEPiHD_fit_pi_phi","BToMuEPiHD_hnl_charge","ProbeTracks_charge[BToMuEPiHD_sel_e_idx]","BToMuEPiHD_hnl_eta","BToMuEPiHD_hnl_phi","BToMuEPiHD_dilepton_mass","BToMuEPiHD_dilepton_pt","Muon_softId[BToMuEPiHD_trg_mu_idx]","ProbeTracks_ndof[BToMuEPiHD_pi_idx]"#,"BToMuEPiHD_hnl_cos2D_star"
	]
}))
 #	rta_charge.append(root2array(inputFiles,"Events",branches= in_branches,selection= "BToMuEPiHD_hnl_charge ==2 && (BToMuEPiHD_hnl_cos2D>0.99 ||BToMuEPiHD_hnl_cos2D<-0.99) "))

#rta.append(rta_muSel)
#rta.append(rta_pfSel)
#rta.append(rta_lowPtSel)
#rta.append(rta_ttSel)
#num_cores = multiprocessing.cpu_count()

#HistPath = "../hists/"
#iprint "printing array..."
#print rta_mu[0]
print len(rta)
print len(rta[0])
print len(rta[0][0])
print len(rta[0][0][0])
#print len(rta_pfSel[0])
#print len(rta_lowPtSel[0])
#print len(rta_ttSel[0])
#print len(rta_mu[0][0][0])
#print len(rta_mu[0][0][0][0])
#histMu = []
histMuSel = []
#histPF = []
histPFSel = []
#histLowpT = []
histLowpTSel = []
#histTT = []
histTTSel = []
#
##inputs = tqdm(rta_mu)
#
##processed_list = Parallel(n_jobs=num_cores)(delayed(fillHistos)(i) for i in rta_mu)
#
#
#
#histMu = createHistos(1, "BToMuMuPi")
#histMuSel= createHistos(0,"BToMuMuPi_sel")
#histPF = createHistos( 1, "BToMuPFPi")
#histPFSel = createHistos( 0,"BToMuPFPi_sel")
#histLowpT=createHistos(1,"BToMuLowPtPi")
#histLowpTSel= createHistos( 0,"BToMuLowPtPi_sel")
#histTT=createHistos(1,"BToMuTPi")
#histTTSel= createHistos( 0,"BToMuTPi_sel")
#print len(histMu)
#tree = ROOT.TTree("Events", "Events")

#lenght = []
#cat = [lenght, lenght, lenght, lenght, lenght, lenght,lenght]
#cat =  ROOT.std.vector(ROOT.std.vector('float'))()
#cat.resize(len(in_branches)) #branches definition for output tree
#for i in range (0,len(in_branches)):
#	temp = ROOT.std.vector('float')()
#	cat.push_back(temp)
#	tree._v = cat.at(i)
#	tree.Branch(in_branches[i],cat.at(i))
gROOT.LoadMacro("/cmshome/ratramon/Analysis/macros/BinnedLikelihood.C")
rfileSel = root_open(sys.argv[2], 'w')
treeSel = Tree("Events", model=Event)
dphi_cut = [0.03,0.08,0.15,0.2]
for j in range(0,len(rta[0])):
    for SigIdx in range (0,4):


#trees after selection
      if ( len(rta[SigIdx][j][0])!=0):
		
        sortBranches(rta[SigIdx][j],0, sig[SigIdx])
	if (rta[SigIdx][j][9][0] < 1):
		treeSel.LxyBin = 0
	elif (rta[SigIdx][j][9][0] > 1 and rta[SigIdx][j][9][0] < 5):
		treeSel.LxyBin = 1
	elif (rta[SigIdx][j][9][0] > 5):# and rta[SigIdx][j][9][0] < 20):
		treeSel.LxyBin = 2
#	elif (rta[SigIdx][j][9][0] > 20):
#		treeSel.LxyBin = 3
	if ( not isMC):
		if ((rta[SigIdx][j][14][0] > 0.05 and rta[SigIdx][j][10][0]> 0.996)):

	
			if(isQCD):
				treeSel.ABCDReg = 0
			else:
				treeSel.ABCDReg = -1
		elif (rta[SigIdx][j][14][0] > 0.05 and rta[SigIdx][j][10][0]< 0.996):
			treeSel.ABCDReg = 1 
		elif (rta[SigIdx][j][14][0] < 0.05 and rta[SigIdx][j][10][0]< 0.996):
			treeSel.ABCDReg = 2 
		elif (rta[SigIdx][j][14][0] < 0.05 and  rta[SigIdx][j][10][0]> 0.996):
			treeSel.ABCDReg = 3
	else:
	
		treeSel.ABCDReg = -1
#	if(abs(dPhi(rta[SigIdx][j][33][0],rta[SigIdx][j][39][0]) )> dphi_cut[treeSel.LxyBin]):
#		continue
	
#	if ( rta[SigIdx][j][0][0]>5 and rta[SigIdx][j][0][0]<5.5 ): # control region-------UPPER SIDEBAND
#	if (rta[SigIdx][j][40][0]==0 ): # control region q ! = 0 
#		if (SigIdx == 1 and rta[SigIdx][j][29][0] >-2 ):
#	     		continue
#		elif (SigIdx !=1 ):
#		continue
	#print treeSel.LxyBin	
#print("for event %d and displacement %d, signal category %d",k,j,SigIdx)
        treeSel.evtIdx =j 
        treeSel.nCand=1
        treeSel.Type=SigIdx
        treeSel.LepQProd= rta[SigIdx][j][41][0]*rta[SigIdx][j][42][0]                                                       
        treeSel.B_eta= rta[SigIdx][j][2][0]                                                       
        treeSel.B_phi= rta[SigIdx][j][3][0]                                                       
        treeSel.TrgMu_pt= rta[SigIdx][j][4][0]                                                    
        treeSel.TrgMu_eta= rta[SigIdx][j][5][0]                                                   
        treeSel.TrgMu_phi= rta[SigIdx][j][6][0]                                                   
        treeSel.dr_trgMu_lep= dR(rta[SigIdx][j][6][0],rta[SigIdx][j][24][0],rta[SigIdx][j][5][0],rta[SigIdx][j][23][0])                       
        treeSel.dr_trgMu_pi= dR(rta[SigIdx][j][6][0],rta[SigIdx][j][33][0],rta[SigIdx][j][5][0],rta[SigIdx][j][32][0]) 
        treeSel.hnl_mass= rta[SigIdx][j][7][0]                                                    
        treeSel.hnl_pt= rta[SigIdx][j][8][0]      
        treeSel.hnl_eta= rta[SigIdx][j][43][0]      
        treeSel.hnl_phi= rta[SigIdx][j][44][0]      
        treeSel.hnl_charge = rta[SigIdx][j][40][0]                                         
        treeSel.hnl_lxy= rta[SigIdx][j][9][0]                                                     
        treeSel.hnl_cos2D= rta[SigIdx][j][10][0]                                                  
        treeSel.hnl_lxy_sig = rta[SigIdx][j][11][0]                                               
        treeSel.hnl_drLepPi = rta[SigIdx][j][12][0]                                               
        treeSel.dr_trgmu_hnl = rta[SigIdx][j][13][0]                                              
        treeSel.hnl_vtxProb= rta[SigIdx][j][14][0]                                                
        treeSel.hnl_vtxChi2= rta[SigIdx][j][15][0]                                                
        treeSel.hnl_vx =  rta[SigIdx][j][16][0]                                                   
        treeSel.hnl_ex =  rta[SigIdx][j][17][0]                                                   
        treeSel.hnl_vy =  rta[SigIdx][j][18][0]                                                   
        treeSel.hnl_ey =  rta[SigIdx][j][19][0]                                                   
        treeSel.hnl_vz = rta[SigIdx][j][20][0]                                                    
        treeSel.hnl_ez = rta[SigIdx][j][21][0]                                                    
        treeSel.hnl_l_pt = rta[SigIdx][j][22][0]                                                  
        treeSel.hnl_l_eta =rta[SigIdx][j][23][0]                                                  
        treeSel.hnl_l_phi =rta[SigIdx][j][24][0]                                                  
        treeSel.hnl_l_dz = rta[SigIdx][j][25][0]   
        if (SigIdx==1):
        					                                              
          treeSel.hnl_l_dzS =rta[SigIdx][j][25][0]/rta[SigIdx][j][26][0]
        else:                                             
          treeSel.hnl_l_dzS =rta[SigIdx][j][26][0]
        treeSel.hnl_l_dxy =rta[SigIdx][j][27][0]                                                  
        if (SigIdx==1):
        					                                               
          treeSel.hnl_l_dxyS =rta[SigIdx][j][27][0]/rta[SigIdx][j][28][0]
        else:                                            
          treeSel.hnl_l_dxyS =rta[SigIdx][j][28][0]
        
        if (SigIdx <3):
          treeSel.hnl_l_DCAS = -99
        else:  
          treeSel.hnl_l_DCAS =rta[SigIdx][j][29][0] 

        if (SigIdx ==1):
         treeSel.hnl_l_mvaId= rta[SigIdx][j][29][0]    
         treeSel.hnl_l_nHits= -99 
        elif (SigIdx ==2):
         treeSel.hnl_l_mvaId= rta[SigIdx][j][30][0]  
         treeSel.hnl_l_nHits= -99 
	elif (SigIdx ==0): 
         treeSel.hnl_l_mvaId= rta[SigIdx][j][30][0]  
         treeSel.hnl_l_nHits= -99 
        else:
         treeSel.hnl_l_mvaId= -99
         treeSel.hnl_l_nHits= rta[SigIdx][j][30][0] 
        treeSel.hnl_pi_pt =rta[SigIdx][j][31][0]                                                  
        treeSel.hnl_pi_eta= rta[SigIdx][j][32][0]                                                 
        treeSel.hnl_pi_phi = rta[SigIdx][j][33][0]                                                
        treeSel.hnl_pi_dz =rta[SigIdx][j][34][0]                                                  
        treeSel.hnl_pi_dzS = rta[SigIdx][j][35][0]                                                
        treeSel.hnl_pi_dxy= rta[SigIdx][j][36][0]                                                 
        treeSel.hnl_pi_dxyS = rta[SigIdx][j][37][0]                                      
        treeSel.hnl_pi_DCAS = rta[SigIdx][j][38][0]                                      
        treeSel.hnl_pi_nHits = rta[SigIdx][j][48][0]                                      
        treeSel.dphi_pi_fitpi = dPhi(rta[SigIdx][j][33][0],rta[SigIdx][j][39][0])  
   #     treeSel.cosTheta_star = rta[SigIdx][j][48][0]                                      
 
	#print(B_mass)
        if (SigIdx <3):
		
	        treeSel.B_pt= rta[SigIdx][j][1][0]                                                        
	        treeSel.B_mass= rta[SigIdx][j][0][0]                                                      
	        treeSel.dilepton_mass = rta[SigIdx][j][45][0]                                      
	        treeSel.dilepton_pt = rta[SigIdx][j][46][0]     
	else:
		mass_hnl=TwoBodyMass(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.hnl_pi_pt,treeSel.hnl_pi_eta,treeSel.hnl_pi_phi,PI_MASS)
		B_mass = TwoBodyMass(SigIdx,treeSel.hnl_pt,treeSel.hnl_eta,treeSel.hnl_phi,mass_hnl,treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                         
		B_pt = TwoBodyPt(SigIdx,treeSel.hnl_pt,treeSel.hnl_eta,treeSel.hnl_phi,mass_hnl,treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                         
        	treeSel.B_pt= B_pt                                                   
        	treeSel.B_mass= B_mass                                                  
		treeSel.dilepton_mass = TwoBodyMass(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                                      
       		treeSel.dilepton_pt = TwoBodyPt(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)    
		                                 
	point  = np.array( [np.float(treeSel.hnl_cos2D),np.float(max(abs(treeSel.hnl_l_dxy),abs(treeSel.hnl_pi_dxy))),np.float(max(abs(treeSel.hnl_l_dxyS),abs(treeSel.hnl_pi_dxyS))),np.float(treeSel.dr_trgmu_hnl),np.float(treeSel.dilepton_mass),np.float(treeSel.hnl_vtxChi2)])
	likelihood = ROOT.BinnedLikelihood(point)
        treeSel.likelihood =  likelihood if not math.isnan(likelihood ) else -99                           
        treeSel.fill()
	break
      else:
	continue
 
##saveHistos(rta_mu, histMu,1, "BToMuMuPi")
#saveHistos(rta_muSel, histMuSel,0,"BToMuMuPi_sel")
#saveHistos(rta_pf, histPF,1, "BToMuPFPi")
#saveHistos(rta_pfSel, histPFSel,0,"BToMuPFPi_sel")
#saveHistos(rta_lowPt, histLowpT,1,"BToMuLowPtPi")
#saveHistos(rta_lowPtSel, histLowpTSel,0,"BToMuLowPtPi_sel")
#saveHistos(rta_tt, histTT,1,"BToMuTPi")
#saveHistos(rta_ttSel, histTTSel,1,"BToMuTPi")
#ree.Write()
treeSel.write()
rfileSel.close()

#array2root(rta_mu[0],"test.root","recreate")

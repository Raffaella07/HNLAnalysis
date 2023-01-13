from root_numpy import root2array, tree2array, fill_hist,array2tree
from root_numpy import testdata
import sys
import array
sys.path.append("/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/python/2.7.14-omkpbe4/lib/python2.7/site-packages/")
from rootpy.tree import Tree, TreeModel
from rootpy.tree import UIntCol,IntCol, FloatCol, FloatArrayCol, CharCol, CharArrayCol
from rootpy.io import root_open
from ctypes import *
import ROOT
import numpy as np
import glob
import math
from ROOT import TH1D,TLorentzVector,gROOT
from operator import itemgetter
from multiprocessing import Pool,Process
from sklearn.externals.joblib import Parallel, delayed
from tqdm import tqdm
from NN_validationPlots import NNinput 
from sklearn.preprocessing import RobustScaler

#define output tree structure 

class Event(TreeModel):
    event = UIntCol()
    PV_npvs = IntCol()
    PV_npvsGood = IntCol()
    evtIdx = IntCol()
    sigFlag = IntCol()
    nCand= IntCol()
    Type = IntCol() # idx for channel category : pi-mu :0 , pi-PF :1, pi-lowpT: 2, track track 3 
    LxyBin = IntCol()	# idx for displacement bin : lxy< 3 : 0, lxy <10 : 1 , lxy < 20:2, lxy >20:3
    ABCDReg = IntCol()
    QCDweight = FloatCol()
    LepQProd = IntCol()
    B_mass= FloatCol()
    B_pt= FloatCol()
    B_eta= FloatCol()
    B_phi= FloatCol()
    TrgMu_pt= FloatCol()
    TrgMu_eta= FloatCol()
    TrgMu_phi= FloatCol()
  #  B_l_dxy= FloatCol()
  #  B_l_dxyS= FloatCol()
  #  B_l_dz= FloatCol()
    B_l_dzS= FloatCol()
    dr_trgMu_lep= FloatCol()
    dz_trgMu_lep= FloatCol()
    dr_trgMu_pi= FloatCol()
    hnl_mass= FloatCol()
    hnl_pt= FloatCol()
    hnl_eta= FloatCol()
    hnl_phi= FloatCol()
    hnl_charge= FloatCol()
    hnl_ct = FloatCol()
    hnl_lxy = FloatCol()
    hnl_cos2D= FloatCol()
    hnl_lxy_sig= FloatCol()
    hnl_drLepPi = FloatCol()
    hnl_to_trgmu = IntCol()
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
    BlepPi_mass = FloatCol()
    dr_Blep_pi = FloatCol()
    dilepton_pt = FloatCol()
    BlepPi_pt = FloatCol()
    likelihood = FloatCol()
    nn_score = FloatCol()
    nn_score_single = FloatCol()
    bdt = FloatCol()
    toEle = IntCol()
    toMu = IntCol()
    PU_weight = FloatCol()
    Trg_weight = FloatCol()
    muId_weight = FloatCol()
#   genMatch_muTrg = IntCol()	
#   genMatch_l = IntCol()	
#   genMatch_pi = IntCol()
#   pdgId_muTrg = IntCol()
#   pdgId_l = IntCol()
#   pdgId_pi = IntCol()
#   genMatch_cat = IntCol() 	
#   genPdgId_muTrg = IntCol()	
#   genPdgId_l = IntCol()	
#   genPdgId_pi = IntCol()	
#   genMotherPdgId_muTrg = IntCol()	
#   genMotherPdgId_l = IntCol()	
#   genMotherPdgId_pi = IntCol()	
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
	if (len(l[0])!=0):
 	   
         test = []
 	 test = l
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



if __name__ == '__main__':


	n_branches = 8
	binMin = [0,0,0,0,0,0,0] 
	binMax = [7,20,100,1,10,50,50] 
	Nbin = [100,100,100,100,100,100,100]
	sig = ["BToMuMuPi","BToMuPFPi","BToMuLowPtPi","BToMutt"]
	
	
	
	lxyMu_bins = ["BToMuMuPi_sv_lxy<3","(BToMuMuPi_sv_lxy>3 && BToMuMuPi_sv_lxy<10)","(BToMuMuPi_sv_lxy>10 && BToMuMuPi_sv_lxy<20)","BToMuMuPi_sv_lxy>20"]
	lxyEle_bins = ["BToMuEPi_sv_lxy<3","BToMuEPi_sv_lxy>3 && BToMuEPi_sv_lxy<10","BToMuEPi_sv_lxy>10 && BToMuEPi_sv_lxy<20","BToMuEPi_sv_lxy>20"]
	lxytt_bins = ["BToMuEPiHD_sv_lxy<3","BToMuEPiHD_sv_lxy>3 && BToMuEPiHD_sv_lxy<10","BToMuEPiHD_sv_lxy>10 && BToMuEPiHD_sv_lxy<20","BToMuEPiHD_sv_lxy>20","ProbeTracks_dz[BToMuLPiHD_pi_idx]"]
	lxy_lables = ["lxyUnder3","lxyUnder10","lxyUnder20","lxyOver20"]
	
	InfoSamples = np.genfromtxt('/cmshome/ratramon/Analysis/data/QCD_datasets_InclusiveFilterEff.csv',delimiter = ' ',dtype=["U200","float64","float64","float64"])
	leaf_lables = ["Hnl_mass","Hnl_pt","lxy","cos2D","Bmass","Bpt","lxy_sig"]
	
	#branches to keep for each channel 
	
	scalar_branches = ["event","PV_npvs","PV_npvsGood"]	
	mu_branches = ["BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_eta","BToMuMuPi_phi",
	               "Muon_pt[BToMuMuPi_trg_mu_idx]","Muon_eta[BToMuMuPi_trg_mu_idx]","Muon_phi[BToMuMuPi_trg_mu_idx]",
	               "BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_sv_lxy_sig","BToMuMuPi_dr_mu_pi",
	                "BToMuMuPi_dimu_vzdiff","BToMuMuPi_sv_prob","BToMuMuPi_sv_chi2","BToMuMuPi_sv_x","BToMuMuPi_sv_xe","BToMuMuPi_sv_y","BToMuMuPi_sv_ye","BToMuMuPi_sv_z","BToMuMuPi_sv_ze",
	                "Muon_pt[BToMuMuPi_sel_mu_idx]","Muon_eta[BToMuMuPi_sel_mu_idx]","Muon_phi[BToMuMuPi_sel_mu_idx]","Muon_dz[BToMuMuPi_sel_mu_idx]","Muon_dzS[BToMuMuPi_sel_mu_idx]",
	                "Muon_dxy[BToMuMuPi_sel_mu_idx]","Muon_dxyS[BToMuMuPi_sel_mu_idx]","Muon_softId[BToMuMuPi_trg_mu_idx]","Muon_looseId[BToMuMuPi_sel_mu_idx]",
	                "BToMuMuPi_pi_pt","BToMuMuPi_pi_eta","BToMuMuPi_pi_phi","BToMuMuPi_pi_dz","BToMuMuPi_pi_dzS",
	                "BToMuMuPi_pi_dxy","BToMuMuPi_pi_dxyS","BToMuMuPi_pi_DCASig","BToMuMuPi_fit_pi_phi","BToMuMuPi_hnl_charge","Muon_charge[BToMuMuPi_sel_mu_idx]","Muon_charge[BToMuMuPi_trg_mu_idx]","BToMuMuPi_hnl_eta","BToMuMuPi_hnl_phi",
	       		"BToMuMuPi_trgmu_mu_mass","BToMuMuPi_trgmu_mu_pt","Muon_mediumId[BToMuMuPi_trg_mu_idx]","BToMuMuPi_pi_ndof","BToMuMuPi_fit_mu_eta","BToMuMuPi_fit_mu_phi","BToMuMuPi_fit_mu_mass","Muon_isPF[BToMuMuPi_trg_mu_idx]","Muon_vx[BToMuMuPi_sel_mu_idx]","Muon_vy[BToMuMuPi_sel_mu_idx]","BToMuMuPi_hnl_ct",#,"GenPart_pdgId[Muon_genPartIdx[BToMuMuPi_trg_mu_idx]]","GenPart_pdgId[Muon_genPartIdx[BToMuMuPi_sel_mu_idx]]","GenPart_pdgId[ProbeTracks_genPartIdx[BToMuMuPi_pi_idx]]","GenPart_pdgId[GenPart_genPartIdxMother[Muon_genPartIdx[BToMuMuPi_trg_mu_idx]]]","GenPart_pdgId[GenPart_genPartIdxMother[Muon_genPartIdx[BToMuMuPi_sel_mu_idx]]]","GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuMuPi_pi_idx]]]"#, "BToMuMuPi_hnl_cos2D_star"
	               "Muon_dz[BToMuMuPi_trg_mu_idx]","Muon_dxy[BToMuMuPi_trg_mu_idx]","Muon_dzS[BToMuMuPi_trg_mu_idx]","Muon_dxyS[BToMuMuPi_trg_mu_idx]","Muon_looseId[BToMuMuPi_trg_mu_idx]","Muon_softId[BToMuMuPi_sel_mu_idx]"
			]
	
	
	
	ele_branches = ["BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_eta","BToMuEPi_phi", 
	               "Muon_pt[BToMuEPi_trg_mu_idx]","Muon_eta[BToMuEPi_trg_mu_idx]","Muon_phi[BToMuEPi_trg_mu_idx]",
	               "BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_sv_lxy_sig"," BToMuEPi_dilep_vydiff",
	               "BToMuEPi_dilep_vzdiff","BToMuEPi_sv_prob","BToMuEPi_sv_chi2","BToMuEPi_sv_x","BToMuEPi_sv_xe","BToMuEPi_sv_y","BToMuEPi_sv_ye","BToMuEPi_sv_z","BToMuEPi_sv_ze",
	               "Electron_pt[BToMuEPi_sel_ele_idx]","Electron_eta[BToMuEPi_sel_ele_idx]","Electron_phi[BToMuEPi_sel_ele_idx]","Electron_dz[BToMuEPi_sel_ele_idx]","Electron_dzErr[BToMuEPi_sel_ele_idx]",
	               "Electron_dxy[BToMuEPi_sel_ele_idx]","Electron_dxyErr[BToMuEPi_sel_ele_idx]","Electron_pfmvaId[BToMuEPi_sel_ele_idx]","Electron_mvaId[BToMuEPi_sel_ele_idx]",
	               "BToMuEPi_fit_pi_pt","BToMuEPi_fit_pi_eta","BToMuEPi_fit_pi_mass","BToMuEPi_pi_dz","BToMuEPi_pi_dzS",
	               "BToMuEPi_pi_dxy","BToMuEPi_pi_dxyS","BToMuEPi_pi_DCASig","BToMuEPi_fit_pi_phi","BToMuEPi_hnl_charge","Electron_charge[BToMuEPi_sel_ele_idx]","Muon_charge[BToMuEPi_trg_mu_idx]","BToMuEPi_hnl_eta","BToMuEPi_hnl_phi"
	,"BToMuEPi_dilepton_mass","BToMuEPi_dilepton_pt","Muon_softId[BToMuEPi_trg_mu_idx]","BToMuEPi_pi_ndof","BToMuEPi_fit_lep_mass","BToMuEPi_fit_lep_eta","BToMuEPi_fit_lep_phi","Electron_dzTrg[BToMuEPi_sel_ele_idx]","Electron_vx[BToMuEPi_sel_ele_idx]","Electron_vy[BToMuEPi_sel_ele_idx]","Muon_dz[BToMuEPi_trg_mu_idx]","Muon_dxy[BToMuEPi_trg_mu_idx]","Muon_dzS[BToMuEPi_trg_mu_idx]","Muon_dxyS[BToMuEPi_trg_mu_idx]","Muon_looseId[BToMuEPi_trg_mu_idx]","Muon_isTriggeringBPark[BToMuEPi_trg_mu_idx]","BToMuEPi_hnl_to_trgmu","BToMuEPi_hnl_ct"#,"BToMuEPi_matching_trg_mu_genIdx","BToMuEPi_matching_sel_lep_genIdx","BToMuEPi_matching_pi_genIdx","BToMuEPi_matching_trg_mu_motherPdgId","BToMuEPi_matching_sel_ele_motherPdgId","BToMuEPi_matching_pi_motherPdgId"
	#,"BToMuEPi_hnl_cos2D_star"
	        ]
	
	print(len(ele_branches))	
	
	tt_branches = ["BToMuEPiHD_mass","BToMuEPiHD_pt","BToMuEPiHD_eta","BToMuEPiHD_phi",
	               "Muon_pt[BToMuEPiHD_trg_mu_idx]","Muon_eta[BToMuEPiHD_trg_mu_idx]","Muon_phi[BToMuEPiHD_trg_mu_idx]",
	               "BToMuEPiHD_hnl_mass","BToMuEPiHD_hnl_pt","BToMuEPiHD_sv_lxy","BToMuEPiHD_hnl_cos2D","BToMuEPiHD_sv_lxy_sig","BToMuEPiHD_dr_lep_pi",
	               "BToMuEPiHD_dr_trgmu_hnl","BToMuEPiHD_sv_prob","BToMuEPiHD_sv_chi2","BToMuEPiHD_sv_x","BToMuEPiHD_sv_xe","BToMuEPiHD_sv_y","BToMuEPiHD_sv_ye","BToMuEPiHD_sv_z","BToMuEPiHD_sv_ze",
	               "ProbeTracks_pt[BToMuEPiHD_sel_e_idx]","ProbeTracks_eta[BToMuEPiHD_sel_e_idx]","ProbeTracks_phi[BToMuEPiHD_sel_e_idx]","ProbeTracks_dz[BToMuEPiHD_sel_e_idx]","ProbeTracks_dzS[BToMuEPiHD_sel_e_idx]",
	               "ProbeTracks_dxy[BToMuEPiHD_sel_e_idx]","ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx]","ProbeTracks_DCASig[BToMuEPiHD_sel_e_idx]","ProbeTracks_ndof[BToMuEPiHD_sel_e_idx]",
	               "ProbeTracks_pt[BToMuEPiHD_pi_idx]","ProbeTracks_eta[BToMuEPiHD_pi_idx]","ProbeTracks_phi[BToMuEPiHD_pi_idx]","ProbeTracks_dz[BToMuEPiHD_pi_idx]","ProbeTracks_dzS[BToMuEPiHD_pi_idx]",
	               "ProbeTracks_dxy[BToMuEPiHD_pi_idx]","ProbeTracks_dxyS[BToMuEPiHD_pi_idx]","ProbeTracks_DCASig[BToMuEPiHD_pi_idx]","BToMuEPiHD_fit_pi_phi","BToMuEPiHD_hnl_charge", "ProbeTracks_charge[BToMuEPiHD_sel_e_idx]","Muon_charge[BToMuEPiHD_trg_mu_idx]","BToMuEPiHD_hnl_eta","BToMuEPiHD_hnl_phi"
	,"BToMuEPiHD_dilepton_mass","BToMuEPiHD_dilepton_pt","Muon_softId[BToMuEPiHD_trg_mu_idx]","ProbeTracks_ndof[BToMuEPiHD_pi_idx]","BToMuEPiHD_fit_pi_eta","BToMuEPiHD_fit_pi_pt","ProbeTracks_dzTrg[BToMuEPiHD_sel_e_idx]","BToMuEPiHD_Tomu","BToMuEPiHD_Toele"
	#,"BToMuEPiHD_hnl_cos2D_star"
			]
	#tt_branches = ["BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_sv_lxy_sig"]
	in_branches = ["BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_sv_lxy_sig"]
	print str(sys.argv)
	
	
	#run options
	
	
	print "initialize list of input"
	isQCD = int(sys.argv[1]) #to run on QCD (does not blind signal region)
	isMC = int(sys.argv[5])#MC or not
	do_orthogonal = False
	print(isMC)
	datatype = isMC
	n = len(sys.argv)
	filenames = sys.argv[6:n]
	
	if isQCD:
		for s in InfoSamples:
			if s[0] in sys.argv[6]:
				QCD_weight = s[1]*s[2]/s[3]
				print(QCD_weight)
	else:
		QCD_weight=1
	#selections to keep only matched events in MC
	
	if datatype == 1:
		mu_sig= "BToMuMuPi_isMatched==1 && "# &&  BToMuMuPi_hnl_charge ==0 "
		pf_sig= "BToMuEPi_isMatched ==1 &&  BToMuEPi_hnl_charge ==0 && Electron_isPF[BToMuEPi_sel_ele_idx] &&"
		lowpT_sig= "BToMuEPi_isMatched==1 &&  BToMuEPi_hnl_charge ==0 && Electron_isLowPt[BToMuEPi_sel_ele_idx] && "
		tt_sig= "BToMuEPiHD_hnl_charge ==0 && BToMuEPiHD_IsMatcheD && "
	#	BToMuEPiHD_IsMatcheD 
	
	#simple categorization for non mc events
	
	elif datatype == 0:
	
		mu_sig= "  " 
		pf_sig=" " 
		lowpT_sig="  Electron_isLowPt[BToMuEPi_sel_ele_idx]" 
		tt_sig= "  "
	#	sel= "(BToMuEPi_mass<2.2 || BToMuEPi_mass>7) && BToMuEPi_hnl_charge ==0"
	#	filenames = glob.glob('/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/ParkingBPH1/crab_data_Run2018A_part1/210729_073310/0000*/*121.root')#/bparknano_56.root'
	print "list of input completed"
	#intree = filenames.Get('Events')
	
	#selection definition
	
	#trigger 
	trg_BPark_mu = "(Muon_fired_HLT_Mu7_IP4[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP3[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP5[BToMuMuPi_trg_mu_idx]==2 || Muon_fired_HLT_Mu9_IP6[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8p5_IP3p5[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP4[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP5[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP6[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu10p5_IP3p5[BToMuMuPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu12_IP6[BToMuMuPi_trg_mu_idx]==1 ) && "
	
	trg_BPark_ele =  " Muon_isTriggeringBPark[BToMuEPi_trg_mu_idx]==1 && "
# || Muon_fired_HLT_Mu8_IP3[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP6[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8p5_IP3p5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP4[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP6[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu10p5_IP3p5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu12_IP6[BToMuEPi_trg_mu_idx]==1 ) "
	
	trg_BPark_trk = "(Muon_fired_HLT_Mu7_IP4[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP3[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP6[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu8p5_IP3p5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP4[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP6[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu10p5_IP3p5[BToMuEPiHD_trg_mu_idx]==1 || Muon_fired_HLT_Mu12_IP6[BToMuEPiHD_trg_mu_idx]==1 ) && "
	#mu_sel = " Muon_pt[BToMuMuPi_trg_mu_idx]>7  &&  Muon_softId[BToMuMuPi_trg_mu_idx]==1  && Muon_looseId[BToMuMuPi_sel_mu_idx]==1  && fabs(Muon_eta[BToMuMuPi_trg_mu_idx])<1.5 &&  Muon_pt[BToMuMuPi_sel_mu_idx]>7 &&  fabs(Muon_eta[BToMuMuPi_sel_mu_idx])<2 &&  ProbeTracks_pt[BToMuMuPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuMuPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuMuPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuMuPi_pi_idx])>0.005  && fabs(ProbeTracks_dzS[BToMuMuPi_pi_idx])>1.5 && fabs(ProbeTracks_dxyS[BToMuMuPi_pi_idx])>3 &&  fabs(ProbeTracks_DCASig[BToMuMuPi_pi_idx])>5   &&  BToMuMuPi_sv_prob>0.001  &&  BToMuMuPi_hnl_cos2D>0.99 &&  BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye>20  &&  fabs(Muon_dxy[BToMuMuPi_sel_mu_idx])>0.001 &&  fabs(Muon_dz[BToMuMuPi_sel_mu_idx])>0.0015 &&  fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx])>1.5 &&  fabs(Muon_dzS[BToMuMuPi_sel_mu_idx])>1";
	
	
	#preselection
	
 	mu_sel = "Muon_pt[BToMuMuPi_trg_mu_idx]>7\
		&& fabs(Muon_eta[BToMuMuPi_trg_mu_idx])<1.5\
 		&&( ( BToMuMuPi_hnl_to_trgmu==0 && Muon_looseId[BToMuMuPi_sel_mu_idx]==1\
					       && Muon_softId[BToMuMuPi_trg_mu_idx]==1\
					       && fabs(Muon_dz[BToMuMuPi_sel_mu_idx])>0.0015\
					       && fabs(Muon_dxy[BToMuMuPi_sel_mu_idx])>0.001\
					       && fabs(Muon_dzS[BToMuMuPi_sel_mu_idx])>1\
					       && fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx])>1.5\
						)\
	        ||(  BToMuMuPi_hnl_to_trgmu==1 && Muon_looseId[BToMuMuPi_trg_mu_idx]==1\
					       && Muon_softId[BToMuMuPi_sel_mu_idx]==1\
					       && fabs(Muon_dz[BToMuMuPi_trg_mu_idx])>0.0015\
					       && fabs(Muon_dxy[BToMuMuPi_trg_mu_idx])>0.001\
					       && fabs(Muon_dzS[BToMuMuPi_trg_mu_idx])>1\
					       && fabs(Muon_dxyS[BToMuMuPi_trg_mu_idx])>1.5\
						))\
 		&& Muon_pt[BToMuMuPi_sel_mu_idx]>1.5\
 		&& fabs(Muon_eta[BToMuMuPi_sel_mu_idx])<2 \
 		&& BToMuMuPi_pi_pt>0.7 \
 		&& fabs(BToMuMuPi_pi_eta)<2 \
 		&& BToMuMuPi_sv_prob>0.0001 \
 		&& BToMuMuPi_hnl_cos2D>0.99\
 		&& fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>20 \
 		&& fabs(BToMuMuPi_pi_dz)>0.005 \
 		&& fabs(BToMuMuPi_pi_dxy)>0.005\
 		&& fabs(BToMuMuPi_pi_dzS)>1.5 \
 		&& fabs(BToMuMuPi_pi_dxyS)>3 \
 		&& fabs(BToMuMuPi_pi_DCASig)>5" 
	mu_sel= mu_sel.replace("	","")
 		#&& Electron_pfmvaId[BToMuEPi_sel_ele_idx]>-3\
 	pf_sel = "Muon_pt[BToMuEPi_trg_mu_idx]>7\
		&& fabs(Muon_eta[BToMuEPi_trg_mu_idx])<1.5\
 		&& (( BToMuEPi_hnl_to_trgmu==0  && Muon_softId[BToMuEPi_trg_mu_idx]==1\
					       && fabs(Electron_dz[BToMuEPi_sel_ele_idx])>0.0015\
					       && fabs(Electron_dxy[BToMuEPi_sel_ele_idx])>0.001\
					       && fabs(Electron_dz[BToMuEPi_sel_ele_idx]/Electron_dzErr[BToMuEPi_sel_ele_idx])>1\
					       && fabs(Electron_dxy[BToMuEPi_sel_ele_idx]/Electron_dxyErr[BToMuEPi_sel_ele_idx])>1.5\
						)\
	        ||(  BToMuEPi_hnl_to_trgmu==1 && Muon_looseId[BToMuEPi_trg_mu_idx]==1\
					       && fabs(Muon_dz[BToMuEPi_trg_mu_idx])>0.0015\
					       && fabs(Muon_dxy[BToMuEPi_trg_mu_idx])>0.001\
					       && fabs(Muon_dzS[BToMuEPi_trg_mu_idx])>1\
					       && fabs(Muon_dxyS[BToMuEPi_trg_mu_idx])>1.5\
						))\
 		&& Electron_pt[BToMuEPi_sel_ele_idx]>1.5\
 		&& fabs(Electron_eta[BToMuEPi_sel_ele_idx])<2 \
 		&& BToMuEPi_pi_pt>0.7 \
 		&& fabs(BToMuEPi_pi_eta)<2 \
 		&& BToMuEPi_sv_prob>0.0001 \
 		&& BToMuEPi_hnl_cos2D>0.99\
 		&& fabs(BToMuEPi_sv_lxy/BToMuEPi_sv_lxye)>20 \
 		&& fabs(BToMuEPi_pi_dz)>0.005 \
 		&& fabs(BToMuEPi_pi_dxy)>0.005\
 		&& fabs(BToMuEPi_pi_dzS)>1.5 \
 		&& fabs(BToMuEPi_pi_dxyS)>3 \
 		&& fabs(BToMuEPi_pi_DCASig)>5" 
	pf_sel= pf_sel.replace("	","")
#	mu_sel = " Muon_pt[BToMuMuPi_trg_mu_idx]>7   && fabs(Muon_eta[BToMuMuPi_trg_mu_idx])<1.5 && Muon_looseId[BToMuMuPi_sel_mu_idx]==1 && Muon_pt[BToMuMuPi_sel_mu_idx]>1.5 &&  fabs(Muon_eta[BToMuMuPi_sel_mu_idx])<2 &&  BToMuMuPi_pi_pt>0.7  &&  fabs(BToMuMuPi_pi_eta)<2 && BToMuMuPi_sv_prob>0.001  &&  BToMuMuPi_hnl_cos2D>0.99 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>20 &&  fabs(BToMuMuPi_pi_dz)>0.005  &&  fabs(BToMuMuPi_pi_dxy)>0.005  && fabs(BToMuMuPi_pi_dzS)>1.5 && fabs(BToMuMuPi_pi_dxyS)>3 &&  fabs(BToMuMuPi_pi_DCASig)>5  "
	#mu_perCat_sel= " && BToMuMuPi_hnl_charge==0 && ((fabs(BToMuMuPi_sv_lxy)<1 && ProbeTracks_pt[BToMuMuPi_pi_idx]>1.1 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>30 && fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx])>5 &&  fabs(ProbeTracks_dxyS[BToMuMuPi_pi_idx])>10  ) ||  (fabs(BToMuMuPi_sv_lxy)>1  && fabs(BToMuMuPi_sv_lxy)<5 && ProbeTracks_pt[BToMuMuPi_pi_idx]>1.2 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>100 && fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx])>12 &&  fabs(ProbeTracks_dxyS[BToMuMuPi_pi_idx])>25  ) || (fabs(BToMuMuPi_sv_lxy)>5 && ProbeTracks_pt[BToMuMuPi_pi_idx]>1.3 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>100 && fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx])>15 &&  fabs(ProbeTracks_dxyS[BToMuMuPi_pi_idx])>20  ))"
	

	#mu_sel = " Muon_pt[BToMuMuPi_trg_mu_idx]>7   && fabs(Muon_eta[BToMuMuPi_trg_mu_idx])<1.5 && Muon_looseId[BToMuMuPi_sel_mu_idx]==1 && Muon_pt[BToMuMuPi_sel_mu_idx]>1.5 &&  fabs(Muon_eta[BToMuMuPi_sel_mu_idx])<2 &&  BToMuMuPi_pi_pt>0.7  &&  fabs(BToMuMuPi_pi_eta)<2 && BToMuMuPi_sv_prob>0.001  &&  BToMuMuPi_hnl_cos2D>0.99 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>20 &&  fabs(BToMuMuPi_pi_dz)>0.005  &&  fabs(BToMuMuPi_pi_dxy)>0.005  && fabs(BToMuMuPi_pi_dzS)>1.5 && fabs(BToMuMuPi_pi_dxyS)>3 &&  fabs(BToMuMuPi_pi_DCASig)>5  "
	mu_lowMass_OS = "( BToMuMuPi_hnl_mass<2 && Muon_charge[BToMuMuPi_trg_mu_idx]*Muon_charge[BToMuMuPi_sel_mu_idx]<0\
							 && ( (fabs(BToMuMuPi_sv_lxy)<1 &&  BToMuMuPi_pi_pt>0.7 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>30 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>0 &&  fabs(1-BToMuMuPi_hnl_cos2D)<2e-03  )\
							 ||  (fabs(BToMuMuPi_sv_lxy)>1  && fabs(BToMuMuPi_sv_lxy)<5 && BToMuMuPi_pi_pt>0.7 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>150 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>10 &&  fabs(1-BToMuMuPi_hnl_cos2D)<1e-04  )\
							 || (fabs(BToMuMuPi_sv_lxy)>5 && BToMuMuPi_pi_pt>2 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>150 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>40 &&  fabs(1-BToMuMuPi_hnl_cos2D)<3e-05  )\
							))" 

	mu_lowMass_SS = " (BToMuMuPi_hnl_mass<2 && Muon_charge[BToMuMuPi_trg_mu_idx]*Muon_charge[BToMuMuPi_sel_mu_idx]>0\
							 && ( (fabs(BToMuMuPi_sv_lxy)<1 && BToMuMuPi_pi_pt>0.7 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>50 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>0 &&  fabs(1-BToMuMuPi_hnl_cos2D)<2e-03  )\
							 ||  (fabs(BToMuMuPi_sv_lxy)>1  && fabs(BToMuMuPi_sv_lxy)<5 && BToMuMuPi_pi_pt>0.7 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>150 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>10 &&  fabs(1-BToMuMuPi_hnl_cos2D)<3e-04  )\
							 || (fabs(BToMuMuPi_sv_lxy)>5 && BToMuMuPi_pi_pt>2 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>150 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>40 &&  fabs(1-BToMuMuPi_hnl_cos2D)<3e-05  )\
							))" 
	mu_mediumMass_OS = " (BToMuMuPi_hnl_mass>2 && BToMuMuPi_hnl_mass<4.5 && Muon_charge[BToMuMuPi_trg_mu_idx]*Muon_charge[BToMuMuPi_sel_mu_idx]<0\
						  	 && ( (fabs(BToMuMuPi_sv_lxy)<1 && BToMuMuPi_pi_pt>0.8 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>90 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>5 &&  fabs(1-BToMuMuPi_hnl_cos2D)<2e-03  )\
							 ||  (fabs(BToMuMuPi_sv_lxy)>1  && fabs(BToMuMuPi_sv_lxy)<5 && BToMuMuPi_pi_pt>1.5 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>300 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>40 &&  fabs(1-BToMuMuPi_hnl_cos2D)<2e-04  )\
							 || (fabs(BToMuMuPi_sv_lxy)>5 && BToMuMuPi_pi_pt>3.5 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>300 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>40 &&  fabs(1-BToMuMuPi_hnl_cos2D)<2e-05  )\
							))" 

	mu_mediumMass_SS = "  (BToMuMuPi_hnl_mass<2 && BToMuMuPi_hnl_mass<4.5 && Muon_charge[BToMuMuPi_trg_mu_idx]*Muon_charge[BToMuMuPi_sel_mu_idx]>0\
							 &&( (fabs(BToMuMuPi_sv_lxy)<1 && BToMuMuPi_pi_pt>0.8 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>100 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>5 &&  fabs(1-BToMuMuPi_hnl_cos2D)<2e-03  )\
							 ||  (fabs(BToMuMuPi_sv_lxy)>1  && fabs(BToMuMuPi_sv_lxy)<5 && BToMuMuPi_pi_pt>1.5 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>300 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>60 &&  fabs(1-BToMuMuPi_hnl_cos2D)<2e-04  )\
							 || (fabs(BToMuMuPi_sv_lxy)>5 && BToMuMuPi_pi_pt>3.5 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>300 && max(fabs(BToMuMuPi_pi_dxyS),fabs(Muon_dxyS[BToMuMuPi_sel_mu_idx]))>40 &&  fabs(1-BToMuMuPi_hnl_cos2D)<3e-05  )\
							))" 
	mu_highMass_OS = " (BToMuMuPi_hnl_mass>4.5 && Muon_charge[BToMuMuPi_trg_mu_idx]*Muon_charge[BToMuMuPi_sel_mu_idx]<0\
							 &&( (fabs(BToMuMuPi_sv_lxy)<1 && BToMuMuPi_pi_pt>2 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>100 )\
							 ||  (fabs(BToMuMuPi_sv_lxy)>1  && fabs(BToMuMuPi_sv_lxy)<5 && BToMuMuPi_pi_pt>4 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>100  )\
							 || (fabs(BToMuMuPi_sv_lxy)>5 && BToMuMuPi_pi_pt>5 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>200  )\
							))" 

	mu_highMass_SS = " (BToMuMuPi_hnl_mass>4.5 && Muon_charge[BToMuMuPi_trg_mu_idx]*Muon_charge[BToMuMuPi_sel_mu_idx]>0\
							 && ( (fabs(BToMuMuPi_sv_lxy)<1 &&  BToMuMuPi_pi_pt>2 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>70 )\
							 ||  (fabs(BToMuMuPi_sv_lxy)>1  && fabs(BToMuMuPi_sv_lxy)<5 && BToMuMuPi_pi_pt>4 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>100  )\
							 || (fabs(BToMuMuPi_sv_lxy)>5 && BToMuMuPi_pi_pt>5 &&  fabs(BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye)>200  )\
							))" 
	mu_perCat_sel= "&& BToMuMuPi_hnl_charge==0 && ("+mu_lowMass_OS +"||"+mu_lowMass_SS +"||"+mu_mediumMass_OS +"||"+mu_mediumMass_SS +"||"+mu_highMass_OS +"||"+mu_highMass_SS +")" 
	mu_perCat_sel= mu_perCat_sel.replace("	","")
#	pf_sel = "  Muon_pt[BToMuEPi_trg_mu_idx]>7  && Muon_softId[BToMuEPi_trg_mu_idx]==1 && fabs(Muon_eta[BToMuEPi_trg_mu_idx])<1.5 && Electron_pfmvaId[BToMuEPi_sel_ele_idx]> -3  && BToMuEPi_fit_pi_pt>0.7  &&  fabs(BToMuEPi_fit_pi_eta)<2  &&  fabs(BToMuEPi_pi_dz)>0.005  &&  fabs(BToMuEPi_pi_dxy)>0.005  &&  fabs(BToMuEPi_pi_dxyS)>3  &&  fabs(BToMuEPi_pi_dzS)>1.5   &&  BToMuEPi_sv_prob>0.001 &&  BToMuEPi_hnl_cos2D>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>3  &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx]/Electron_dxyErr[BToMuEPi_sel_ele_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx]/Electron_dzErr[BToMuEPi_sel_ele_idx])>1 "# &&  deltaR_trgMuPi<2;
	#lowpt missing cut on lowptmvaId           ____________________________&&  fabs(BToMuEPi_pi_DCASig)>5
	
	
	
	lowPt_sel = "   Muon_pt[BToMuEPi_trg_mu_idx]>7&&  Muon_softId[BToMuEPi_trg_mu_idx]==1   &&  fabs(Muon_eta[BToMuEPi_trg_mu_idx])<1.5 && Electron_mvaId[BToMuEPi_sel_ele_idx]> -3  && ProbeTracks_pt[BToMuEPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuEPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPi_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPi_pi_idx])>5  &&  BToMuEPi_sv_prob>0.001 &&  BToMuEPi_hnl_cos2D>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>3  &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx]/Electron_dxyErr[BToMuEPi_sel_ele_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx]/Electron_dzErr[BToMuEPi_sel_ele_idx])>1 "# &&  deltaR_trgMuPi<2;
	
	lowPtAsTT_sel = "  Muon_pt[BToMuEPi_trg_mu_idx]>7 &&  Muon_softId[BToMuEPi_trg_mu_idx]==1 &&  fabs(Muon_eta[BToMuEPi_trg_mu_idx])<1.5 &&  fabs(ProbeTracks_eta[BToMuEPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPi_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPi_pi_idx])>5  &&  BToMuEPi_sv_prob>0.001 && BToMuEPi_hnl_cos2D>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>20  &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx]/Electron_dxyErr[BToMuEPi_sel_ele_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx]/Electron_dzErr[BToMuEPi_sel_ele_idx])>1"# &&  deltaR_trgMuPi<2;
	
	
	tt_sel = "  Muon_pt[BToMuEPiHD_trg_mu_idx]>7 &&  Muon_softId[BToMuEPiHD_trg_mu_idx]==1 &&  fabs(Muon_eta[BToMuEPiHD_trg_mu_idx])<1.5 &&  fabs(ProbeTracks_eta[BToMuEPiHD_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPiHD_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPiHD_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPiHD_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPiHD_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPiHD_pi_idx])>5  &&  BToMuEPiHD_sv_prob>0.001 && BToMuEPiHD_hnl_cos2D>0.99  &&  BToMuEPiHD_sv_lxy/BToMuEPiHD_sv_lxye>20  &&  fabs(ProbeTracks_dxy[BToMuEPiHD_sel_e_idx])>0.008 &&  fabs(ProbeTracks_dz[BToMuEPiHD_sel_e_idx])>0.008 &&  fabs(ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx])>1.5 &&  fabs(ProbeTracks_dzS[BToMuEPiHD_sel_e_idx])>1"# &&  deltaR_trgMuPi<2;


	trgMu_match = "(BToMuEPi_matching_trg_mu_genIdx !=-1)"	
	lep_match = "(BToMuEPi_matching_sel_lep_genIdx !=-1)"	
	pi_match = "(BToMuEPi_matching_pi_genIdx !=-1)"
	
	GenMatch = " && "+trgMu_match+" && "+lep_match+" && "+pi_match
	noGenMatch = " && !"+trgMu_match+" && !"+lep_match+" && !"+pi_match
 	only_trgMuMatch = ""+trgMu_match+" && !"+lep_match+" && !"+pi_match+" && " 	
 	only_eleMatch = "&& !"+trgMu_match+" && "+lep_match+" && !"+pi_match 	
 	only_piMatch = "&& !"+trgMu_match+" && !"+lep_match+" && "+pi_match 	
 	mu_piMatch = "&& "+trgMu_match+" && !"+lep_match+" && "+pi_match 	
 	ele_piMatch = "&& !"+trgMu_match+" && "+lep_match+" && "+pi_match 	
 	ele_muMatch = "&& "+trgMu_match+" && "+lep_match+" && !"+pi_match 	
	
	#have a nice look @ inputs while running
	inputFiles = tqdm(filenames)
#	ROOT.EnableImplicitMT(); # enable multi-threading
	
	# and convert the TTree into an array
	
	
	#use root to array to convert tree in array, while performning a columnar analysis ( if one of the elements of the matrix does not pass a selection, all the row - meaning all the other variables of the instance not passing a selection - is removed)
	#"GenPart_pdgId[BToMuEPi_matching_trg_mu_genIdx]","GenPart_pdgId[BToMuEPi_matching_sel_lep_genIdx]","GenPart_pdgId[BToMuEPi_matching_pi_genIdx]",
	rta = []
	singleB_mu = []
	rta.append(root2array(inputFiles,"Events",branches= mu_branches,object_selection={mu_sig+trg_BPark_mu+mu_sel:mu_branches}))
#	rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel:ele_branches}))
	#rta.append(root2array(inputFiles,"Events",branches= ele_branches+["GenPart_pdgId[BToMuEPi_matching_trg_mu_genIdx]"],object_selection={only_trgMuMatch+pf_sig+trg_BPark_ele+pf_sel:ele_branches+["GenPart_pdgId[BToMuEPi_matching_trg_mu_genIdx]"]}))
	#rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel+only_eleMatch:ele_branches+["GenPart_pdgId[BToMuEPi_matching_sel_lep_genIdx]"]}))
	#rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel+only_piMatch:ele_branches+["GenPart_pdgId[BToMuEPi_matching_pi_genIdx]"]}))
	#rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel+mu_piMatch:ele_branches+["GenPart_pdgId[BToMuEPi_matching_trg_mu_genIdx]","GenPart_pdgId[BToMuEPi_matching_sel_lep_genIdx]"]}))
	#rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel+ele_piMatch:ele_branches+["GenPart_pdgId[BToMuEPi_matching_sel_lep_genIdx]","GenPart_pdgId[BToMuEPi_matching_pi_genIdx]"]}))
	#rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel+ele_muMatch:ele_branches+["GenPart_pdgId[BToMuEPi_matching_trg_mu_genIdx]","GenPart_pdgId[BToMuEPi_matching_pi_genIdx]"]}))
	#rta.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+trg_BPark_ele+pf_sel+GenMatch:ele_branches+["GenPart_pdgId[BToMuEPi_matching_trg_mu_genIdx]","GenPart_pdgId[BToMuEPi_matching_sel_lep_genIdx]","GenPart_pdgId[BToMuEPi_matching_pi_genIdx]"]}))
	singleB = root2array(inputFiles,"Events",branches= scalar_branches)
	#genParts = root2array(inputFiles,"Events",branches= ["GenPart_pdgId"])
	
	
	#print rta_mu[0]
	print len(rta)
	print len(rta[0])
	print len(rta[0][0])
	print len(rta[0][0][0])
	singleB = np.array(singleB.tolist(),dtype="float32")
	#initialize machines for likelihood and BDT scoring
	
	gROOT.LoadMacro("/cmshome/ratramon/Analysis/macros/BinnedLikelihood.C")
	rfileSel = root_open(sys.argv[2], 'w')
	treeSel = Tree("Events", model=Event)
	dphi_cut = [0.03,0.08,0.15,0.2]
	#reader = ROOT.TMVA.Reader( )
	#
	#B_pt = array.array('f',[0])
	#minPt = array.array('f',[0])
	#maxPt = array.array('f',[0])
	#vtxProb = array.array('f',[0])
	#MaxDxyS = array.array('f',[0])
	#hnl_l_dxyS = array.array('f',[0])
	#hnl_l_dxy = array.array('f',[0])
	#hnl_l_dz = array.array('f',[0])
	#hnl_pi_dxy = array.array('f',[0])
	#dilepton_mass = array.array('f',[0])
	#hnl_cos2D = array.array('f',[0])
	#hnl_l_mvaId = array.array('f',[0])
	#hnl
#Mass = array.array('f',[0])
	#hnlLxy = array.array('f',[0])
	#reader.AddVariable( "B_pt",B_pt );
	#reader.AddVariable( "min(hnl_l_pt,hnl_pi_pt)",minPt );
	#reader.AddVariable( "max(hnl_l_pt,hnl_pi_pt)",maxPt );
	#reader.AddVariable( "hnl_vtxProb", vtxProb);
	#reader.AddVariable( "max(abs(hnl_l_dxyS),abs(hnl_pi_dxyS))", MaxDxyS );
	#reader.AddVariable( "abs(hnl_l_dxyS)", hnl_l_dxyS );
	#reader.AddVariable( "abs(hnl_l_dxy)", hnl_l_dxy );
	#reader.AddVariable( "abs(hnl_l_dz)", hnl_l_dz );
	#reader.AddVariable( "abs(hnl_pi_dxy)", hnl_pi_dxy );
	#reader.AddVariable( "dilepton_mass", dilepton_mass );
	#reader.AddVariable( "hnl_cos2D", hnl_cos2D );
	#reader.AddVariable( "hnl_l_mvaId", hnl_l_mvaId );
	#reader.AddSpectator( "hnl_mass", hnlMass );
	#reader.AddSpectator("abs(hnl_lxy)", hnlLxy );
	#reader.BookMVA("BDT method","../macros/dataset/weights//MultipleSignal_BDT.xml")
	entries = np.arange(0,len(rta[0]))

	#pile up weight input 
	histFile= ROOT.TFile.Open("../data/pileup_weight_datatot_sigAug21.root")

	PU_histo = histFile.Get("hist_weight")

	#Mu trg SF input 
	hFile= ROOT.TFile.Open("../data/trgMu_scale_factors_D1.root")

	trgMu_histo = hFile.Get("hist_scale_factor")

	#Mu ID SF input 
	muIdFile= ROOT.TFile.Open("../data/RunABCD_SF_MuonID_2018.root")

	softMu_histo = muIdFile.Get("NUM_SoftID_DEN_genTracks_pt_abseta") #soft muId applied on if muon from B vtx
	looseMu_histo = muIdFile.Get("NUM_LooseID_DEN_genTracks_pt_abseta") #loose muID applied if muon from HNL vtx


	# pNN score
#	paths =  glob.glob("model/trainNN_*_21Sep2022_17h50m22s/")
#	print paths
#	#h5 = "saved-model-0045_val_loss_0.1131_val_acc_0.9614.h5"
#	modelSpecs = []
#	lxy_sel =["hnl_lxy<1","hnl_lxy>1 && hnl_lxy<5","hnl_lxy>5"]
#	lxy_cat =["Under1","Over1Under5","Over5"]
#	sign_cat = ["LepQProd<0","LepQProd>0"]	
#	sign_sel = ["OS","SS"]	
#	for p in paths:
#		splitPath = p.split("_")
#		cat = splitPath[1]+"_"+splitPath[2]
#		modelSpecs.append([cat,NNinput(p)])

	#single mass training
#	path = "./model/trainNN_5vars_singleMass_3p0_ctau100p0/"
#	model_single,qt_single,features_single = NNinput(path) 

#	print(singleB[0:30,3])
#	singleB_mu = singleB_mu.tolist()
#	print(singleB)
	
#	print(len(singleB))
#	print(len(singleB[:,0]))
	for entry in entries:
	    for SigIdx in range (0,1):
	
	      rows = []	
	#trees after selection
	      if ( len(rta[SigIdx][entry][0])!=0):
		treeSel.event = singleB[entry][0]
		treeSel.PV_npvs = singleB[entry][1]
		treeSel.PV_npvsGood = singleB[entry][2]
		treeSel.sigFlag = isMC
#	

		
	        sortBranches(rta[SigIdx][entry],0, sig[SigIdx]) #sort candidates in each event wrt to a chosen variable
		# define useful categorization flags
		if (rta[SigIdx][entry][9][0] < 1): 
			treeSel.LxyBin = 0
		elif (rta[SigIdx][entry][9][0] > 1 and rta[SigIdx][entry][9][0] < 5):
			treeSel.LxyBin = 1
		elif (rta[SigIdx][entry][9][0] > 5):# and rta[SigIdx][entry][9][0] < 20):
			treeSel.LxyBin = 2
	#	elif (rta[SigIdx][entry][9][0] > 20):
	#		treeSel.LxyBin = 3
#		if ( not isMC):
#			if ((rta[SigIdx][entry][14][0] > 0.05 and rta[SigIdx][entry][10][0]> 0.996)):
#	
#		
#				if(isQCD):
#					treeSel.ABCDReg = 0
#				else:
#					treeSel.ABCDReg = -1
#			elif (rta[SigIdx][entry][14][0] > 0.05 and rta[SigIdx][entry][10][0]< 0.996):
#				treeSel.ABCDReg = 1 
#			elif (rta[SigIdx][entry][14][0] < 0.05 and rta[SigIdx][entry][10][0]< 0.996):
#				treeSel.ABCDReg = 2 
#			elif (rta[SigIdx][entry][14][0] < 0.05 and  rta[SigIdx][entry][10][0]> 0.996):
#				treeSel.ABCDReg = 3
#		else:
		
		treeSel.ABCDReg = -1
		treeSel.QCDweight = QCD_weight
########	if (SigIdx==0):
########		if (len(rta[SigIdx+1][entry][0])!=0):
########			treeSel.PFMuOverlap = 1
########		else:
########			treeSel.PFMuOverlap = 0
########	
########	else:
########		if (len(rta[SigIdx-1][entry][0])!=0 ):

########			treeSel.PFMuOverlap = 1
########		else:
########			treeSel.PFMuOverlap = 0
	#	if(abs(dPhi(rta[SigIdx][entry][33][0],rta[SigIdx][entry][39][0]) )> dphi_cut[treeSel.LxyBin]):
	#		continue
		
	#	if ( rta[SigIdx][entry][0][0]>5 and rta[SigIdx][entry][0][0]<5.5 ): # control region-------UPPER SIDEBAND
	#	if (rta[SigIdx][entry][40][0]==0 ): # control region q ! = 0 
	#		if (SigIdx == 1 and rta[SigIdx][entry][29][0] >-2 ):
	#	     		continue
	#		elif (SigIdx !=1 ):
	#		continue
		#print treeSel.LxyBin	
	#print("for event %d and displacement %d, signal category %d",k,entry,SigIdx)
	
		#fill the output tree	
	
	        treeSel.evtIdx =entry 
	        treeSel.nCand=1
	        treeSel.Type=0
 	        treeSel.LepQProd= rta[SigIdx][entry][41][0]*rta[SigIdx][entry][42][0]                                                       
                treeSel.B_eta= rta[SigIdx][entry][2][0]                                                       
                treeSel.B_phi= rta[SigIdx][entry][3][0]                                                       
                treeSel.TrgMu_pt= rta[SigIdx][entry][4][0]                                                    
                treeSel.TrgMu_eta= rta[SigIdx][entry][5][0]                                                   
                treeSel.TrgMu_phi= rta[SigIdx][entry][6][0]                                                   
		if (isMC ==0 and isQCD==0):
			treeSel.PU_weight = 1 
	               	treeSel.Trg_weight = 1
		else: 
			treeSel.PU_weight = PU_histo.GetBinContent(PU_histo.GetXaxis().FindBin(treeSel.PV_npvs)) 
	               	treeSel.Trg_weight = trgMu_histo.GetBinContent(trgMu_histo.GetXaxis().FindBin(treeSel.TrgMu_pt),trgMu_histo.GetYaxis().FindBin(abs(rta[SigIdx][entry][58][0])))
                treeSel.dr_trgMu_lep= dR(rta[SigIdx][entry][6][0],rta[SigIdx][entry][24][0],rta[SigIdx][entry][5][0],rta[SigIdx][entry][23][0])                       
 	        treeSel.dz_trgMu_lep=  rta[SigIdx][entry][51][0]            
 	        treeSel.hnl_mass= rta[SigIdx][entry][7][0]                                                    
 	        treeSel.hnl_pt= rta[SigIdx][entry][8][0]      
  	        treeSel.hnl_eta= rta[SigIdx][entry][43][0]      
  	        treeSel.hnl_phi= rta[SigIdx][entry][44][0]      
  	        treeSel.hnl_charge = rta[SigIdx][entry][40][0]                                         
  	        treeSel.hnl_lxy= rta[SigIdx][entry][9][0]                                                     
  	        treeSel.hnl_ct= rta[SigIdx][entry][55][0]*10                                                     
  	        treeSel.hnl_to_trgmu= rta[SigIdx][entry][61][0]                                                     
  	        treeSel.hnl_cos2D= rta[SigIdx][entry][10][0]                                                  
  	        treeSel.hnl_lxy_sig = rta[SigIdx][entry][11][0]                                               
  	        treeSel.hnl_drLepPi = rta[SigIdx][entry][12][0]                                               
  	        treeSel.dr_trgmu_hnl = rta[SigIdx][entry][13][0]                                              
  	        treeSel.hnl_vtxProb= rta[SigIdx][entry][14][0]                                                
  	        treeSel.hnl_vtxChi2= rta[SigIdx][entry][15][0]                                                
  	        treeSel.hnl_vx =  rta[SigIdx][entry][16][0]                                                   
  	        treeSel.hnl_ex =  rta[SigIdx][entry][17][0]                                                   
  	        treeSel.hnl_vy =  rta[SigIdx][entry][18][0]                                                   
  	        treeSel.hnl_ey =  rta[SigIdx][entry][19][0]                                                   
  	        treeSel.hnl_vz = rta[SigIdx][entry][20][0]                                                    
  	        treeSel.hnl_ez = rta[SigIdx][entry][21][0]                                                    
  	        treeSel.hnl_l_pt = rta[SigIdx][entry][22][0]                                                  
  	        treeSel.hnl_l_eta =rta[SigIdx][entry][23][0]                                                  
  	        treeSel.hnl_l_phi =rta[SigIdx][entry][24][0]                                                  
  	        treeSel.hnl_l_dz = rta[SigIdx][entry][25][0]   
  	        if (SigIdx==0 or SigIdx==1):
  	        					                                              
  	          treeSel.hnl_l_dzS =rta[SigIdx][entry][25][0]/rta[SigIdx][entry][26][0]
  	        else:                                             
  	          treeSel.hnl_l_dzS =rta[SigIdx][entry][26][0]
  	          treeSel.hnl_l_dxy =rta[SigIdx][entry][27][0]                                                  
  	        if (SigIdx==0 or SigIdx==1):
  	        					                                               
  	          treeSel.hnl_l_dxyS =rta[SigIdx][entry][27][0]/rta[SigIdx][entry][28][0]
  	        else:                                            
  	          treeSel.hnl_l_dxyS =rta[SigIdx][entry][28][0]
  	        
  	        if (SigIdx <3):
  	          treeSel.hnl_l_DCAS = -99
  	        else:  
  	          treeSel.hnl_l_DCAS =rta[SigIdx][entry][29][0] 
  	
  	        if (SigIdx ==0):
	#	 print(rta[SigIdx][entry][29][0])
  	         treeSel.hnl_l_mvaId= rta[SigIdx][entry][29][0]    
  	         treeSel.hnl_l_nHits= -99 
  	        elif (SigIdx ==1):
  	         treeSel.hnl_l_mvaId= rta[SigIdx][entry][30][0]  
  	         treeSel.hnl_l_nHits= -99 
  	         treeSel.hnl_l_nHits= rta[SigIdx][entry][30][0] 
  	        treeSel.hnl_pi_pt =rta[SigIdx][entry][31][0]                                                  
  	        treeSel.hnl_pi_eta= rta[SigIdx][entry][32][0]                                                 
  	        treeSel.hnl_pi_phi = rta[SigIdx][entry][33][0]                                                
  	        treeSel.hnl_pi_dz =rta[SigIdx][entry][34][0]                                                  
  	        treeSel.hnl_pi_dzS = rta[SigIdx][entry][35][0]                                                
  	        treeSel.hnl_pi_dxy= rta[SigIdx][entry][36][0]                                                 
  	        treeSel.hnl_pi_dxyS = rta[SigIdx][entry][37][0]                                      
  	        treeSel.hnl_pi_DCAS = rta[SigIdx][entry][38][0]                                      
  	        treeSel.hnl_pi_nHits = rta[SigIdx][entry][48][0]                                      
  	        treeSel.dphi_pi_fitpi = dPhi(rta[SigIdx][entry][33][0],rta[SigIdx][entry][39][0])  
  	   #    treeSel.genMatch_muTrg = rta[SigIdx][entry][63][0]
  	   #    treeSel.genMatch_l = rta[SigIdx][entry][64][0]
  	   #    treeSel.genMatch_pi = rta[SigIdx][entry][65][0]
	   #    if treeSel.genMatch_muTrg != -1:
  	   #   		 treeSel.pdgId_muTrg = genParts[entry][0][treeSel.genMatch_muTrg]
	   #    else:
  	   #   		 treeSel.pdgId_muTrg = -1
	   #    if treeSel.genMatch_l != -1:
  	   #   		 treeSel.pdgId_l = genParts[entry][0][treeSel.genMatch_l]
	   #    else:
  	   #   		 treeSel.pdgId_l = -1
	   #    if treeSel.genMatch_pi != -1:
	   #            treeSel.pdgId_pi = genParts[entry][0][treeSel.genMatch_pi]
  	   #    else:
	   #    	treeSel.pdgId_pi = -1
	   #    treeSel.genMotherPdgId_muTrg = rta[SigIdx][entry][66][0]
  	   #    treeSel.genMotherPdgId_l = rta[SigIdx][entry][67][0]
  	   #    treeSel.genMotherPdgId_pi = rta[SigIdx][entry][68][0]
  	#	if treeSel.genMatch_muTrg == -1 and treeSel.genMatch_l == -1 and treeSel.genMatch_pi == -1:	
  	#		cat = 0			
  	#	if treeSel.genMatch_muTrg != -1 and treeSel.genMatch_l == -1 and treeSel.genMatch_pi == -1:	
  	#		cat = 1			
  	#	if treeSel.genMatch_muTrg == -1 and treeSel.genMatch_l != -1 and treeSel.genMatch_pi == -1:	
  	#		cat = 2			
  	#	if treeSel.genMatch_muTrg == -1 and treeSel.genMatch_l == -1 and treeSel.genMatch_pi != -1:	
  	#		cat = 3			
  	#	if treeSel.genMatch_muTrg != -1 and treeSel.genMatch_l != -1 and treeSel.genMatch_pi == -1:	
  	#		cat = 4			
  	#	if treeSel.genMatch_muTrg != -1 and treeSel.genMatch_l == -1 and treeSel.genMatch_pi != -1:	
  	#		cat = 5			
  	#	if treeSel.genMatch_muTrg == -1 and treeSel.genMatch_l != -1 and treeSel.genMatch_pi != -1:	
  	#		cat = 6		
  	#	if treeSel.genMatch_muTrg != -1 and treeSel.genMatch_l != -1 and treeSel.genMatch_pi != -1:	
  	#		cat = 7			
  	#	treeSel.genMatch_cat = cat
  	   #     treeSel.cosTheta_star = rta[SigIdx][entry][48][0]                                      
  		treeSel.bdt = -99 
  		#print(B_mass)
  	        if SigIdx ==-1:
  			
  			treeSel.B_pt= rta[SigIdx][entry][1][0]                                                        
  			treeSel.B_mass= rta[SigIdx][entry][0][0]                                                      
  			treeSel.dilepton_mass = rta[SigIdx][entry][45][0]                                      
  			treeSel.dilepton_pt = rta[SigIdx][entry][46][0]    
  		#	treeSel.dilepton_mass = TwoBodyMass(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                                  
  		#	treeSel.dilepton_pt = TwoBodyPt(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)    
  			treeSel.likelihood =  -99
  			treeSel.toEle = -1                 
  			treeSel.toMu =  -1                 


  		elif SigIdx==0 :

			#define vector 
				
  			treeSel.B_pt= rta[SigIdx][entry][1][0]                                                        
  			treeSel.B_mass= rta[SigIdx][entry][0][0]                                                      
  			treeSel.toMu =  -1                 
			
			
  			treeSel.B_pt= rta[SigIdx][entry][1][0]                                                        
  			treeSel.B_mass= rta[SigIdx][entry][0][0]                                                      
  			treeSel.dilepton_mass = rta[SigIdx][entry][45][0]                                      
  			treeSel.dilepton_pt = rta[SigIdx][entry][46][0]  
			 
  			treeSel.dilepton_mass = TwoBodyMass(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                                  
	  		
  			treeSel.BlepPi_mass = TwoBodyMass(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                                  
  			treeSel.BlepPi_pt = TwoBodyPt(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)    
			
			if (treeSel.hnl_to_trgmu==0): 
				
  	          		treeSel.hnl_l_dzS =rta[SigIdx][entry][25][0]/rta[SigIdx][entry][26][0]
  	          		treeSel.hnl_l_dz =rta[SigIdx][entry][25][0]
  	          		treeSel.hnl_l_dxyS =rta[SigIdx][entry][27][0]/rta[SigIdx][entry][28][0]
  	          		treeSel.hnl_l_dxy =rta[SigIdx][entry][27][0]
  	          		treeSel.muId_weight = softMu_histo.GetBinContent(softMu_histo.GetXaxis().FindBin(treeSel.TrgMu_pt),softMu_histo.GetYaxis().FindBin(abs(treeSel.TrgMu_eta)))
	  			point_l  = np.array( [np.float(treeSel.hnl_cos2D),np.float(abs(treeSel.hnl_lxy_sig)),np.float(treeSel.hnl_vtxProb),np.float(treeSel.hnl_pi_pt)])
  				treeSel.BlepPi_mass = TwoBodyMass(SigIdx,treeSel.hnl_pi_pt,treeSel.hnl_pi_eta,treeSel.hnl_pi_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                                  
  				treeSel.BlepPi_pt = TwoBodyPt(SigIdx,treeSel.hnl_pi_pt,treeSel.hnl_pi_eta,treeSel.hnl_pi_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)   
 	        		treeSel.dr_Blep_pi= dR(treeSel.hnl_pi_phi,treeSel.TrgMu_phi,treeSel.hnl_pi_eta,treeSel.TrgMu_eta)  #primary lepton is muon
				 
			elif (treeSel.hnl_to_trgmu==1): 
  	          		treeSel.hnl_l_dzS =rta[SigIdx][entry][57][0]
  	          		treeSel.hnl_l_dz =rta[SigIdx][entry][55][0]
  	          		treeSel.hnl_l_dxyS =rta[SigIdx][entry][58][0]
  	          		treeSel.hnl_l_dxy =rta[SigIdx][entry][56][0]
  	          		treeSel.muId_weight = looseMu_histo.GetBinContent(looseMu_histo.GetXaxis().FindBin(treeSel.TrgMu_pt),looseMu_histo.GetYaxis().FindBin(abs(treeSel.TrgMu_eta)))
	  			point_l  = np.array( [np.float(treeSel.hnl_cos2D),np.float(abs(treeSel.hnl_lxy_sig)),np.float(treeSel.hnl_vtxProb),np.float(treeSel.hnl_pi_pt)])
  				treeSel.BlepPi_mass = TwoBodyMass(SigIdx,treeSel.hnl_pi_pt,treeSel.hnl_pi_eta,treeSel.hnl_pi_phi,mass[SigIdx],treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,ELECTRON_MASS)                                  
  				treeSel.BlepPi_pt = TwoBodyPt(SigIdx,treeSel.hnl_pi_pt,treeSel.hnl_pi_eta,treeSel.hnl_pi_phi,mass[SigIdx],treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,ELECTRON_MASS)    
 	        		treeSel.dr_Blep_pi= dR(treeSel.hnl_pi_phi,treeSel.hnl_l_phi,treeSel.hnl_pi_eta,treeSel.hnl_l_eta) #primary lepton is electron
	  			#point  = np.array( [np.float(treeSel.hnl_cos2D),np.float(max(np.abs(treeSel.hnl_l_dxyS),np.abs(treeSel.hnl_pi_dxyS))),np.float(treeSel.hnl_vtxProb),np.float(treeSel.hnl_pi_pt)])
  			likelihood = ROOT.BinnedLikelihood(SigIdx+1,int(sys.argv[3]),int(sys.argv[4]),point_l)
  			treeSel.likelihood =  likelihood if not math.isnan(likelihood ) else -99                           
  			treeSel.toEle = -1                 
  			treeSel.toMu =  -1          
  	#		B_pt[0] = treeSel.B_pt
  	#                minPt[0] = min(treeSel.hnl_pi_pt,treeSel.hnl_l_pt)
  	#                maxPt[0] = max(treeSel.hnl_pi_pt,treeSel.hnl_l_pt)
  	#                hnlMass[0] = treeSel.hnl_mass
  	#                hnlLxy[0] = abs(treeSel.hnl_lxy)
  	
  	
  		#	treeSel.bdt = reader.EvaluateMVA("BDT method");
  		elif SigIdx==1:
  			treeSel.B_pt= rta[SigIdx][entry][1][0]                                                        
  			B_mass = TwoBodyMass(SigIdx,treeSel.hnl_pt,treeSel.hnl_eta,treeSel.hnl_phi,mass_hnl,treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                         
  			B_mass = TwoBodyMass(SigIdx,treeSel.hnl_pt,treeSel.hnl_eta,treeSel.hnl_phi,mass_hnl,treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                         
  			B_pt = TwoBodyPt(SigIdx,treeSel.hnl_pt,treeSel.hnl_eta,treeSel.hnl_phi,mass_hnl,treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                         
  	        	treeSel.B_pt= B_pt                                                   
  	        	treeSel.B_mass= B_mass                                                  
  			treeSel.dilepton_mass = TwoBodyMass(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)                                      
  	       		treeSel.dilepton_pt = TwoBodyPt(SigIdx,treeSel.hnl_l_pt,treeSel.hnl_l_eta,treeSel.hnl_l_phi,mass[SigIdx],treeSel.TrgMu_pt,treeSel.TrgMu_eta,treeSel.TrgMu_phi,MUON_MASS)    
  			                                 
  			point  = np.array( [np.float(treeSel.hnl_cos2D),np.float(max(abs(treeSel.hnl_l_dzS),abs(treeSel.hnl_pi_dzS))),np.float(max(abs(treeSel.hnl_l_dxyS),abs(treeSel.hnl_pi_dxyS))),np.float(min(treeSel.hnl_l_nHits,treeSel.hnl_pi_nHits)),np.float(treeSel.hnl_vtxProb)])
  			likelihood = ROOT.BinnedLikelihood(SigIdx,int(sys.argv[3]),int(sys.argv[4]),point)
  			treeSel.likelihood =  likelihood if not math.isnan(likelihood ) else -99   
  			treeSel.toEle =-1# rta[SigIdx][entry][52][0]                        
  			treeSel.toMu = -1# rta[SigIdx][entry][53][0]                       
	        treeSel.fill()
		if (do_orthogonal):
			break
	      else:
		continue
	
			
	treeSel.Write() 
	#loop on array entries after preselection
	#genMatch_titles = [" #mu_{trg} l #pi - 0","#mu_{trg} - 1, l #pi - 0","l - 1, #mu_{trg}, #pi - 0","#pi, #mu_{trg} - 1, l - 0","#mu_{trg}, l - 1, #pi - 0","#mu_{trg}, #pi - 1,  l - 0 ","pi, l - 1 , #mu_{trg} - 0","#mu_{trg}, l, pi 1"]
#        pool = Pool(processes=20)
#        pool.map(fillTreeEntry, entries)
#        pool.close()
#        pool.join()
##saveHistos(rta_mu, histMu,1, "BToMuMuPi")
#saveHistos(rta_muSel, histMuSel,0,"BToMuMuPi_sel")
#saveHistos(rta_pf, histPF,1, "BToMuPFPi")
#saveHistos(rta_pfSel, histPFSel,0,"BToMuPFPi_sel")
#saveHistos(rta_lowPt, histLowpT,1,"BToMuLowPtPi")
#saveHistos(rta_lowPtSel, histLowpTSel,0,"BToMuLowPtPi_sel")
#saveHistos(rta_tt, histTT,1,"BToMuTPi")
#saveHistos(rta_ttSel, histTTSel,1,"BToMuTPi")
#ree.Write()

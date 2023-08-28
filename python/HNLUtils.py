import os
import json 
os.sys.path.insert(1,"~/Analysis/python/")
from decays import HNLDecays,Decays

import numpy as np
import math
from scipy.stats import expon  



# constants 
const_pi = math.pi
const_hbar = 6.582119569e-22 * 1e-03 # GeV s  # from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
const_c = 299792458. # m / s                  # from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf


def BToMuX(mass=-99,vv =1, B_flav=''):

    dec    = Decays(mass=mass, mixing_angle_square=vv)
    br = dec.BR_tot_mu
    if 'Bu' in B_flav :
    	br = dec.BR_B_mu
    elif 'Bd' in B_flav :
    	br = dec.BR_B0_mu
    elif 'Bs' in B_flav :
    	br = dec.BR_Bs_mu
    elif 'Bc' in B_flav :
    	br = dec.BR_tot_mu_Bc
	
    return br

def BToEX(mass=-99,vv =1,B_flav= ''):

    dec    = Decays(mass=mass, mixing_angle_square=vv)
    br = dec.BR_tot_e
    if 'Bu' in B_flav :
    	br = dec.BR_B_e
    elif 'Bd' in B_flav :
    	br = dec.BR_B0_e
    elif 'Bs' in B_flav :
    	br = dec.BR_Bs_e
    elif 'Bc' in B_flav :
    	br = dec.BR_tot_e_Bc

    return br

def BcToMuX(mass=-99,vv =1):

    dec    = Decays(mass=mass, mixing_angle_square=vv)
    br = dec.BR_tot_mu_Bc
    return br

def BcToEX(mass=-99,vv =1):

    dec    = Decays(mass=mass, mixing_angle_square=vv)
    br = dec.BR_tot_e_Bc
    return br

def BToMuX_partialBR(mass=-99,vv =1):

    dec    = Decays(mass=mass, mixing_angle_square=vv)
    br = [0]*4
 #  br[0] = dec.B_to_uHNL.BR *1000
 #  br[1] = dec.B_to_D0uHNL.BR*1000
 #  br[2] = dec.B_to_D0staruHNL.BR*1000
 #  br[3] = dec.B_to_pi0uHNL.BR*1000
 # #br[4] = dec.B_to_rho0uHNL.BR*1000
   #br[0] = dec.B0_to_piuHNL.BR *1000
   #br[1] = dec.B0_to_DuHNL.BR*1000
   #br[2] = dec.B0_to_DstaruHNL.BR*1000
   #br[3] = dec.B0_to_rhouHNL.BR*1000
    br[0] = dec.Bs_to_DsuHNL.BR*1000
    br[1] = dec.Bs_to_DsstaruHNL.BR*1000
    br[2] = dec.Bs_to_KuHNL.BR *1000
    br[3] = dec.Bs_to_KstaruHNL.BR*1000
    return br

def BToEX_partialBR(mass=-99,vv =1):

    dec    = Decays(mass=mass, mixing_angle_square=vv)
    br = [0]*4
#   br[0] = dec.B_to_eHNL.BR*1000
#   br[1] = dec.B_to_D0eHNL.BR*1000
#   br[2] = dec.B_to_D0stareHNL.BR*1000
#   br[3] = dec.B_to_pi0eHNL.BR*1000
#   br[4] = dec.B_to_rho0eHNL.BR*1000
    br[0] = dec.Bs_to_DseHNL.BR*1000
    br[1] = dec.Bs_to_DsstareHNL.BR*1000
    br[2] = dec.Bs_to_KeHNL.BR *1000
    br[3] = dec.Bs_to_KstareHNL.BR*1000
#    for b in br:
#	print(b)
	
    return br


def ctau_from_gamma(gamma):
    tau_natural = 1. / gamma                  # 1/GeV
    tau = tau_natural * const_hbar            # s
    ctau = tau * const_c * 1000               # mm
    return ctau

def gamma_total(mass,vv):
    '''
    Total width for N (Dirac)
    '''
    gamma_total =  HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['tot']   # GeV
    return gamma_total

def gamma_partial(mass,vv,f):
    '''
    Partial width for N->mupi (Dirac)
    '''
    gamma_partial_mu = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['mupi'] # GeV
    gamma_partial_e = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['epi'] # GeV
    gamma_partial_tau = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['taupi'] # GeV
    print "mu",(gamma_partial_mu)
    print "e",(gamma_partial_e)
    print "tau",(gamma_partial_tau)

    if f == 'mu':
  	  gamma_partial =  gamma_partial_mu
    else:
  	  gamma_partial =  gamma_partial_e
	
    return gamma_partial

def BR_HNLmupion(mass): # vv is irrelevant, as it cancels out in the ratio
    print gamma_partial(mass=mass,vv=1.,f='mu')#/gamma_total(mass=mass,vv=1.)
    return gamma_partial(mass=mass,vv=1.,f='mu')/gamma_total(mass=mass,vv=1.)

def BR_HNLelepion(mass): # vv is irrelevant, as it cancels out in the ratio
    #print gamma_partial(mass=mass,vv=1.)/gamma_total(mass=mass,vv=1.)
    return gamma_partial(mass=mass,vv=1.,f='ele')/gamma_total(mass=mass,vv=1.)

def getVV(mass=-99.,ctau=-99.,ismaj=True):
    '''
    Helper function to go from ctau,m -> vv
    '''
    mult = 2. if ismaj else 1.
    ref_m = 1. # GeV
    ref_vv = 1. 
    ref_ctau = ctau_from_gamma(mult*gamma_total(mass=ref_m,vv=ref_vv))
#    print (ref_ctau)

    k = ref_ctau * np.power(ref_m,5) * ref_vv

    return k/(pow(mass, 5) * ctau)

def getCtau(mass=-99.,v2=-99.,ismaj=True):
    '''
    Helper function to go from vv,m -> ctau
    '''
    mult = 2. if ismaj else 1.
    ref_m = 1. # GeV
    ref_vv = 1. 
    ref_ctau = ctau_from_gamma(mult*gamma_total(mass=ref_m,vv=ref_vv))
#    print (ref_ctau)

    k = ref_ctau * np.power(ref_m,5) * ref_vv

    return k/(pow(mass, 5) * v2)


def SigYield(acceptance,lumi,HNLMass,HNLctau,filter_eff):
	
    Xsec_b = 472.8 * pow(10,9) #fb
    f_u = 0.4
    FullLumi = 41.6 #fb -1
    ProcessedLumi = lumi #fb -1
    br_MuNuInclusive = BToMuX(HNLMass,1)
    br_ENuInclusive = BToEX(HNLMass,1)
    v2 =getVV(HNLMass,HNLctau, True) 
    print(br_MuNuInclusive)
    v2 = getVV(HNLMass,HNLctau,True)
    print(v2)
    br_NtoMuPi =BR_HNLmupion(HNLMass) 
    br_NtoEPi =BR_HNLelepion(HNLMass) 
    Nsig = acceptance*2 * ProcessedLumi * Xsec_b/f_u*(br_MuNuInclusive *   v2   *  br_NtoMuPi)  * filter_eff
    print(Nsig)

def sigPars(mass,cat):

    with open('/cmshome/ratramon/Analysis/data/sig_par.txt','r') as f:
	pars = f.read()

    pars = json.loads(pars)

    mean = pars[cat]['mean'][1]*mass + pars[cat]['mean'][0]
    sigma = pars[cat]['sigma'][1]*mass + pars[cat]['sigma'][0]
    alpha1 = 1.2
    n1 = 2.5
    alpha2 = 1.5
    n2 = 4.0
    print mean, sigma
    return mean,sigma,alpha1,n1,alpha2,n2

def coupling_scaler(f_mu_ref, f_e_ref, f_tau_ref,f_mu,f_e,f_tau,ch,mass):

        vv = 1.0
    
	gamma_partial_mu = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['mupi'] # GeV
	gamma_partial_e = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['epi'] # GeV
	gamma_partial_tau = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['taupi'] # GeV
	gamma_factor =( gamma_partial_mu * f_mu_ref + gamma_partial_e * f_e_ref + gamma_partial_tau * f_tau_ref )/( gamma_partial_mu * f_mu + gamma_partial_e * f_e + gamma_partial_tau * f_tau)


	print gamma_factor
	if ch == 'mumu':
		return (f_mu * gamma_factor)**2
	if ch == 'mue':
		return f_mu*f_e*( gamma_factor)**2


vv = getVV(1,300,True)
ctau = getCtau(4.5,0.001457,True)

print 'ctau',ctau
#sigPars(2.0,'lxysig0to50_SS')
#gamma = gamma_total(2,vv)
#print 'total gamma', gamma 
#gamma = gamma_partial(4.5,vv,'mu')

f_mu = coupling_scaler(1,0,0,0.33,0.33,0.33,'mumu',4.5)
f_e = coupling_scaler(0.5,0.5,0,0.33,0.33,0.33,'mue',4.5)


print 'corrections',f_mu, f_e
#br_mu= BR_HNLmupion(4.5)
#br_e= BR_HNLelepion(4.5)
#print br_mu, br_e
#print("BR electrons %f,BR muons %f"%(br_e,br_mu))
#print(BToMuX(3))
print(BToEX(3))
print(BToEX(4))
#print(BToMuX_partialBR(2))
#print(BToEX_partialBR(2))
#SigYield(0.044,41.6,3,100,0.193)

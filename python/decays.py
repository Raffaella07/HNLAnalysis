'''
This script collects the particle definition and creates the Decays class, that contains
the B->HNL branching ratios of interest, as a function of HNL mass and the mixing angle
'''
import os
#os.sys.path.insert(1,"~/Analysis/python")
from objects import Particle, Decay, HNLDecay
from QCD_corr_table import QCD_corr_table
import math


# use approximations of Shuve Peskin
global PESKIN
PESKIN = False

## CONSTANTS ##

# masses [in GeV] 
# mesons
m_B_pdg       = 5.27934 
m_B0_pdg      = 5.27965
m_B_sub_c_pdg = 6.2749
m_B_sub_s_pdg = 5.36688
m_D_pdg       = 1.86965  
m_D0_pdg      = 1.86483
m_D_sub_s_pdg = 1.96834
m_D0star_pdg  = 2.00685 # check
m_Dstar_pdg   = 2.001021
m_D_sub_sstar_pdg = 2.1122
m_K_pdg       = 0.493677
m_K0_pdg      = 0.497611
m_Kstar_pdg   = 0.89166 
m_K0star_pdg  = 0.89555 #0.824 #1.425
m_rho0_pdg    = 0.77526 #0.769 # 0.770 # 0.77526 # make sure this is the correct meson to consider # in GeV
m_rho_pdg     = 0.77511
m_pi0_pdg     = 0.1349768
m_pi_pdg      = 0.13957039 if not PESKIN else 0

# leptons
m_el_pdg  = 0.510999 * 1e-3 if not PESKIN else 0
m_mu_pdg  = 0.105658        if not PESKIN else 0
m_tau_pdg = 1.77686         if not PESKIN else 0

# quarks
m_uq_pdg = 2.16 * 1e-3 if not PESKIN else 0
m_dq_pdg = 4.67 * 1e-3 if not PESKIN else 0
m_sq_pdg = 93   * 1e-3 if not PESKIN else 0
m_cq_pdg = 1.27        if not PESKIN else 0
m_bq_pdg = 4.18
m_tq_pdg = 172.76  # here for completeness

# meson decay constant --> to be checked
# in [GeV]
dc_B       = 0.1871 # Table8 Bondarenko et al. 0.19   # 176 MeV according to pdg
dc_B_sub_c = 0.434  # Table8 Bondarenko et al. # 0.48
dc_D       = 0.212  # Table8 Bondarenko et al. # 0.2226   
dc_D_sub_s = 0.249  # Table8 Bondarenko et al

# lifetime in [GeV-1]
sToGeVconv       = 1.5198e24
lifetime_B       = 2.49e12 
lifetime_B0      = 1519e-15  * sToGeVconv 
lifetime_B_sub_c = 0.51e-12  * sToGeVconv
lifetime_B_sub_s = 1.527e-12 * sToGeVconv
lifetime_D       = 1.04e-12  * sToGeVconv
lifetime_D_sub_s = 5.04e-13  * sToGeVconv
lifetime_D0      = 4.101e-13 * sToGeVconv

# B-fractions
fraction_B       = 0.4
fraction_B0      = 0.4
fraction_B_sub_s = 0.1
fraction_B_sub_c = 0.00261

# CKM matrix elements
Vud_pdg = 0.97417 if not PESKIN else 1.
Vus_pdg = 0.2248  if not PESKIN else 0.
Vub_pdg = 0.00409 if not PESKIN else 0.
Vcb_pdg = 40.5e-3 if not PESKIN else 0.
Vcd_pdg = 0.220   if not PESKIN else 0.
Vcs_pdg = 0.995   if not PESKIN else 1.
Vtd_pdg = 0.0080  if not PESKIN else 0. 
Vts_pdg = 0.00388 if not PESKIN else 0.
Vtb_pdg = 1.01    if not PESKIN else 1.

## PARTICLES ##

# beauty mesons
B_meson           = Particle('B_meson'          , 'meson', m_B_pdg          , dc_B      , lifetime_B               , fraction=fraction_B)
B0_meson          = Particle('B0_meson'         , 'meson', m_B0_pdg                     , lifetime=lifetime_B0     , fraction=fraction_B0)
B_sub_c_meson     = Particle('B_sub_c_meson'    , 'meson', m_B_sub_c_pdg    , dc_B_sub_c, lifetime_B_sub_c         , fraction=fraction_B_sub_c)
B_sub_s_meson     = Particle('B_sub_s_meson'    , 'meson', m_B_sub_s_pdg                , lifetime=lifetime_B_sub_s, fraction=fraction_B_sub_s)

# charm mesons
D_meson           = Particle('D_meson'          , 'meson', m_D_pdg          , dc_D      , lifetime_D)
D0_meson          = Particle('D0_meson'         , 'meson', m_D0_pdg         , dc_D      , lifetime_D0)
D_sub_s_meson     = Particle('D_sub_s_meson'    , 'meson', m_D_sub_s_pdg    , dc_D_sub_s, lifetime_D_sub_s)
D0star_meson      = Particle('D0star_meson'     , 'meson', m_D0star_pdg)
Dstar_meson       = Particle('Dstar_meson'      , 'meson', m_Dstar_pdg)
D_sub_sstar_meson = Particle('D_sub_sstar_meson', 'meson', m_D_sub_sstar_pdg)

# strange mesons
K_meson           = Particle('K_meson'          , 'meson', m_K_pdg)
K0_meson          = Particle('K0_meson'         , 'meson', m_K0_pdg)
Kstar_meson       = Particle('Kstar_meson'      , 'meson', m_Kstar_pdg)
K0star_meson      = Particle('K0star_meson'     , 'meson', m_K0star_pdg)

# light mesons
rho0_meson        = Particle('rho0_meson'       , 'meson', m_rho0_pdg)
rho_meson         = Particle('rho_meson'        , 'meson', m_rho_pdg)
pi0_meson         = Particle('pi0_meson'        , 'meson', m_pi0_pdg) 
pi_meson          = Particle('pi_meson'         , 'meson', m_pi_pdg) 

# leptons
el                = Particle('el'               , 'lepton', m_el_pdg) 
mu                = Particle('mu'               , 'lepton', m_mu_pdg) 
tau               = Particle('tau'              , 'lepton', m_tau_pdg) 

# neutrinos
nu_el             = Particle('nu_el'            , 'neutrino', 0.)
nu_mu             = Particle('nu_mu'            , 'neutrino', 0.)
nu_tau            = Particle('nu_tau'           , 'neutrino', 0.)

# quarks
uq                = Particle('uq'              , 'quark',  m_uq_pdg)
dq                = Particle('dq'              , 'quark',  m_dq_pdg)
cq                = Particle('cq'              , 'quark',  m_cq_pdg)
sq                = Particle('sq'              , 'quark',  m_sq_pdg)
bq                = Particle('bq'              , 'quark',  m_bq_pdg)
tq                = Particle('tq'              , 'quark',  m_tq_pdg)

## DECAYS ##

class Decays(object):
  def __init__(self, mass, mixing_angle_square): # add model label? 
    self.mass = mass
    self.mixing_angle_square = mixing_angle_square

    # define the HNL
    hnl = Particle('hnl', 'lepton', self.mass)
    
    # get the model
    V_el_square = self.mixing_angle_square
    V_mu_square = self.mixing_angle_square
    V_tau_square = self.mixing_angle_square

    # list of the decays of interest
    # leptonic
    self.B_to_eHNL  = Decay(B_meson,       el, hnl, V_el_square, Vub_pdg, 'leptonic') 
    self.Bc_to_eHNL = Decay(B_sub_c_meson, el, hnl, V_el_square, Vcb_pdg, 'leptonic') 
    self.D_to_eHNL  = Decay(D_meson,       el, hnl, V_el_square, Vcd_pdg, 'leptonic') 
    self.Ds_to_eHNL = Decay(D_sub_s_meson, el, hnl, V_el_square, Vcs_pdg, 'leptonic')
    
    self.B_to_uHNL  = Decay(B_meson,       mu, hnl, V_mu_square, Vub_pdg, 'leptonic') 
    self.Bc_to_uHNL = Decay(B_sub_c_meson, mu, hnl, V_mu_square, Vcb_pdg, 'leptonic') 
    self.D_to_uHNL  = Decay(D_meson,       mu, hnl, V_mu_square, Vcd_pdg, 'leptonic') 
    self.Ds_to_uHNL = Decay(D_sub_s_meson, mu, hnl, V_mu_square, Vcs_pdg, 'leptonic') 
    
    self.B_to_tHNL  = Decay(B_meson,       tau, hnl, V_tau_square, Vub_pdg, 'leptonic') 
    self.Bc_to_tHNL = Decay(B_sub_c_meson, tau, hnl, V_tau_square, Vcb_pdg, 'leptonic')  
    self.D_to_tHNL  = Decay(D_meson,       tau, hnl, V_tau_square, Vcd_pdg, 'leptonic') 
    self.Ds_to_tHNL = Decay(D_sub_s_meson, tau, hnl, V_tau_square, Vcs_pdg, 'leptonic') 
    
    # semileptonic into pseudoscalar meson
    self.B_to_D0eHNL   = Decay(B_meson      , [D0_meson, el]     , hnl, V_el_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.B_to_pi0eHNL  = Decay(B_meson      , [pi0_meson, el]    , hnl, V_el_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi') 
    self.B0_to_pieHNL  = Decay(B0_meson     , [pi_meson, el]     , hnl, V_el_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi') 
    self.B0_to_DeHNL   = Decay(B0_meson     , [D_meson, el]      , hnl, V_el_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.Bs_to_KeHNL   = Decay(B_sub_s_meson, [K_meson, el]      , hnl, V_el_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='Bs_to_K') 
    self.Bs_to_DseHNL  = Decay(B_sub_s_meson, [D_sub_s_meson, el], hnl, V_el_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.D_to_K0eHNL   = Decay(D_meson      , [K0_meson, el]     , hnl, V_el_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ) 
    self.D_to_pi0eHNL  = Decay(D_meson      , [pi0_meson, el]    , hnl, V_el_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi') 
    self.D0_to_pieHNL  = Decay(D0_meson     , [pi_meson, el]     , hnl, V_el_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi')
    self.D0_to_KeHNL   = Decay(D0_meson     , [K_meson, el]      , hnl, V_el_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ) 

    self.B_to_D0uHNL   = Decay(B_meson      , [D0_meson, mu]     , hnl, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.B_to_pi0uHNL  = Decay(B_meson      , [pi0_meson, mu]    , hnl, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi') 
    self.B0_to_piuHNL  = Decay(B0_meson     , [pi_meson, mu]     , hnl, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi') 
    self.B0_to_DuHNL   = Decay(B0_meson     , [D_meson, mu]      , hnl, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.Bs_to_KuHNL   = Decay(B_sub_s_meson, [K_meson, mu]      , hnl, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='Bs_to_K') 
    self.Bs_to_DsuHNL  = Decay(B_sub_s_meson, [D_sub_s_meson, mu], hnl, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.D_to_K0uHNL   = Decay(D_meson      , [K0_meson, mu]     , hnl, V_mu_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ) 
    self.D_to_pi0uHNL  = Decay(D_meson      , [pi0_meson, mu]    , hnl, V_mu_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi') 
    self.D0_to_piuHNL  = Decay(D0_meson     , [pi_meson, mu]     , hnl, V_mu_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi')
    self.D0_to_KuHNL   = Decay(D0_meson     , [K_meson, mu]      , hnl, V_mu_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ) 

    self.B_to_D0tHNL   = Decay(B_meson      , [D0_meson, tau]     , hnl, V_tau_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.B_to_pi0tHNL  = Decay(B_meson      , [pi0_meson, tau]    , hnl, V_tau_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi')
    self.B0_to_pitHNL  = Decay(B0_meson     , [pi_meson, tau]     , hnl, V_tau_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi')
    self.B0_to_DtHNL   = Decay(B0_meson     , [D_meson, tau]      , hnl, V_tau_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ) 
    self.Bs_to_KtHNL   = Decay(B_sub_s_meson, [K_meson, tau]      , hnl, V_tau_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='Bs_to_K') 
    self.Bs_to_DstHNL  = Decay(B_sub_s_meson, [D_sub_s_meson, tau], hnl, V_tau_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' )
    
    
    # semileptonic into vector meson
    self.B_to_rho0eHNL    = Decay(B_meson      , [rho0_meson, el]       , hnl, V_el_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ) 
    self.B_to_D0stareHNL  = Decay(B_meson      , [D0star_meson, el]     , hnl, V_el_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  )
    self.B0_to_DstareHNL  = Decay(B0_meson     , [Dstar_meson, el]      , hnl, V_el_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  )
    self.B0_to_rhoeHNL    = Decay(B0_meson     , [rho_meson, el]        , hnl, V_el_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    )
    self.Bs_to_DsstareHNL = Decay(B_sub_s_meson, [D_sub_sstar_meson, el], hnl, V_el_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Dsstar')  
    self.Bs_to_KstareHNL  = Decay(B_sub_s_meson, [Kstar_meson, el]      , hnl, V_el_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Kstar' ) 
    self.D_to_K0stareHNL  = Decay(D_meson      , [K0star_meson, el]     , hnl, V_el_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  ) 
    self.D0_to_KstareHNL  = Decay(D0_meson     , [Kstar_meson, el]      , hnl, V_el_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  )
    
    self.B_to_rho0uHNL    = Decay(B_meson      , [rho0_meson, mu]       , hnl, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ) 
    self.B_to_D0staruHNL  = Decay(B_meson      , [D0star_meson, mu]     , hnl, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ) 
    self.B0_to_DstaruHNL  = Decay(B0_meson     , [Dstar_meson, mu]      , hnl, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ) 
    self.B0_to_rhouHNL    = Decay(B0_meson     , [rho_meson, mu]        , hnl, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ) 
    self.Bs_to_DsstaruHNL = Decay(B_sub_s_meson, [D_sub_sstar_meson, mu], hnl, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Dsstar')  
    self.Bs_to_KstaruHNL  = Decay(B_sub_s_meson, [Kstar_meson, mu]      , hnl, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Kstar' ) 
    self.D_to_K0staruHNL  = Decay(D_meson      , [K0star_meson, mu]     , hnl, V_mu_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  )
    self.D0_to_KstaruHNL  = Decay(D0_meson     , [Kstar_meson, mu]      , hnl, V_mu_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  )

    self.B_to_rho0tHNL    = Decay(B_meson      , [rho0_meson, tau]       , hnl, V_tau_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    )
    self.B_to_D0startHNL  = Decay(B_meson      , [D0star_meson, tau]     , hnl, V_tau_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ) 
    self.B0_to_DstartHNL  = Decay(B0_meson     , [Dstar_meson, tau]      , hnl, V_tau_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ) 
    self.B0_to_rhotHNL    = Decay(B0_meson     , [rho_meson, tau]        , hnl, V_tau_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ) 
    
    self.Bs_to_KstartHNL  = Decay(B_sub_s_meson, [Kstar_meson, tau]      , hnl, V_tau_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Kstar' )
  
    self.BR_tot_mu = (B_meson.fraction*(self.B_to_uHNL.BR + self.B_to_D0uHNL.BR+self.B_to_pi0uHNL.BR+self.B_to_rho0uHNL.BR+self.B_to_D0staruHNL.BR )+
	             B0_meson.fraction*(self.B0_to_piuHNL.BR + self.B0_to_DuHNL.BR+self.B0_to_DstaruHNL.BR+self.B0_to_rhouHNL.BR )+			   
	      B_sub_s_meson.fraction*(self.Bs_to_KuHNL.BR + self.Bs_to_DsuHNL.BR+self.Bs_to_DsstaruHNL.BR+self.Bs_to_KstaruHNL.BR ))		
   
    self.BR_tot_e = B_meson.fraction       * (self.B_to_eHNL.BR    + self.B_to_D0eHNL.BR  + self.B_to_pi0eHNL.BR     + self.B_to_rho0eHNL.BR   + self.B_to_D0stareHNL.BR) \
                   + B0_meson.fraction      * (self.B0_to_pieHNL.BR + self.B0_to_DeHNL.BR  + self.B0_to_DstareHNL.BR  + self.B0_to_rhoeHNL.BR                            ) \
                   + B_sub_s_meson.fraction * (self.Bs_to_KeHNL.BR  + self.Bs_to_DseHNL.BR + self.Bs_to_DsstareHNL.BR + self.Bs_to_KstareHNL.BR                          )

    self.BR_B_e =  (self.B_to_eHNL.BR    + self.B_to_D0eHNL.BR  + self.B_to_pi0eHNL.BR     + self.B_to_rho0eHNL.BR   + self.B_to_D0stareHNL.BR) 
    self.BR_B0_e =   (self.B0_to_pieHNL.BR + self.B0_to_DeHNL.BR  + self.B0_to_DstareHNL.BR  + self.B0_to_rhoeHNL.BR ) 
    self.BR_Bs_e =   (self.Bs_to_KeHNL.BR  + self.Bs_to_DseHNL.BR + self.Bs_to_DsstareHNL.BR + self.Bs_to_KstareHNL.BR )

  
    self.BR_B_mu = (self.B_to_uHNL.BR    + self.B_to_D0uHNL.BR  + self.B_to_pi0uHNL.BR     + self.B_to_rho0uHNL.BR   + self.B_to_D0staruHNL.BR) 
    self.BR_B0_mu =   (self.B0_to_piuHNL.BR + self.B0_to_DuHNL.BR  + self.B0_to_DstaruHNL.BR  + self.B0_to_rhouHNL.BR ) 
    self.BR_Bs_mu =  (self.Bs_to_KuHNL.BR  + self.Bs_to_DsuHNL.BR + self.Bs_to_DsstaruHNL.BR + self.Bs_to_KstaruHNL.BR )
 
    self.BR_tot_mu_Bc = B_sub_c_meson.fraction * self.Bc_to_uHNL.BR
    self.BR_tot_e_Bc = B_sub_c_meson.fraction * self.Bc_to_eHNL.BR
   
#    self.BR_tot_mu_Bc = self.Bc_to_uHNL.BR
 #   self.BR_tot_e_Bc =  self.Bc_to_eHNL.BR
## HNL DECAYS ##
class HNLDecays(object):
  def __init__(self, mass, mixing_angle_square): 
    self.mass = mass
    self.mixing_angle_square = mixing_angle_square

    # define the HNL
    hnl = Particle('hnl', 'lepton', self.mass)
    
    # get the model
    V_mu_square =  self.mixing_angle_square*1
    V_tau_square = self.mixing_angle_square*0#self.mixing_angle_square # uncomment for figure 13 right
    V_el_square =  self.mixing_angle_square*0#self.mixing_angle_square # uncomment for figure 13 right 
  #  V_tau_square = 0.#self.mixing_angle_square # uncomment for figure 13 right
   # V_el_square =  0.#self.mixing_angle_square # uncomment for figure 13 right 
    QCD_corr = QCD_corr_table[mass] if not PESKIN else 0.
    #QCD_corr = 0.18 if not PESKIN else 0.
    special_V_mu_square = V_mu_square if not PESKIN else 0.

    # list of the decays of interest
    decay_rates = {}
    self.decay_rate = {}

    self.decay_rate['mupi'] = HNLDecay(hnl, [mu,pi_meson], V_mu_square, Vud_pdg, 'mupi').decay_rate
    self.decay_rate['epi'] = HNLDecay(hnl, [el,pi_meson], V_mu_square, Vud_pdg, 'epi').decay_rate
    self.decay_rate['taupi'] = HNLDecay(hnl, [tau,pi_meson], V_mu_square, Vud_pdg, 'epi').decay_rate
                            
    decay_rates['cc_lep'] = [ 
                              HNLDecay(hnl, [mu,el,nu_el],   V_mu_square, 1, 'cc_lep').decay_rate, 
                              HNLDecay(hnl, [mu,tau,nu_tau], special_V_mu_square, 1, 'cc_lep').decay_rate, 
                              # Vel
                              HNLDecay(hnl, [el,mu,nu_mu],   V_el_square, 1, 'cc_lep').decay_rate, 
                              HNLDecay(hnl, [el,tau,nu_tau], V_el_square, 1, 'cc_lep').decay_rate, 
                              # Vtau
                              HNLDecay(hnl, [tau,el,nu_el],  V_tau_square,1, 'cc_lep').decay_rate, 
                              HNLDecay(hnl, [tau,mu,nu_mu],  V_tau_square,1, 'cc_lep').decay_rate, 

                              # these contributions should already be included in nc_lep
                              #HNLDecay(hnl, [mu,mu,nu_mu],   V_mu_square, 1, 'cc_lep').decay_rate,
                              #HNLDecay(hnl, [el,el,nu_el],   V_el_square, 1, 'cc_lep').decay_rate,
                              #HNLDecay(hnl, [tau,tau,nu_tau],V_tau_square, 1, 'cc_lep').decay_rate, 
                              
                            ]
    self.decay_rate['cc_lep'] = sum(decay_rates['cc_lep'])

    decay_rates['cc_had'] = [ 
                              HNLDecay(hnl, [mu,uq,dq], V_mu_square, Vud_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [mu,uq,sq], V_mu_square, Vus_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [mu,uq,bq], V_mu_square, Vub_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [mu,cq,dq], V_mu_square, Vcd_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [mu,cq,sq], V_mu_square, Vcs_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [mu,cq,bq], V_mu_square, Vcb_pdg, 'cc_had').decay_rate,
                              # decays to top not possible for hNL mass < 6 GeV
                              # Vel
                              HNLDecay(hnl, [el,uq,dq], V_el_square, Vud_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [el,uq,sq], V_el_square, Vus_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [el,uq,bq], V_el_square, Vub_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [el,cq,dq], V_el_square, Vcd_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [el,cq,sq], V_el_square, Vcs_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [el,cq,bq], V_el_square, Vcb_pdg, 'cc_had').decay_rate,
                              # Vtau
                              HNLDecay(hnl, [tau,uq,dq], V_tau_square, Vud_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [tau,uq,sq], V_tau_square, Vus_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [tau,uq,bq], V_tau_square, Vub_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [tau,cq,dq], V_tau_square, Vcd_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [tau,cq,sq], V_tau_square, Vcs_pdg, 'cc_had').decay_rate,
                              HNLDecay(hnl, [tau,cq,bq], V_tau_square, Vcb_pdg, 'cc_had').decay_rate,
                            ]  
    self.decay_rate['cc_had'] = sum(decay_rates['cc_had']) * (1 + QCD_corr)  # (three-loop correction for quark->hadrons, with alphas_s(m_tau=1.8 GeV) )    

    decay_rates['nc_lep'] = [ 
                              HNLDecay(hnl, [nu_mu,el,el],   V_mu_square, 1, 'nc_lep').decay_rate,
                              HNLDecay(hnl, [nu_mu,mu,mu],   V_mu_square, 1, 'nc_lep').decay_rate,
                              HNLDecay(hnl, [nu_mu,tau,tau], special_V_mu_square, 1, 'nc_lep').decay_rate,
                              # Vel
                              HNLDecay(hnl, [nu_el,el,el],   V_el_square, 1, 'nc_lep').decay_rate,
                              HNLDecay(hnl, [nu_el,mu,mu],   V_el_square, 1, 'nc_lep').decay_rate,
                              HNLDecay(hnl, [nu_el,tau,tau], V_el_square, 1, 'nc_lep').decay_rate,
                              # Vtau
                              HNLDecay(hnl, [nu_tau,el,el],   V_tau_square, 1, 'nc_lep').decay_rate,
                              HNLDecay(hnl, [nu_tau,mu,mu],   V_tau_square, 1, 'nc_lep').decay_rate,
                              HNLDecay(hnl, [nu_tau,tau,tau], V_tau_square, 1, 'nc_lep').decay_rate,
                            ]
    self.decay_rate['nc_lep'] = sum(decay_rates['nc_lep'])

    decay_rates['nc_had'] = [ 
                              HNLDecay(hnl, [nu_mu,uq,uq],      V_mu_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_mu,dq,dq],      V_mu_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_mu,cq,cq],      V_mu_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_mu,sq,sq],      V_mu_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_mu,bq,bq],      V_mu_square, 1, 'nc_had').decay_rate, 
                              # decays to bbbar and ttbar not possible for HNL mass < 6 GeV
                              # Vel
                              HNLDecay(hnl, [nu_el,uq,uq],      V_el_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_el,dq,dq],      V_el_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_el,cq,cq],      V_el_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_el,sq,sq],      V_el_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_el,bq,bq],      V_el_square, 1, 'nc_had').decay_rate, 
                              # Vtau
                              HNLDecay(hnl, [nu_tau,uq,uq],     V_tau_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_tau,dq,dq],     V_tau_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_tau,cq,cq],     V_tau_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_tau,sq,sq],     V_tau_square, 1, 'nc_had').decay_rate,
                              HNLDecay(hnl, [nu_tau,bq,bq],     V_tau_square, 1, 'nc_had').decay_rate, 

                            ]
    self.decay_rate['nc_had'] = sum(decay_rates['nc_had']) * (1 + QCD_corr)  # (three-loop correction for quark->hadrons, with alphas_s(m_tau=1.8 GeV) )    

    decay_rates['nc_neu'] = [
                              HNLDecay(hnl, [nu_mu,nu_el,nu_el],   V_mu_square, 1, 'nc_neu').decay_rate,
                              HNLDecay(hnl, [nu_mu,nu_mu,nu_mu],   V_mu_square, 1, 'nc_neu').decay_rate,
                              HNLDecay(hnl, [nu_mu,nu_tau,nu_tau], V_mu_square, 1, 'nc_neu').decay_rate,
                              # Vel
                              HNLDecay(hnl, [nu_el,nu_el,nu_el],   V_el_square, 1, 'nc_neu').decay_rate,
                              HNLDecay(hnl, [nu_el,nu_mu,nu_mu],   V_el_square, 1, 'nc_neu').decay_rate,
                              HNLDecay(hnl, [nu_el,nu_tau,nu_tau], V_el_square, 1, 'nc_neu').decay_rate,
                              # Vtau
                              HNLDecay(hnl, [nu_tau,nu_el,nu_el],   V_tau_square, 1, 'nc_neu').decay_rate,
                              HNLDecay(hnl, [nu_tau,nu_mu,nu_mu],   V_tau_square, 1, 'nc_neu').decay_rate,
                              HNLDecay(hnl, [nu_tau,nu_tau,nu_tau], V_tau_square, 1, 'nc_neu').decay_rate,

                            ]
    self.decay_rate['nc_neu'] = sum(decay_rates['nc_neu'])

    # partial sums
    self.decay_rate['tot_lep'] = self.decay_rate['cc_lep'] + self.decay_rate['nc_lep']
    self.decay_rate['tot_neu'] = self.decay_rate['nc_neu']
    self.decay_rate['tot_had'] = self.decay_rate['cc_had'] + self.decay_rate['nc_had']

    self.decay_rate['tot'] = self.decay_rate['tot_lep'] + self.decay_rate['tot_neu'] + self.decay_rate['tot_had']
    #self.ctau = ctau_from_gamma(gamma=self.decay_rate)

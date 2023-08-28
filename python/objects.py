'''
This srcripts collects the necessary classes and functions for the computation of the 
HNL production rates 

'''

from scipy.integrate import quad as integrate
import math 

const_GF =  1.1663787e-05
const_pi = math.pi
const_sinThW_square =  0.23121 # from https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf 

global debug
debug = False
  
class Particle(object):
  def __init__(self, name, particle_type, mass=None, decay_constant=None, lifetime=None, fraction=None):
    self.name           = name
    self.particle_type  = particle_type
    self.mass           = mass
    self.decay_constant = decay_constant
    self.lifetime       = lifetime
    self.fraction       = fraction
    
    if self.particle_type not in ['meson', 'lepton', 'neutrino', 'quark']: 
      raise RuntimeError('You entered an unkown type of particles, please check')

class FormFactor(object):
  def __init__(self, label, chi, m1, m2):
    self.label = label
    self.chi   = chi
    self.m1    = m1
    self.m2    = m2

    # different types of decays require different expressions of the form factors
    if self.label in ['B_to_D', 'B_to_pi', 'Bs_to_K', 'D_to_K', 'D_to_pi']: # complete! 
      self.decay_type = 'semileptonic_pseudoscalar'

    elif self.label in ['B_to_rho', 'B_to_Dstar', 'Bs_to_Dsstar', 'Bs_to_Kstar', 'D_to_Kstar']: # complete!
      self.decay_type = 'semileptonic_vector'

    else: raise RuntimeError('The requested form factor has not been registered. Please check')

    if self.decay_type == 'semileptonic_pseudoscalar':
      if self.label == 'B_to_D':
        Mpole_plus = -99    # case of infinite mass pole
        a0_plus = 0.909
        a1_plus = -7.11
        a2_plus = 66
        Mpole_null = -99    # case of infinite mass pole
        a0_null = 0.794
        a1_null = -2.45
        a2_null = 33

      elif self.label == 'B_to_pi':
        Mpole_plus = 5.325
        a0_plus = 0.404
        a1_plus = -0.68
        a2_plus = -0.86 
        Mpole_null = -99 # 5.65
        a0_null = 0.49
        a1_null = -1.61
        a2_null = 0.93

      elif self.label == 'Bs_to_K':
        Mpole_plus = 5.325
        a0_plus = 0.36
        a1_plus = -0.828
        a2_plus = 1.1
        Mpole_null = 5.65
        a0_null = 0.233
        a1_null = 0.197
        a2_null = 0.18

      elif self.label == 'D_to_K':
        f0_plus = 0.7647
        c_plus  = 0.066
        P_plus  = 0.224
        f0_null = 0.7647
        c_null  = 2.084
        P_null  = 0.

      elif self.label == 'D_to_pi':
        f0_plus = 0.6117
        c_plus  = 1.985
        P_plus  = 0.1314
        f0_null = 0.6117
        c_null  = 1.188
        P_null  = 0.0342

      else: raise RuntimeError('It seems that the constants for the form factor computation have not been given. Please check')

      tplus  = (self.m1 + self.m2) * (self.m1 + self.m2)
      t0     = (self.m1 + self.m2) * (math.pow(self.m1, 0.5) - math.pow(self.m2, 0.5)) * (math.pow(self.m1, 0.5) - math.pow(self.m2, 0.5))
      z      = (math.sqrt(tplus - self.chi) - math.sqrt(tplus - t0)) / (math.sqrt(tplus - self.chi) + math.sqrt(tplus - t0))
      z0     = (math.sqrt(tplus) - math.sqrt(tplus - t0)) / (math.sqrt(tplus) + math.sqrt(tplus - t0))

      if 'B_' in self.label or 'Bs_' in self.label:
        prefactor_plus = 1.0 if Mpole_plus == -99 else 1 / (1 - self.chi / (Mpole_plus * Mpole_plus))
        prefactor_null = 1.0 if Mpole_null == -99 else 1 / (1 - self.chi / (Mpole_null * Mpole_null))

        self.expr_fplus = prefactor_plus * (a0_plus + a1_plus * (z - (1/3.) * math.pow(z, 3)) + a2_plus * (z*z + (2/3.) * math.pow(z, 3)))
        self.expr_fnull = prefactor_null * (a0_null + a1_null * (z - (1/3.) * math.pow(z, 3)) + a2_null * (z*z + (2/3.) * math.pow(z, 3)))

      elif 'D_' in self.label:
        self.expr_fplus = (f0_plus - (c_plus * (z - z0) * (1 + ((z + z0) / 2)))) / (1 - P_plus * self.chi)
        self.expr_fnull = (f0_null - (c_null * (z - z0) * (1 + ((z + z0) / 2)))) / (1 - P_null * self.chi)

      else: raise RuntimeError('It seems that no expression for the form factor of that transition has been given. Please check')
    
    elif self.decay_type == 'semileptonic_vector':
      if self.label == 'B_to_rho':
        fVhhp = 0.295
        fA0hhp = 0.231
        fA1hhp = 0.269
        fA2hhp = 0.282
        sigmaVhhp = 0.875
        sigmaA0hhp = 0.796
        sigmaA1hhp = 0.54
        sigmaA2hhp = 1.34
        chiVhhp = 0
        chiA0hhp = 0.055
        chiA1hhp = 0
        chiA2hhp = -0.21
        MPh = 5.279
        MVh = 5.325

      elif self.label == 'B_to_Dstar':
        fVhhp = 0.76 
        fA0hhp = 0.69
        fA1hhp = 0.66
        fA2hhp = 0.62
        sigmaVhhp = 0.57
        sigmaA0hhp = 0.59
        sigmaA1hhp = 0.78
        sigmaA2hhp = 1.4
        chiVhhp = 0
        chiA0hhp = 0
        chiA1hhp = 0
        chiA2hhp = 0.41
        MPh = 6.275
        MVh = 6.331

      elif self.label == 'Bs_to_Dsstar':
        fVhhp = 0.95  
        fA0hhp = 0.67
        fA1hhp = 0.7
        fA2hhp = 0.75
        sigmaVhhp = 0.372
        sigmaA0hhp = 0.350
        sigmaA1hhp = 0.463
        sigmaA2hhp = 1.04
        chiVhhp = 0.561
        chiA0hhp = 0.6
        chiA1hhp = 0.51
        chiA2hhp = 0.07
        MPh = 6.275
        MVh = 6.331
        
      elif self.label == 'Bs_to_Kstar':
        fVhhp = 0.291   
        fA0hhp = 0.289
        fA1hhp = 0.287
        fA2hhp = 0.286
        sigmaVhhp = -0.516
        sigmaA0hhp = -0.383
        sigmaA1hhp = 0.
        sigmaA2hhp = 1.05
        chiVhhp = 2.1
        chiA0hhp = 1.58
        chiA1hhp = 1.06
        chiA2hhp = -0.074
        MPh = 5.367
        MVh = 5.415

      elif self.label == 'D_to_Kstar':
        fVhhp = 1.03    
        fA0hhp = 0.76
        fA1hhp = 0.66
        fA2hhp = 0.49
        sigmaVhhp = 0.27
        sigmaA0hhp = 0.17
        sigmaA1hhp = 0.3
        sigmaA2hhp = 0.67
        chiVhhp = 0.
        chiA0hhp = 0.
        chiA1hhp = 0.2
        chiA2hhp = 0.16
        MPh = 1.969
        MVh = 2.112

      else: raise RuntimeError('It seems that the constants for the form factor computation have not been given. Please check')

      ratioV = chi /(MVh**2)
      ratioP = chi / (MPh**2)

      # form factor g
      Vhhp_num  = fVhhp
      Vhhp_deno = (1-ratioV) * (1 - sigmaVhhp * ratioV - chiVhhp * ratioV**2)
      Vhhp = Vhhp_num / Vhhp_deno
      self.expr_g = Vhhp / (m1 + m2)

      # form factor f
      A1hhp_num  = fA1hhp
      A1hhp_deno = 1 - sigmaA1hhp * ratioV - chiA1hhp * ratioV**2
      A1hhp = A1hhp_num / A1hhp_deno
      self.expr_f = A1hhp * (m1 + m2)

      # form factor a+
      A2hhp_num  = fA2hhp
      A2hhp_deno = 1 - sigmaA2hhp * ratioV - chiA2hhp * ratioV**2
      A2hhp = A2hhp_num / A2hhp_deno
      self.expr_aplus = (-1) * A2hhp / (m1 + m2)
      
      # form factor a-
      A0hhp_num = fA0hhp
      A0hhp_deno = (1 - ratioP) * (1- sigmaA0hhp * ratioV - chiA0hhp * ratioV**2)
      A0hhp = A0hhp_num / A0hhp_deno
      self.expr_aminus = 1/chi * (2 * m2 * A0hhp - self.expr_f - ((m1**2 - m2**2) * self.expr_aplus)) 



# defining the elements needed for the computation of the semileptonic decay rates (Teixeira et al. (2014) & Bondarenko et al. (2018))
def Lambda(a, b, c):
  return a*a + b*b + c*c -2*a*b -2*a*c -2*b*c


def IntGammas(qsquare, m, m1, m2, m3, decay_type, formFactorLabel):

  if decay_type == 'semileptonic_pseudoscalar':
    formFactor_fplus = FormFactor(formFactorLabel, qsquare, m, m3).expr_fplus 
    formFactor_fnull = FormFactor(formFactorLabel, qsquare, m, m3).expr_fnull 
    deltamsquare = m**2 - m3**2

    integral1 = (1/3.) * formFactor_fplus**2 * math.pow(Lambda(qsquare, m**2, m3**2), 1.5) * math.pow(Lambda(qsquare, m1**2, m2**2), 1.5) * math.pow(qsquare, -3)
    integral2 = 0.5 * formFactor_fplus**2 * math.pow(Lambda(qsquare, m**2, m3**2), 1.5) * math.sqrt(Lambda(qsquare, m1**2, m2**2)) * math.pow(qsquare, -3) * (qsquare*(m1**2 + m2**2) - ((m1**2 - m2**2)**2))
    integral3 = 0.5 * formFactor_fnull**2 * (deltamsquare / qsquare)**2 * math.sqrt(Lambda(qsquare, m**2, m3**2)) * (1/qsquare) * math.sqrt(Lambda(qsquare, m1**2, m2**2)) * (qsquare * (m1**2 + m2**2) - ((m1**2 - m2**2)**2))

    return [integral1, integral2, integral3]
  
  elif decay_type == 'semileptonic_vector':

    formFactor_f = FormFactor(formFactorLabel, qsquare*m*m, m, m3).expr_f 
    formFactor_g = FormFactor(formFactorLabel, qsquare*m*m, m, m3).expr_g
    formFactor_aplus = FormFactor(formFactorLabel, qsquare*m*m, m, m3).expr_aplus
    formFactor_aminus = FormFactor(formFactorLabel, qsquare*m*m, m, m3).expr_aminus

    yl = m1 / m
    yN = m2 / m
    yhprime = m3 / m

    bigLambda = math.sqrt(Lambda(1, yhprime**2, qsquare)) * math.sqrt(Lambda(qsquare, yN**2, yl**2))
    bigF = (1 - qsquare)**2 - 2 * yhprime**2 * (1 + qsquare) + math.pow(yhprime, 4)
    Gminus = qsquare * (yN**2 + yl**2) - (yN**2 - yl**2)**2
    Gplus = qsquare * (yN**2 + yl**2) + (yN**2 - yl**2)**2
    
    integral4  = (1/3.) * m**2 * yhprime**2 * math.pow(qsquare, -2) * formFactor_g**2 * bigLambda * bigF * (2*qsquare**2 - Gplus)
    integral5  = (1/24.) * math.pow(m, -2) * math.pow(qsquare, -3) * formFactor_f**2 * bigLambda * (3 * bigF * (qsquare**2 - (yl**2 - yN**2)**2) - bigLambda**2 + 12 * yhprime**2 * qsquare * (2*qsquare**2 - Gplus))
    integral6  = (1/24.) * m**2 * math.pow(qsquare, -3) * formFactor_aplus**2 * bigLambda * bigF * (bigF * (2*qsquare**2 - Gplus) + 3*Gminus*(1-yhprime**2)**2)
    integral7  = (1/8.) * m**2 * math.pow(qsquare, -1) * formFactor_aminus**2 * bigLambda * bigF * Gminus
    integral8  = (1/12.) * math.pow(qsquare, -3) * formFactor_f * formFactor_aplus * bigLambda * (3*qsquare * bigF * Gminus + (1-qsquare-yhprime**2) * (3*bigF * (qsquare**2 - (yl**2 - yN**2)**2) - bigLambda**2))
    integral9  = (1/4.) * math.pow(qsquare, -2) * formFactor_f * formFactor_aminus * bigLambda * bigF * Gminus
    integral10 = (1/4.) * m**2 * math.pow(qsquare, -2) * formFactor_aplus * formFactor_aminus * bigLambda * bigF * Gminus * (1 - yhprime**2)

    return [integral4, integral5, integral6, integral7, integral8, integral9, integral10]


class Decay(object):
  def __init__(self, mother, daughters, hnl, mixing_angle_square, ckm_coupling, decay_type, formFactorLabel=None):
    self.mother = mother
    self.daughters = daughters
    self.hnl = hnl
    self.mixing_angle_square = mixing_angle_square
    self.ckm_coupling = ckm_coupling
    self.decay_type = decay_type
    self.formFactorLabel = formFactorLabel

    if self.decay_type == 'leptonic':
      lep = self.daughters
      yl = lep.mass / self.mother.mass
      yN = self.hnl.mass / self.mother.mass 

      if self.hnl.mass < (self.mother.mass - lep.mass): 
        expr = const_GF**2 * self.mother.decay_constant**2 * (1/(8*const_pi)) \
                          * math.pow(self.mother.mass, 3) * self.ckm_coupling**2 * self.mixing_angle_square \
                          *(yN**2 + yl**2 - (yN**2 - yl**2)**2) \
                          * math.sqrt(Lambda(1, yN**2, yl**2))
        self.decay_rate = expr if expr > 0 else 0
      else:
        #print 'outside kinematic range'
        self.decay_rate = 0

    elif self.decay_type in ['semileptonic_pseudoscalar', 'semileptonic_vector']:                    
      
      # fetch the daughters
      if self.daughters[0].particle_type == 'meson' and self.daughters[1].particle_type != 'lepton' \
        or self.daughters[0].particle_type == 'lepton' and self.daughters[1].particle_type != 'meson':
        raise RuntimeError('Final state not valid, please check')

      daughter_meson = self.daughters[0] if self.daughters[0].particle_type == 'meson' else self.daughters[1]
      daughter_lepton = self.daughters[0] if self.daughters[0].particle_type == 'lepton' else self.daughters[1]

      # definition 
      m  = self.mother.mass
      m1 = daughter_lepton.mass
      m2 = self.hnl.mass
      m3 = daughter_meson.mass

      if self.decay_type == 'semileptonic_pseudoscalar':

        if m2 < (m - m1 - m3): 
          # we define the different components of the expression according to Teixeira et al. (2014)
          # integral boundaries:
          lower_bound = (m1+m2)**2
          upper_bound = (m-m3)**2
          
          # compute the integrals 
          int1 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[0], lower_bound, upper_bound)[0]
          int2 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[1], lower_bound, upper_bound)[0]
          int3 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[2], lower_bound, upper_bound)[0]
          
          # clebsh-gordan coefficient
          CK = 1 / math.sqrt(2) if daughter_meson.name == 'pi0_meson' else 1
          
          # compute the decay rate
          expr = const_GF**2 * CK**2 * self.ckm_coupling**2 * self.mixing_angle_square / (64 * math.pow(const_pi, 3) * math.pow(m, 3)) * (int1 + int2 + int3)
          self.decay_rate = expr if expr > 0 else 0

        else:
          #print 'outside kinematic range'
          self.decay_rate = 0


      elif self.decay_type == 'semileptonic_vector':
        
        if m2 < (m - m1 - m3): 
          # definition
          yl = m1 / m
          yN = m2 / m
          yhprime = m3 / m

          # integral boundaries:
          lower_bound = (yl + yN)**2
          upper_bound = (1 - yhprime)**2

          # compute the integrals 
          int1 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[0], lower_bound, upper_bound)[0]
          int2 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[1], lower_bound, upper_bound)[0]
          int3 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[2], lower_bound, upper_bound)[0]
          int4 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[3], lower_bound, upper_bound)[0]
          int5 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[4], lower_bound, upper_bound)[0]
          int6 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[5], lower_bound, upper_bound)[0]
          int7 = integrate(lambda x: IntGammas(x, m, m1, m2, m3, self.decay_type, self.formFactorLabel)[6], lower_bound, upper_bound)[0]

          # clebsh-gordan coefficient
          CK = 1 / math.sqrt(2) if daughter_meson.name == 'rho0_meson' else 1
          
          # compute the decay rate
          expr = const_GF**2 * math.pow(m, 7) * CK**2 * self.ckm_coupling**2 * self.mixing_angle_square / (64 * math.pow(const_pi, 3) * m3**2) * (int1 + int2 + int3 + int4 + int5 + int6 + int7)
          self.decay_rate = expr if expr > 0 else 0
        
        else:
          #print 'outside kinematic range'
          self.decay_rate = 0
        
    else: raise RuntimeError("You have entered an unkown decay type. Please choose among ['leptonic', 'semileptonic_pseudoscalar', 'semileptonic_vector']")

    # compute the branching ratio
    self.BR = self.decay_rate * self.mother.lifetime


#def loopCorr_quarkToHadrons(scale):
#  alphaS_over_pi = 
#  return (1 + alphaS_over_pi + 5.2 * alphaS_over_pi**2 + 26.4 * alphaS_over_pi**3 )
def L(x):
  num = 1. - 3 * x**2 - (1. - x**2) * math.sqrt(1 - 4. * x**2)
  den = x**2 * (1. + math.sqrt(1 - 4. * x**2) ) 
  #print('x={}, num={}, den={}, num/den={}'.format(x,num,den,num/den))
  #print('l(x) = {}'.format(math.log( num / den )))
  if num==0 or den==0 or x<0.005: return 0
  else: return math.log( num / den )


class HNLDecay(object):
  def __init__(self, hnl, daughters, mixing_angle_square, ckm_coupling, decay_type):
    self.hnl = hnl
    self.daughters = daughters
    self.mixing_angle_square = mixing_angle_square
    self.ckm_coupling = ckm_coupling
    self.decay_type = decay_type

    if debug: 
      print('\n******************************************')
      print('Decaytype={}, daughters={}'.format(self.decay_type,[x.name for x in self.daughters]))

    allowed_decay_types = ['cc_lep','cc_had', 'nc_lep', 'nc_had', 'nc_neu', 'mupi',"epi"]
    if self.decay_type not in allowed_decay_types: raise RuntimeError('Decay type %s not allowed, please check' % self.decay_type)

    # the decay for the analysis
    if self.decay_type == 'mupi' or self.decay_type == 'epi':
      if len(self.daughters)!=2: raise RuntimeError('Expected exactly two particles for this decay')
      if self.daughters[0].particle_type != 'lepton': raise RuntimeError('Expected a lepton as first daughter')
      if self.daughters[1].particle_type != 'meson' : raise RuntimeError('Expected a meson as second daughter')
      mass_check = ( (self.daughters[0].mass + self.daughters[1].mass) < self.hnl.mass )
      if not mass_check:
        self.decay_rate = 0.
      else:
        xl = self.daughters[0].mass / self.hnl.mass
        xh = self.daughters[1].mass / self.hnl.mass
        fh = 130.2 * 1e-03 # in GeV, table 8 
        if debug: print('xl={}, xh={}'.format(xl,xh))
        self.decay_rate = const_GF**2 * fh**2 * self.ckm_coupling**2 * self.mixing_angle_square * self.hnl.mass**3 / (16. * const_pi) \
                          * ( (1-xl**2)**2 - xh**2*(1+xl**2) ) * math.sqrt(Lambda(1,xh**2,xl**2))
    # all other decays
    else: 

      # mass check
      if len(self.daughters)!=3: raise RuntimeError('Expected exactly three particles for this decay')
      mass_check = ( (self.daughters[0].mass + self.daughters[1].mass + self.daughters[2].mass) < self.hnl.mass  )
  
      if not mass_check:
        self.decay_rate = 0.
        if debug:
          print('mass check not passed')
          print('daugters mass = {} , hnl mass = {}'.format((self.daughters[0].mass + self.daughters[1].mass + self.daughters[2].mass), self.hnl.mass))
          print(self.daughters[0].__dict__)
          print(self.daughters[1].__dict__)
          print(self.daughters[2].__dict__)
      else:
        ### charged-current leptonic and quarks
        if self.decay_type == 'cc_lep' or self.decay_type == 'cc_had':
          
          lep = self.daughters[0]
          U = self.daughters[1]
          D = self.daughters[2]
  
          if lep.particle_type != 'lepton': raise RuntimeError('For this decay first daughter must be a charged lepton, please check')
          # other checks
  
          if self.decay_type == 'cc_lep':
            NW = 1.
          else: 
            NW = 3. * self.ckm_coupling**2
  
          xl =    float(lep.mass) / float(self.hnl.mass)
          xU =    float(U.mass)   / float(self.hnl.mass)
          xD =    float(D.mass    / float(self.hnl.mass))
  
          # integration
          upper_bound = (1-xU)**2
          lower_bound = (xD+xl)**2
          integrand = lambda x : 12. * (x-xl**2-xD**2) / x * (1 + xU**2 - x) * math.sqrt( Lambda(x,xl**2,xD**2) * Lambda(1,x,xU**2) )
          (integral,error) = integrate(integrand, lower_bound, upper_bound)
          if error/integral > 0.1: print('Warning: error on integral higher than 10\% ') 
  
          self.decay_rate =  NW * const_GF**2 * self.hnl.mass**5 * self.mixing_angle_square / (192 * const_pi**3) * integral
  
        ### neutral-currents (leptons and quarks)
        elif self.decay_type == 'nc_lep' or self.decay_type == 'nc_had':
  
          nu_alpha = self.daughters[0]
          f    = self.daughters[1]
          fbar = self.daughters[2]
  
          if nu_alpha.particle_type != 'neutrino': raise RuntimeError('For this decay first daughter must be a neutrino, please check')
  
          if self.decay_type == 'nc_lep':
            if f.particle_type == 'neutrino' : raise RuntimeError('should not be a neutrino')
            NZ = 1.
            if debug: print('f.name in nu_alpha.name = {}'.format((f.name in nu_alpha.name)))
            if f.name in nu_alpha.name:
             C1f = 0.25 * ( 1 + 4 * const_sinThW_square + 8 * const_sinThW_square**2 )
             C2f = 0.5 * const_sinThW_square * ( 2 * const_sinThW_square + 1 )
            else:
             C1f = 0.25 * ( 1 - 4 * const_sinThW_square + 8 * const_sinThW_square**2 )
             C2f = 0.5 * const_sinThW_square * ( 2 * const_sinThW_square - 1 )
          else:
            NZ = 3.
            if f.name == 'uq' or f.name == 'cq' or f.name == 'tq':
              C1f = 0.25 * ( 1 - 8./3. * const_sinThW_square + 32./9. * const_sinThW_square**2 ) 
              C2f = 1./3. * const_sinThW_square * (4./3. * const_sinThW_square - 1)
            elif f.name == 'dq' or f.name == 'sq' or f.name == 'bq':
              C1f = 0.25 * ( 1 - 4./3. * const_sinThW_square + 8./9. * const_sinThW_square**2 )
              C2f = 1./6. * const_sinThW_square * (2./3. * const_sinThW_square - 1)
  
          x = f.mass / self.hnl.mass
          par1 = (1. - 14. * x**2 - 2 * x**4 - 12. * x**6) * math.sqrt(1. - 4. * x**2) + 12. * x**4 * (x**4 - 1.) * L(x)
          par2 = x**2 * (2. + 10. * x**2 - 12 * x**4) * math.sqrt(1 - 4. * x**2) + 6 * x**4 * (1 - 2. * x**2 + 2 * x**4) * L(x) 
          bigpar = C1f * par1 + 4. * C2f * par2
          if debug: print('C1f={}, par1={}, C2f={}, par2={}'.format(C1f,par1,C2f,par2))
  
          self.decay_rate = NZ * const_GF**2 * self.hnl.mass**5 * self.mixing_angle_square / (192 * const_pi**3) * bigpar
  
        ### neutral-current (all neutrinos)
        elif self.decay_type == 'nc_neu':
          
          nu_alpha = self.daughters[0]
          nu_beta = self.daughters[1]
          nubar_beta = self.daughters[2]
  
          if nu_alpha.particle_type != 'neutrino': raise RuntimeError('For this decay first daughter must be a neutrino, please check')
  
          delta_alpha_beta = 1. if nu_alpha.name == nu_beta.name else 0.
          self.decay_rate = (1. + delta_alpha_beta) * const_GF**2 * self.hnl.mass**5 * self.mixing_angle_square / (768 * const_pi**3) 
  
    if debug: print('decay_rate={}'.format(self.decay_rate))

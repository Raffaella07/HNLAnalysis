from ROOT import TGraphAsymmErrors, TGraph
from ROOT import TH1D, TH2D
from ROOT import gROOT
import ROOT
import argparse
from CMSlumi import CMSlumi
from HNLUtils import *
import numpy as np
import glob
import os
import math
import array
import sys
import scipy as sp
from numpy import ma
#import skimage
#from skimage import measure
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm

class signal(object):

   def __init__(self,tag):

        self.file = tag

        if os.path.exists(self.file) == False:
                print "limit file does not exist, will skip point"
                self.mass = -1
        #extract mass from path 
        splitter = tag.split("/")
        for chunk in splitter:
                if 'Mass' in chunk:
                        m_string = chunk.split('_')
                        m_string = m_string[0].replace('Mass','')
                        mass = float(m_string.replace('p','.'))

        self.mass = mass

class datacard(object):

   def __init__(self,tag):

        self.file = tag
        features = self.file.split('_')
        self.mass = -1
        self.ctau = -1
        self.v2 = -1
        print features
        features
        if 'm' in features:

                self.mass = features[features.index('m') + 1]
        if 'ctau' in features:
                self.ctau = features[features.index('ctau') + 1]
                self.v2 = getVV(float(self.mass.replace('p','.')),float(self.ctau.replace('p','.')),1)

        if 'v2' in features:
                self.v2 = features[features.index('v2') + 1]
                self.ctau = getCtau(float(self.mass.replace('p','.')),float((self.v2.replace('p','.')).replace("em","E-")),1)

        def isMuE(self):
                if 'electron' in self.file:
                        return True
                else:
                        return False

        def isMuMu(self):
                if 'muon' in self.file:
                        return True
                else:
                        return False


class f_couplings(object):

  def __init__(self,f_mu,f_e):

        self.mu = f_mu
        self.e = f_e
        print f_mu, f_e
        self.MuMu_scale = f_mu * f_mu
        self.MuE_scale = f_mu * f_e

  def getMuMuTag(self):

        return '_MuMu{:.2f}'.format(self.mu).replace('.',"p")

  def getMuETag(self):

        return '_MuE{:.2f}'.format(self.e).replace('.',"p")


def prepareOutput(combDir):
        print combDir
        os.mkdir(combDir)
        os.mkdir(combDir+'/CategoriesCombination')
        os.mkdir(combDir+'/ScaledDatacards')
        os.mkdir(combDir+'/MuMu_MuE_Datacards')
        os.mkdir(combDir+'/MuMu_MuE_Limits')
        os.mkdir(combDir+'/plots')

def limSurface(tag,doDirac):

        gROOT.LoadMacro('/cmshome/ratramon/Analysis/test/LimitsPlotter.C')

	filename = "LimitsCombined_param.txt" #.txt format: v2 expected -1sigma -2sigma +1sigma +2sigma observed
	label = "Majorana"

	if doDirac:

		filename = "LimitsCombined_param_dirac.txt"
		label = "quasi-Dirac"

	print '/cmshome/ratramon/Analysis/output/Mass*'+tag+'/LimitsCombined_param.txt'
        lims = glob.glob('/cmshome/ratramon/Analysis/output/Mass*'+tag+'/'+filename)

	massLims =[]
        for lim in lims:
		 print lim	
                 massLims.append(signal(lim))

	masses = []	
	v2 = []
	lims = []
	fig, ax = plt.subplots()
	counter = 0 
	graph = []
        for t in range(0,6):
        	for idx,lim in enumerate(massLims):
			if (lim.mass>5.2):
				continue	
	                if(t==0):
				 ROOT.limPlotter(lim.file,lim.mass,0.00001, 100, 0.0001, 100, tag+"_", "PF",1)

			#do not compute limits in vetoed region
	     
	                if (  lim.mass not in [  1.0,1.02, 1.04, 1.06,1.08,1.1,1.12,1.14,1.16,1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.4,1.42, 1.44, 1.46,1.48, 1.5, 1.53 , 1.56, 1.59, 1.62, 1.65, 1.68, 1.71,1.74, 1.83, 1.89, 1.92, 1.95, 2.0,2.05, 2.1, 2.15,2.2,2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9,2.95,3.0, 3.2, 3.25, 3.3, 3.35,3.40,3.45,3.50,3.55,3.60,3.65]):
	                        continue
#			ef = 1
				
			#prepare limits arrays
	                with open(lim.file) as f:

				#read lines from file - format: v2 expected -1sigma -2sigma +1sigma +2sigma observed
	                        limits  = f.readlines()

	
				#remove 0 limit rows
	                        limits = [line.split() for line in limits]
	                        limits = np.asarray(limits,dtype=np.float32)
				zeroes = np.where(limits[:,1] ==0)
				for i in zeroes:
#					print limits[i,:]
					limits = np.delete(limits,i,0)
				limits = limits[limits[:, 0].argsort()]
#				print limits
				if limits[0,1] ==0:
					continue
	                        band_excl = np.zeros(6)
			
                                if t ==0 or t ==5:
                                         points =np.c_[limits[:,0],limits[:,t+1]] #central expected + observed 
                                elif t==1 or t==3:
                                         points = np.c_[limits[:,0],np.subtract(limits[:,1],limits[:,t+1])] #build central - 1(2)sigma
                                elif t==2 or t==4:
                                         points = np.c_[limits[:,0],np.add(limits[:,1],limits[:,t+1])] #build central+ 1(2)sigma
			
			
				X = np.full(len(points),lim.mass) #masses
				Y = points[:,0] #v2			
				Z = points[:,1] # 95% CL limit
#			print X
#			print Z
#			X,Y,Z = np.meshgrid(X,Y,Z) 

			#append respectively the three arrays for 3D surface plotting
				if (lim.mass ==1.0):
					masses = X
					v2 = Y
					lims = Z
				#	print 'lenght', len(Z),len(Y), len(masses),len(lims),lim.mass
				else:
					masses = np.concatenate([masses,X],axis=0)
					v2 = np.concatenate([v2,Y],axis=0)
					lims = np.hstack((lims,Z))
				#	print 'lenght', len(Z),len(Y), len(masses),len(lims),lim.mass
			#	print 'len mass',len(masses)
	
				print masses.shape,v2.shape,lims.shape #check if shape is consistent - needs to be the same
	
		#3D limits surface 
		masses = masses[lims<10]	
		v2 = v2[lims<10]	
		lims = lims[lims<10]

		#for standard contour =1 extraction build TGraph2D
		surf = ROOT.TGraph2D('surf','surf',len(masses.astype(float)),masses.astype(float),v2.astype(float),lims.astype(float))
	
		#for multiple contour plotting - just a nice crosscheck
		surf_log = ROOT.TGraph2D('surf','surf',len(masses.astype(float)),masses.astype(float),np.log10(v2.astype(float)),np.log10(lims.astype(float)))
			
		
		#extract DIRECTLY FROM TGRAPH (no interpolation on surface) the z = 1 contour
		graph.append(surf.GetContourList(1))
		canva = ROOT.TCanvas('canva')
		surf.Draw('TRI1') #for 3D surface plotting use Delaunay triangles interpolation - http://www.cs.cornell.edu/info/people/chew/Delaunay.html
		
		
		canva.SetTheta(3)
		canva.SetPhi(120)
		canva.SetBottomMargin(0.20)
		canva.SetLeftMargin(0.20)
		canva.SetRightMargin(0.20)
		surf.SetMaximum(10)
		surf.SetMinimum(0.001)
		surf.GetXaxis().SetLimits(0.7,3.6)
		surf.GetXaxis().SetTitle('m_{N}(GeV)')
		surf.GetXaxis().SetTitleOffset(2.15)
		surf.GetYaxis().SetLimits(0.000001,0.01)
		surf.GetYaxis().SetTitleOffset(2.15)
		surf.GetYaxis().SetTitle("#frac{|V_{#mu} V_{e}|^{2}}{|V_{#mu}|^{2}+|V_{e}|^{2}}")
		surf.GetZaxis().SetTitle("95% C.L limit on #mu")
		surf.GetZaxis().SetTitleOffset(1.15)
		ROOT.gStyle.SetOptStat(0)
		canva.SetLogz()
		canva.SetLogy()
		canva.SaveAs('smoothened_surface_'+str(t)+'.png')

		#several contour plotting - just a cross check: N.B. here interpolation is used 
		cc = ROOT.TCanvas('log10')
		surf_log.Draw('CONT1Z') #for 3D surface plotting use Delaunay triangles interpolation - http://www.cs.cornell.edu/info/people/chew/Delaunay.html
		
		
		cc.SetTheta(3)
		cc.SetPhi(120)
		cc.SetBottomMargin(0.20)
		cc.SetLeftMargin(0.20)
		cc.SetRightMargin(0.20)
		surf_log.SetMaximum(2)
		surf_log.SetMinimum(-2)
		surf_log.GetXaxis().SetLimits(0.7,3.6)
		surf_log.GetXaxis().SetTitle('m_{N}(GeV)')
		surf_log.GetXaxis().SetTitleOffset(2.15)
		surf_log.GetYaxis().SetLimits(-7,0)
		surf_log.GetYaxis().SetTitleOffset(2.15)
		surf_log.GetYaxis().SetTitle("log_{10}#frac{|V_{#mu} V_{e}|^{2}}{|V_{#mu}|^{2}+|V_{e}|^{2}}")
		surf_log.GetZaxis().SetTitle("log_{10}95% C.L limit on #mu")
		surf_log.GetZaxis().SetTitleOffset(1.15)
		ROOT.gStyle.SetOptStat(0)
		cc.SaveAs('log_contour_'+str(t)+'.png')
	
	# z = 1 contours gymnastics for easy plotting
	plt.clf()
        ax = plt.axes()
	contour_x = []
	contour_y = []

	#graph index runs on the 'type' of limits: 0 -> expected, 1-> -1sigma, 2-> +1sigma, 3-> -2 sigma, 4-> +2sigma, 5 -> observed
	for ig,g in enumerate(graph):
		x = []
		y = []
		gidx =0
		# checks whether for each graph element there is more than 1 TGraph built and appends all the points to the same np.arrays
		while g.At(gidx):
			print gidx
			if len(x) ==0:
				x = np.ndarray(g.At(gidx).GetN(), 'd', g.At(gidx).GetX())
				y = np.ndarray(g.At(gidx).GetN(), 'd', g.At(gidx).GetY())
				if len(x) != 0:
					print x[-1]
					print y[-1]
					 
					while x[len(x)-1]<3.0 and y[len(y)-1]>0.0001:
						x = np.delete(x,-1)
						y = np.delete(y,-1)
			else:
				add_x = np.ndarray(g.At(gidx).GetN(), 'd', g.At(gidx).GetX()) 
				add_y = np.ndarray(g.At(gidx).GetN(), 'd', g.At(gidx).GetY()) 
				# second crossing values for m<3.0 are a consequence of the TGraph2D construction - no physics meaning --> we cut them out
				#N.B. this might need further refinement for higher masses 
				add_y = add_y[add_x>=3.0]
				add_x = add_x[add_x>=3.0]
				#depending on contour extraction, reverse array ordering for upper part of nose might be needed
			#	if ig ==5:
			#		add_y = add_y[::-1]
			#		add_x = add_x[::-1]
				x = np.concatenate([x,add_x],axis = 0)
				y = np.concatenate([y,add_y], axis =0)
			
			gidx = gidx + 1 
		contour_x.append(x)
		contour_y.append(y)
	
	#contours plotting 
	
	plt.fill(np.append(contour_x[3],contour_x[4][::-1]), np.append(contour_y[3],contour_y[4][::-1]), color='gold', label=r'$\pm 2 \sigma$')
	plt.fill(np.append(contour_x[1],contour_x[2][::-1]), np.append(contour_y[1],contour_y[2][::-1]), color='forestgreen', label=r'$\pm 1 \sigma$')
        plt.plot(contour_x[0], contour_y[0], color='red', label='central expected', linewidth=2)	
        plt.plot(contour_x[5], contour_y[5], color='black', label='observed', linewidth=2)
	print contour_x[5]	
	
	#legend + canvas plotting	
        plt.legend()
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
        plt.yscale('log')
        plt.ylim(0.0000001,1)
        plt.savefig('contour_1.png')

def exclusionPlot(tag,doDirac,coupling_scale,f_mu_ref, f_e_ref, f_tau_ref,f_mu,f_e,f_tau,ch):

        m_mc = [1.0,1.2,1.3,1.4,1.5,1.8,2.0,2.2,2.3,2.4,2.8]#,3.0,3.1,3.2,3.4]
        enhanceFactors = [0.001,0.01,1,1,100,100]
	filename = "LimitsCombined_param.txt"
	label = "Majorana"
	if doDirac:

		filename = "LimitsCombined_param_dirac.txt"
		label = "quasi-Dirac"
	print '/cmshome/ratramon/Analysis/output/Mass*'+tag+'/LimitsCombined_param.txt'
        lims = glob.glob('/cmshome/ratramon/Analysis/output/Mass*'+tag+'/'+filename)
       # lims = glob.glob('/cmshome/ratramon/Analysis/scripts/combination/output/LimitsCombined_mumu_m_*_0.50_0.50_0.00_coup.txt')
	print lims
        massLims = []
	

        gROOT.LoadMacro('/cmshome/ratramon/Analysis/test/LimitsPlotter.C')
        for lim in lims:
		 print lim	
                 massLims.append(signal(lim))
        masses = []
        exclusions = []
        for idx,lim in enumerate(massLims):
		if (lim.mass>5.2):
			continue	
                ROOT.limPlotter(lim.file,lim.mass,0.00001, 100, 0.0001, 100, tag+"_", "PF",1)
                if (  lim.mass not in [  1.0,1.02, 1.04, 1.06, 1.08,1.1,1.12,1.14,1.16,1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.4,1.42, 1.44,1.46,1.48, 1.5, 1.53 , 1.56, 1.59, 1.62, 1.65, 1.68, 1.71,1.74,1.77, 1.83, 1.89, 1.92, 1.95, 1.98, 2.0,2.05, 2.1, 2.15,2.2,2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9,2.95,3.0, 3.2, 3.25, 3.3, 3.35,3.40,3.45,3.50,3.55,3.60,3.80]):
                #if (  lim.mass not in [1.0, 1.3, 1.4, 1.5, 1.95,2.0,2.05, 2.1, 2.15, 2.2, 2.4, 2.45, 2.55, 2.65, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35,3.40,3.50,3.55,3.6,3.70]):
               # if (  lim.mass>3.3  or  lim.mass==1.2 or lim.mass==3.25 or lim.mass==1.62  or (lim.mass>1.79 and lim.mass<1.99) or (lim.mass>3.047 and lim.mass<3.147) or (lim.mass>1.814 and lim.mass<1.914) or (lim.mass>3.636 and lim.mass<3.736)):
                        continue
                print lim.mass
		ef = 1
			
             #   for im, m in enumerate(m_mc):

              #          if lim.mass == m:
               #                 ef = 1# enhanceFactors[im] * 1/np.sqrt(41.6/4.91)
#               print lim.mass, ef
                with open(lim.file) as f:
                        limits  = f.readlines()
	#		if coupling_scale:
	#			ef = 1.0/coupling_scaler(f_mu_ref, f_e_ref, f_tau_ref,f_mu,f_e,f_tau,ch,lim.mass)
                        limits = [line.split() for line in limits]
                        limits = np.asarray(limits,dtype=np.float32)
			limits[limits[:, 0].argsort()]
			print '___________________________________________________',np.where(limits[:,1] ==0)
			zeroes = np.where(limits[:,1] ==0)
			for i in zeroes:
				print limits[i,:]
				limits = np.delete(limits,i,0)
			print limits
			if limits[0,1] ==0:
				continue
                        band_excl = np.zeros(6)
                        m = 0
                	masses.append(lim.mass)
                        for i in range(0,6):
                                if i ==0 or i ==5:
                                         points =np.c_[limits[:,0],limits[:,i+1]]
                                         m,band_excl[i] = interpolator(points,1,ef,m)
                                elif i==1 or i==3:
                                     #    print m
                                         points = np.c_[limits[:,0],np.subtract(limits[:,1],limits[:,i+1])]
                                         m_temp,band_excl[i] = interpolator(points,1,ef,m)
                                elif i==2 or i==4:
                                         points = np.c_[limits[:,0],np.add(limits[:,1],limits[:,i+1])]
                                         m_temp,band_excl[i] = interpolator(points,1,ef,m)


                               # print i+1, band_excl[i]
                        exclusions.append(band_excl)
        exclusions =np.asarray(exclusions, dtype= np.float32)
	print exclusions
        masses =np.asarray(masses, dtype= np.float32)


        arr1inds = masses.argsort()
        masses = masses[arr1inds[::]]
        masses = masses.reshape(len(masses))
      #  print masses.shape
        exclusions[:,0] = exclusions[:,0][arr1inds[::]]
        exclusions[:,1] = exclusions[:,1][arr1inds[::]]
        exclusions[:,2] = exclusions[:,2][arr1inds[::]]
        exclusions[:,3] = exclusions[:,3][arr1inds[::]]
        exclusions[:,4] = exclusions[:,4][arr1inds[::]]
        exclusions[:,5] = exclusions[:,5][arr1inds[::]]


        exclusions[:,1] = [abs(exclusions[i,0]-exclusions[i,1]) for i in range(0,len(exclusions)) ]
        exclusions[:,2] = [abs(exclusions[i,2]-exclusions[i,0]) for i in range(0,len(exclusions)) ]
        exclusions[:,3] = [abs(exclusions[i,0]-exclusions[i,3]) for i in range(0,len(exclusions)) ]
        exclusions[:,4] = [abs(exclusions[i,4]-exclusions[i,0]) for i in range(0,len(exclusions)) ]


	if coupling_scale:
		ex_file = "exclusions"+tag+"mumu_{:.2f}_{:.2f}_{:.2f}.txt".format(f_mu,f_e,f_tau)
	else:
		ex_file = "exclusions"+tag+"mumu_{:.2f}_{:.2f}_{:.2f}.txt".format(f_mu,f_e,f_tau)

        with open(ex_file, "a") as exc:

		for i in range(0,len(masses)):

                	exc.write("%e %e %e %e %e %e \n"%(masses[i],exclusions[i,0],exclusions[i,1],exclusions[i,2],exclusions[i,3],exclusions[i,4]))	

	print masses,exclusions[:,0]
        g_exp = TGraph(len(masses),masses.astype('float'),exclusions[:,0].astype('float'))
        g_obs = TGraph(len(masses),masses.astype('float'),exclusions[:,5].astype('float'))
	#vetoed regions with tBox
       # print g.GetX()[1], g.GetY()[1]
        g_exp.SetLineWidth(2)
        g_exp.SetLineColor(ROOT.kRed)
       # g_exp.SetLineStyle(7)
  #      g_exp.SetLineStyle(7)
        g_exp.SetTitle("expected limit")
        g_obs.SetLineWidth(2)
        g_obs.SetTitle("observed limit")

        g_1sigma = TGraphAsymmErrors(len(masses),masses.astype('float'),exclusions[:,0].astype('float'),np.zeros(len(masses)),np.zeros(len(masses)),exclusions[:,1].astype('float'),exclusions[:,2].astype('float'))
        g_1sigma.SetFillColor(ROOT.kGreen+2)
        g_1sigma.SetTitle("68% expected")

        g_2sigma = TGraphAsymmErrors(len(masses),masses.astype('float'),exclusions[:,0].astype('float'),np.zeros(len(masses)),np.zeros(len(masses)),exclusions[:,3].astype('float'),exclusions[:,4].astype('float'))
        g_2sigma.SetFillColor(ROOT.kOrange-2)
        g_2sigma.SetTitle("95% expected")

        canva = ROOT.TCanvas()
	canva.SetLeftMargin(0.25)
        g_2sigma.Draw("A3")
        g_1sigma.Draw("same3")
        g_exp.Draw("sameL")
        g_obs.Draw("sameL")
         
        g_2sigma.GetYaxis().SetRangeUser(0.000001,0.09)
        g_2sigma.GetXaxis().SetTitle("m_{N} (GeV)")
        g_2sigma.GetYaxis().SetTitle("#frac{|V_{#mu} V_{e}|^{2}}{|V_{#mu}|^{2}+|V_{e}|^{2}}")
        g_2sigma.GetYaxis().SetTitleOffset(1.8)
      #  g_exp.GetYaxis().SetRangeUser(0.000001,0.09)
      #  g_exp.GetXaxis().SetTitle("m_{N} (GeV)")
      #  g_exp.GetYaxis().SetTitle("#frac{|V_{#mu} V_{e}|^{2}}{|V|^{2}}")
        l= canva.BuildLegend(0.6,0.65,0.9,0.85)
	l.SetHeader(label)
	veto_d0 = ROOT.TBox(1.814,0.000001,1.914,0.09)
	veto_d0_white = ROOT.TBox(1.814,0.000001,1.914,0.09)
	veto_jpsi = ROOT.TBox(3.047,0.000001,3.147,0.009)
	veto_jpsi_white = ROOT.TBox(3.05,0.000001,3.15,0.09)
	veto_psi2s = ROOT.TBox(3.65,0.000001,3.75,0.009)
	veto_psi2s_white = ROOT.TBox(3.65,0.000001,3.75,0.09)
	veto_d0.SetFillStyle(3013)
	veto_d0.SetLineWidth(0)
	veto_d0_white.SetFillColor(ROOT.kWhite)
	veto_d0.SetFillColor(ROOT.kGray)
	veto_jpsi_white.SetFillColor(ROOT.kWhite)
	veto_jpsi.SetFillStyle(3013)
	veto_jpsi.SetFillColor(ROOT.kGray)
	veto_d0.SetLineWidth(0)
 	veto_d0_white.Draw('same')
 	veto_d0.Draw('same')
	veto_jpsi_white.Draw('same')
	veto_jpsi.Draw('same')
	veto_psi2s_white.Draw('same')
	veto_psi2s.Draw('same')
        canva.SetLogy()
	canva.SetGrid(1,1)
	canva.Update()
	canva.RedrawAxis()
	l.Draw()

	CMSlumi(canva,
        	iPeriod = 5,
       		iPosX = 11,
   	        lumiText = '',
        	extraText="Preliminary") 


	if coupling_scale:
	
        	canva.SaveAs("ExcludedV_mass"+tag+"_unsmoothed_"+label+"_{:.2f}_{:.2f}_{:.2f}.pdf".format(f_mu,f_e,f_tau))
	else:
        	canva.SaveAs("ExcludedV_mass"+tag+"_unsmoothed_"+label+".pdf")

	print masses,exclusions[:,0]
	print min(exclusions[:,5]), max(exclusions[:,5])
	return g_obs,g_exp,g_1sigma,g_2sigma

def exclusionUnblind(g_exp,g_1sigma,g_2sigma,g_obsi,tag):

        g_exp.SetLineWidth(2)
     #   g_exp.SetLineStyle(7)
        g_exp.SetLineStyle(9)
        g_exp.SetTitle("expected limit, 5.302 fb^{-1} projected to 41.6 fb^{-1}")
        g_obs.SetLineWidth(2)
        g_obs.SetTitle("observed limit 41.6 fb^{-1}")

        g_1sigma.SetFillColor(ROOT.kGreen+2)
        g_1sigma.SetTitle("68% expected")

        g_2sigma.SetFillColor(ROOT.kOrange-2)
        g_2sigma.SetTitle("95% expected")
     
        canva = ROOT.TCanvas()
        g_2sigma.Draw("A3")
        g_1sigma.Draw("same3")
        g_exp.Draw("sameL")
        g_obs.Draw("sameL")
         
        g_2sigma.GetYaxis().SetRangeUser(0.000001,0.009)
        g_2sigma.GetXaxis().SetRangeUser(1,3)
        g_2sigma.GetXaxis().SetTitle("m_{N} (GeV)")
        g_2sigma.GetYaxis().SetTitle("#frac{|V_{#mu} V_{e}|^{2}}{|V_{#mu}|^{2}+|V_{e}|^{2}}")
        canva.SetLogy()
	canva.SetGrid(1,1)
	canva.Update()
        canva.BuildLegend(0.25,0.17,0.7,0.4)

#	CMSlumi(canva,
	sigma3_line=ROOT.TLine(0.75,0.0011,4.0,0.0011) 
	sigma2_line=ROOT.TLine(0.75,0.0202,4.0,0.0202) 
	sigma1_line=ROOT.TLine(0.75,0.1469,4.0,0.1469) 
	sigma3_line.SetLineColor(ROOT.kRed)
	sigma2_line.SetLineColor(ROOT.kRed)
	sigma3_line=ROOT.TLine(0.75,0.0011,4.0,0.0011) 
	sigma2_line=ROOT.TLine(0.75,0.0202,4.0,0.0202) 
	sigma1_line=ROOT.TLine(0.75,0.1469,4.0,0.1469) 
	sigma3_line.SetLineColor(ROOT.kRed)
	sigma2_line.SetLineColor(ROOT.kRed)
	sigma1_line.SetLineColor(ROOT.kRed)
	sigma3_line.SetLineWidth(2)
	sigma2_line.SetLineWidth(2)
	sigma1_line.SetLineWidth(2)
	#veto_d0.SetLineWidth(0)
#	veto_d0_white.SetFillColor(ROOT.kWhite)
#	veto_d0.SetFillColor(ROOT.kGray)
	veto_jpsi_white.SetFillColor(ROOT.kWhite)
	veto_jpsi.SetFillStyle(3013)
	veto_jpsi.SetFillColor(ROOT.kGray)
#	veto_d0.SetLineWidth(0)
 	#veto_d0_white.Draw('same')
 #	veto_d0.Draw('same')
	veto_jpsi_white.Draw('same')
	veto_jpsi.Draw('same')
	sigma3_line.Draw('same')
	sigma2_line.Draw('same')
	sigma1_line.Draw('same')
	canva.Update()
	canva.RedrawAxis()
	CMSlumi(canva,
	 	iPeriod = 5,
  		iPosX = 11,
	        lumiText = '',
	 	extraText="Preliminary") 
		
		
	canva.SaveAs("pval_observed_"+tag+".pdf")


def plot_pval(tag):

        g_pval = TGraph('../test/pValue_vs_mass_07_08_dirac.txt',"%lg %lg")
        g_pval.SetLineColor(ROOT.kBlack)
    	
	ROOT.gStyle.SetLegendBorderSize(0) 
	ROOT.gStyle.SetLegendFillColor(0) 
        g_pval.SetName('Observed p-value') 
        g_pval.SetTitle('Observed p-value') 
        canva = ROOT.TCanvas()
	canva.SetLeftMargin(0.15)
        g_pval.Draw("AL")
#        g_1sigma.Draw("same3")
#        g_exp.Draw("sameL")
#        g_obs.Draw("sameL")
        g_pval.GetYaxis().SetRangeUser(0.0001,5)
        g_pval.GetXaxis().SetLimits(0.75,4.0)
        g_pval.GetXaxis().SetTitle("m_{N} (GeV)")
        g_pval.GetYaxis().SetTitle("#frac{|V_{#mu} V_{e}|^{2}}{|V_{#mu}|^{2}+|V_{e}|^{2}}")
        g_pval.GetYaxis().SetTitleOffset(1.8)
        canva.SetLogy()
	canva.SetGrid(1,1)
	canva.Update()
        canva.BuildLegend(0.25,0.1,0.5,0.3)

#	CMSlumi(canva,
	sigma3_line=ROOT.TLine(0.75,0.0011,4.0,0.0011) 
	sigma2_line=ROOT.TLine(0.75,0.0202,4.0,0.0202) 
	sigma1_line=ROOT.TLine(0.75,0.1469,4.0,0.1469) 
	sigma3_line.SetLineColor(ROOT.kRed)
	sigma2_line.SetLineColor(ROOT.kRed)
	sigma3_line=ROOT.TLine(0.75,0.0011,4.0,0.0011) 
	sigma2_line=ROOT.TLine(0.75,0.0202,4.0,0.0202) 
	sigma1_line=ROOT.TLine(0.75,0.1469,4.0,0.1469) 
	sigma3_line.SetLineColor(ROOT.kRed)
	sigma2_line.SetLineColor(ROOT.kRed)
	sigma1_line.SetLineColor(ROOT.kRed)
	sigma3_line.SetLineWidth(2)
	sigma2_line.SetLineWidth(2)
	sigma1_line.SetLineWidth(2)
	#veto_d0.SetLineWidth(0)
#	veto_d0_white.SetFillColor(ROOT.kWhite)
#	veto_d0.SetFillColor(ROOT.kGray)
#	veto_jpsi_white.SetFillColor(ROOT.kWhite)
#	veto_jpsi.SetFillStyle(3013)
#	veto_jpsi.SetFillColor(ROOT.kGray)
#	veto_d0.SetLineWidth(0)
 	#veto_d0_white.Draw('same')
 #	veto_d0.Draw('same')
#	veto_jpsi_white.Draw('same')
#	veto_jpsi.Draw('same')
	sigma3_line.Draw('same')
	sigma2_line.Draw('same')
	sigma1_line.Draw('same')
	canva.Update()
	canva.RedrawAxis()
	CMSlumi(canva,
	 	iPeriod = 5,
  		iPosX = 11,
	        lumiText = '',
	 	extraText="Preliminary") 
		
		
	canva.SaveAs("pval_observed_"+tag+"_dirac.pdf")

		

limSurface("pNN0p990_"+sys.argv[1],int(sys.argv[2]))
#g_1D_obs,g_1D_exp,g_1D_1sigma,g_1D_2sigma =  exclusionPlot("pNN0p990_"+sys.argv[1],int(sys.argv[2]),False,1.0,0,0,0.60,0.30,0.10,'mumu')
#plot_pval("pNN0p990_"+sys.argv[1])
#g_obs,g_exp,g_1sigma,g_2sigma =exclusionPlot("pNN0p990_"+sys.argv[2])
#exclusionUnblind(g_1D_exp,g_1D_1sigma,g_1D_2sigma,g_obs,"noVetos")

import os
import sys 
from numpy import loadtxt
import numpy as np
from datetime import datetime
from os import path
import glob
import ROOT
import multiprocessing
from multiprocessing import Process,Value,Array, Pool
from multiprocessing.pool import ThreadPool as threadPool

ROOT.gROOT.LoadMacro('../macros/setStyle.C')

def pseudoexp(datacard,nToys,rmin,rmax,inject):
	
	log = datacard.strip('.txt')+'.log'
	toysFile = datacard.strip('.txt')+'_'+nToys+'toys.root'
	splitter = datacard.split('/')	
	card = splitter[len(splitter)-1]
	out = card.strip('.txt')+'_Inject_{}'.format('{:.1f}'.format(inject).replace('.','p'))
	
	splitter = datacard.split('_')
#	category = splitter[len(splitter)-2]+'_'+splitter[len(splitter)-1].strip('.txt')
	categories = 'pdfindex_lxysig50to150_OS_2018_13TeV,pdfindex_lxysig50to150_SS_2018_13TeV,pdfindex_lxysig0to50_OS_2018_13TeV,pdfindex_lxysig0to50_SS_2018_13TeV,pdfindex_lxysiggt150_OS_2018_13TeV,pdfindex_lxysiggt150_SS_2018_13TeV'

	print('combine -M FitDiagnostics {}  -t {} --rMin {} --rMax {} --freezeParameters {} --cminDefaultMinimizerStrategy=0 --expectSignal {} --name {} --out ~/Analysis/python/  >> {}'.format(datacard,nToys,rmin,rmax,categories,inject,out,log))
	os.system('combine -M FitDiagnostics {} --toysFrequentist  -t {} --rMin {} --rMax {} --freezeParameters {} --cminDefaultMinimizerStrategy=0 --expectSignal {} --name {} >> {}'.format(datacard,nToys,rmin,rmax,categories,inject,out,log))

	out = 'fitDiagnostics'+out+'.root'
#	os.system('mv fitDiagnosticsTest.root '+out)
    	return out 

def plotAndFit(inject,label,out):


	ROOT.setStyle()
	#out = 'fitDiagnosticsTest.root'
	mc_file = ROOT.TFile.Open(out)

	histo = ROOT.TH1D('histo','histo',25,-5,5)
	
	sb_tree =mc_file.Get('tree_fit_sb')
	sb_tree.Draw('(r-'+np.str(inject)+')/(0.5*(rHiErr+rLoErr))>>histo','fit_status != -1')
	sb_tree.Draw('r>>h')
	h = ROOT.gDirectory.Get("h");
	mean = h.GetMean()
	sigma = h.GetStdDev()
	

	fit = histo.Fit('gaus','s','F',-4,4)

	c = ROOT.TCanvas('c','c',800,600)
	histo.SetFillColorAlpha(9, 0.571)
	histo.SetLineColor(9)
	histo.SetLineWidth(2)
	ROOT.gROOT.GetFunction('gaus').SetLineColor(ROOT.kRed)
	ROOT.gROOT.GetFunction('gaus').SetLineWidth(2)
	ROOT.gROOT.GetFunction('gaus').SetRange(-4,4)
	histo.Draw('hist')
	ROOT.gROOT.GetFunction('gaus').Draw('same')
	
	l = ROOT.TLatex()

	l.SetTextSize(0.035)
	l.SetTextAlign(13)
	l.DrawLatex(-4.5,histo.GetMaximum(),'injected r = {}, nToys = {}'.format(inject,nToys))
	l.DrawLatex(-4.5,histo.GetMaximum()*0.9,'fit mean = {:.3E}'.format(fit.Parameter(1)))
	l.DrawLatex(-4.5,histo.GetMaximum()*0.8,'fit sigma = {:.3E} '.format(fit.Parameter(2)))

	histo.GetXaxis().SetTitle('(S_{observed}-S_{injected})/#sigma(S_{observed})')
	
	splitter = out.split('/')	
	target = splitter[len(splitter)-1]
	
	c.SaveAs('../plots/SigInjection/'+label+'_'+target.strip('.root')+'_toys_expSignal_{}.pdf'.format('{:.1f}'.format(inject).replace('.','p')))
	c.SaveAs('../plots/SigInjection/'+label+'_'+target.strip('.root')+'_toys_expSignal_{}.png'.format('{:.1f}'.format(inject).replace('.','p')))
	
	return mean,sigma

def throwAndFit(inject):
	
	out = pseudoexp(datacard,nToys,rmin,rmax,inject)
	

	mean = plotAndFit(inject,label,out)
	
	return mean

def plotInjections(mean,injections):

	ROOT.setStyle()
	c = ROOT.TCanvas('c','c',800,600)
	ex = np.zeros(len(injections))
	print mean[:,0].astype('float')
	g = ROOT.TGraphErrors(len(injections),np.array(injections),mean[:,0].astype('float'))#,ex,mean[:,1].astype('float'))
	g.SetMarkerStyle(8)
	g.SetMarkerColor(ROOT.kBlue)
	line = ROOT.TF1('bi','x',-0.5,5)
	line.SetLineColor(13)
	line.SetLineStyle(7)
	line.SetLineWidth(1)
	g.GetXaxis().SetTitle('injected signal strenght')
	g.GetXaxis().SetLimits(-0.5,5)
	g.GetYaxis().SetRangeUser(-0.5,5)
	g.GetYaxis().SetTitle('extracted signal strenght')
	g.Draw('AP')
	line.Draw('same')
	
	c.SaveAs('../plots/SigInjection/Injections'+label+'.pdf')

if __name__ == '__main__':


	datacard= sys.argv[1]
	rmin = -20
	rmax = 20
	nToys = sys.argv[2]
	inject = sys.argv[3]
	label = sys.argv[4]


	injections = [0,0.5,1,2,4]

#out = pseudoexp(datacard,nToys,rmin,rmax,inject)
#mean[0]= plotAndFit(inject,label,out)

 	Tpool = Pool(len(injections))
#	mean[0] = throwAndFit(injections[0])
  	mean = Tpool.map(throwAndFit,injections) 
#	singleTraining(specs[0])
  #	Tpool.termminate()
 	Tpool.close()
 	Tpool.join()
	mean = np.asarray(mean,dtype=np.float32)
	print mean[:,1]
	
	plotInjections(mean,injections)


import os
import sys
import time
import re
import glob
import numpy as np
import ROOT

if (len(sys.argv) != 6):
    print "usage limScan_plotter.py outdirName signalTag minNN maxNN NNstep"
    sys.exit(1)
NN_grains = np.arange(float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]))
sigTag = sys.argv[2]
Mass = sigTag.split("_")[1]
print NN_grains

categories = ["lxysig0to50_OS","lxysig0to50_SS","lxysig50to150_OS","lxysig50to150_SS","lxysiggt150_OS","lxysiggt150_SS"]
 
lims = np.zeros((len(NN_grains),len(categories)))
cuts = np.zeros(len(NN_grains))
temp = np.zeros(len(NN_grains))
cut = np.zeros(len(NN_grains))

for j,nn_cut in enumerate(NN_grains):

	dirname =  "../output/Mass"+Mass+("_{:.2f}".format(nn_cut)).replace(".","p")+"_"+sys.argv[1]+"/"+sigTag+"/"
	print dirname
	with open(dirname+'Limits.txt') as f:
		lines = f.readlines()
		lines = [line.split(" ") for line in lines]

	for i,cat in enumerate(categories):
	
		for line in lines:
			
			if line[0] == cat:
#				print line
				lims[j][i] = line[1]	
			        print cat,nn_cut, line[1]	
gs = []
mg = ROOT.TMultiGraph("mg","mg")
canva = ROOT.TCanvas("c","c",600,800)
for i,cat in enumerate(categories):
	npoints = 0
	for j,nn_cut in enumerate(NN_grains):
		if lims[j][i]!=0:
			temp[j] = lims[j][i]
			cut[j] = nn_cut
			npoints = npoints+1
			print nn_cut, lims[j][i],temp[j]
	print temp
	g = ROOT.TGraph(npoints,cut,temp)
	g.SetMarkerStyle(8)
	g.SetName(cat)
	g.SetTitle(cat)
	gs.append(g)
	mg.Add(gs[i],"PL")		
	print ""

mg.Draw("APMC PLC")

mg.GetXaxis().SetTitle("NN score cut")
mg.GetYaxis().SetTitle("expected 95% limit")
mg.GetYaxis().SetRangeUser(0.001,100)
mg.Draw("APMC PLC")
#canva.BuildLegend()

canva.SetLogy()
canva.SaveAs("../plots/"+sigTag+"_limScan.pdf")

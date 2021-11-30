from ROOT import TH1D,TH2D, TLegend,TCanvas,TChain,TStyle, gStyle, gPad
import numpy as np
import os 
import ROOT


fbranches ="branches.txt"

chLable = ["Muon","PF","LowPt","Tracks"]
lxyLable = ["Under3","Over3Under10","Over10Under20","Over20"]

with open("branches.txt") as fbkg:
   histos =[ l.split(" ") for l in fbkg.readlines()]
#histos = np.genfromtxt(fbranches,dtype=['U35','<f8','<f8','<f8','<f8','U80']) # array contains 4 columns:  signal channel, displacement category,MC signal yield,MC yield error,MC mean, MC sigma
print(histos)
bkg = []
with open("Bkgdatasets.txt") as fbkg:
   bkg =[ l.split(" ") for l in fbkg.readlines()]
print(bkg[0][0])
print(bkg[1][0])
with open("Sigdatasets.txt") as fsig:
   sig =[ l.split(" ") for l in fsig.readlines()]
print(sig)
nSig =len(sig)
nBkg =len(bkg)
print(nBkg) 
nHistos = len(histos)
#print(nDataset)
print(nHistos)
nCategories = 4 
nLxyBin = 4
os.system("mkdir MCVsQCD")
os.system("cd MCVsQCD")
h = []
gStyle.SetPalette(91)
gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

for i in range(0,nHistos):
	 
	#exp = ROOT.RDataFrame("Events","BToMuLPiSmall_ChunkMCPiAsEle.root",{"B_pt"})
	for cat in range(2,4):
		if cat == 2:
			continue
	#	h[i][cat]= []
		for l in range(0,1):
	#		h[i][cat][l] = []
			os.system("mkdir "+np.str("MCVsQCD/lxy_"+lxyLable[l]))
			hSig= []
			hBkg= []
			maxi = -99
			for j in range(0,nSig):
				c = TChain("Events")
				c.Add(np.str(sig[j][0]))
			#	hSig.append(TH1D(histos[i][0]+"_"+chLable[cat]+"_"+lxyLable[l],histos[i][0]+"_"+chLable[cat]+"_"+lxyLable[l],int(histos[i][1]),float(histos[i][2]),float(histos[i][3])))
				hSig.append(TH1D(np.str(i)+"_"+chLable[cat]+lxyLable[l],np.str(sig[j][1]),int(histos[i][1]),float(histos[i][2]),float(histos[i][3])))
			#	h[i][cat][l][j]  exp.Filter("Type =="+np.str(cat)+" && LxyBin=="+np.str(l)).Histo1D({"h","h",np.uint(int(histos[i][1])),float(histos[i][2]),float(histos[i][3])},np.str(histos[i][0])).GetValue() 
#				c.Draw("fabs(max(hnl_l_dxyS,hnl_pi_dxyS))>>"+hSig[j].GetName(),"Type =="+np.str(cat))
				if(cat==3 ):
					massCut =" && hnl_mass>"+ np.str(np.int(sig[j][2])-0.3)+ " && hnl_mass<"+np.str(np.int(sig[j][2])+0.3)
					print(massCut)
					c.Draw(np.str(histos[i][0])+">>"+hSig[j].GetName(),"Type =="+np.str(cat)+" && hnl_charge==0 "+massCut)
				else:
					c.Draw(np.str(histos[i][0])+">>"+hSig[j].GetName(),"Type =="+np.str(cat)+" && hnl_charge==0 && dilepton_mass<5.3")
				hSig[j].SetBinContent(int(histos[i][1]),hSig[j].GetBinContent(int(histos[i][1]))+hSig[j].GetBinContent(int(histos[i][1])+1))
				hSig[j].Scale(1/hSig[j].Integral())
				if maxi < hSig[j].GetMaximum():
					maxi =hSig[j].GetMaximum() 
				hSig[j].SetLineWidth(3)
				hSig[j].GetXaxis().SetTitle(np.str(histos[i][5]))
		#		hSig[j].SetLineColor(j+2)
				
				print(histos[i][0])
			for j in range(0,nBkg):
				c = TChain("Events")
				c.Add(np.str(bkg[j][0]))
			#	hBkg.append( TH1D(histos[i][0]+"_"+chLable[cat]+"_"+lxyLable[l],histos[i][0]+"_"+chLable[cat]+"_"+lxyLable[l],int(histos[i][1]),float(histos[i][2]),float(histos[i][3])))
				hBkg.append( TH1D(np.str(i)+"_"+chLable[cat],np.str(bkg[j][2]),int(histos[i][1]),float(histos[i][2]),float(histos[i][3])))

				if(cat==3):
					massCut =" && hnl_mass>"+ np.str(np.int(sig[0][2])-0.3)+ " && hnl_mass<"+np.str(np.int(sig[0][2])+0.3)
					c.Draw(np.str(histos[i][0])+">>"+hBkg[j].GetName(),"Type =="+np.str(cat)+" && hnl_charge==0  && B_mass<7 && dilepton_mass<5.3 && hnl_l_pt!=hnl_pi_pt")
				else:
					c.Draw(np.str(histos[i][0])+">>"+hBkg[j].GetName(),"Type =="+np.str(cat)+np.str(l)+" && hnl_charge==0 ")
				hBkg[j].SetBinContent(int(histos[i][1]),hBkg[j].GetBinContent(int(histos[i][1]))+hBkg[j].GetBinContent(int(histos[i][1])+1))
				hBkg[j].Scale(np.float(bkg[j][1]))
				if not j== 0:
				
					hBkg[0].Add(hBkg[j]);
				
			hBkg[0].Scale(1/hBkg[0].Integral())
			if maxi < hBkg[0].GetMaximum():
				maxi =hBkg[0].GetMaximum() 
			hBkg[0].SetLineWidth(3)
			hBkg[0].GetXaxis().SetTitle(np.str(histos[i][5]))
			hBkg[0].GetYaxis().SetRangeUser(0,maxi*1.2)
		#	hBkg[0].SetLineColor(8)
			
			canva = TCanvas("c","c",1000,700)		
			#leg = TLegend(0.6,0.7,0.9,0.9)
			plotter = TH1D
			
			hBkg[0].SetFillStyle(3003)
			hBkg[0].Draw("hist PFC")
			for j in range(0,nSig):
				print( "Signal efficiency %f"%(1.0*hSig[j].Integral(hSig[j].GetXaxis().FindBin(float(histos[i][4])),int(histos[i][1])) /hSig[j].Integral()))
        			hSig[j].Draw("hist same PLC")
				hSig[j].SetName("sig");
#				hSig[j].SaveAs(np.str("MCVsQCD/pdf/"+histos[i][6]+"_"+chLable[cat]+"_"+np.str(sig[j][1])+".root"))
		#		leg.AddEntry(hSig[j],np.str(sig[j][1]),"l")
		#	leg.AddEntry(hBkg[0],np.str(bkg[j][2]),"l")
			
			print( "Bkg efficiency %f"%(1.0*hBkg[0].Integral(hBkg[0].GetXaxis().FindBin(float(histos[i][4])),int(histos[i][1]))/hBkg[0].Integral() ))
		#	leg.Draw("PLC")
			gPad.BuildLegend(0.4,0.2,0.4,0.2,chLable[cat]+"Channel","l")
				
			canva.SaveAs(np.str("MCVsQCD/lxy_inclusive/"+histos[i][6]+"_"+chLable[cat]+".pdf"))
			canva.SaveAs(np.str("MCVsQCD/lxy_inclusive/"+histos[i][6]+"_"+chLable[cat]+".C"))
			hBkg[0].SetName("bkg");
			hBkg[0].SetTitle("bkg");
#			hBkg[0].SaveAs(np.str("MCVsQCD/pdf/"+histos[i][6]+"_"+chLable[cat]+"_bkg.root"))
				
		
os.system("cd ..")

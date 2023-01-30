import os 
import numpy as np
import ROOT
import glob
from multiprocessing import Process,Value,Array, Pool
from multiprocessing.pool import ThreadPool as threadPool
#from HNLUtils import BR_HNLmupion,BR_HNLelepion,getVV,BToMuX,BToEX

def vectorizeList(inputs):

	Vec = ROOT.std.vector(ROOT.std.string)()
	for file_ in inputs:
		Vec.push_back(file_)

	return Vec


dataIn ="/cmshome/ratramon/Analysis/data/HNLFlatTuples/Parking1A0_newSig_4varsL_section0/*.root"
QCDIn = glob.glob("/cmshome/ratramon/Analysis/data/HNLFlatTuples/Pt-*_newSig/*.root")
print QCDIn
electron_id_sf_file = ROOT.TFile.Open('EGamma_stdEta_SF.root')
electron_id_sf= electron_id_sf_file.Get('EGamma_SF2D')
ROOT.gInterpreter.ProcessLine("auto h_sf = EGamma_SF2D;")
getSF ="  {double id_sf;\
	  id_sf = h_sf->GetBinContent(h_sf->GetXaxis()->FindBin(hnl_l_eta),h_sf->GetYaxis()->FindBin(hnl_l_pt));\
          return id_sf;}"
data = ROOT.RDataFrame("Events",dataIn)
qcd = ROOT.RDataFrame("Events",vectorizeList(QCDIn))


qcd = qcd.Define("id_sf",getSF)
qcd = qcd.Define("tot_weight","QCDweight*PU_weight")
qcd_sf = qcd.Histo1D("id_sf")

toPlot = [
          [ROOT.RDF.TH1DModel("name", "title", 50, 0, 100),"hnl_pi_DCAS"," #pi DCA Sig"],
          [ROOT.RDF.TH1DModel("name", "title", 30, -2.5, 2.5),"hnl_l_eta","#eta"],
          [ROOT.RDF.TH1DModel("name", "title", 20, -8, 8),"hnl_l_mvaId","PF mvaId"],
       #   [ROOT.RDF.TH1DModel("name", "title", 100, 0, 20),"hnl_l_dxyS","Signficance of the xy impact parameter (cm)"],
	]

for branch in toPlot:

	c = ROOT.TCanvas()
	print branch
	dataHist = data.Histo1D(branch[0],str(branch[1]))	
	qcdHist = qcd.Histo1D(branch[0],str(branch[1]),"tot_weight")
	dataHist.SetMarkerStyle(8)
	qcdHist.SetFillStyle(1001)
	qcdHist.SetFillColor(ROOT.kViolet-9)
	qcdHist.SetLineColor(ROOT.kViolet-9)
	qcdHist.SetLineColor(ROOT.kViolet-9)
	qcdHist.SetTitle("MC QCD p_{T} 15-300 GeV;"+branch[2]+";a.u.")
	dataHist.SetTitle("data;"+branch[2]+";a.u.")
	h_data = dataHist.GetPtr()
	h_data.Scale(1./h_data.Integral())
	h_qcd = qcdHist.GetPtr()
	h_qcd.Sumw2(ROOT.kTRUE)
	h_qcd.Scale(1./h_qcd.GetSumOfWeights())
	rp = ROOT.TRatioPlot(h_qcd,h_data)
	#qcdHist.DrawNormalized("hist")
	#dataHist.DrawNormalized("sameEP")
	rp.Draw()
	rp.GetLowerRefYaxis().SetRangeUser(0.0,2)	
	rp.GetLowerRefGraph().SetMarkerStyle(8)
	c.Update()
#	c.BuildLegend()
	c.SaveAs("../plots/QCD_data_"+branch[1]+"_withEB_EE_SF.pdf")

import os
import sys
import ROOT
from NN_validationPlots import humanLabel
from HNLUtils import BR_HNLmupion,BR_HNLelepion,getVV,BToMuX,BToEX,BcToMuX,BcToEX
#from collections import OrderedDict
import CMSlumi

#from tools import Tools
#sys.path.append('../objects')
#from categories import categories

#tools = Tools()

output_label = 'V12_30Jul23' 
tag = 'pNN0p990_pNN_hits' 

tag_0 = 'pNN0p000_presel' 

masses = [1.0, 1.5, 2.0, 3.0,4.5]
sigma = [ 0.0014 +0.0074 * m for m in masses]
sigma_old = [ 0.010,0.013,0.017,0.025,0.035]
scales = [1000, 100, 1, 1, 0.001]
lumi = 41.6 
processedLumi = 4.91
scales = [s*41.6/4.91 for s in scales]
ctaus = [0.1,10,1000.0,10000.]
#categories = ['LxySUnder50_OS','LxySOver50Under150_OS','LxySOver150_OS','LxySUnder50_SS','LxySOver50Under150_SS','LxySOver150_SS']
categories = ["lxysig0to50_OS","lxysig50to150_OS","lxysiggt150_OS","lxysig0to50_SS","lxysig50to150_SS","lxysiggt150_SS"]

def getCouplingLabel( v2):
  coupling = "{:e}".format(v2)
  part1 = coupling[:coupling.find('e')]
  part1 = str(round(float(part1), 1))
  part2 = coupling[coupling.find('e'):]
  return (part1+part2)

def getSignalYields(line,scale):
  print line
  yields = float(line.split(' ')[3].strip('\t'))
  print line[3].strip('\t')
  #yields = float(line[5:line.rfind('1.0')-1])*scale
  return yields
#  return '{:.2e}'.format(yields)


def getBackgroundYields(workspace,scale):
#  print 'hereee_____________________________'
#  print workspace.data("data_obs").sumEntries()
  yields = workspace.data("data_obs").sumEntries()
#  print  "yields _____________________",yields
  yields = yields * 0.2 *scale#(2 sigma / 10 sigma)
#  return '{:.2e}'.format(yields)
  return yields

def getBackgroundYields_pre(workspace,scale):
#  print 'hereee_____________________________'
#  print workspace.data("data_obs").sumEntries()
  yields = workspace.data("data_obs").sumEntries()
#  print  "yields _____________________",yields
  yields = yields * 0.2 *scale#(2 sigma / 10 sigma)
#  return '{:.2e}'.format(yields)
  return yields

def plotSuperposition(path,presel,pNN,category):

	
	binWidth = pNN.var('hnl_mass').getBinWidth(0,'plotter')
#	min = presel.var('hnl_mass').getMin()
#	max = presel.var('hnl_mass').getMax()
#	data_pre = presel.data("data").createHistogram('',presel.var('hnl_mass'),ROOT.RooFit.Binning(60))
#	data_pNN = pNN.data("data").createHistogram('',pNN.var('hnl_mass'),ROOT.RooFit.Binning(60))
	data_pre = presel.data("data_obs")
	data_pNN = pNN.data("data_obs")
	if not data_pre:
		print 'no plot'
		return 0 
#	print data_pre.GetEntries(), data_pNN.GetEntries()
 #	data_pre.SetName("data - preselection+baseline")
 #	data_pNN.SetName("data - preselection+baseline+pNN score>0.95")
 #	data_pre.SetTitle("data - preselection+baseline;HNL mass(GeV);normalized to unity")
 #	data_pNN.SetTitle("data - preselection+baseline+pNN score>0.95")
	ROOT.gStyle.SetOptStat(0)	
	ROOT.gStyle.SetPadLeftMargin(0.16)

	c = ROOT.TCanvas("c","c",800,800)
	p_pnn = pNN.var('hnl_mass').frame()
	p = presel.var('hnl_mass').frame()
	presel.var('hnl_mass').setBins(60,'plotter')
	pNN.var('hnl_mass').setBins(60,'plotter_pnn')
	binWidth = presel.var('hnl_mass').getBinning('plotter').binWidth(0)
	binWidth_pnn = pNN.var('hnl_mass').getBinning('plotter_pnn').binWidth(0)
	
	plot_presel=data_pre.plotOn(p,ROOT.RooFit.Name('data'),ROOT.RooFit.Binning('plotter'),ROOT.RooFit.Rescale(1./(data_pre.sumEntries())),ROOT.RooFit.DrawOption('B'),ROOT.RooFit.FillColor(ROOT.kBlue-3),ROOT.RooFit.LineColor(ROOT.kBlue-3),ROOT.RooFit.LineWidth(2),ROOT.RooFit.FillStyle(3344))
	plot_pnn= data_pNN.plotOn(p_pnn,ROOT.RooFit.Name('data'),ROOT.RooFit.Binning('plotter_pnn'),ROOT.RooFit.Rescale(1./(data_pNN.sumEntries())),ROOT.RooFit.DrawOption('EP'),ROOT.RooFit.MarkerStyle(8),ROOT.RooFit.MarkerColor(ROOT.kBlue-3),ROOT.RooFit.FillColor(ROOT.kBlue-3),ROOT.RooFit.LineColor(ROOT.kBlue-3))
#	data_pre.plotOn(p,ROOT.RooFit.Name('data - preselection+baseline'),ROOT.RooFit.FillStyle(3144),ROOT.RooFit.FillColor(ROOT.kBlue-3),ROOT.RooFit.LineColor(ROOT.kBlue-3),ROOT.RooFit.DrawOption('LF'),ROOT.RooFit.VLines()
#	data_pNN.plotOn(p_pnn,ROOT.RooFit.Name('data - preselection+baseline+pNN score>0.95'),ROOT.RooFit.MarkerStyle(8),ROOT.RooFit.FillColor(ROOT.kBlue-3),ROOT.RooFit.LineColor(ROOT.kBlue-3),ROOT.RooFit.DrawOption('EP1

	p.SetMaximum(max(plot_presel.GetMaximum()/data_pre.sumEntries(),plot_pnn.GetMaximum()/data_pNN.sumEntries())*1.5)
	p.SetTitle("")
	p.GetYaxis().SetTitle('normalized to unity')
	p.GetXaxis().SetTitle('m_{N}(GeV)')
	p.Draw("LF")
	p_pnn.Draw("same")

 	l = ROOT.TLegend(0.2,0.55,0.75,0.65)
 	l.AddEntry(p.findObject("data"),'data - no pNN requirement','F')
 	l.AddEntry(p_pnn.findObject("data"),'data - pNN score > 0.99','P')
 	l.SetBorderSize(0)
 	l.Draw()
	CMSlumi.CMSlumi(c,iPeriod=0,iPosX=11,lumiText='',extraText='Preliminary '+humanLabel(category))
#	c.BuildLegend()
#	data_pre.SetFillStyle(3144)
#	data_pre.SetLineColor(ROOT.kBlue-3)
#	data_pNN.SetMarkerColor(ROOT.kBlue-3)
#	data_pNN.SetLineColor(ROOT.kBlue-3)
#	data_pNN.SetMarkerStyle(8)
#	data_pre.Scale(1.0/60.)
#	data_pre.Draw("HIST F")
#	data_pNN.Scale(1.0/60.)
#	data_pNN.Draw("EPsame")
	print path+'/superposition_pNN0p0_'+category+'.pdf'
	c.SaveAs(path+'/superposition_pNN0p0_'+category+'.pdf')

def getCategoryTitle(category):
  title = category
  title = title.replace('lxysig50to150', '50$<$L$_\mathrm{xy}/\sigma<$150')
  title = title.replace('lxysig0to50', '$L$_\mathrm{xy}/\sigma<$50')
  title = title.replace('lxysiggt150', '$L$_\mathrm{xy}/\sigma>$150')
  title = title.replace('lxysig', 'L$_\mathrm{xy}/\sigma$')
  title = title.replace('#sigma_{xy}', '$\\sigma$')
  title = title.replace('Under', '$\leq$')
 # title = title.replace('Under', '$<$')
  title = title.replace('Over', '$>$')
  title = title.replace('$$', '')
  title = title.replace('_OS', ', OS')
  title = title.replace('_SS', ', SS')
  return title

table_yields = open('table_yields.txt', 'w+')

for imass,mass in enumerate(masses):
  coupling_line = '{} & &'.format(mass)
  ctau_line = '& & & '
  for ictau, ctau in enumerate(ctaus):
    v2 = getVV(mass, ctau)
    coupling = getCouplingLabel(v2)
    if ictau != len(ctaus)-1:
      coupling_line += ' & {}'.format(coupling)
      ctau_line = ctau_line + '{:.1f} mm &'.format(ctau)

    else:
      coupling_line += ' & {} \\\ '.format(coupling)
      ctau_line = ctau_line + '{:.1f} mm \\'.format(ctau)
  if imass ==0:
	  table_yields.write('\hline')
	  table_yields.write('\n'+ctau_line)
  table_yields.write('\hline')
  table_yields.write('\n' + coupling_line)
  table_yields.write('\n' + '\hline')

  for icat, category in enumerate(categories): 
    signal_yields = []
    signal_yields_pre = []
    signal_efficiency = []
 #   print mass, sigma[imass]
   # path_0 = '../output/Mass{}_{}/Background/DiscreteProfiling/window_m{:.2f}_s{:.3f}_ns10/ws/'.format(str(mass).replace('.', 'p'), tag_0,mass,sigma_old[imass])
    #path_0 = '../output/Mass{}_{}/Background/DiscreteProfiling/window_m{:.2f}_s{:.3f}_ns10/ws/'.format(str(mass).replace('.', 'p'), tag_0,mass,sigma[imass])
    path_0 = '../output/Mass{}_{}/Background/DiscreteProfiling/window_m{:.2f}_s{:.3f}_ns10/ws/'.format(('{:.2f}'.format(mass)).replace('.', 'p'), tag_0,mass,sigma[imass])
    path = '../output/Mass{}_{}/Background/DiscreteProfiling/window_m{:.2f}_s{:.3f}_ns10/ws/'.format(('{:.2f}'.format(mass)).replace('.', 'p'), tag,mass,sigma[imass])
    #data_obs_name = 'CMS-BHNL_multipdf_{}.root'.format(category)
    #data_obs_name_0= 'workspace_multipdf_bhnl_m_{}_cat_{}.root'.format(('{:.1f}'.format(mass)).replace('.', 'p'),category)
    data_obs_name_0= 'workspace_multipdf_bhnl_m_{}_cat_{}.root'.format(('{:.2f}'.format(mass)).replace('.', 'p'),category)
    data_obs_name= 'workspace_multipdf_bhnl_m_{}_cat_{}.root'.format(('{:.2f}'.format(mass)).replace('.', 'p'),category)
    try:
 #     print(path+data_obs_name)
      data_obs_file = ROOT.TFile.Open(path+data_obs_name, 'READ')
      workspace = data_obs_file.Get("multipdf")
    #  print category,workspace.data('data_obs').sumEntries()
      background_yields = getBackgroundYields(workspace,1.0)
    except:
      background_yields = '-'
    try:
      print(path_0+data_obs_name_0)
      data_obs_file_pre = ROOT.TFile.Open(path_0+data_obs_name_0, 'READ')
      workspace_pre = data_obs_file_pre.Get("multipdf")
      background_yields_pre = getBackgroundYields_pre(workspace_pre,1.0)
      print category, workspace_pre.data('data_obs').sumEntries()
    except:
      background_yields_pre = '-'
 #   print background_yields
    print background_yields,background_yields_pre
    if background_yields != '-' and background_yields_pre != '-':
#	 print background_yields,background_yields_pre
   	 table_entry = ' & {} & {:.2f} \%'.format(getCategoryTitle(category), float(background_yields)/float(background_yields_pre)*100)
   	# table_entry = ' & {} & {:.2E} '.format(getCategoryTitle(category), float(background_yields))
    else: 
   	 table_entry = ' & {} & {}'.format(getCategoryTitle(category),background_yields)
    plotSuperposition(path+'../../.',workspace_pre,workspace,category)
    for ictau, ctau in enumerate(ctaus):
      v2 = getVV(mass, ctau)
      coupling = getCouplingLabel(v2)

     # path_0 = '../output/Mass{}_{}/HNL_{}_ctau{}_PF/Datacards'.format(str(mass).replace('.', 'p'), tag_0,str(mass).replace('.', 'p'),'{:.1f}'.format(ctau).replace('.', 'p'))
     # path_0 = '../output/Mass{}_{}/HNL_Mass{}_ctau{}_PF/Datacards'.format(str(mass).replace('.', 'p'), tag_0,str(mass).replace('.', 'p'),'{:.3f}'.format(ctau).replace('.', 'p'))
      path_0 = '../output/Mass{}_{}/HNL_Mass{}_ctau{}_PF/Datacards'.format('{:.2f}'.format(mass).replace('.', 'p'),tag_0,'{:.2f}'.format(mass).replace('.', 'p'),'{:.3f}'.format(ctau).replace('.', 'p'))
      path = '../output/Mass{}_{}/HNL_Mass{}_ctau{}_PF/Datacards'.format('{:.2f}'.format(mass).replace('.', 'p'),tag,'{:.2f}'.format(mass).replace('.', 'p'),'{:.3f}'.format(ctau).replace('.', 'p'))
  #    datacard_name_0 = 'HNL_M{:.1f}_ctau{:.1f}_PF_{}.txt'.format(mass, ctau, category)
      #datacard_name_0 = 'HNL_m_{}_ctau_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'),'{:.3f}'.format(ctau).replace('.', 'p'), category)
      #datacard_name_0 = 'HNL_m_{}_ctau_{}_cat_{}.txt'.format('{:.2f}'.format(mass).replace('.', 'p'),'{:.3f}'.format(ctau).replace('.', 'p'), category)
      datacard_name_0 = 'HNL_m_{}_ctau_{}_cat_{}.txt'.format('{:.2f}'.format(mass).replace('.', 'p'),'{:.3f}'.format(ctau).replace('.', 'p'), category)
      datacard_name = 'HNL_m_{}_ctau_{}_cat_{}.txt'.format('{:.2f}'.format(mass).replace('.', 'p'),'{:.3f}'.format(ctau).replace('.', 'p'), category)
      try:
      	print('{}/{}'.format(path_0, datacard_name_0))
        card = open('{}/{}'.format(path_0, datacard_name_0))
        lines = card.readlines()
        for line in lines:
          if 'rate' not in line: continue
          signal_yields_pre.append(getSignalYields(line,1.0))
          print line
         # print ctau
      except:
          signal_yields_pre.append('-')
      try:
     # 	print('{}/{}'.format(path, datacard_name))
        card = open('{}/{}'.format(path, datacard_name))
        lines = card.readlines()
        for line in lines:
          if 'rate' not in line: continue
          signal_yields.append(getSignalYields(line,1.0))
#	  print line
#	  print getSignalYields(line,1.0)
#	  print ctau
      except:
          signal_yields.append('-')
    print 'signal ',signal_yields, signal_yields_pre
    for ictau, ctau in enumerate(ctaus):
      if ictau != len(ctaus)-1: 
#	print "_________________________________________________________here!"
#	print  signal_yields[ictau],ctau
	if signal_yields[ictau] != '-' and signal_yields_pre[ictau] != '-' :
	  #print ctau, signal_yields[ictau],signal_yields_pre[ictau]
      	  table_entry += ' &  {:.2f} \%'.format(float(signal_yields[ictau])/float(signal_yields_pre[ictau])*100)
#      	  table_entry += ' &  {:.2E} '.format(float(signal_yields[ictau]))
	else:
      	  table_entry += ' &  {} '.format(signal_yields[ictau])
      else:
	if signal_yields[ictau] != '-' and signal_yields_pre[ictau] != '-' :
	  print ctau, signal_yields[ictau],signal_yields_pre[ictau]
      	  table_entry += ' &  {:.2f} \% \\ \\'.format(float(signal_yields[ictau])/float(signal_yields_pre[ictau])*100)
      #	  table_entry += ' &  {:.2E}  \\ \\'.format(float(signal_yields[ictau]))
	else:
      	  table_entry += ' &  {} \\ \\'.format(signal_yields[ictau])
    table_yields.write('\n' + table_entry)
    if icat == len(categories)-1:
      table_yields.write('\n' + '\hline')


table_yields.close()
print '-> table_yields.txt created'
os.system('cat table_yields.txt')

    



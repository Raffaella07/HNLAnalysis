import os
import sys


class window(object):
  def __init__(self,m,s,ns):
    self.m = m
    self.s = s
    self.ns = ns



if __name__ == "__main__":
 
  #########################
  ## Options
  ######################## 
  tag = sys.argv[4]
  mass_in = float(sys.argv[1])
  cat_in = sys.argv[2]
  n_sigma = int(sys.argv[5])
  doBc = int(sys.argv[6])
  print mass_in  
  doCopyEOS = False
  doCleanEOS = False
  doRun = True
  doBlind = True
  windows = [window(mass_in, 0.0014 +0.0074 *mass_in,n_sigma) ]
# windows = [
#window(1.0, 0.0014 +0.0074 *1.00,10),
#window(1.2,0.01,10),
#window(1.5, 0.0014 +0.0074 *1.5,10),
#window(2.0, 0.0014 +0.0074 *2.0,10),
#window(2.5, 0.0014 +0.0074 *1.005),
#window(2.75 0.0014 +0.0074 *1.005,10),
#window(3.0, 0.0014 +0.0074 *3.0,10),
####window(3 0.0014 +0.0074 *1.00025,5),
# window(3.3 0.0014 +0.0074 *1.005,10),
#window(4.0, 0.0014 +0.0074 *1.00,10),
#window(4.5, 0.0014 +0.0074 *4.5,10),
#window(5.0,0.035,10),
# window(5.5,0.035,10),
# ]

  cats = [
 "lxysig0to50_OS", 
 "lxysig50to150_OS", 
 "lxysiggt150_OS", 
 "lxysig0to50_SS", 
 "lxysig50to150_SS", 
 "lxysiggt150_SS", 
  ]
  infiles = [
 "Data_PFele_0_OS.root",
 "Data_PFele_1_OS.root",
 "Data_PFele_2_OS.root",
 "Data_PFele_0_SS.root",
 "Data_PFele_1_SS.root",
 "Data_PFele_2_SS.root"
  ]
  #######################
  #######################
  B_c = ''
  #doBlind = False
  if doCopyEOS and doCleanEOS:
    # clean existing dir in eos
    os.system("rm -rf /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/{tag}/".format(tag=tag))
    os.system("mkdir /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/{tag}/".format(tag=tag))
    os.system("cp HTACCESS /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/{tag}/.htaccess".format(tag=tag))
 
  
  for w in windows:
    wtag = "m{:.2f}_s{:.3f}_ns{}".format(w.m,w.s,w.ns)
    if mass_in != w.m:
	print type(mass_in),type(w.m)
	continue
    for i,cat in enumerate(cats):
      if cat_in != cat: 
	continue
      inpath = sys.argv[3]
      if doBc:
	print infiles[i]
	infiles[i] = infiles[i].replace('.root','_Bc.root')
	print infiles[i]
        cat = cat+'Bc_'
      inws = inpath + infiles[i]
      outws = inpath+"/{tag}/window_{wtag}/ws/workspace_multipdf_bhnl_m_{mass}_cat_{cat}.root".format(tag=tag,wtag=wtag,mass='{:.2f}'.format(w.m).replace('.','p'),cat=cat)
      outdir = inpath+"/{tag}/window_{wtag}/".format(tag=tag,wtag=wtag)
      command = "/cmshome/ratramon/CMSSW_10_2_13/src/flashggFinalFit/Background/bin/fTest -i {inws} -v  --saveMultiPdf {outws}  -D {outdir} -f {cat} --mN {m} --sigma {s} --nsigma {ns}  \
                   --isData 1 --year 2018 --catOffset {catoffset} {doblind}".format(
                     inws=inws,
                     outws=outws,
                     outdir=outdir,
                     cat=cat,
                     m=w.m,
                     s=w.s,
                     ns=w.ns,
                     catoffset=i,
                     doblind="--blind" if doBlind else ''
                   )
	

      if not os.path.exists(outdir+'ws'): 
        os.makedirs(outdir+'ws')

      print "\n\n\n\n"
      print "*********************************"
      print "===> F-test analysis:"
      print "         mass   : {:.2f}".format(w.m)
      print "         sigma  : {:.2f}".format(w.s)
      print "         nsigma : {:.2f}".format(w.ns)
      print "         cat    : {}".format(cat)
      print "*********************************"
      print command
      if doRun:
        os.system(command) 

    if doCopyEOS:
      print "===> Copying results to EOS" 
      os.system("cp -r {outdir} /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/{outdir}".format(outdir=outdir))
      os.system("cp HTACCESS_SUB /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/{outdir}/.htaccess".format(outdir=outdir))
      print "*********************************"






# Config file: options for signal fitting

backgroundScriptCfg = {
  
  # Setup
  'inputWSDir':'/work/mratti/BParkingNano/Oct_CMSSW_10_2_15/CMSSW_10_2_15/src/PhysicsTools/BHNLAnalysis/exampleWS/wsAnalysis/', # location of 'allData.root' file
  #'inputWSDir':'/work/anlyon/flashgg/CMSSW_10_2_13/src/flashggFinalFit/Background/test/mine/data_cats_sel_large/', # location of 'allData.root' file
  #'inputWSDir':'/work/anlyon/flashgg/CMSSW_10_2_13/src/flashggFinalFit/Background/test/mine/data_cats_sel/', # location of 'allData.root' file
  'cats':'auto', # auto: automatically inferred from input ws
  'catOffset':0, # add offset to category numbers (useful for categories from different allData.root files)  
  'ext':'testMGRealData', # extension to add to output directory
  'year':'2018', # Use combined when merging all years in category (for plots)

  # Job submission options
  'batch':'local', # [condor,SGE,IC,local]
  'queue':'hep.q' # for condor e.g. microcentury
  
}

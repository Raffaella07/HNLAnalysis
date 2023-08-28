#! /usr/bin/env python
import os
import sys
import time
import re
import glob
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################

if (len(sys.argv)!=6):
	print "usage: python SkimNano_submit doMCSignal doMC_QCD doBParking_1A doBc"
	print "all arguments are boolean"
        sys.exit(1)	
doMCSignal =int( sys.argv[1] )
doMC_QCD = int(sys.argv[2] )
doBParking = int(sys.argv[3])
doBc= int(sys.argv[4])
doPrivate= int(sys.argv[5])

QCDTuplesPath = '/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/' 
DataTuplesPath = '/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/ParkingBPH' 
MCTuplesPath = '/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/' 

MCPrivate_datasets = ["1p2_ctau10p0_aa","1p2_ctau10p0_ab","1p2_ctau10p0_ac","1p2_ctau10p0_ad",
		     "1p3_ctau10p0_aa","1p3_ctau10p0_ab","1p3_ctau10p0_ac","1p3_ctau10p0_ad",
		     "1p4_ctau10p0_aa","1p4_ctau10p0_ab","1p4_ctau10p0_ac","1p4_ctau10p0_ad",
		     "1p8_ctau10p0_aa","1p8_ctau10p0_ab","1p8_ctau10p0_ac",
		     "2p2_ctau10p0","2p3_ctau10p0","2p4_ctau10p0",
		     "2p5_ctau10p0_aa","2p5_ctau10p0_ab","2p5_ctau10p0_ac","2p5_ctau10p0_ad",
		     "2p8_ctau10p0","2p9_ctau10p0_aa","2p9_ctau10p0_ab",
		     "3p1_ctau1p0","3p2_ctau1p0","3p4_ctau1p0","3p7_ctau1p0",
		]
MCSignal_datasets = ["1p0_ctau1000p0","1p0_ctau100p0","1p0_ctau10p0",
                     "1p5_ctau1000p0","1p5_ctau100p0","1p5_ctau10p0",
		     "2p0_ctau1000p0","2p0_ctau100p0","2p0_ctau10p0",
		     "3p0_ctau1000p0","3p0_ctau100p0","3p0_ctau10p0","3p0_ctau1p0",
		     "4p5_ctau100p0","4p5_ctau10p0","4p5_ctau1p0","4p5_ctau0p1"
		]

MCSignal_Bc_datasets = ["3p0_ctau1000p0","3p0_ctau100p0","3p0_ctau10p0","3p0_ctau1p0",
		        "4p5_ctau10p0","4p5_ctau1p0","4p5_ctau0p1",
		        "5p5_ctau10p0","5p5_ctau1p0","5p5_ctau0p1","5p5_ctau0p01",
			]
QCDlables = [ "Pt-15to20_newSig",
              "Pt-20to30_newSig",
              "Pt-30to50_newSig",
              "Pt-50to80_newSig",
              "Pt-80to120_newSig",
              "Pt-80to120_ext_newSig",
              "Pt-120to170_newSig",
              "Pt-120to170_ext_newSig",
              "Pt-170to300_newSig"
	     ]
MC_QCD_datasets =["QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-15to20/220819_**/0000/",
		  "QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-20to30/220819_**/0000/",
		  "QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-30to50/220819_**/0000/",
		  "QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-50to80/220819_**/0000/",
		  "QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-80to120/220819_**/0000/",
		  "QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-80to120_ext/220819_**/0000/", 
		  "QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-120to170_ext/220819_**/0000/",
		  "QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-120to170/220819_**/0000/",
		  "QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-170to300/220819_**/0000/"
		]
BParking_datasets  = [
 			['B','1','230313'],
 			['B','2','230314'],
 			['B','3','230315'],
 			['B','4','230315'],
 			['B','5','230315'],
 			['B','6','230315'],
 			['C','1','230312'],
 			['C','2','230312'],
 			['C','3','230312'],
 			['C','4','230312'],
 			['C','5','230312'],
 			['D','2','230315'],
 			['D','3','230312'],
 			['D','4','230312'],
 			['D','5','230315'],
#			['D','1','230310'],
			['A','1','230319'],
	 		['A','2','230319'],
	 		['A','3','230319'],
	 		['A','4','230321'],
	 		['A','5','230319'],
	 		['A','6','230319'],
			]



if (doMCSignal and not doBc and not doPrivate):
	for name in MCSignal_datasets:
		path="/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples//BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL"+name+"mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass"+name+"/220810_**/0000/"
		
		nfiles = len(glob.glob(path+"*.root"))
		
		os.system( "python sendOnBatch.py "+ "NewSig_hits_"+name+" "+path+" "+str(nfiles));
if (doMCSignal and not doBc and doPrivate):
	path="/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/private_B_mass4p6*/BuToHNLToLPi_2023Aug*/*/0000/"
	MCPrivate_datasets = glob.glob(path)
	for name in MCPrivate_datasets:
		print name
		
		sig_tag = (name.split('/'))[9]
		print sig_tag
		sig_tag = sig_tag.split('_')
		if len(sig_tag) == 5:
			s_tag = sig_tag[2].strip('mass')+'_'+sig_tag[3]+'_'+sig_tag[4]
		else:
			s_tag = sig_tag[2].strip('mass')+'_'+sig_tag[3]
#		print s_tag
		nfiles = len(glob.glob(name+"*.root"))
		
		print( "python sendOnBatch.py "+ "NewSig_ordered_"+s_tag+" "+name+" "+str(nfiles));
		os.system( "python sendOnBatch.py "+ "NewSig_flav_"+s_tag+" "+name+" "+str(nfiles));
if (doMCSignal and doBc and doPrivate):
	path="/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/private_B_Bc_mass3p2_*/BuToHNLToLPi_2023Aug11/*/0000/"
	MCPrivate_datasets = glob.glob(path)
	for name in MCPrivate_datasets:
		print name
		
		sig_tag = (name.split('/'))[9]
		sig_tag = sig_tag.split('_')
		s_tag = sig_tag[2].strip('mass')+'_'+sig_tag[3]
#		print s_tag
		nfiles = len(glob.glob(name+"*.root"))
		
	#	print( "python sendOnBatch.py "+ "NewSig_private_"+s_tag+" "+name+" "+str(nfiles));
		os.system( "python sendOnBatch.py "+ "NewSig_private_"+s_tag+" "+name+" "+str(nfiles));
if (doMCSignal and doBc and not doPrivate):
	for name in MCSignal_Bc_datasets:
		path="/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/BcToNMuX_NToEMuPi_SoftQCD_b_mN"+name+"mm_TuneCP5_13TeV-pythia8-evtgen/crab_Mass"+name+"_Bc/221018_**/0000/"
		
		nfiles = len(glob.glob(path+"*.root"))
		print nfiles
		os.system( "python sendOnBatch.py "+ "NewSig_hits_"+name+"_Bc "+path+" "+str(nfiles));
if (doMC_QCD):
	for i in range(0,len(MC_QCD_datasets)):
		#retrieve njobs such that each job has max ~0.5G input -> homogeneus processing 
		#print(TuplesPath+MC_QCD_datasets[i])
		dataset_size = 1.0*sum([os.path.getsize(f) for f in glob.glob((QCDTuplesPath+MC_QCD_datasets[i]+"*.root"))])/1000000 #getsize output is in byte, /1000 in MB
		nfiles = len(glob.glob(QCDTuplesPath+MC_QCD_datasets[i]+"*.root"))
		nfilesXjob =int(300*nfiles/dataset_size)

		os.system("mkdir /cmshome/ratramon/Analysis/data/HNLFlatTuples/"+QCDlables[i])
		print ("Dataset size %f, nFiles in dataset %d: will process %d files x job to have ~1G jobs in batch \n"%(dataset_size,nfiles,nfilesXjob))
		os.system( " python sendOnBatch.py  "+QCDlables[i]+"_Mu12 "+QCDTuplesPath+MC_QCD_datasets[i]+" "+str(nfilesXjob));
if (doBParking):
	for tag in BParking_datasets:
			datapath = DataTuplesPath+tag[1]+'/crab_data_Run2018'+tag[0]+'_part'+tag[1]+'_section*'
			section_paths = glob.glob(datapath)
			for isec,sec in enumerate(section_paths):
				subpath = sec+'/'+tag[2]+'_*/*/'
				subdirs = glob.glob(subpath)
				print subdirs	
				for isub,sub in enumerate(subdirs):
					
				#		print( "python sendOnBatch.py Parking_{}_{}_{}_JanProd {} 50".format(tag[1]+tag[0],str(isec),str(isub),sub));
						os.system( "python sendOnBatch.py Parking_hits_{}_{}_{} {} 50".format(tag[1]+tag[0],str(isec),str(isub),sub));

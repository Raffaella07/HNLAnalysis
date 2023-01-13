# first neural network with keras tutorial

import os 
from numpy import loadtxt
import numpy as np
from datetime import datetime
from os import path
import random
from random import uniform
import glob

import tensorflow as tf
import keras
from keras.models import Sequential,Model
from keras.layers import Dense, Input
from keras.constraints import unit_norm
from keras.callbacks import EarlyStopping, Callback, ReduceLROnPlateau, ModelCheckpoint
import sklearn as sk
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.preprocessing import RobustScaler
import pickle
import pandas as pd

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from itertools import product
from root_numpy import root2array, tree2array, fill_hist,array2tree
import numpy.lib.recfunctions as rfn

import ROOT
import seaborn as sns

from NN_trainer import saveFig
import multiprocessing
from multiprocessing import Process,Value,Array, Pool
from multiprocessing.pool import ThreadPool as threadPool

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 20,
        }

def NNinput(outdir):



        model_paths = glob.glob(outdir+'*.h5')
	#print outdir+'*.h5'
	#print model_paths
	model_paths_matrix = [path.split("_") for path in model_paths]
	accuracies =[ float(model_paths_matrix[i][10].strip(".h5")) for i in range (0,len(model_paths_matrix))]
        best_model_idx = accuracies.index(max(accuracies))
#	print best_model_idx,model_paths[best_model_idx]
        model = tf.keras.models.load_model(model_paths[best_model_idx])
        scaler_filename = '/'.join([outdir, 'input_tranformation_weighted.pck'])
        qt = pickle.load(open(scaler_filename, 'rb'))

        features_filename = '/'.join([outdir, 'input_features.pck'])
        features = pickle.load(open(features_filename, 'rb'))

      #  print features

        features = features# +["hnl_mass/7."]

	return model,qt,features

def inputHandler(selection,List,sig,features):
   
	'''
	Function setting signal and background flags + preparing paramMass variable when going for pNN
	''' 

	'''
	 Function setting signal and background flags + preparing paramMass variable when going for pNN
	''' 

  	#hnl mass for mass discretization - signal
#  	MCMass = [1.0,1.5,2.0,3.0,4.5]
  	MCMass = [1.0,1.5,2.0,3.0,4.5]
  	#hnl mass for mass discretization - background
  	#massBin=[0.0,1.3,1.8,2.5,3.5,5.5,7]
	sigma = [0.010,0.011,0.025,0.025,0.035,0.035]
	nSigmas = 10
	
	arrays = []

	print "check list____",List
 	max = 0 
	if sig:
 # 		for i,sample in enumerate(List):
		#build a dataset for each signal hypothesis for either signal or background
		# events in the dataset are selected in a nSignma window around the signal hypothesis
  		arrays.append(root2array(List,"Events",branches= features, selection=selection))

	else:
		array_size = -1
		Multiplier =1 
	  	for i,mass in enumerate(MCMass):
  			arrays.append(root2array(List,"Events",branches= features, selection=selection+" && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i]) ))
			

#  	if len(mcList) ==1:
# 		print "\n"
# 		print "********* PREPARING INPUTS ************"
# 		print "     training on single sample"
# 		print "***************************************"
# 		print "\n"
# 	else:
# 		print "\n"
# 		print "*************************** PREPARING INPUTS **************************************"
# 		print "                training on multiple samples: pNN training"
# 		print "             will build a new feature 'massParam' with values:"
# 		print "       MC truth for signal: sharp value in ", MCMass
# 		print "  Mean over an interval for bkg: discretization bins sizes: ", sigma
# 		print "***********************************************************************************"
# 		print "\n"
# 

 
	A = np.recarray(arrays[0].shape,dtype = arrays[0].dtype.descr)
  	for i,array in enumerate(arrays):
		#dataset processing - add flag + mass parameter
		if not sig:
#		#	print int(array_size)
  	#		array.resize((int(array_size),))
			
#		#	print  len(array)
			array = array[0:int(array_size)]
#		#	print  len(array)
  		new_dt = np.dtype(array.dtype.descr)
		if sig:
  			array_flag = np.ones(len(array))
  			array_flag = np.ones(len(array))
		else:
  			array_flag = np.zeros(len(array))
  			array_flag = np.zeros(len(array))

  		discrete_mass = np.full(len(array),-1.)
  		MCdiscrete = -1 
		if sig:
  		    for mass in MCMass:
  		      	if "_"+str(mass).replace(".","p") in List:
  		    			
					MCdiscrete = mass
  		    			break
		else:
		     MCdiscrete = MCMass[i]
	  	discrete_mass = np.full(len(array),MCdiscrete)
		array = rfn.append_fields(array,['massParam','flag'],[discrete_mass,array_flag])
		A = array if i == 0 else rfn.stack_arrays((A,array),usemask=False)
        		
	return A


def splitSamples(splitArray,samples):

	splitted=[]
	for i in splitArray:
		temp = []
		for path in samples:
		
			if "_"+str(i).replace(".","p") in path:

				temp.append(path) 
	        splitted.append(temp)	

	return splitted

def predictScore(outdir,features,model,qt,sample,sig,selection):
  '''
    Return score with scaled input features
  '''
  print "in predict score",sig
  print(selection)
  x = inputHandler(selection,sample,sig,features)
  x = rfn.stack_arrays((x), usemask = False ) 
  x = np.array(x.tolist(),dtype="float64")
  print x
  if (len(x))==0:
	return 0 
  x = x[x[:,len(features)] != -1.0]
  #root2array(sample,"Events",branches = features)
  #x= np.array(x.tolist(),dtype="float64")
  # apply the scaler
  #print x[0:10]
  xx = qt.transform(x[:,0:(len(features)+1)])

  # predict
  score = model.predict(xx)

  return score

def Add_WPs(model,qt,sig,bkg, outdir, lable,selection):

  MCMass = [1.0,1.5,2.0,3.0,4.5]
  #hnl mass for mass discretization - background
  #massBin=[0.0,1.3,1.8,2.5,3.5,5.5,7]
  sigma = [0.010,0.011,0.025,0.025,0.035,0.035]
  nSigmas = 10

  s_array = []
  b_array = []

  #cut based implementations
  lowMassOS = "(hnl_mass<2 && LepQProd<0 && (\
                                                   (hnl_lxy<1 && hnl_pi_pt>0.7 && hnl_lxy_sig>30 && hnl_cos2D>1-2e-03)\
                                                || ( hnl_lxy>1 && hnl_lxy<5 && hnl_pi_pt>1.5 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-05 )\
                                                || ( hnl_lxy>5 && hnl_pi_pt>0.7 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-04 ))\
                                                ) "
  lowMassSS = "(hnl_mass<2 && LepQProd>0 && (\
                                                   ( hnl_lxy<1 && hnl_pi_pt>0.7 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-04 )\
                                                || ( hnl_lxy>1 && hnl_lxy<5 && hnl_pi_pt>0.7 && hnl_lxy_sig>50 && (1-hnl_cos2D)<2e-03 )\
                                                || ( hnl_lxy>5 && hnl_pi_pt>1.5 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-05 ))\
                                                )"	
  mediumMassOS = "(hnl_mass>2 && hnl_mass<4.5 && LepQProd<0 && (\
                                                    ( hnl_lxy<1 && abs(B_mass-5.3)<.15 && hnl_pi_pt>0.8 && hnl_lxy_sig>90 && (1-hnl_cos2D)<2e-03 )\
                                                 || ( hnl_lxy>1 && hnl_lxy<5 && abs(B_mass-5.3)<.35 && hnl_pi_pt>1.5 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-05 )\
                                                 || ( hnl_lxy>5 && abs(B_mass-5.3)<.2  && hnl_pi_pt>1.5 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-04 ))) "
  mediumMassSS = "(hnl_mass>2 && hnl_mass<4.5 && LepQProd>0 && (\
                                                   ( hnl_lxy<1 && fabs(B_mass-5.3)<.2 && hnl_pi_pt>0.8 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-04 )\
                                                || ( hnl_lxy>1 && hnl_lxy<5  && fabs(B_mass-5.3)<.15 && hnl_pi_pt>0.8 && hnl_lxy_sig>100 && (1-hnl_cos2D)<2e-03 )\
                                                || ( hnl_lxy>5 && fabs(B_mass-5.3)<.35 && hnl_pi_pt>2.5 && hnl_lxy_sig>150 && (1-hnl_cos2D)<2e-05 ))) "
  highMassOS = "( hnl_mass>=4.5 && LepQProd<0 && (\
                                                   ( hnl_lxy<1 && fabs(B_mass-5.3)<.15 && hnl_pi_pt>2.0 && hnl_lxy_sig>30 && (1-hnl_cos2D)<2e-03 )\
                                                || ( hnl_lxy>1 && hnl_lxy<5  && fabs(B_mass-5.3)<.35 && hnl_pi_pt>3.5 && hnl_lxy_sig>50 && (1-hnl_cos2D)<2e-05 )\
                                                || ( hnl_lxy>5 && fabs(B_mass-5.3)<.2  && hnl_pi_pt>3.0 && hnl_lxy_sig>50 && (1-hnl_cos2D)<2e-04 ))) "
  highMassSS = "( hnl_mass>=4.5 && LepQProd>0 && (\
                                                   ( hnl_lxy<1 && fabs(B_mass-5.3)<.2 && hnl_pi_pt>3.0 && hnl_lxy_sig>50 && (1-hnl_cos2D)<2e-04 )\
                                                || ( hnl_lxy>1 && hnl_lxy<5  && fabs(B_mass-5.3)<.15 && hnl_pi_pt>2.0 && hnl_lxy_sig>30 && (1-hnl_cos2D)<2e-03 )\
                                                || ( hnl_lxy>5 && fabs(B_mass-5.3)<.35 &&  hnl_pi_pt>3.5 && hnl_lxy_sig>50 && (1-hnl_cos2D)<2e-05 ))) "

#  selection_cutbased = mediumMassOS
  selection_cutbased = "("+ lowMassSS +"|| "+lowMassSS  +" || "+mediumMassOS +" || "+mediumMassSS +" || "+highMassOS +" || "+highMassSS +")"  

  s_array = root2array(sig,"Events",branches= ["likelihood","B_mass","hnl_pi_pt","hnl_lxy_sig","hnl_cos2D"], selection=selection )
  s_array_cut = root2array(sig,"Events",branches= ["likelihood","B_mass","hnl_pi_pt","hnl_lxy_sig","hnl_cos2D"], selection=selection+" && "+selection_cutbased )
  sig_eff = len(s_array[s_array['likelihood']>0.8])*1.0/len(s_array)
  sig_eff_cutBased = len(s_array_cut)*1.0/len(s_array)

  for i,mass in enumerate(MCMass):
 	if "_"+str(mass).replace(".","p") in sig:
	  	b_array = root2array(bkg,"Events",branches= ["likelihood","B_mass-5.3","hnl_pi_pt","hnl_lxy_sig","1-hnl_cos2D"], selection=selection+" && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i]) )
	  	b_array_cut = root2array(bkg,"Events",branches= ["likelihood","B_mass-5.3","hnl_pi_pt","hnl_lxy_sig","1-hnl_cos2D"], selection=selection+" && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i])+" && "+selection_cutbased )
		bkg_rate = (len(b_array[b_array['likelihood']>0.8])*1.0/len(b_array))
		bkg_rate_cutBased = len(b_array_cut)*1.0/len(b_array)
  print sig_eff,sig_eff_cutBased
  print bkg_rate,bkg_rate_cutBased,len(b_array[b_array['likelihood']>0.8]),len(b_array)
	
  return sig_eff, sig_eff_cutBased, bkg_rate,bkg_rate_cutBased

def superScore(model,qt,sigs,bkg, features,outdir,lable,selection):
  
  temp = [sigs[i].split("/") for  i in range (0, len(sigs))]	   
  lables= [temp[i][3].replace("NewSignature_","") for  i in range (0, len(temp))]
  s_colors = [["coral","salmon","red","darkred"],
              ["darkorange","orange","gold"],				
              ["greenyellow","lawngreen","green"],				
              ["turquoise","aqua","deepskyblue","dodgerblue"],				
              ["blueviolet","mediumorchid","magenta","hotpink"]				
		]
  masses = [1.0,1.5,2.0,3.0,4.5]
  colors = []
  temp = [sigs[i].split("/") for  i in range (0, len(sigs))]	
  lables= [temp[i][3].replace("NewSignature_","") for  i in range (0, len(temp))]
  print lables
  splitter = lable.split("_")
  print splitter
  category = splitter[2]+" "+splitter[3]	
  for ic,mass in enumerate(masses):
	if str(mass).replace(".","p") in lables[0]:
		colors = s_colors[ic]
		print(colors)
 # sig_score = []
  fig = plt.figure(figsize=(8,8),dpi=300)
  ax = fig.add_subplot(111)
  for j,sig in enumerate(sigs): 
 	 sig_score = predictScore(outdir,features,model,qt,sig,True,selection)
	 print (colors[j])
  	 ax.hist(sig_score, bins=np.arange(0,1.025,0.025), alpha=1, label=lables[j], histtype='step', linewidth=2, normed=True, color=colors[j] )
  
  bkg_score =predictScore(outdir,features,model,qt, bkg,False,selection)
  bkg_name=['data-driven background']
  ax.hist(bkg_score, bins=np.arange(0,1.025,0.050), stacked=True, alpha=0.5, label=bkg_name, normed=True, color = "cornflowerblue")
  ax.legend(loc='upper left',prop={'size': 12})
 # ax.set_title("Score distribution of signal and background for testing set", fontsize=20)
  ax.set_xlabel('Score',fontsize=18)
  #plt.ylim(0.1,100)
  #plt.yscale('log')
#  fig.xlabel('False Positive Rate',fontdict=font)
#  fig.yscale('linear')
#  fig.xscale('log')
 # plt.ylabel('True Positive Rate',fontdict=font)
  plt.legend(loc = 'upper left',fontsize='medium',frameon=False)
  axes = fig.gca()
  axes.tick_params(axis='x', labelsize=14)
  axes.tick_params(axis='y', labelsize=14)
  xmin,xmax = axes.get_xlim()
  ymin,ymax = axes.get_ylim()
  plt.text(xmin,ymax*1.02,"CMS ",fontsize=20,fontweight='demibold')
  plt.text(xmin+0.15,ymax*1.02,"Preliminary",fontsize=18, fontstyle='italic')
  plt.text(0.5,ymax*0.95,category,fontsize=15)
  splitter = lable.split("_")
  category = splitter[2]+" "+splitter[3]
#  plt.text(0.1,1.05,category)
  saveFig(outdir,fig, 'score'+lable)
  fig.clear()

def superROC(model,qt,sigs,bkg, features, outdir, lable,selection):
	
	#predict signal outputs
  	b_array = inputHandler(selection,bkg,False,features)
  	#b_array = root2array(bkg,"Events",branches = features)
  	#b_array = np.array(b_array.tolist(),dtype="float64")
        b_array = rfn.stack_arrays((b_array), usemask = False )
        b_array = np.array(b_array.tolist(),dtype="float64")
  	b_array = b_array[b_array[:,len(features)] != -1.0]
  	b_std = qt.transform(b_array[:,0:len(features)+1])

	train_pred_bkg = model.predict(b_std)
	train_pred_bkg = [i[0] for i in train_pred_bkg]	
	y_b = np.zeros(len(train_pred_bkg)) 	

	temp = [sigs[i].split("/") for  i in range (0, len(sigs))]	
	lables= [temp[i][3].replace("NewSignature_","") for  i in range (0, len(temp))]
#	print lables
	splitter = lable.split("_")
#	print splitter
        category = splitter[2]+" "+splitter[3]	
	plt.figure(figsize=(7,7),dpi=300)
	for idx,s in enumerate(sigs):
		
  		array = inputHandler(selection,s,True,features)
        	array = rfn.stack_arrays((array), usemask = False )
       	 	array = np.array(array.tolist(),dtype="float64")
		if len(array)==0:
			continue
  		array = array[array[:,len(features)] != -1.0]
  	#	array = root2array(s,"Events",branches = features)
  	#        array = np.array(array.tolist(),dtype="float64")
  		s_std = qt.transform(array[:,0:len(features)+1])

		train_pred = model.predict(s_std)
		train_pred = [i[0] for i in train_pred]	
		y_s = np.ones(len(train_pred)) 	
		y = np.concatenate((y_s,y_b),axis=0 )
		train_pred = train_pred + train_pred_bkg
	# 	print idx 	
               # l_sigEff,cut_sigEff,l_bkgRate, cut_bkgRate = Add_WPs(model,qt,s,bkg, outdir, lable,selection)
  		fpr, tpr, wps = roc_curve(y, train_pred) 
		print lables[int(idx)]
  		roc = plt.plot(fpr, tpr, label='ROC for signal'+lables[int(idx)], linewidth=2.0)
		#plt.plot(l_bkgRate,l_sigEff,label='L>0.8 signal:'+lables[int(idx)],marker="^", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
		#plt.plot(cut_bkgRate,cut_sigEff,label='cutBased signal:'+lables[int(idx)],marker="o", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
 # 		print("AUC train",sk.metrics.auc(fpr,tpr))
 	plt.xlabel('False Positive Rate',fontdict=font)
 	plt.yscale('linear')
 	plt.xscale('log')
 	plt.ylabel('True Positive Rate',fontdict=font)
 	plt.legend(loc = 'lower right',fontsize='medium',frameon=False)
 	axes = plt.gca()
	axes.tick_params(axis='x', labelsize=14)
	axes.tick_params(axis='y', labelsize=14)
 	xmin,xmax = axes.get_xlim()
 	print xmin
 	plt.text(xmin,1.01,"CMS ",fontsize=20,fontweight='demibold')
 	plt.text(xmin+xmin*5,1.01,"Preliminary",fontsize=18, fontstyle='italic')
 	plt.text(xmin+0.1*xmin,0.9,category,fontsize=15)
			
  	saveFig(outdir,plt,'SuperROC'+lable)
	plt.clf()

def plotter(spec):

        #outdir = model.replace("X",specs[0])
	print spec[0]
        outdir = glob.glob(spec[0])
	print outdir
#	print len(outdir)
	model,qt,features = NNinput(outdir[0])

	split_byMass =splitSamples(masses,samples)
	split_byCtau =splitSamples(ctau,samples)

#	print split_byCtau
#	print split_byMass

	lables = spec[0].split("_")

	for i in range(0,len(split_byMass)):

		superScore(model,qt,split_byMass[i],data, features, outdir[0],"_"+str(masses[i])+"_"+lables[1]+"_"+lables[2],spec[1])
#		superROC(model,qt,split_byMass[i],data, features, outdir[0],"_"+str(masses[i])+"_"+lables[1]+"_"+lables[2],spec[1])
	
#	for i in range(0,len(split_byCtau)):
	#	superScore(model,qt,split_byCtau[i],data, features, outdir[0],"_"+ctau[i]+"_"+lables[1]+"_"+lables[2],spec[1])
#		superROC(model,qt,split_byCtau[i],data, features, outdir[0],"_"+ctau[i]+"_"+lables[1]+"_"+lables[2],spec[1])

if __name__ == '__main__':


	#outdir ="./model/trainNN_ctauParam_19Sep2022_17h20m21s/" 
	#model_h5 = "saved-model-0040_val_loss_0.0287_val_acc_0.9906.h5"
########model = tf.keras.models.load_model(model_path)
########scaler_filename = '/'.join([outdir, 'input_tranformation_weighted.pck'])
########qt = pickle.load(open(scaler_filename, 'rb'))

########features_filename = '/'.join([outdir, 'input_features.pck'])
########features = pickle.load(open(features_filename, 'rb'))

########print features	
	
	model = "model/trainNN_cat_12Oct2022_19h27m11s/"
	specs = [
 		 	["LxySUnder50_OS","hnl_charge==0 && hnl_lxy_sig<50 && LepQProd<0"],
         	 	["LxySUnder50_SS","hnl_charge==0 && hnl_lxy_sig<50 && LepQProd>0"],
         		["LxySOver50Under150_OS","hnl_charge==0 && hnl_lxy_sig>50 && hnl_lxy_sig<150 && LepQProd<0"],
         		["LxySOver50Under150_SS","hnl_charge==0 && hnl_lxy_sig>50 && hnl_lxy_sig<150 && LepQProd>0"],
        		["LxySOver150_OS","hnl_charge==0 && hnl_lxy_sig>150 && LepQProd<0"],
         		["LxySOver150_SS","hnl_charge==0 && hnl_lxy_sig>150 && LepQProd>0"]
			]
	for i,spec in enumerate(specs):
		specs[i] = [model.replace("cat",spec[0]),spec[1]]
	print specs 
#	features = features# +["hnl_mass"]

	masses = [1.0,1.5,2.0,3.0,4.5]
#	masses = [4.5]
	ctau = ["ctau10p0","ctau100p0"]

	samples =[
        	"../data/HNLFlatTuples/NewSignature_1p0_ctau10p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_1p0_ctau100p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_1p0_ctau1000p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_1p5_ctau10p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_1p5_ctau100p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_1p5_ctau1000p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_2p0_ctau10p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_2p0_ctau100p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_2p0_ctau1000p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_3p0_ctau1p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_3p0_ctau10p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_3p0_ctau100p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_3p0_ctau1000p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_4p5_ctau0p1/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_4p5_ctau1p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_4p5_ctau10p0/HNLFlat_0.root",
        	"../data/HNLFlatTuples/NewSignature_4p5_ctau100p0/HNLFlat_0.root"
		]

	data = [
		"../data/HNLFlatTuples/Parking1A0_newSig_4varsL_section0/HNLFlat_*.root",
		"../data/HNLFlatTuples/Parking1A0_newSig_4varsL_section1/HNLFlat_*.root",
		"../data/HNLFlatTuples/Parking1A0_newSig_4varsL_section2/HNLFlat_*.root",
		"../data/HNLFlatTuples/Parking1A0_newSig_4varsL_section3/HNLFlat_*.root",
		]


   	Tpool = Pool(len(specs))
 #   	Tpool = Pool(1)
   	Tpool.map(plotter,specs) 
  	Tpool.close()
  	Tpool.join()
#	plotter(specs[0])

#	split_byMass =splitSamples(masses,samples)
#	split_byCtau =splitSamples(ctau,samples)
#
#	print split_byCtau
#	print split_byMass
#
#
#	for i in range(0,len(split_byMass)):
#
#		#superScore(model,qt,split_byMass[i],data, features, outdir,"_"+str(masses[i]),"hnl_charge==0")
#		superROC(model,qt,split_byMass[i],data, features, outdir,"_"+str(masses[i]),"hnl_charge==0")
#	
#	for i in range(0,len(split_byCtau)):
#		#superScore(model,qt,split_byCtau[i],data, features, outdir,"_"+ctau[i],"hnl_charge==0")
#		superROC(model,qt,split_byCtau[i],data, features, outdir,"_"+ctau[i],"hnl_charge==0")
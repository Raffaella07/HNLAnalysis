# first neural network with keras tutorial

import os
import sys 
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

Bc = 0

def humanLabel(label):
	
	sign = label.split('_')[1]	
	
	if ('gt' in label):

		return 'L_{xy}/\sigma_{xy} > 150, '+sign
	elif ('sig0to50' in label):
		
		return 'L_{xy}/\sigma_{xy} < 50, '+sign
	elif ('sig50to150' in label):
		
		return '50 < L_{xy}/\sigma_{xy} < 150, '+sign
		
def NNinput(outdir):

	Bc = False
	print type(outdir)
	print outdir
	if 'Bc' in outdir:
		Bc = True
        model_paths = glob.glob(outdir+'*.h5')
	outdir = str(outdir)
	print outdir+'*.h5'
	print model_paths
	model_paths_matrix = [path.split("_") for path in model_paths]
	idx = 0 
	if Bc:
		idx =11
	else:
		idx =10
	if len( model_paths_matrix)>1:
		accuracies =[ float(model_paths_matrix[i][idx].strip(".h5")) for i in range (0,len(model_paths_matrix))]
	        best_model_idx = accuracies.index(max(accuracies))
	else:
		best_model_idx = 0
	print best_model_idx,model_paths[best_model_idx]
        model = tf.keras.models.load_model(model_paths[best_model_idx])
        scaler_filename = '/'.join([outdir, 'input_tranformation_weighted.pck'])
        qt = pickle.load(open(scaler_filename, 'rb'))

        features_filename = '/'.join([outdir, 'input_features.pck'])
        features = pickle.load(open(features_filename, 'rb'))

      #  print features

        features = features# +["hnl_mass/7."]

	return model,qt,features

def inputHandler(selection,List,sig,features,do_singlesig):
   
	'''
	Function setting signal and background flags + preparing paramMass variable when going for pNN
	''' 

	'''
	 Function setting signal and background flags + preparing paramMass variable when going for pNN
	''' 

  	#hnl mass for mass discretization - signal
#  	MCMass = [1.0,1.5,2.0,3.0,4.5]
	if Bc:
  		MCMass = [3.0,4.5,5.5]
	else:
  		MCMass = [1.0,1.5,2.0,3.0,4.5]
  	#hnl mass for mass discretization - background
  	#massBin=[0.0,1.3,1.8,2.5,3.5,5.5,7]
	sigma = [  0.0014 +0.0074 * mass for mass in MCMass]
	nSigmas = 10
	
	arrays = []

	print "check list____",List
 	max = 0 
	print 
	if sig:
 # 		for i,sample in enumerate(List):
		#build a dataset for each signal hypothesis for either signal or background
		# events in the dataset are selected in a nSignma window around the signal hypothesis
  		arrays.append(root2array(List,"Events",branches= features, selection=selection))
	else:
		array_size = -1
		Multiplier =1 
	  	for i,mass in enumerate(MCMass):
			print List
			if do_singlesig:
				print 'check ', mass, do_singlesig
				if mass == do_singlesig:
					print selection+" && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i])
  					arrays.append(root2array(List,"Events",branches= features, selection=selection+" && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i]) ))
			else:
				print selection+" && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i])
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
		     if do_singlesig:
		     	MCdiscrete = do_singlesig
		     else:
		     	MCdiscrete = MCMass[i]
			
	  	discrete_mass = np.full(len(array),MCdiscrete)
		array = rfn.append_fields(array,['massParam','flag'],[discrete_mass,array_flag])
		print array[0:10]
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
  x = inputHandler(selection,sample,sig,features,0)
  x = rfn.stack_arrays((x), usemask = False ) 
  x = np.array(x.tolist(),dtype="float64")
  print x
  if (len(x))==0:
	x = [0]
	return x
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
  lables= [temp[i][3].replace("NewSig_ordered","") for  i in range (0, len(temp))]
  s_colors = [["coral","salmon","red","darkred"],
              ["darkorange","orange","gold"],				
              ["greenyellow","lawngreen","green","forestgreen"],				
              ["turquoise","aqua","deepskyblue","dodgerblue"],				
              ["blueviolet","mediumorchid","magenta","hotpink"]				
		]
  if Bc:
	masses = [3.0,4.5,5.5]
  else:
	masses = [1.0,1.5,2.0,3.0,4.5]
  colors = []
  temp = [sigs[i].split("/") for  i in range (0, len(sigs))]	
  lables= [temp[i][3].replace("NewSig_","") for  i in range (0, len(temp))]
  lables= [temp[i][3].replace("ordered","") for  i in range (0, len(lables))]
  print 'labels       ',lables
  splitter = lable.split("_")
 # print splitter
  category = splitter[2]+"_"+splitter[3]	
  for ic,mass in enumerate(masses):
	if str(mass).replace(".","p") in lables[0]:
		colors = s_colors[ic]
		print(colors)
 # sig_score = []
  fig = plt.figure(figsize=(8,8),dpi=300)
  ax = fig.add_subplot(111)
  lables = [l.split('_') for l in lables]
  print 'labels       ',lables
  lables = [l[2].replace('p','.')+' GeV, lifetime '+l[3].strip('ctau').replace('p','.')+' mm'for l in lables]
  print 'labels       ',lables
  for j,sig in enumerate(sigs): 
 	sig_score = predictScore(outdir,features,model,qt,sig,True,selection)
	print len(sig_score)
	print (colors[j])
  	ax.hist(sig_score, bins=np.arange(0,1.025,0.025), alpha=1, label=lables[j], histtype='step', linewidth=2, weights=np.full(len(sig_score),1./len(sig_score)), color=colors[j] )
  
  bkg_score =predictScore(outdir,features,model,qt, bkg,False,selection)
  bkg_name=['data-driven background']
  ax.hist(bkg_score, bins=np.arange(0,1.025,0.025), stacked=True, alpha=0.5, label=bkg_name, weights=np.full(len(bkg_score),1./len(bkg_score)), color = "cornflowerblue")
  ax.legend(loc='upper left',prop={'size': 12})
 # ax.set_title("Score distribution of signal and background for testing set", fontsize=20)
  ax.set_xlabel('Score',fontsize=18)
  #plt.ylim(0.1,100)
  plt.yscale('log')
#  fig.xlabel('False Positive Rate',fontdict=font)
#  fig.yscale('linear')
#  fig.xscale('log')
 # plt.ylabel('True Positive Rate',fontdict=font)
  plt.legend(loc = 'upper center',fontsize='medium',frameon=False)
  axes = fig.gca()
  axes.tick_params(axis='x', labelsize=14)
  axes.tick_params(axis='y', labelsize=14)
  plt.ylim(0.001,1)
  xmin,xmax = axes.get_xlim()
  ymin,ymax = axes.get_ylim()
  plt.text(xmin,ymax*1.1,"CMS ",fontsize=20,fontweight='demibold')
  plt.text(xmin+0.15,ymax*1.1,"Preliminary",fontsize=18,fontstyle='italic')
  plt.text(0.65*xmax,ymax*1.1,humanLabel(category),fontsize=18)
  splitter = lable.split("_")
 # plt.rc('text', usetex=True)
 # plt.rc('font', family='serif')
  category = splitter[2]+" "+splitter[3]
#  plt.text(0.1,1.05,category)
  saveFig(outdir,fig, 'score'+lable)
  fig.clear()

def superROC(model,qt,sigs,bkg, features, outdir, lable,selection):
	
	#predict signal outputs
  	b_array = inputHandler(selection,bkg,False,features,0)
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
#	lables= [temp[i][3].replace("NewSig_withIso","") for  i in range (0, len(temp))]
#	lables = [l.split('_')[0].replace('p','.')+', '+l.split('_')[1].replace('p','.') for l in lables]
  	lables= [temp[i][3].replace("NewSig_","") for  i in range (0, len(temp))]
	lables= [temp[i][3].replace("ordered","") for  i in range (0, len(lables))]
	print lables
	splitter = lable.split("_")
#	print splitter
        category = splitter[2]+"_"+splitter[3]	
	lables = [l.split('_') for l in lables]
	plt.figure(figsize=(7,7),dpi=300)
	lables = [l[2].replace('p','.')+' GeV, lifetime '+l[3].strip('ctau').replace('p','.')+' mm'for l in lables]
	for idx,s in enumerate(sigs):
		
  		array = inputHandler(selection,s,True,features,0)
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
  		roc = plt.plot(fpr, tpr, label='ROC for signal '+lables[int(idx)], linewidth=2.0)
		#plt.plot(l_bkgRate,l_sigEff,label='L>0.8 signal:'+lables[int(idx)],marker="^", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
		#plt.plot(cut_bkgRate,cut_sigEff,label='cutBased signal:'+lables[int(idx)],marker="o", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
 # 		print("AUC train",sk.metrics.auc(fpr,tpr))
 	plt.xlabel('False Positive Rate',fontdict=font)
 	plt.yscale('linear')
 	plt.xscale('log')
 	plt.ylabel('True Positive Rate',fontdict=font)
 	plt.legend(loc = 'upper left',fontsize='medium',frameon=False)
 	axes = plt.gca()
	axes.tick_params(axis='x', labelsize=14)
	axes.tick_params(axis='y', labelsize=14)
 	xmin,xmax = axes.get_xlim()
 	print xmin
 	plt.text(xmin,1.01,"CMS ",fontsize=20,fontweight='demibold')
 	plt.text(xmin+xmin*5,1.01,"Preliminary",fontsize=18, fontstyle='italic')
#  	plt.text(xmin+xmin*5,1.02,"$Preliminary}$",fontsize=18)# fontstyle='italic')
#	plt.rc('text', usetex=True)
#	plt.rc('font', family='serif')
 	plt.text(xmax*0.015,1.02,humanLabel(category),fontsize=18)
			
  	saveFig(outdir,plt,'SuperROC'+lable)
	plt.clf()

def compareROCs(model,model1,qt_,qt1,sig,bkg, feat,feat1, outdir, lable,selection):
	
	models = []
	qt = []
	features = []
	wpidx = []
	models.append(model1)
	models.append(model)
	qt.append(qt1)
	qt.append(qt_)
	features.append(feat1)
	features.append(feat)
	print qt
	temp = sig.split("/") 
#	lables= [temp[i][3].replace("NewSig_withIso","") for  i in range (0, len(temp))]
#	lables = [l.split('_')[0].replace('p','.')+', '+l.split('_')[1].replace('p','.') for l in lables]
  	lables= temp[3].replace("NewSig_","")
	lables= temp[3].replace("ordered","")
	print lables
	lables= lables.split('_')
	mass = lables[2]
	mod_labels = ['reference training','updated training']
	#predict signal outputs
      	
	plt.figure(figsize=(7,7),dpi=300)
	for imod, mod in enumerate(models):
		print float(mass.replace('p','.'))
		wpidx.append(-1)
	  	b_array = inputHandler(selection,bkg,False,features[imod],float(mass.replace('p','.')))
  	#b_array = root2array(bkg,"Events",branches = features)
  	#b_array = np.array(b_array.tolist(),dtype="float64")
		print b_array[0:10]
		print ''
       		b_array = rfn.stack_arrays((b_array), usemask = False )
		print b_array[0:10]
		b_array = np.array(b_array.tolist(),dtype="float64")
 	 	b_array = b_array[b_array[:,len(features[imod])] != -1.0]
  		b_std = qt[imod].transform(b_array[:,0:len(features[imod])+1])

		print b_std[0:10]
		train_pred_bkg = mod.predict(b_std)
		print train_pred_bkg[0:10]
		train_pred_bkg = [i[0] for i in train_pred_bkg]	
		y_b = np.zeros(len(train_pred_bkg)) 	

		splitter = lable.split("_")
#	print splitter
       		category = splitter[2]+"_"+splitter[3]	
  		array = inputHandler(selection,sig,True,features[imod],float(mass.replace('p','.')))
        	array = rfn.stack_arrays((array), usemask = False )
  	 	array = np.array(array.tolist(),dtype="float64")
		if len(array)==0:
			continue
		array = array[array[:,len(features[imod])] != -1.0]
  #	array = root2array(s,"Events",branches = features)
  #	 array = np.array(array.tolist(),dtype="float64")
  		s_std = qt[imod].transform(array[:,0:len(features[imod])+1])

		print s_std[0:10]
		train_pred = mod.predict(s_std)
		print train_pred[0:10]
		train_pred = [i[0] for i in train_pred]	
		y_s = np.ones(len(train_pred)) 	
		y = np.concatenate((y_s,y_b),axis=0 )
		train_pred = train_pred + train_pred_bkg
#	print idx 	
                l_sigEff,cut_sigEff,l_bkgRate, cut_bkgRate = Add_WPs(mod,qt,sig,bkg, outdir, lable,selection)
  		fpr, tpr, wps = roc_curve(y, train_pred) 
  		roc = plt.plot(fpr, tpr, label=mod_labels[imod], linewidth=2.0)
		print wps[0:10]
		print fpr[0:10]
		print tpr[0:10]
		for idx, thrs in enumerate(wps):
			if thrs >= 0.99: wpidx[imod] = idx
       			if thrs < 0.99: break
		print wpidx[imod],fpr[wpidx[imod]], tpr[wpidx[imod]]
		plt.plot(fpr[wpidx[imod]], tpr[wpidx[imod]], '*', label='score > 0.99 '+ mod_labels[imod],markersize=10)
		#plt.plot(l_bkgRate,l_sigEff,label='L>0.8 signal:'+lables[int(idx)],marker="^", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
		#plt.plot(cut_bkgRate,cut_sigEff,label='cutBased signal:'+lables[int(idx)],marker="o", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
 # 		print("AUC train",sk.metrics.auc(fpr,tpr))
 	plt.xlabel('False Positive Rate',fontdict=font)
 	plt.yscale('linear')
 	plt.xscale('log')
	plt.xlim(0.00001,1)
 	plt.ylabel('True Positive Rate',fontdict=font)
 	plt.legend(loc = 'upper left',fontsize='medium',frameon=False)
 	axes = plt.gca()
	axes.tick_params(axis='x', labelsize=14)
	axes.tick_params(axis='y', labelsize=14)
 	xmin,xmax = axes.get_xlim()
 	print xmin
 	plt.text(xmin,1.01,"CMS ",fontsize=20,fontweight='demibold')
 	plt.text(xmin+xmin*5,1.01,"Preliminary",fontsize=18, fontstyle='italic')
#  	plt.text(xmin+xmin*5,1.02,"$Preliminary}$",fontsize=18)# fontstyle='italic')
#	plt.rc('text', usetex=True)
#	plt.rc('font', family='serif')
 	plt.text(xmax*0.015,1.02,humanLabel(category),fontsize=18)
			
  	saveFig(outdir,plt,'ModelROC_'+mass+'_'+lable)
	plt.clf()

def compareNets(model,model_single,qt_,qt_single,feat,feat_single,sigs,bkg,outdir, lable,selection):

	#predict signal outputs
	auc = np.zeros((2,len(sigs)))
	models = []
	qt = []
	features = []
	models.append(model_single)
	models.append(model)
	qt.append(qt_single)
	qt.append(qt_)
	features.append(feat_single)
	features.append(feat)
	for imod, mod in enumerate(models):
  		b_array = inputHandler(selection,bkg,False,features[imod],0)
  	#b_array = root2array(bkg,"Events",branches = features)
  	#b_array = np.array(b_array.tolist(),dtype="float64")
		b_array = rfn.stack_arrays((b_array), usemask = False )
		b_array = np.array(b_array.tolist(),dtype="float64")
	  	b_array = b_array[b_array[:,len(features[imod])] != -1.0]
  		b_std = qt[imod].transform(b_array[:,0:len(features[imod])+1])

		train_pred_bkg = models[imod].predict(b_std)
		train_pred_bkg = [i[0] for i in train_pred_bkg]	
		y_b = np.zeros(len(train_pred_bkg)) 	

		temp = [sigs[i].split("/") for  i in range (0, len(sigs))]	
		lables= [temp[i][3].replace("NewSig_ordered","") for  i in range (0, len(temp))]
		lables = [l.split('_')for l in lables]
		masses = [float(l[2].replace('p','.')) for l in lables]
		print masses
		splitter = lable.split("_")
		print 'splitter',splitter
		#masses.append(splitter[])
		category = splitter[2]+"_"+splitter[3]	
		plt.figure(figsize=(7,7),dpi=300)
		lables = [l[2].replace('p','.')+' GeV, lifetime '+l[3].strip('ctau').replace('p','.')+' mm'for l in lables]
		for idx,s in enumerate(sigs):
  			array = inputHandler(selection,s,True,features[imod],0)
        		array = rfn.stack_arrays((array), usemask = False )
       	 		array = np.array(array.tolist(),dtype="float64")
			if len(array)==0:
				continue
	  		array = array[array[:,len(features[imod])] != -1.0]
	  	#	array = root2array(s,"Events",branches = features)
	  	#        array = np.array(array.tolist(),dtype="float64")
	  		s_std = qt[imod].transform(array[:,0:len(features[imod])+1])
	
			train_pred = models[imod].predict(s_std)
			train_pred = [i[0] for i in train_pred]	
			y_s = np.ones(len(train_pred)) 	
			y = np.concatenate((y_s,y_b),axis=0 )
			train_pred = train_pred + train_pred_bkg
		# 	print idx 	
	               # l_sigEff,cut_sigEff,l_bkgRate, cut_bkgRate = Add_WPs(model,qt,s,bkg, outdir, lable,selection)
	  		fpr, tpr, wps = roc_curve(y, train_pred)
			auc[imod][idx] = sk.metrics.auc(fpr,tpr) 
			print imod, idx,auc[imod][idx]
#			print lables[int(idx)]
			print 'ROC for signal'+lables[int(idx)]+str(imod)
	  		roc = plt.plot(fpr, tpr, label='ROC for signal'+lables[int(idx)]+str(imod), linewidth=2.0)
			#plt.plot(l_bkgRate,l_sigEff,label='L>0.8 signal:'+lables[int(idx)],marker="^", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
			#plt.plot(cut_bkgRate,cut_sigEff,label='cutBased signal:'+lables[int(idx)],marker="o", markersize=7, markeredgecolor=roc[0].get_color(), markerfacecolor=roc[0].get_color())
	 # 		print("AUC train",sk.metrics.auc(fpr,tpr))
 	plt.xlabel('False Positive Rate',fontdict=font)
 	plt.yscale('linear')
 	plt.xscale('log')
 	plt.ylabel('True Positive Rate',fontdict=font)
 	plt.legend(loc = 'upper left',fontsize='medium',frameon=False)
 	axes = plt.gca()
 	axes.tick_params(axis='x', labelsize=14)
 	axes.tick_params(axis='y', labelsize=14)
 	xmin,xmax = axes.get_xlim()
 	print xmin
 	plt.text(xmin,1.01,"CMS ",fontsize=20,fontweight='demibold')
 	plt.text(xmin+xmin*5,1.01,"Preliminary",fontsize=18, fontstyle='italic')
#  	plt.text(xmin+xmin*5,1.02,"$Preliminary}$",fontsize=18)# fontstyle='italic')
#	plt.rc('text', usetex=True)
#	plt.rc('font', family='serif')
 	plt.text(xmax*0.015,1.02,humanLabel(category),fontsize=18)
			
#  	saveFig(outdir,plt,'ModelROC'+lable)
			
	aucs_NN = ROOT.TGraph(len(sigs),np.asarray(masses).astype('float'),np.asarray(auc[0]).astype('float'))
	aucs_pNN = ROOT.TGraph(len(sigs),np.asarray(masses).astype('float'),np.asarray(auc[1]).astype('float'))
	aucs_NNtrain = ROOT.TGraph(1,np.asarray([3]).astype('float'),np.asarray(auc[0][3]).astype('float'))

	aucs_NN.SetTitle('pNN - hits info')
	aucs_pNN.SetTitle('pNN - pi Hits info only')
	aucs_NNtrain.SetTitle('NN - training signal sample')


	aucs_NN.SetMarkerStyle(8)
	aucs_NNtrain.SetMarkerStyle(8)
	aucs_pNN.SetMarkerStyle(8)
	aucs_NN.SetMarkerSize(1.3)
	aucs_NNtrain.SetMarkerSize(1.3)
	aucs_pNN.SetMarkerSize(1.3)
	aucs_NN.SetLineStyle(8)
	aucs_pNN.SetLineStyle(8)
	aucs_NN.SetMarkerColor(9)
	aucs_NNtrain.SetMarkerColor(9)
	aucs_pNN.SetMarkerColor(46)
	aucs_NN.SetLineColor(9)
	aucs_pNN.SetLineColor(46)
	aucs_NN.SetLineWidth(2)
	aucs_pNN.SetLineWidth(2)

	multi = ROOT.TMultiGraph()
	multi.Add(aucs_NN,"PL")
	multi.Add(aucs_pNN,"PL")
#	multi.Add(aucs_NNtrain,"PL")
	ROOT.gStyle.SetLegendBorderSize(0)
	ROOT.gStyle.SetLegendFillColor(0)
	c = ROOT.TCanvas('c','c',600,600)
	multi.Draw("APL")
	multi.GetXaxis().SetTitle('HNL mass (GeV)')
	multi.GetYaxis().SetTitle('AUC')
	multi.GetYaxis().SetRangeUser(0.8,1.0)
	c.BuildLegend(0.15,0.2,0.75,0.45,humanLabel(category))
	cms = ROOT.TLatex()
	cms.DrawLatex(1,1.03,'CMS Preliminary')
	c.SaveAs(outdir+'AUC_'+lable+'.png')
	c.SaveAs(outdir+'AUC_'+lable+'.pdf')

#	superROC(model,qt,split_byCtau[i],data, features, outdir[0],"_"+ctau[i]+"_"+lables[1]+"_"+lables[2],spec[1])
#	plt.xlabel('False Positive Rate',fontdict=font)
# 	plt.yscale('linear')
# 	plt.xscale('log')
# 	plt.ylabel('True Positive Rate',fontdict=font)
# 	plt.legend(loc = 'lower right',fontsize='medium',frameon=False)
# 	axes = plt.gca()
#	axes.tick_params(axis='x', labelsize=14)
#	axes.tick_params(axis='y', labelsize=14)
# 	xmin,xmax = axes.get_xlim()
# 	print xmin
# 	plt.text(xmin,1.01,"CMS ",fontsize=20,fontweight='demibold')
# 	plt.text(xmin+xmin*5,1.01,"Preliminary",fontsize=18, fontstyle='italic')
# 	plt.text(xmin+0.1*xmin,0.9,category,fontsize=15)
			
 # 	saveFig(outdir,plt,'SuperROC'+lable)
#	plt.clf()
	

def plotDistr(model,qt,feat,sig,bkg,outdir,mass, lable,selection):
		print feat

		feat.append('hnl_pi_PixelLayers')	
  		b_array = inputHandler(selection,bkg,False,feat,mass)
		b_array = rfn.stack_arrays((b_array), usemask = False )
		b_array = np.array(b_array.tolist(),dtype="float64")
	  	b_topred = b_array[b_array[:,len(feat)-1] != -1.0]
  		b_std = qt.transform(b_topred[:,0:len(feat)])

		train_pred_bkg = model.predict(b_std)
		train_pred_bkg = [i[0] for i in train_pred_bkg]	
  		s_array = inputHandler(selection,sig,False,feat,mass)
		s_array = rfn.stack_arrays((s_array), usemask = False )
		s_array = np.array(s_array.tolist(),dtype="float64")
	  	s_topred = s_array[s_array[:,len(feat)-1] != -1.0]
  		s_std = qt.transform(s_topred[:,0:len(feat)])

		train_pred= model.predict(s_std)
		train_pred = [i[0] for i in train_pred]
		print 'dimensions  ',len(train_pred),len(s_std)	
		s_array = np.column_stack((s_array,train_pred))
		b_array = np.column_stack((b_array,train_pred_bkg))
		
		s_sel = s_array[s_array[:,len(s_array[0])-1]>0.99]
		b_sel = b_array[b_array[:,len(b_array[0])-1]>0.99]

		print s_sel[0:10]
		print b_sel[0:10]
		print s_sel[0:10,len(s_array[0])-4]

  		fig = plt.figure(figsize=(8,8),dpi=300)
		ax = fig.add_subplot(111)
		#px_sig
		
		ax.hist(b_sel[len(s_array[0])-4], bins=np.arange(0,7,1), stacked=True, alpha=0.4, label='10 #sigma window background', weights=np.full(len(b_sel[len(s_array[0])-4]),1./len(b_sel[len(s_array[0])-4])), color = "cornflowerblue")	
		ax.hist(s_sel[len(s_array[0])-4], bins=np.arange(0,7,1), stacked=True, alpha=1, label=lable,histtype='step', weights=np.full(len(s_sel[len(s_array[0])-4]),1./len(s_sel[len(s_array[0])-4])), linewidth = 2,color = "magenta")	

  		ax.legend(loc='upper left',prop={'size': 12})
 # ax.set_title("Score distribution of signal and background for testing set", fontsize=20)
 		ax.set_xlabel('pi nPixelLayers',fontsize=18)
  #plt.ylim(0.1,100)
#  plt.yscale('log')
#  fig.xlabel('False Positive Rate',fontdict=font)
#  fig.yscale('linear')
#  fig.xscale('log')
 # plt.ylabel('True Positive Rate',fontdict=font)
		plt.legend(loc = 'upper center',fontsize='medium',frameon=False)
 		axes = fig.gca()
		axes.tick_params(axis='x', labelsize=14)
		axes.tick_params(axis='y', labelsize=14)
		plt.ylim(0.001,1)
		xmin,xmax = axes.get_xlim()
		ymin,ymax = axes.get_ylim()
		splitter = lable.split("_")
		category = splitter[2]+"_"+splitter[3]
		plt.text(xmin,ymax*1.1,"CMS ",fontsize=20,fontweight='demibold')
		plt.text(xmin+0.15,ymax*1.1,"Preliminary",fontsize=18,fontstyle='italic')
		plt.text(0.65*xmax,ymax*1.1,humanLabel(category),fontsize=18)
 # plt.rc('text', usetex=True)
 # plt.rc('font', family='serif')
#  plt.text(0.1,1.05,category)
 		saveFig(outdir[0],fig, 'pixel'+lable)

def plotter(spec):

        #outdir = model.replace("X",specs[0])
	print spec[0]
	label = spec[0].split('_')[1]+'_'+spec[0].split('_')[2]
	print label
	outdir_single =''
	for s in single:
		if label in s[0]:
		
        		outdir_single = s[0]
        		sel_single = s[1]
		
#	print single[0]
        outdir = glob.glob(spec[0])
#	print outdir

	print 'outdir single', outdir_single
	print outdir
	model,qt,features = NNinput(outdir[0])
	model_single,qt_single,features_single = NNinput(outdir_single)
	if Bc:
		split_byMass =splitSamples(masses,samples_Bc)
		split_byCtau =splitSamples(ctau,samples_Bc)

	else:
		split_byMass =splitSamples(masses,samples)
		split_byCtau =splitSamples(ctau,samples)

	print split_byCtau
#	print split_byMass

	lables = spec[0].split("_")
#	for i in range(0,len(split_byMass)):
#
#		print masses[i]
#		superScore(model,qt,split_byMass[i],data, features, outdir[0],"_"+str(masses[i])+"_"+lables[1]+"_"+lables[2],spec[1])
#		superROC(model,qt,split_byMass[i],data, features, outdir[0],"_"+str(masses[i])+"_"+lables[1]+"_"+lables[2],spec[1])

#	compareNets(model,model_single,qt,qt_single,features,features_single,split_byCtau[1],data,outdir_single,"_"+ctau[1]+"_"+lables[1]+"_"+lables[2],sel_single)
	for sample in comparison_samples:
		sig_lable =sample.split("/") 
		sig_lable[3].replace('NewSig_hits_','')
	        plotDistr(model,qt,features,sample,data,outdir,4.5,sig_lable[3],sel_single) 
#		sig_lable[3].replace('NewSig_hits_','')
#		compareROCs(model,model_single,qt,qt_single,sample,data,features,features_single,outdir_single,"_"+sig_lable[3]+"_"+lables[1]+"_"+lables[2],sel_single)
			
#	for i in range(0,len(split_byCtau)):
	#	superScore(model,qt,split_byCtau[i],data, features, outdir[0],"_"+ctau[i]+"_"+lables[1]+"_"+lables[2],spec[1])
#		superROC(model,qt,split_byCtau[i],data, features, outdir[0],"_"+ctau[i]+"_"+lables[1]+"_"+lables[2],spec[1])

if __name__ == '__main__':

#	features = ["hnl_pi_pt","hnl_l_pt","TrgMu_pt","B_mass","hnl_cos2D","hnl_lxy_sig","hnl_vtxProb","hnl_vtxChi2","B_pt","dilepton_mass","BlepPi_mass","dr_trgMu_lep","dr_Blep_pi","TrgMu_relIso","hnl_l_relIso"]

	#outdir ="./model/trainNN_ctauParam_19Sep2022_17h20m21s/" 
	#model_h5 = "saved-model-0040_val_loss_0.0287_val_acc_0.9906.h5"
########model = tf.keras.models.load_model(model_path)
########scaler_filename = '/'.join([outdir, 'input_tranformation_weighted.pck'])
########qt = pickle.load(open(scaler_filename, 'rb'))

########features_filename = '/'.join([outdir, 'input_features.pck'])
########features = pickle.load(open(features_filename, 'rb'))

########print features
	Bc = np.float(sys.argv[1])
#	model = "model/trainNN_cat_24Jan2023_17h26m30s/"
#	model = "model/trainNN_cat_19Jul2023_23h13m32s/" new Bc
#	model = "model/trainNN_cat_19Jul2023_18h24m56s/"
#	model = "model/trainNN_cat_26Jul2023_17h22m57s/" # iso + hits - no chi2 
#	model = "model/trainNN_cat_21Jul2023_09h06m30s/"
	model = "model/trainNN_cat_12Oct2022_15h25m05s/"
#	model_single = "model/trainNN_cat_30Jun2023_02h56m27s/"
#	model_single = "model/trainNN_cat_21Jul2023_11h03m06s/"
#	model_single = "model/trainNN_cat_28Jul2023_09h12m33s/" #no iso
	model_single = "model/trainNN_cat_30Jul2023_12h16m53s/" #no iso, ele lost hits
#	model_single = "model/trainNN_cat_12Oct2022_15h25m05s/"
#	model_single = "model/trainNN_cat_29Jun2023_18h30m17s/"
#	model_single = "model/trainNN_cat_26Jan2023_00h25m43s/" old Bc
	specs = [
 		 	["lxysig0to50_OS","hnl_charge==0 && hnl_lxy_sig<50 && LepQProd<0 && fabs(hnl_l_eta)<1.442 && hnl_l_mvaId>-3 && dr_Blep_pi>0.01"],
         	 	["lxysig0to50_SS","hnl_charge==0 && hnl_lxy_sig<50 && LepQProd>0 && fabs(hnl_l_eta)<1.442 && hnl_l_mvaId>-3 && dr_Blep_pi>0.01"],
         		["lxysig50to150_OS","hnl_charge==0 && hnl_lxy_sig>50 && hnl_lxy_sig<150 && LepQProd<0 && fabs(hnl_l_eta)<1.442 && hnl_l_mvaId>-3 && dr_Blep_pi>0.01"],
         		["lxysig50to150_SS","hnl_charge==0 && hnl_lxy_sig>50 && hnl_lxy_sig<150 && LepQProd>0 && fabs(hnl_l_eta)<1.442 && hnl_l_mvaId>-3 && dr_Blep_pi>0.01"],
        		["lxysiggt150_OS","hnl_charge==0 && hnl_lxy_sig>150 && LepQProd<0 && fabs(hnl_l_eta)<1.442 && hnl_l_mvaId>-3 "],
         		["lxysiggt150_SS","hnl_charge==0 && hnl_lxy_sig>150 && LepQProd>0 && fabs(hnl_l_eta)<1.442 && hnl_l_mvaId>-3 && dr_Blep_pi>0.01"]
			]
	
  	Bc_sel = "&& B_mass<5.7"
  	Bc_string = ""
	single = []
  	if Bc:
		Bc_sel=" && B_mass>5.7" 
		Bc_string="_Bc" 
	for i,spec in enumerate(specs):
		
		specs[i] = [model.replace("cat",spec[0]+Bc_string),spec[1]+Bc_sel]
		single.append([model_single.replace("cat",spec[0]+Bc_string),spec[1]+Bc_sel])
	print specs 
#	features = features# +["hnl_mass"]

	masses = []
	if Bc:
		masses = [3.0,4.5,5.5]
	else:
		masses = [1.0,1.5,2.0,3.0,4.5]
#	masses = [4.5]
	ctau = ["ctau10p0","ctau100p0"]

	samples = [
		"../data/HNLFlatTuples/NewSig_hits_4p5_ctau0p1/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_4p5_ctau1p0/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_3p0_ctau1p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_1p0_ctau10p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_1p5_ctau10p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_2p0_ctau10p0/HNLFlat_0.root",
	        "../data/HNLFlatTuples/NewSig_hits_3p0_ctau10p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_4p5_ctau10p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_1p0_ctau100p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_1p5_ctau100p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_2p0_ctau100p0/HNLFlat_0.root",
	        "../data/HNLFlatTuples/NewSig_hits_3p0_ctau100p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_4p5_ctau100p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_1p0_ctau1000p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_1p5_ctau1000p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_2p0_ctau1000p0/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_3p0_ctau1000p0/HNLFlat_0.root"
		]
	samples_Bc = [
		"../data/HNLFlatTuples/NewSig_hits_5p5_ctau0p01_Bc/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_5p5_ctau0p1_Bc/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_4p5_ctau0p1_Bc/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_5p5_ctau1p0_Bc/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_4p5_ctau1p0_Bc/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_3p0_ctau1p0_Bc/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_5p5_ctau10p0_Bc/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_4p5_ctau10p0_Bc/HNLFlat_0.root",
	        "../data/HNLFlatTuples/NewSig_hits_3p0_ctau10p0_Bc/HNLFlat_0.root",
	        "../data/HNLFlatTuples/NewSig_hits_3p0_ctau100p0_Bc/HNLFlat_0.root",
	        "../data/HNLFlatTuples/NewSig_hits_3p0_ctau1000p0_Bc/HNLFlat_0.root"
	]
	comparison_samples = [
	#	"../data/HNLFlatTuples/NewSig_hits_1p0_ctau1000p0/HNLFlat_0.root",
#	 	"../data/HNLFlatTuples/NewSig_hits_1p5_ctau1000p0/HNLFlat_0.root",
	# 	"../data/HNLFlatTuples/NewSig_hits_2p0_ctau100p0/HNLFlat_0.root",
	#        "../data/HNLFlatTuples/NewSig_hits_3p0_ctau10p0/HNLFlat_0.root",
	 	"../data/HNLFlatTuples/NewSig_hits_4p5_ctau0p1/HNLFlat_0.root",
	]
	data = [
	"../data/HNLFlatTuples/Parking_hits_1D_0_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_0_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_0_2/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_1_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_1_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_13_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_13_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_14_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_14_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_14_2/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_2_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_2_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_3_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_3_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_4_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_4_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_5_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_5_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_6_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_6_1/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_7_0/HNL*.root",
	"../data/HNLFlatTuples/Parking_hits_1D_7_1/HNL*.root",
#	"../data/HNLFlatTuples/Parking_1D_1_1_Bmass/HNLFlat_*1.root",
#	"../data/HNLFlatTuples/Parking_1D_2_0_Bmass/HNLFlat_*1.root",
#	"../data/HNLFlatTuples/Parking_1D_2_1_Bmass/HNLFlat_*1.root",
#	"../data/HNLFlatTuples/Parking_1D_3_0_Bmass/HNLFlat_*1.root",
		]
	plotter(specs[5])

#   	Tpool = Pool(len(specs))
    	#Tpool = Pool(1)
  # 	Tpool.map(plotter,specs) 
 #	Tpool.close()
  #	Tpool.join()

#	split_byMass =splitSamples(masses,samples)
#	split_byCtau =splitSamples(ctau,samples)
#
#	print split_byCtau
#	print split_byMass
#
#
#	for i in range(0,len(split_byMass)):
#
		#superScore(model,qt,split_byMass[i],data, features, outdir,"_"+str(masses[i]),"hnl_charge==0")
#		superROC(model,qt,split_byMass[i],data, features, outdir,"_"+str(masses[i]),"hnl_charge==0")
#	
#	for i in range(0,len(split_byCtau)):

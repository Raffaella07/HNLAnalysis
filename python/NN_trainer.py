# first neural network with keras tutorial

import os
import random
import sys
from random import uniform
from numpy import loadtxt
import numpy as np
from datetime import datetime
from os import path

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
from keras.optimizers import SGD, Adam

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
import multiprocessing
from multiprocessing import Process,Value,Array, Pool
from multiprocessing.pool import ThreadPool as threadPool


class Trainer(object):

  def __init__(self, features,n_features, epochs, batch_size,signal_size,bkg_size,learning_rate, scaler_type, do_early_stopping, do_reduce_lr, dirname, baseline_selection,train_samples):
    self.features = features
    self.n_features = n_features
    self.epochs = epochs
    self.batch_size = batch_size
    self.signal_size = signal_size
    self.bkg_size = bkg_size
    self.learning_rate = learning_rate
    self.min_mass_size = bkg_size
    self.training_size = bkg_size
    self.valid_size = -1
    self.scaler_type = scaler_type
    self.do_early_stopping = do_early_stopping
    self.do_reduce_lr = do_reduce_lr
    self.dirname = dirname + '_' + datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss')

    self.baseline_selection = baseline_selection
    self.train_samples = train_samples

  #  self.tag = tag

    self.target_branch = 'flag'

  def createOutDir(self):
     '''
       This function creates the output directory
     '''
     outdir = './model/{}'.format(self.dirname)
     if not path.exists(outdir):
       os.system('mkdir -p {}'.format(outdir))
     return outdir

  def inputHandler(self,mcList,bkgList,sig):
     
	'''
	 Function setting signal and background flags + preparing paramMass variable when going for pNN
	''' 

  	#hnl mass for mass discretization - signal
#  	MCMass = [1.0,1.5,2.0,3.0,4.5]
	MCMass = []
	if Bc:
  		MCMass = [3.0,4.5,5.5]
	else:
  		MCMass = [1.0,1.5,2.0,3.0,4.5]
  	#	MCMass = [3.0]
  	#hnl mass for mass discretization - background
  	#massBin=[0.0,1.3,1.8,2.5,3.5,5.5,7]
	h_sig = [ROOT.TH1D('sig'+str(m),'sig'+str(m),100,0.6,6) for m in MCMass]
	print h_sig	
	h_bkg = [ROOT.TH1D('bkg'+str(m),'bkg'+str(m),100,0.6,6) for m in MCMass]	
	full = ROOT.TH1D('bkg','bkg',100,0.6,6)
	sigma = [0.0014 +0.0074 *mass for mass in MCMass]
	nSigmas = 10
	feat_hist_sig = []
	feat_hist_bkg = []
	n_ctau = np.zeros(len(MCMass))
	sig_len = np.zeros(len(MCMass))
	per_mass = np.zeros(len(MCMass))
	temp = [mcList[i].split("/") for  i in range (0, len(mcList))]	
	self.train_samples= [temp[i][3].replace("NewSignature_","") for  i in range (0, len(temp))]
	#	print self.train_samples
		 
	arrays = []
	sig_mass =[[]]*5
	sigA =[] 
	sigMass =[] 
	bkg_mass = []
	bkg_mass_full = []
	if sig:
		print mcList
		for mIdx,mass in enumerate(MCMass):
			for i,sample in enumerate(mcList):
  		         	if "_"+str(mass).replace(".","p") in sample and sig:
		#build a dataset for each signal hypothesis for either signal or background
		# events in the dataset are selected in a nSignma window around the signal hypothesis
		#build a dataset for each signal hypothesis for either signal or background
		# events in the dataset are selected in a nSignma window around the signal hypothesis
  					arrays.append(root2array(sample,"Events",branches= self.features, selection=self.baseline_selection))
  					sigMass.append(root2array(sample,"Events",branches=["hnl_mass"], selection=self.baseline_selection))
					print sample, len(arrays),mass, self.signal_size
					print ''
					if mass == 4.5:
						print mIdx, i, len(arrays)
						print arrays[len(arrays)-1][0:10]["B_mass"] 
			#		if len(arrays[i]) < self.signal_size:
			sig_mass[mIdx] = sigMass	
					
		print(len(arrays))	
  		for idx,sample in enumerate(mcList):
       		 	for j,m in enumerate(MCMass):
				print 'in shape counter'
#				print sample, m 
        			if "_"+str(m).replace(".","p")+"_" in sample:
        				sig_len[j] = sig_len[j]+len(arrays[idx])
					n_ctau[j] +=1
       		 		print "before reshaping",sample,arrays[idx].shape

      		print sig_len,n_ctau
		self.min_mass_size = min(sig_len)
		for idx in range(0,len(arrays)):
        		for j,m in enumerate(MCMass):
        			if "_"+str(m).replace(".","p")+"_" in sample:
	 		 		if len(arrays[idx]) >self.min_mass_size/n_ctau[j]:
        					arrays[idx].resize((int(self.min_mass_size/n_ctau[j]),), refcheck=True)
					for i_m in range(0,len(sig_mass[j])):		
						
	 		 			if len(sig_mass[j][i_m]) >self.min_mass_size/n_ctau[j]:
        						sig_mass[j][i_m].resize((int(self.min_mass_size/n_ctau[j]),), refcheck=True)
						#	print "__________________",sample,len(sig_mass[j][i_m])
      		#print arrays.fields()
 #     		print arrays[0].fields()
  #  		print arrays[0][0].fields()
		count = []
		for mIdx,mass in enumerate(MCMass):
			temp = []
			temp_mass= []
			for i,sample in enumerate(mcList):
  		         	if "_"+str(mass).replace(".","p") in sample and sig:
						count.append([0])
						print mIdx,i, len(count)
						if len(temp) == 0:
							temp = arrays[len(count)-1]	
							temp_mass = sigMass[len(count)-1]	
						else:
							temp = rfn.stack_arrays((temp,arrays[len(count)-1]),usemask=False)
							temp_mass = rfn.stack_arrays((temp_mass,sigMass[len(count)-1]),usemask=False)
					#	break
						#np.concatenate(temp, arrays[mIdx])	
			#arr = np.recarray()
	#		if mIdx==0:
	#			arr = sig_mass[mIdx]
	#		else:
			sig_feat = []
			for i_f,feat in enumerate(features):
				hsig = ROOT.TH1D('sig '+feat,'sig '+feat,feat_binning[i_f][0],feat_binning[i_f][1],feat_binning[i_f][2])
#				print '______________', mass
#				print '______________', temp[feat][0:10], feat
							#arr = rfn.stack_arrays((temp[:][feat][:]),usemask=False)
							#arr = rfn.stack_arrays((arrays[mIdx][:][feat]),usemask=False)
				fill_hist(hsig, temp[feat], weights=None, return_indices=False)
				sig_feat.append(hsig)
			feat_hist_sig.append(sig_feat)
					#	print arr[0:10]
					#	print ""
			sigA.append(temp_mass)


		per_mass = np.zeros(len(MCMass))
	  	for i,sample in enumerate(mcList):
        		for j,m in enumerate(MCMass):
        			if "_"+str(m).replace(".","p")+"_" in sample:
        				per_mass[j] = per_mass[j]+len(arrays[i])

		print "after reshaping for similar per mass signal input ", per_mass
		self.signal_size = per_mass.sum()	
		#N.B  signal.size parameter is set in signal input and used in background - homogeneus sized background wrt to signal all over the mass hypothesis
	else:
		List = bkgList 
		array_size = -1
		Multiplier =1
		print "reshape background" 
		print self.signal_size
		self.bkg_size = 0
  		bkg_mass_full = (root2array(bkgList,"Events",branches=['hnl_mass'], selection=self.baseline_selection+" && hnl_mass>0 && hnl_mass<6"))
	  	for i,mass in enumerate(MCMass):
  			arrays.append(root2array(bkgList,"Events",branches= self.features, selection=self.baseline_selection+" && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i])))
			feat_bkg = []
  			bkg_mass.append(root2array(bkgList,"Events",branches=['hnl_mass'], selection=self.baseline_selection+"  && hnl_mass> "+str(mass-nSigmas*sigma[i])+" && hnl_mass<"+str(mass+nSigmas*sigma[i])))
			print arrays[i].shape
			if (i==4 and self.signal_size < len(arrays[i])): 
				scale = self.min_mass_size/len(arrays[i])
			else: 
				scale =1
			if (self.signal_size < len(arrays[i])):
				arrays[i].resize((int(self.min_mass_size)), refcheck=True)
				bkg_mass[i].resize((int(self.min_mass_size)), refcheck=True)
			print arrays[i].shape
			for i_f,feat in enumerate(features):
				hbkg = ROOT.TH1D('bkg_'+feat,'bkg_'+feat,feat_binning[i_f][0],feat_binning[i_f][1],feat_binning[i_f][2])
				arr = rfn.stack_arrays((arrays[i][:][feat]),usemask=False)
				fill_hist(hbkg, arr, weights=None, return_indices=False)
				feat_bkg.append(hbkg)
			feat_hist_bkg.append(feat_bkg)	
			self.bkg_size = self.bkg_size + len(arrays[i])
	print self.bkg_size
	#dataset processing - add class flag + mass parameter
	A = np.recarray(arrays[0].shape,dtype = arrays[0].dtype.descr)
	if sig:
		for i,a in enumerate(sigA):
#			print a
			fill_hist(h_sig[i], a, weights=None, return_indices=False)
			print h_sig[i].GetMean()
			h_sig[i].SetTitle('mass parameter {:.1f}'.format(MCMass[i]))
			h_sig[i].SetName('mass_parameter'.format(MCMass[i]))
		c_sig = ROOT.TCanvas('c','c',800,600)
		ROOT.gStyle.SetPalette(91)
		ROOT.gStyle.SetLegendBorderSize(0)
		ROOT.gStyle.SetLegendFillColor(0)
		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetOptTitle(0)
		h_sig[0].SaveAs('testHist.root')
		h_sig[0].GetXaxis().SetTitle("m_{N}(GeV)")
		h_sig[0].GetYaxis().SetTitle("entries")
		h_sig[0].SetLineWidth(2)
		h_sig[0].SetMarkerSize(0)
		h_sig[0].Draw("hist plc pfc")
		for i in range(1,len(MCMass)):
		#	if i==0:
#			else:
		
			h_sig[i].SetLineWidth(2)
			h_sig[i].SetMarkerSize(0)
			h_sig[i].Draw("same hist plc pfc")
		c_sig.BuildLegend(0.58,0.58,0.88,0.88,"","F")
		c_sig.RedrawAxis()
		c_sig.SaveAs("sig_par.pdf")		
	else:
		c_bkg = ROOT.TCanvas('c','c',800,600)
		fill_hist(full, bkg_mass_full, weights=None, return_indices=False)
		full.Scale(scale)
		full.GetXaxis().SetTitle("m_{N}(GeV)")
		full.SetTitle('not used in training')
		full.SetName('not used in training')
		full.GetXaxis().SetTitle("m_{N}(GeV)")
		full.GetYaxis().SetRangeUser(0.0,full.GetMaximum()*1.1)
		full.SetLineColor(ROOT.kGray+1)
		full.SetFillColor(ROOT.kGray)
		full.SetFillStyle(3004)
		full.SetMarkerSize(0)
		full.GetYaxis().SetTitle("entries")
		full.Draw("hist")
	
		for i,a in enumerate(bkg_mass):
			
			fill_hist(h_bkg[i], a, weights=None, return_indices=False)

			h_bkg[i].SetTitle('mass parameter {:.1f}'.format(MCMass[i]))
			h_bkg[i].SetName('mass parameter {:.1f}'.format(MCMass[i]))
			h_bkg[i].SetMarkerSize(0)
			h_bkg[i].Draw('same hist PLC pfc')
	
		c_bkg.BuildLegend(0.58,0.58,0.88,0.88,"","F")
		c_bkg.SaveAs("sig_par_bkg.pdf")		
	
	idx_discrete = []

	if sig:
	       for mass in MCMass:
			for mc in mcList:
	        	 	if "_"+str(mass).replace(".","p") in mc:
					print 'matcher___________________', mc, mass
	       				idx_discrete.append(mass)
	       				
	print idx_discrete
  	for i,array in enumerate(arrays):
		if not sig:
			print int(array_size)
  	#		array.resize((int(array_size),))
			
			print  len(array)
			array = array[0:int(array_size)]
			print  len(array)
  		new_dt = np.dtype(array.dtype.descr)
		if sig:
  			array_flag = np.ones(len(array))
  			array_flag = np.ones(len(array))
		else:
  			array_flag = np.zeros(len(array))
  			array_flag = np.zeros(len(array))

  		discrete_mass = np.full(len(array),-1.)
		if (len(mcList)!=1):
  		   MCdiscrete = -1 
		   if sig:
  		      # for mass in MCMass:
		#		for mc in mcList:
  		 #       	 	if "_"+str(mass).replace(".","p") in mc:
		#				print 'matcher___________________', mc, mass
  		       MCdiscrete = idx_discrete[i]
  		   
		   else:
			MCdiscrete = MCMass[i]
		   print MCdiscrete
	  	   discrete_mass = np.full(len(array),MCdiscrete)
		   print discrete_mass[0:10]
  		else: 		
  			new_dt = np.dtype(array.dtype.descr + [('flag', 'int64')])
		print 'setted',MCdiscrete
	        array = rfn.append_fields(array,['massParam','flag'],[discrete_mass,array_flag])
		print discrete_mass[0:10] 
		
		print '_________________________________________',MCdiscrete,array[0:10]
		A = array if i == 0 else rfn.stack_arrays((A,array),usemask=False)
        	
	if sig:	
		return A, feat_hist_sig
	else: 
		return A, feat_hist_bkg
	
  def doScaling(self, X):
      '''
        Normalise the input features with a keras scaler 
      '''
      if self.scaler_type == 'robust':
        qt = RobustScaler()
      elif self.scaler_type == 'standard':
        qt = StandardScaler()
      else:
        raise RuntimeError('Unknown scaler "{}" - Aborting...'.format(self.scaler_type))

      qt.fit(X)
      xx = qt.transform(X)
 #     print xx.shape

      return xx, qt

  def preprocessing(self, mc_df, bkg_df):
    '''
      Preprocessing of data before training/testing the NN
      This includes:
        - building the main_df
        - building the scaler
        - get the scaled features xx
        - get the target Y
    '''

    #merge signal and background dataset + shuffle
    print "training - mc ",len(rfn.stack_arrays((mc_df), usemask = False ))
    print "training - data ",len(rfn.stack_arrays((bkg_df), usemask = False ))
    dataset = rfn.stack_arrays((mc_df,bkg_df), usemask = False )
     
    print dataset[0:20]
  #  print dataset[0:10]
    dataset = np.array(dataset.tolist(),dtype="float64")
    np.random.shuffle(dataset)
    
    self.n_features+=1

    
    #check obtained structure

   
    X = dataset[:,0:self.n_features]
    print X[0:20]
  #  print X.shape
    Y  = dataset[:,self.n_features]


    # scale the features
    # this is an important step!
    xx, qt = self.doScaling(X)
  #  print xx.shape

    # and save the scaler, which will have to be used throughout the full process, even at evaluation time
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted.pck'])
    pickle.dump(qt,open(scaler_filename, 'wb'))
    print ' --> {} created'.format(scaler_filename)

    # save the exact list of features
    features_filename = '/'.join([self.outdir, 'input_features.pck'])
    pickle.dump(self.features, open(features_filename, 'wb' ))
    print ' --> {} created'.format(features_filename)

    return dataset, qt, xx, Y

  def defineModel(self):
    '''
      Define the NN
    '''
    #NOTE for the moment, everything is hardcoded

    activation = 'relu'

    # define the net
    input  = Input((self.n_features,))
    layer1  = Dense(64, activation=activation, name='dense1', kernel_constraint=unit_norm())(input)
    output = Dense(1 , activation='sigmoid' , name='output', )(layer1)

    # Define outputs of your model
    model = Model(input, output)
    optimiser = Adam(lr=self.learning_rate)
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['mae', 'acc'])

    print(model.summary())

    return model

  def defineCallbacks(self):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    statistics = (self.signal_size+self.bkg_size)
	
    if statistics < 25000:
        patience_es = 5
        patience_lr = 3
    else:
        patience_es = 10
        patience_lr = 5
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=patience_es)
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=patience_lr, min_lr=0.00001, cooldown=10, verbose=True)
   

    # reduce learning rate when at plateau, fine search the minimum

    # save the model every now and then
    filepath = '/'.join([self.outdir, 'saved-model-{epoch:04d}_val_loss_{val_loss:.4f}_val_acc_{val_acc:.4f}.h5'])
    save_model = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=False, mode='auto', period=1)

    callbacks = [save_model]
    if self.do_early_stopping:
      callbacks.append(es)

    if self.do_reduce_lr:
      callbacks.append(reduce_lr)

    return callbacks

  def splitAndTransform(self,X,y,qt):

	# split train and test
	x_train, x_val, y_train, y_val = train_test_split(X, y, test_size=0.2, shuffle= True)
	# split train and val
  #      x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.2, shuffle= True)
	
	 # scale the features
      #  x_train = qt.transform(x_train)
       # x_val = qt.transform(x_val)

	self.training_size = len(x_train)
	self.valid_size = len(x_val)


	return x_train,x_val,y_train,y_val


  def train(self, model, x_train, y_train, x_val, y_val, callbacks):
    '''
      Perform the training
    '''
  #  print 'train sample size ', len(x_train)
    history = model.fit(x_train, y_train, validation_data=(x_val, y_val), epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True)  
    return history


def saveFig(outdir,plt, name):
    '''w
      Save python figure
    '''
    plt.savefig('{}/{}.pdf'.format(outdir, name))    
    plt.savefig('{}/{}.png'.format(outdir, name))    
    print ' --> {}/{}.png created'.format(outdir, name)

def plotLoss(epochs,history,outdir):
   '''
     Plot the loss for training and validation sets
   '''
   plt.clf()
   loss_train = history.history['loss']
   loss_val = history.history['val_loss']
   epochs_range = range(1, len(loss_train)+1)
   epochs = epochs_range
   plt.plot(epochs, loss_train, 'g', label='Training loss')
   plt.plot(epochs, loss_val, 'b', label='Validation loss')
   plt.title('Training and Validation Loss')
   plt.xlabel('Epochs')
   plt.ylabel('Loss')
   plt.legend()
   saveFig(outdir,plt, 'loss')


def plotAccuracy(epochs,history,outdir):
  '''
    Plot the accuracy for training and validation sets
  '''
  plt.clf()
  acc_train = history.history['acc']
  acc_val = history.history['val_acc']
  epochs_range = range(1, len(acc_train)+1)
  epochs = epochs_range
  plt.plot(epochs, acc_train, 'g', label='Training accuracy')
  plt.plot(epochs, acc_val, 'b', label='Validation accuracy')
  plt.title('Training and Validation Accuracy')
  plt.xlabel('Epochs')
  plt.ylabel('Accuracy')
  plt.legend()
  saveFig(outdir,plt, 'accuracy')

def plotROC(model,x_train_std,x_test_std,y_train,y_test,outdir):

  train_pred = model.predict(x_train_std)
  train_pred = [i[0] for i in train_pred]
  
  # let sklearn do the heavy lifting and compute the ROC curves for you
  #print("train",train_pred)
  fpr, tpr, wps = roc_curve(y_train, train_pred) 
  plt.plot(fpr, tpr, label='train ROC')
  print("AUC train",sk.metrics.auc(fpr,tpr))
  plt.xlabel('False Positive Rate')
  plt.ylabel('True Positive Rate')
  
  
  test_pred = model.predict(x_test_std)
  test_pred = [i[0] for i in test_pred]
  fpr, tpr, wps = roc_curve(y_test, test_pred)
  test = np.ma.array(test_pred,mask=y_test)
  #print(test[test.mask])
  plt.plot(fpr, tpr, label='test ROC')
  print("AUC test",sk.metrics.auc(fpr,tpr))
  
  xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
  plt.plot(xy, xy, color='grey', linestyle='--')
  plt.yscale('linear')
  plt.xscale('log')
  plt.legend()
  saveFig(outdir,plt,'roc_weighted')
  print("RoC cuerve saved")


def predictScore(outdir,features,model, sample):
  '''
    Return score with scaled input features
  '''
  x = pd.DataFrame(sample)
  # apply the scaler
  scaler_filename = '/'.join([outdir, 'input_tranformation_weighted.pck'])
 # scaler_filename = '/'.join(['model_tmvaVars', 'input_tranformation_weighted.pck'])
  qt = pickle.load(open(scaler_filename, 'rb'))
  
  xx = qt.transform(x)

  # predict
  score = model.predict(xx)

  return score

def plotScore( outdir, features,model, mc, bkg,lable,train_samples):
  '''
    Plot the score distributions for signal and background
  '''
  # get the score for the test dataframe
  sig_score =predictScore(outdir,features,model, mc)
  
  bkg_score =predictScore(outdir,features,model, bkg)

  # add the score to the dataframes
#  mc_test_df['score'] = sig_score
 # data_test_df['score'] = bkg_score

  # plot the score distributions
  fig = plt.figure(figsize=(20,10))
  ax = fig.add_subplot(111)
 # bkg_score = [data_test_df['score']]
  bkg_name=['data-driven background']
  ax.hist(bkg_score, bins=np.arange(0,1.025,0.025), stacked=True, alpha=0.5, label=bkg_name, normed=True)
  ax.hist(sig_score, bins=np.arange(0,1.025,0.025), alpha=1, label='signal', histtype='step', linewidth=2,normed=True)
  ax.legend(loc='upper left',prop={'size': 12})
  ax.set_title("Score distribution of signal and background for testing set", fontsize=20)
  ax.set_xlabel('Score',fontsize=18)
  ymin,ymax = ax.get_ylim()
  for i, sample in enumerate(train_samples):	
	ax.text(0.75,ymax*0.95, "Training samples:")
	ax.text(0.75,ymax*(0.9-i*0.05), sample)
  #fig.savefig('outputs/score.pdf')
  #fig.savefig('outputs/score.png')
  saveFig(outdir,fig, 'score'+lable)


# load the dataset
#inputSig = ["../data/HNLFlatTuples/NewSignature_1p0_ctau10p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_3p0_ctau10p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_4p5_ctau10p0/HNLFlat_0.root"]
#otherSig = ["../data/HNLFlatTuples/NewSignature_2p0_ctau100p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_3p0_ctau10p0/HNLFlat_0.root"]
#inputBkg =[
#	"../data/HNLFlatTuples/Parking1A0_newSig_4varsL_section0/HNLFlat_*11.root",
##	"../data/HNLFlatTuples/Pt-15to20_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-20to30_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-30to50_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-50to80_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-80to120_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-80to120_ext_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-120to170_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-120to170_ext_newSig/HNLFlat_*.root",
##	"../data/HNLFlatTuples/Pt-170to300_newSig/HNLFlat_*.root",
#
#	]
##features
##ele_branches = ["hnl_cos2D","abs(hnl_pi_DCAS)","hnl_pi_pt","B_mass","QCDweight"]
#MCMass = [1.0,1.5,2.0,3.0,4.5]
##hnl mass integrals for mass discretization
#massBin=[0.0,1.3,1.8,2.5,4,5.5,7]
#ele_branches = ["hnl_cos2D","abs(hnl_lxy_sig)","B_mass","max(hnl_l_dxyS,hnl_pi_dxyS)","QCDweight"]
#signals = []
#signal_size = 5000
#Sig = []
##nparray with features for signal
#for i,sample in enumerate(inputSig):
#	signal = root2array(sample,"Events",branches= ele_branches)
#	signal.resize((5000,))
#	MCdiscrete = -1 
#	for mass in MCMass:
#		if str(mass).replace(".","p") in sample:
#			print(mass)	
#			MCdiscrete = mass
#			break
#	discrete_Sigmass = np.full(len(signal),MCdiscrete)
#	print len(discrete_Sigmass)
#	signal_flag = np.ones(len(signal))
#	new_dt = np.dtype(signal.dtype.descr + [('massParam','float64')]+ [('flag', 'int64')])
#	Sig.append(np.zeros(signal.shape,dtype=new_dt))
#	print(Sig[i])
#	for b in ele_branches:
#		Sig[i][b] = signal[b]
#	Sig[i]['massParam'] = discrete_Sigmass
#	Sig[i]['flag'] = signal_flag
#
##sig_preprocessed = rfn.stack_arrays((Sig), usemask = True ) 
##print(sig_preprocessed[0:10])
#signals.append(root2array(otherSig[0],"Events",branches= ele_branches[0:4]+["hnl_mass"]))
#signals.append(root2array(otherSig[1],"Events",branches= ele_branches[0:4]+["hnl_mass"]))
##create output column - 1 for signal
#
##append it to features with value fixed to 1 for signal
#
##same for bkg:nparray with features
#bkg_array = []
#for bkg in inputBkg:
#	bkg = root2array(inputBkg,"Events",branches= ele_branches)
#	bkg_mass = root2array(inputBkg,"Events",branches= ["hnl_mass"])
#	bkg.resize((50000,))
#	bkg_mass.resize((50000,))
#	bkg_mass = np.array(bkg_mass.tolist(),dtype="float64")
##create output column - 0 for bkg
#	bkg_flag = np.zeros(len(bkg))
#	discrete_bkgMass = np.zeros(len(bkg))
#	for i in range(0,len(bkg)):
#		for j in range(0,len(massBin)-1):
#			if (massBin[j] <= bkg_mass[i]) and (bkg_mass[i] <= massBin[j+1]):
#				discrete_bkgMass[i] = bkg_mass[(massBin[j] <= bkg_mass) & (bkg_mass <= massBin[j+1])].mean()
#	new_dt = np.dtype(bkg.dtype.descr  + [("massParam","float64")]+ [('flag', 'int64')])
#	Bkg = np.zeros(bkg.shape,dtype=new_dt)
#	for b in ele_branches:
#		Bkg[b] = bkg[b]
#	Bkg['massParam'] = discrete_bkgMass
#	Bkg['flag'] = bkg_flag
#	print Bkg[0:10]
##	Bkg = np.array(Bkg.tolist(),dtype="float64")
#	bkg_array.append(Bkg)
#	
##merge signal and background dataset + shuffle
#dataset = rfn.stack_arrays((Sig,bkg_array), usemask = False ) 
#dataset = np.array(dataset.tolist(),dtype="float64")
#np.random.shuffle(dataset)
#signals[0] = np.array(signals[0].tolist(),dtype="float64")
#signals[1] = np.array(signals[1].tolist(),dtype="float64")
#print(dataset[0:100])
##define training and score features
#X = dataset[:,0:len(ele_branches)+1]
#print(len(X))
#y = dataset[:,len(ele_branches)+1]
#print y[0:10]
#
##build keras sequential model
#activation = 'relu'
#
## define the net
#inp  = Input((len(ele_branches),))
#layer  = Dense(64, activation=activation, name='dense1', kernel_constraint=unit_norm())(inp)
#output = Dense(1 , activation='sigmoid' , name='output', )(layer)
#
## Define outputs of your model
#model = Model(inp, output)
##model = Sequential()
##model.add(Dense(12, input_shape=(len(ele_branches)-1,),activation='relu'))
##model.add(Dense(8, activation='relu'))
##model.add(Dense(1, activation='sigmoid'))
#epoch= 50
## compile the keras model
#opt = keras.optimizers.Adam(lr=0.001)
#	model.compile(optimizer=opt, loss='binary_crossentropy', metrics=['mae', 'acc'])
##model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
#print(model.summary())
#
##split test,train valid
#	x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.4, shuffle= True)
#x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.3, shuffle= True)
#	weight_train = x_train[:,len(ele_branches)-1]
#	weight_test = x_test[:,len(ele_branches)-1]
#	weight_val = x_val[:,len(ele_branches)-1]
#	print x_train[:,0:len(ele_branches)-1].shape
#	print x_train[:,len(ele_branches)].shape
#	x_train = np.delete(x_train,len(ele_branches)-1,axis=1)
#	x_test = np.delete(x_test,len(ele_branches)-1,axis=1)
#x_val = np.delete(x_val,len(ele_branches)-1,axis=1)
#	print x_train[0:50]
#print len(y_train)
##standardize + normalize features
#	trans = RobustScaler()
#	x_train_std = trans.fit_transform(x_train)
#	x_test_std = trans.transform(x_test)
#x_val_std = trans.transform(x_val)
#	pickle.dump(trans, open( '/'.join(('model_tmvaVars', 'input_tranformation_weighted.pck')), 'wb' ) )
#	print len(x_train_std)
#print len(y_train)
#
## fit the keras model on the dataset
#history = model.fit(x_train_std, y_train, epochs=50, batch_size=50,validation_data=(x_val_std, y_val))
#
#
#####################################
####### ROC CURVE ###################
#####################################
#train_pred = model.predict(x_train_std)
#	train_pred = [i[0] for i in train_pred]
#
## let sklearn do the heavy lifting and compute the ROC curves for you
##print("train",train_pred)
#fpr, tpr, wps = roc_curve(y_train, train_pred) 
#	plt.plot(fpr, tpr, label='train ROC')
#	print("AUC train",sk.metrics.auc(fpr,tpr))
#	plt.xlabel('False Positive Rate')
#	plt.ylabel('True Positive Rate')
#
#
#test_pred = model.predict(x_test_std)
#	test_pred = [i[0] for i in test_pred]
#	fpr, tpr, wps = roc_curve(y_test, test_pred)
#test = np.ma.array(test_pred,mask=y_test)
##print(test[test.mask])
#	plt.plot(fpr, tpr, label='test ROC')
#	print("AUC test",sk.metrics.auc(fpr,tpr))
#
#	xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
#	plt.plot(xy, xy, color='grey', linestyle='--')
#	plt.yscale('linear')
#	plt.xscale('log')
#plt.legend()
#	plt.savefig('/'.join(['model_tmvaVars/roc_weighted.pdf']) )
#	print("RoC cuerve saved")
#
#plt.clf()
#
#	xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
#	plt.plot(xy, xy, color='grey', linestyle='--')
#
#	model.save_weights('net_model_weights.h5')
#
###################################
####### TEST overfitting ##########
###################################
#
##kolmogorov test
#
#print len(x_train_std)
#
##x_train_sigMask =
##x_test_sigMask = np.ma.array(x_test_std,mask= y_test)
#	train_pred_pass = model.predict(x_train_std.compress(y_train.astype(bool),axis=0))
#	len(train_pred_pass)
#test_pred_pass = model.predict(x_test_std.compress(y_test.astype(bool),axis=0))
#
##plot
#	h1 = ROOT.TH1F("","train",30,0,1)
#	h2 = ROOT.TH1F("","test",30,0,1)
#	for t,b in zip(train_pred_pass,test_pred_pass):
#	h1.Fill(t)
#	h2.Fill(b)
#	c1=ROOT.TCanvas()
#	h1.Scale(1/h1.Integral())
#	h2.Scale(1/h2.Integral())
#	c1.Draw()
#	h1.Draw("hist")
#	h2.SetLineColor(ROOT.kRed)
#	h2.SetFillColor(ROOT.kWhite)
#	h1.SetFillColor(ROOT.kWhite)
#	h2.Draw("hist SAME")
#	c1.BuildLegend()
#	ks_score = h1.KolmogorovTest(h2)
#	ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
#	ks_value.AddText('ks score pass = %.4f' %ks_score)
#	ks_value.SetFillColor(0)
#	ks_value.Draw('EP same')
#
#	c1.SaveAs("model_tmvaVars/KS_overfitting_test_pass_"+str(epoch)+".pdf")
#
##ks_fail_rw = h_ratio_passwpass.KolmogorovTest(h_ratio_passpass)
#	print("KS score: ",ks_score, len(train_pred_pass),len(test_pred_pass))
##print("KS fail rw: ",ks_fail_rw)
#
#	train_pred_fail = model.predict(x_train_std.compress(~y_train.astype(bool),axis=0))
#test_pred_fail = model.predict(x_test_std.compress(~y_test.astype(bool),axis=0))
#
##plot
#	h1 = ROOT.TH1F("","train",30,0,1)
#	h2 = ROOT.TH1F("","test",30,0,1)
#	for t,b in zip(train_pred_fail,test_pred_fail):
#	h1.Fill(t)
#	h2.Fill(b)
#	c1=ROOT.TCanvas()
#	h1.Scale(1/h1.Integral())
#	h2.Scale(1/h2.Integral())
#	c1.Draw()
#	h1.Draw("hist")
#	h2.SetLineColor(ROOT.kRed)
#	h2.SetFillColor(ROOT.kWhite)
#	h1.SetFillColor(ROOT.kWhite)
#	h2.Draw("hist SAME")
#	c1.BuildLegend()
#	ks_score = h1.KolmogorovTest(h2)
#	ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
#	ks_value.AddText('ks score fail = %.4f' %ks_score)
#	ks_value.SetFillColor(0)
#	ks_value.Draw('EP same')
#
#	c1.SaveAs("model_tmvaVars/KS_overfitting_test_fail_"+str(epoch)+".pdf")
#
#
###########################################################################################
######   CORRELATION MATRIX SIGNAL
###########################################################################################
## Compute the correlation matrix for the signal
#
##passing_array = qt.transform(passing[features])
##failing_array = qt.transform(failing[features])
##print(passing_array)
##print(len(passing_array))
##print(passing.shape)
##print(model.predict(passing_array))
##print(passing)
##print([i[0] for i in model.predict(passing_array)])
##pas.loc[:,'nn'] = [i[0] for i in model.predict(passing_array)]
##failing.loc[:,'nn'] = [i[0] for i in model.predict(failing_array)]
#	print test_pred_pass[0:10]
#	corr_test = np.append(x_test_std.compress(y_test.astype(bool),axis=0),test_pred_pass,axis=1)
#corr = np.corrcoef(corr_test,rowvar=False)
#	print corr.shape
#
## Set up the matplotlib figure
#f, ax = plt.subplots(figsize=(11, 9))
#
## Generate a custom diverging colormap
#cmap = sns.diverging_palette(220, 10, as_cmap=True)
#
## Draw the heatmap with the mask and correct aspect ratio
#	g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
#			square=True, linewidths=.8, cbar_kws={"shrink": .8})
#
## rotate axis labels
#	g.set_xticklabels(ele_branches[0:len(ele_branches)-1]+['score'], rotation='vertical')
#	g.set_yticklabels(ele_branches[0:len(ele_branches)-1]+['score'], rotation='horizontal')
#
## plt.show()
#	plt.title('linear correlation matrix - pass')
#plt.tight_layout()
#	plt.savefig('model_tmvaVars/corr_pass.png' )
#plt.clf()
#
###########################################################################################
######   CORRELATION MATRIX BACKGROUND
###########################################################################################
## Compute the correlation matrix for the signal
#	print test_pred_fail[0:10]
#	corr_test = np.append(x_test_std.compress(~y_test.astype(bool),axis=0),test_pred_fail,axis=1)
#corr = np.corrcoef(corr_test,rowvar=False)
#	print corr.shape
#
## Set up the matplotlib figure
#f, ax = plt.subplots(figsize=(11, 9))
#
## Generate a custom diverging colormap
#cmap = sns.diverging_palette(220, 10, as_cmap=True)
#
## Draw the heatmap with the mask and correct aspect ratio
#	g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
#			square=True, linewidths=.8, cbar_kws={"shrink": .8})
#
## rotate axis labels
#	g.set_xticklabels(ele_branches[0:len(ele_branches)-1]+['score'], rotation='vertical')
#	g.set_yticklabels(ele_branches[0:len(ele_branches)-1]+['score'], rotation='horizontal')
#
## plt.show()
#	plt.title('linear correlation matrix - fail')
#plt.tight_layout()
#	plt.savefig('model_tmvaVars/corr_fail.png' )
#
#
#
#	plotLoss(epoch,history)
#plotAccuracy(epoch,history)
#	plotScore(ele_branches[0:len(ele_branches)-1]+["hnl_mass"],model,x_test.compress(y_test.astype(bool),axis=0),x_test.compress(~y_test.astype(bool),axis=0) ,"")
#
#	plotScore(ele_branches[0:len(ele_branches)-1]+["hnl_mass"],model,signals[0],x_test.compress(~y_test.astype(bool),axis=0),"mass2p0_ctau10" )
#	plotScore(ele_branches[0:len(ele_branches)-1]+["hnl_mass"],model,signals[1],x_test.compress(~y_test.astype(bool),axis=0),"mass3p0_ctau10" )
#
def feat_plotter(outdir,hist_sig,hist_bkg):

	color = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kOrange+7,ROOT.kMagenta]
  	MCMass = [1.0,1.5,2.0,3.0,4.5]

	print 'here_______________________________',len(hist_sig),len(hist_sig[0])
	print 'here_______________________________',len(hist_bkg),len(hist_bkg[0])
	for i in range(0,len(hist_sig)):
		for j in range(0,len(hist_sig[0])):
				
			c = ROOT.TCanvas()
			if hist_sig[i][j].GetEntries() < 0 or hist_bkg[i][j].GetEntries() < 0:
				continue
			print type(hist_sig[i][j])	
			hist_sig[i][j].SetLineWidth(2)
			hist_sig[i][j].SetLineColor(color[i])
			hist_bkg[i][j].SetLineColor(ROOT.kBlue-7)
			hist_bkg[i][j].SetFillColor(ROOT.kBlue-7)
			hist_bkg[i][j].SetFillStyle(3004)
		#	print 'entries',hist_bkg[i][j].GetEntries()	
			hist_bkg[i][j].GetXaxis().SetTitle(features[j])
			hist_bkg[i][j].GetYaxis().SetTitle('arbitry units')

		#	print 'here',hist_bkg[i][j].GetName()	
		
			hist_bkg[i][j].SetTitle('Background 10 #sigma window')
			hist_sig[i][j].SetTitle('signal mass '+str(MCMass[i])+ ' GeV ')
	
			hist_bkg[i][j].Scale(1./hist_bkg[i][j].Integral())
			hist_sig[i][j].Scale(1./hist_sig[i][j].Integral())
			hist_bkg[i][j].Draw('hist')
			hist_sig[i][j].Draw('same hist')
		
			hist_bkg[i][j].GetYaxis().SetRangeUser(0,1.3*max(hist_bkg[i][j].GetMaximum(),hist_sig[i][j].GetMaximum()))
			c.BuildLegend()
			c.RedrawAxis()
#			hist_bkg[i][j].SaveAs('/'.join([outdir,(features[j]+'_mass'+str(MCMass[i])+'bkg.root')]))
		#	print 'before saving pdf', '/'.join([outdir,(features[j]+'_mass'+str(MCMass[i])+'.pdf')])
			c.SaveAs('/'.join([outdir,(features[j]+'_mass'+str(MCMass[i])+'.pdf')]))
			
		



def singleTraining(specs):



  #NOTE add optimiser, learning rate etc? 

  Bc_string = ""
  if Bc:
	Bc_string="_Bc" 
  trainer = Trainer(
      features = features, 
      n_features = n_features,
      epochs = epochs,
      batch_size = batch_size,
      signal_size = signal_size,
      bkg_size = bkg_size,
      learning_rate = learning_rate,
      scaler_type = scaler_type,
      do_early_stopping = do_early_stopping,
      do_reduce_lr = do_reduce_lr,
      dirname = "trainNN_"+specs[0]+Bc_string,
      baseline_selection = specs[1],
      train_samples = samples
)
  trainer.outdir = trainer.createOutDir()     

  sig_len = np.zeros(5)
  if Bc:
 	 Sig,hist_sig = trainer.inputHandler(inputSig_Bc,inputBkg,True)	
 	 Bkg,hist_bkg = trainer.inputHandler(inputSig_Bc,inputBkg,False)
  else:  

 	 Sig,hist_sig = trainer.inputHandler(inputSig,inputBkg,True)
	 print 'check hits___________________',type(hist_sig[0][0])	
 	 Bkg,hist_bkg = trainer.inputHandler(inputSig,inputBkg,False)
	 print 'check hits___________________',type(hist_bkg[0][0])	
  
  feat_plotter(trainer.outdir,hist_sig,hist_bkg)
# dataset, qt, xx, Y = trainer.preprocessing(Sig,Bkg)

# model = trainer.defineModel() 	
# callbacks = trainer.defineCallbacks() 	
# 
# x_train, x_val,y_train, y_val = trainer.splitAndTransform(xx,Y,qt)
#       
# #print training stats on file
# with open(trainer.outdir+"/stats.txt","w") as file:

#       file.write("aimed stats per mass point %i \n"%trainer.min_mass_size)
#       file.write("full stat: MC %i \n"%trainer.signal_size)
#       file.write("full stat: Data %i \n"%trainer.bkg_size)
#       file.write("full stat: Data+MC %i \n"%(trainer.signal_size+trainer.bkg_size))
#       file.write("training stats: %i \n"%trainer.training_size)
#       file.write("validation stats: %i \n"%trainer.valid_size)


# #print(len(y_val.compress(y_val.astype(bool),axis=0)))
# history = trainer.train(model, x_train, y_train, x_val, y_val, callbacks)
# model.save_weights('net_model_weights.h5')
# print callbacks 
# trainer.epochs = min(callbacks[1].stopped_epoch,trainer.epochs)
# print trainer.epochs
# plotROC(model,x_train,x_val,y_train,y_val,trainer.outdir)
# plotLoss(trainer.epochs,history,trainer.outdir)
# plotAccuracy(trainer.epochs,history,trainer.outdir)
# plotScore(trainer.outdir,trainer.features,model,qt.inverse_transform(x_val).compress(y_val.astype(bool),axis=0),qt.inverse_transform(x_val).compress(~y_val.astype(bool),axis=0),"",trainer.train_samples)
##plotScore(trainer.outdir,trainer.features,model,x_val.compress(y_val.astype(bool),axis=0),x_val.compress(~y_val.astype(bool),axis=0) ,"",)
#
###################################
####### TEST overfitting ##########
###################################
#
##kolmogorov test
#
#
# train_pred_pass = model.predict(x_train.compress(y_train.astype(bool),axis=0))
# len(train_pred_pass)
# test_pred_pass = model.predict(x_val.compress(y_val.astype(bool),axis=0))
# 
# #plot
# h1 = ROOT.TH1F("","train",30,0,1)
# h2 = ROOT.TH1F("","test",30,0,1)
# for t,b in zip(train_pred_pass,test_pred_pass):
#     h1.Fill(t)
#     h2.Fill(b)
# c1=ROOT.TCanvas()
# h1.Scale(1/h1.Integral())
# h2.Scale(1/h2.Integral())
# c1.Draw()
# h1.Draw("hist")
# h2.SetLineColor(ROOT.kRed)
# h2.SetFillColor(ROOT.kWhite)
# h1.SetFillColor(ROOT.kWhite)
# h2.Draw("hist SAME")
# c1.BuildLegend()
# ks_score = h1.KolmogorovTest(h2)
# ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
# ks_value.AddText('ks score pass = %.4f' %ks_score)
# ks_value.SetFillColor(0)
# ks_value.Draw('EP same')
# 
# c1.SaveAs(trainer.outdir+"/KS_overfitting_test_pass_"+str(trainer.epochs)+".pdf")
# 
# #ks_fail_rw = h_ratio_passwpass.KolmogorovTest(h_ratio_passpass)
# print("KS score: ",ks_score, len(train_pred_pass),len(test_pred_pass))
# #print("KS fail rw: ",ks_fail_rw)
# 
# train_pred_fail = model.predict(x_train.compress(~y_train.astype(bool),axis=0))
# test_pred_fail = model.predict(x_val.compress(~y_val.astype(bool),axis=0))
# 
# #plot
# h1 = ROOT.TH1F("","train",30,0,1)
# h2 = ROOT.TH1F("","test",30,0,1)
# for t,b in zip(train_pred_fail,test_pred_fail):
#     h1.Fill(t)
#     h2.Fill(b)
# c1=ROOT.TCanvas()
# h1.Scale(1/h1.Integral())
# h2.Scale(1/h2.Integral())
# c1.Draw()
# h1.Draw("hist")
# h2.SetLineColor(ROOT.kRed)
# h2.SetFillColor(ROOT.kWhite)
# h1.SetFillColor(ROOT.kWhite)
# h2.Draw("hist SAME")
# c1.BuildLegend()
# ks_score = h1.KolmogorovTest(h2)
# ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
# ks_value.AddText('ks score fail = %.4f' %ks_score)
# ks_value.SetFillColor(0)
# ks_value.Draw('EP same')
# 
# c1.SaveAs(trainer.outdir+"/KS_overfitting_test_fail_"+str(trainer.epochs)+".pdf")
# 
# 
# ##########################################################################################
# #####   CORRELATION MATRIX SIGNAL
# ##########################################################################################
# # Compute the correlation matrix for the signal
# 
# #passing_array = qt.transform(passing[features])
# #failing_array = qt.transform(failing[features])
# #print(passing_array)
# #print(len(passing_array))
# #print(passing.shape)
# #print(model.predict(passing_array))
# #print(passing)
# #print([i[0] for i in model.predict(passing_array)])
# #pas.loc[:,'nn'] = [i[0] for i in model.predict(passing_array)]
# #failing.loc[:,'nn'] = [i[0] for i in model.predict(failing_array)]
# #print test_pred_pass[0:10]
# corr_test = np.append(x_val.compress(y_val.astype(bool),axis=0),test_pred_pass,axis=1)
# corr = np.corrcoef(corr_test,rowvar=False)
# print corr.shape
# 
# # Set up the matplotlib figure
# f, ax = plt.subplots(figsize=(11, 9))
# 
# # Generate a custom diverging colormap
# cmap = sns.diverging_palette(220, 10, as_cmap=True)
# 
# # Draw the heatmap with the mask and correct aspect ratio
# g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
#                 square=True, linewidths=.8, cbar_kws={"shrink": .8})
# 
# # rotate axis labels
# g.set_xticklabels(features+['score'], rotation='vertical')
# g.set_yticklabels(features+['score'], rotation='horizontal')
# 
# # plt.show()
# plt.title('linear correlation matrix - pass')
# plt.tight_layout()
# plt.savefig(trainer.outdir+'/corr_pass.png' )
# plt.clf()
   
if __name__ == '__main__':
	Bc = int(sys.argv[1])
  ##########################################################################################
#	features = ["hnl_cos2D","hnl_vtxProb","hnl_pi_pt","B_mass","hnl_vtxChi2","hnl_l_pt","TrgMu_pt","hnl_lxy_sig","hnl_pi_nHits"]
	if Bc:
		features = ["hnl_pi_pt","l0_pt","l_pt","B_mass","hnl_cos2D","hnl_lxy_sig","hnl_pi_DCAS","hnl_vtxProb","B_pt","dilepton_mass","BlepPi_mass","dr_trgMu_lep","dr_Blep_pi","hnl_pi_PixelLayers","hnl_pi_TrackerLayers","TrgMu_PixelLayers","TrgMu_TrackerLayers","hnl_l_lostHits"]
		#features = ["hnl_pi_pt","l0_pt","l_pt","B_mass","hnl_cos2D","hnl_lxy_sig","hnl_pi_DCAS","hnl_vtxProb","B_pt","dilepton_mass","BlepPi_mass","dr_trgMu_lep","dr_Blep_pi","l0_relIso","l_relIso","hnl_pi_PixelLayers","hnl_pi_TrackerLayers","TrgMu_PixelLayers","TrgMu_TrackerLayers"]
	else:
		features = ["hnl_pi_pt","l0_pt","l_pt","B_mass","hnl_cos2D","hnl_lxy_sig","hnl_pi_DCAS","hnl_vtxProb","B_pt","dilepton_mass","BlepPi_mass","dr_trgMu_lep","dr_Blep_pi","hnl_pi_PixelLayers","hnl_pi_TrackerLayers","TrgMu_PixelLayers","TrgMu_TrackerLayers","hnl_l_lostHits"]
	feat_binning = [[50,0,50],[50,0,50],[50,0,50],[50,0,15],[50,0.99,1],[100,0,700],[50,0,100],[50,0,1],[50,0,50],[50,0,7],[50,0,15],[50,0,15],[50,0,15],[6,0,6],[20,0,20],[6,0,6],[20,0,20],[4,0,4]] 
	inputSig = [
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
 	 	"../data/HNLFlatTuples/NewSig_hits_1p0_ctau1000p0/HNLFlat_0.root",
 	 	"../data/HNLFlatTuples/NewSig_hits_1p5_ctau1000p0/HNLFlat_0.root",
 	 	"../data/HNLFlatTuples/NewSig_hits_2p0_ctau1000p0/HNLFlat_0.root",
		"../data/HNLFlatTuples/NewSig_hits_3p0_ctau1000p0/HNLFlat_0.root"
		]#,"../data/HNLFlatTuples/NewSignature_1p5_ctau100p0/HNLFlat_0.root"]#,"../data/HNLFlatTuples/NewSignature_3p0_ctau1p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_2p0_ctau10p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_1p0_ctau10p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_4p5_ctau0p1/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_1p5_ctau10p0/HNLFlat_0.root"
	otherSig = ["../data/HNLFlatTuples/NewSignature_3p0_ctau100p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_3p0_ctau1p0/HNLFlat_0.root"]

	inputSig_Bc = [
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
 	inputBkg =[
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
#	"../data/HNLFlatTuples/Parking_1D_3_0_Bmass/HNLFlat_*1.root",
#	"../data/HNLFlatTuples/Pt-15to20_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-20to30_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-30to50_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-50to80_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-80to120_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-80to120_ext_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-120to170_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-120to170_ext_newSig/HNLFlat_*.root",
#	"../data/HNLFlatTuples/Pt-170to300_newSig/HNLFlat_*.root",

	]

	n_features = len(features)
	epochs = 100
	batch_size = 32
	signal_size = 10000
	bkg_size = 100000
	scaler_type = 'robust'
	learning_rate = 0.01
	do_early_stopping = True
	do_reduce_lr = True
	dirname = 'test'
	if Bc:
		baseline_selection = "hnl_charge==0 && fabs(hnl_l_eta)<1.442 && B_mass>5.7 && B_mass<7 && hnl_l_mvaId>-3 && dr_BlepPi>0.01"
	else:
		baseline_selection = "hnl_charge==0 && fabs(hnl_l_eta)<1.442 && B_mass<5.7 && hnl_l_mvaId>-3 && dr_BlepPi>0.01 "
	samples = []
	tag = "test"
   
     #   ctx = multiprocessing.set_start_method('spawn')
	specs = [
 		 	["lxysig0to50_OS","hnl_charge==0 && hnl_lxy_sig<50 && LepQProd<0 && fabs(hnl_l_eta)<1.442"],
         	 	["lxysig0to50_SS","hnl_charge==0 && hnl_lxy_sig<50 && LepQProd>0 && fabs(hnl_l_eta)<1.442"],
         		["lxysig50to150_OS","hnl_charge==0 && hnl_lxy_sig>50 && hnl_lxy_sig<150 && LepQProd<0 && fabs(hnl_l_eta)<1.442"],
         		["lxysig50to150_SS","hnl_charge==0 && hnl_lxy_sig>50 && hnl_lxy_sig<150 && LepQProd>0 && fabs(hnl_l_eta)<1.442"],
        		["lxysiggt150_OS","hnl_charge==0 && hnl_lxy_sig>150 && LepQProd<0 && fabs(hnl_l_eta)<1.442"],
         		["lxysiggt150_SS","hnl_charge==0 && hnl_lxy_sig>150 && LepQProd>0 && fabs(hnl_l_eta)<1.442"]
			]
 
#	singleTraining(specs[0])
	#output = Parallel(n_jobs=6)(delayed(singleTraining)(i) for i in specs)	 
   	Tpool = Pool(len(specs))
 #    	Tpool = Pool(1)
    	Tpool.map(singleTraining,specs) 
   	Tpool.close()
   	Tpool.join()
#	results = Parallel(n_jobs=6, prefer='threads')(delayed(singleTraining)(i) for i in specs)
  	#results = Tpool.map_async(singleTraining,specs) 
#	singleTraining(specs[0])
  #	Tpool.termminate()
 	#Tpool.close()
 	#Tpool.join()











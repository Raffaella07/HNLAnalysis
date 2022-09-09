# first neural network with keras tutorial

import os 
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


class Trainer(object):

  def __init__(self, features,n_features, epochs, batch_size,signal_size,bkg_size, scaler_type, do_early_stopping, do_reduce_lr, dirname, baseline_selection,train_samples):
    self.features = features
    self.n_features = n_features
    self.epochs = epochs
    self.batch_size = batch_size
    self.signal_size = signal_size
    self.bkg_size = bkg_size
    self.scaler_type = scaler_type
    self.do_early_stopping = do_early_stopping
    self.do_reduce_lr = do_reduce_lr
    self.dirname = dirname + '_' + datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss')

    self.baseline_selection = baseline_selection
    self.train_samples = train_samples

    self.target_branch = 'flag'

  def createOutDir(self):
     '''
       This function creates the output directory
     '''
     outdir = './model/{}'.format(self.dirname)
     if not path.exists(outdir):
       os.system('mkdir -p {}'.format(outdir))
     return outdir

  def inputHandler(self,mcList,bkgList,sig ):
     
	'''
	 Function setting signal and background flags + preparing paramMass variable when going for pNN
	''' 

  	#hnl mass for mass discretization - signal
  	MCMass = [1.0,1.5,2.0,3.0,4.5]
  	#hnl mass for mass discretization - background
  	massBin=[0.0,1.3,1.8,2.5,3.5,5.5,7]
	
	temp = [mcList[i].split("/") for  i in range (0, len(mcList))]	
	self.train_samples= [temp[i][3].replace("NewSignature_","") for  i in range (0, len(temp))]
	print self.train_samples
	 

 	A = [] 
	if sig:
		List = mcList 
		array_size = signal_size
	else:
		List = bkgList 
		array_size = bkg_size

  	if len(mcList) ==1:
  		print "\n"
  		print "********* PREPARING INPUTS ************"
  		print "     training on single sample"
  		print "***************************************"
  		print "\n"
  	else:
  		print "\n"
  		print "*************************** PREPARING INPUTS **************************************"
  		print "                training on multiple samples: pNN training"
  		print "             will build a new feature 'massParam' with values:"
  		print "       MC truth for signal: sharp value in ", MCMass
  		print "  Mean over an interval for bkg: discretization bins: ", massBin
  		print "***********************************************************************************"
  		print "\n"
  
  
  	for i,sample in enumerate(List):
  		array = root2array(sample,"Events",branches= self.features, selection="hnl_charge==0")
		mass = root2array(sample,"Events",branches= ["hnl_mass"])
		mass = np.array(mass.tolist(),dtype="float64")
  		array.resize((int(array_size),))
  		mass.resize((int(array_size),))
  		new_dt = np.dtype(array.dtype.descr)
		if sig:
  			array_flag = np.ones(len(array))
  			array_flag = np.ones(len(array))
		else:
  			array_flag = np.zeros(len(array))
  			array_flag = np.zeros(len(array))

  		discrete_mass = np.full(len(array),-1.)
  		if len(mcList) !=1:
  			new_dt = np.dtype(array.dtype.descr + [('massParam','float64')]+ [('flag', 'int64')])
			if sig:
  				MCdiscrete = -1 
  				for mass in MCMass:
  					if str(mass).replace(".","p") in sample:
  						MCdiscrete = (mass--massBin[0])/(massBin[len(massBin)-1]-massBin[0])
  						break
  				discrete_mass = np.full(len(array),MCdiscrete)
			else:
			 	for k in range(0,len(mass)):
 					for j in range(0,len(massBin)-1):
 						if (massBin[j] <= mass[k]) and (mass[k] <= massBin[j+1]):
 							discrete_mass[k] = (mass[(massBin[j] <= mass) & (mass <= massBin[j+1])].mean()-massBin[0])/(massBin[len(massBin)-1]-massBin[0])
							
  		else: 		
  			new_dt = np.dtype(array.dtype.descr + [('flag', 'int64')])
  		A.append(np.zeros(array.shape,dtype=new_dt))
	 	
		if len(A)>1 or sig: 
			for b in self.features:
				A[i][b] = array[b]
			if len(mcList) !=1:
				A[i]['massParam'] = discrete_mass
				if i ==0:
					self.n_features+=1
			A[i]['flag'] = array_flag
	  		print(A[i][0:3])
	  		print("\n")
		else:
			for b in self.features:
				A[0][b] = array[b]
			if len(mcList) !=1:
				A[0]['massParam'] = discrete_mass
			A[0]['flag'] = array_flag
	  		print(A[0][0:10])
	  		print("\n")
	  		
	return A

	
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

      qt.fit(X[:,0:self.n_features-1])
      xx = qt.transform(X[:,0:self.n_features-1])
      print xx.shape
      print X[:,(self.n_features-1)].shape
      xx = np.insert(xx,(self.n_features-1),X[:,(self.n_features-1)],axis=1)

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
    dataset = rfn.stack_arrays((mc_df,bkg_df), usemask = False ) 
    dataset = np.array(dataset.tolist(),dtype="float64")
    np.random.shuffle(dataset)

    #check obtaine structure
    print(dataset[0:100])

    #define training and score features
    print '____________________________________________________________________________________________________________', self.n_features
    X = dataset[:,0:self.n_features]
    print X.shape
    Y  = dataset[:,self.n_features]


    # scale the features
    # this is an important step!
    xx, qt = self.doScaling(X)
    print xx.shape

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
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['mae', 'acc'])

    print(model.summary())

    return model

  def defineCallbacks(self):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=50)

    # reduce learning rate when at plateau, fine search the minimum
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=5, min_lr=0.00001, cooldown=10, verbose=True)

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
	x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.4, shuffle= True)
	# split train and val
        x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.3, shuffle= True)
	
	 # scale the features
      #  x_train = qt.transform(x_train)
       # x_val = qt.transform(x_val)


	return x_train,x_test, x_val,y_train,y_test, y_val


  def train(self, model, x_train, y_train, x_val, y_val, callbacks):
    '''
      Perform the training
    '''
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
   epochs_range = range(1, epochs+1)
   epochs = epochs_range
   plt.plot(epochs, loss_train, 'g', label='Training loss')
   plt.plot(epochs, loss_val, 'b', label='Validation loss')
   plt.title('Training and Validation Loss')
   plt.xlabel('Epochs')
   plt.ylabel('Loss')
   plt.legend()
   saveFig(outdir,plt, 'loss')


def plotAccuracy(epochs,historyy,outdir):
  '''
    Plot the accuracy for training and validation sets
  '''
  plt.clf()
  acc_train = history.history['acc']
  acc_val = history.history['val_acc']
  epochs_range = range(1, epochs+1)
  epochs = epochs_range
  plt.plot(epochs, acc_train, 'g', label='Training accuracy')
  plt.plot(epochs, acc_val, 'b', label='Validation accuracy')
  plt.title('Training and Validation Accuracy')
  plt.xlabel('Epochs')
  plt.ylabel('Accuracy')
  plt.legend()
  saveFig(outdir,plt, 'accuracy')

def plotROC(model,x_train_std,x_test_std,outdir):

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




if __name__ == '__main__':

  features = ["hnl_cos2D","B_mass","hnl_pi_pt","hnl_pi_DCAS"]
 
  inputSig = ["../data/HNLFlatTuples/NewSignature_3p0_ctau100p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_2p0_ctau100p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_1p0_ctau100p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_4p5_ctau100p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_1p5_ctau100p0/HNLFlat_0.root"]#,"../data/HNLFlatTuples/NewSignature_3p0_ctau1p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_2p0_ctau10p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_1p0_ctau10p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_4p5_ctau0p1/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_1p5_ctau10p0/HNLFlat_0.root"]
  otherSig = ["../data/HNLFlatTuples/NewSignature_3p0_ctau100p0/HNLFlat_0.root","../data/HNLFlatTuples/NewSignature_3p0_ctau1p0/HNLFlat_0.root"]

  inputBkg =[
	"../data/HNLFlatTuples/Parking1A0_newSig_4varsL_section0/HNLFlat_*11.root",
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
  epochs = 50
  batch_size = 32
  signal_size = 3000
  bkg_size = 50000
  scaler_type = 'robust'
  do_early_stopping = False
  do_reduce_lr = False
  dirname = 'test'
  baseline_selection = 'hnl_charge==0'
  samples = []


  #NOTE add optimiser, learning rate etc? 


  trainer = Trainer(
      features = features, 
      n_features = n_features,
      epochs = epochs,
      batch_size = batch_size,
      signal_size = signal_size,
      bkg_size = bkg_size,
      scaler_type = scaler_type,
      do_early_stopping = do_early_stopping,
      do_reduce_lr = do_reduce_lr,
      dirname = "trainNN",
      baseline_selection = baseline_selection,
      train_samples = samples
)
  trainer.outdir = trainer.createOutDir()     
  Sig = trainer.inputHandler(inputSig,inputBkg,True)	
  Bkg = trainer.inputHandler(inputSig,inputBkg,False)
    

 
  dataset, qt, xx, Y = trainer.preprocessing(Sig,Bkg)

  model = trainer.defineModel() 	
  callbacks = trainer.defineCallbacks() 	
  
  x_train, x_test, x_val,y_train, y_test, y_val = trainer.splitAndTransform(xx,Y,qt)

  print(len(y_test.compress(y_test.astype(bool),axis=0)))
  history = trainer.train(model, x_train, y_train, x_val, y_val, callbacks)
   
  model.save_weights('net_model_weights.h5')


  plotROC(model,x_train,x_test,trainer.outdir)
  plotLoss(trainer.epochs,history,trainer.outdir)
  plotAccuracy(trainer.epochs,history,trainer.outdir)
  plotScore(trainer.outdir,trainer.features,model,qt.inverse_transform(x_test).compress(y_test.astype(bool),axis=0),qt.inverse_transform(x_test).compress(~y_test.astype(bool),axis=0),"",trainer.train_samples)
#  plotScore(trainer.outdir,trainer.features,model,x_test.compress(y_test.astype(bool),axis=0),x_test.compress(~y_test.astype(bool),axis=0) ,"")

  ##################################
  ###### TEST overfitting ##########
  ##################################
  
  #kolmogorov test
  
  
  train_pred_pass = model.predict(x_train.compress(y_train.astype(bool),axis=0))
  len(train_pred_pass)
  test_pred_pass = model.predict(x_test.compress(y_test.astype(bool),axis=0))
  
  #plot
  h1 = ROOT.TH1F("","train",30,0,1)
  h2 = ROOT.TH1F("","test",30,0,1)
  for t,b in zip(train_pred_pass,test_pred_pass):
      h1.Fill(t)
      h2.Fill(b)
  c1=ROOT.TCanvas()
  h1.Scale(1/h1.Integral())
  h2.Scale(1/h2.Integral())
  c1.Draw()
  h1.Draw("hist")
  h2.SetLineColor(ROOT.kRed)
  h2.SetFillColor(ROOT.kWhite)
  h1.SetFillColor(ROOT.kWhite)
  h2.Draw("hist SAME")
  c1.BuildLegend()
  ks_score = h1.KolmogorovTest(h2)
  ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
  ks_value.AddText('ks score pass = %.4f' %ks_score)
  ks_value.SetFillColor(0)
  ks_value.Draw('EP same')
  
  c1.SaveAs(trainer.outdir+"/KS_overfitting_test_pass_"+str(trainer.epochs)+".pdf")
  
  #ks_fail_rw = h_ratio_passwpass.KolmogorovTest(h_ratio_passpass)
  print("KS score: ",ks_score, len(train_pred_pass),len(test_pred_pass))
  #print("KS fail rw: ",ks_fail_rw)
  
  train_pred_fail = model.predict(x_train.compress(~y_train.astype(bool),axis=0))
  test_pred_fail = model.predict(x_test.compress(~y_test.astype(bool),axis=0))
  
  #plot
  h1 = ROOT.TH1F("","train",30,0,1)
  h2 = ROOT.TH1F("","test",30,0,1)
  for t,b in zip(train_pred_fail,test_pred_fail):
      h1.Fill(t)
      h2.Fill(b)
  c1=ROOT.TCanvas()
  h1.Scale(1/h1.Integral())
  h2.Scale(1/h2.Integral())
  c1.Draw()
  h1.Draw("hist")
  h2.SetLineColor(ROOT.kRed)
  h2.SetFillColor(ROOT.kWhite)
  h1.SetFillColor(ROOT.kWhite)
  h2.Draw("hist SAME")
  c1.BuildLegend()
  ks_score = h1.KolmogorovTest(h2)
  ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
  ks_value.AddText('ks score fail = %.4f' %ks_score)
  ks_value.SetFillColor(0)
  ks_value.Draw('EP same')
  
  c1.SaveAs(trainer.outdir+"/KS_overfitting_test_fail_"+str(trainer.epochs)+".pdf")
  
  
  ##########################################################################################
  #####   CORRELATION MATRIX SIGNAL
  ##########################################################################################
  # Compute the correlation matrix for the signal
  
  #passing_array = qt.transform(passing[features])
  #failing_array = qt.transform(failing[features])
  #print(passing_array)
  #print(len(passing_array))
  #print(passing.shape)
  #print(model.predict(passing_array))
  #print(passing)
  #print([i[0] for i in model.predict(passing_array)])
  #pas.loc[:,'nn'] = [i[0] for i in model.predict(passing_array)]
  #failing.loc[:,'nn'] = [i[0] for i in model.predict(failing_array)]
  print test_pred_pass[0:10]
  corr_test = np.append(x_test.compress(y_test.astype(bool),axis=0),test_pred_pass,axis=1)
  corr = np.corrcoef(corr_test,rowvar=False)
  print corr.shape
  
  # Set up the matplotlib figure
  f, ax = plt.subplots(figsize=(11, 9))
  
  # Generate a custom diverging colormap
  cmap = sns.diverging_palette(220, 10, as_cmap=True)
  
  # Draw the heatmap with the mask and correct aspect ratio
  g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
                  square=True, linewidths=.8, cbar_kws={"shrink": .8})
  
  # rotate axis labels
  g.set_xticklabels(features+['score'], rotation='vertical')
  g.set_yticklabels(features+['score'], rotation='horizontal')
  
  # plt.show()
  plt.title('linear correlation matrix - pass')
  plt.tight_layout()
  plt.savefig(trainer.outdir+'/corr_pass.png' )
  plt.clf()
  
  ##########################################################################################
  #####   CORRELATION MATRIX BACKGROUND
  ##########################################################################################
  # Compute the correlation matrix for the signal
  print test_pred_fail[0:10]
  corr_test = np.append(x_test.compress(~y_test.astype(bool),axis=0),test_pred_fail,axis=1)
  corr = np.corrcoef(corr_test,rowvar=False)
  print corr.shape
  
  # Set up the matplotlib figure
  f, ax = plt.subplots(figsize=(11, 9))
  
  # Generate a custom diverging colormap
  cmap = sns.diverging_palette(220, 10, as_cmap=True)
  
  # Draw the heatmap with the mask and correct aspect ratio
  g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
                  square=True, linewidths=.8, cbar_kws={"shrink": .8})
  
  # rotate axis labels
  g.set_xticklabels(features+['score'], rotation='vertical')
  g.set_yticklabels(features+['score'], rotation='horizontal')
  
  # plt.show()
  plt.title('linear correlation matrix - fail')
  plt.tight_layout()
  plt.savefig(trainer.outdir+'/corr_fail.png' )
  
  signals = []
  signals.append(root2array(otherSig[0],"Events",branches= features))
  signals.append(root2array(otherSig[1],"Events",branches= features))
  test_bkg = root2array("../data/HNLFlatTuples/Parking1A0_newSig_4varsL_section0/HNLFlat_*12.root","Events",branches= features+["hnl_mass"])
  test_bkg.resize((10000,))


  signals[0] = np.array(signals[0].tolist(),dtype="float64")
  signals[1] = np.array(signals[1].tolist(),dtype="float64")

  plotScore(trainer.outdir,features,model,signals[0],qt.inverse_transform(x_test).compress(~y_test.astype(bool),axis=0),"mass2p0_ctau100",trainer.train_samples )
  plotScore(trainer.outdir,features,model,signals[1],qt.inverse_transform(x_test).compress(~y_test.astype(bool),axis=0),"mass3p0_ctau1",trainer.train_samples )












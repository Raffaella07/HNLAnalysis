import sys 
import glob
import numpy as np
import ROOT
from sklearn.preprocessing import RobustScaler
from NN_validationPlots import NNinput

modelSpecs = []
def LoadModel(p):

	print "_________________________________________________________Loading model" 
	global modelSpecs
	modelSpecs = NNinput(p)
	print "__________________________________________________________model loaded"


def ComputeNNScore(hnl_cos2D,hnl_vtxProb,hnl_pi_pt,B_mass,hnl_vtxChi2,hnl_l_pt,TrgMu_pt, hnl_lxy_sig,hnl_pi_DCAS,massParam):
	
	#modelSpecs = NNinput(p)
	
	point_nn = np.array([hnl_cos2D,hnl_vtxProb,hnl_pi_pt,B_mass,hnl_vtxChi2,hnl_l_pt,TrgMu_pt, hnl_lxy_sig,hnl_pi_DCAS,massParam])
	point_nn= point_nn.reshape(1, -1)
	model = modelSpecs[0]
	qt = modelSpecs[1]
	xx = qt.transform(point_nn)
	nn_score = model.predict(xx)
	print nn_score
	return nn_score

#print modelSpecs
#ComputeNNScore(0.999986,0.697965,5.839892,5.273438,0.150597,2.334133,8.726894,83.153740,7.763214,3.000000)
#print len(sys.argv)
#print ComputeNNScore(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13])

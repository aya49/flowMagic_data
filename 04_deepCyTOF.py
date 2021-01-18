## date: 2021-01-12
## author: deepCyTOF authors, modified by alice yue
## input: 2D matrices (cell x marker) + their clr (cell x cell population true/false matrix)
## ouput: deepCyTOF results


'''
This script will train a DeepCyTOF neural network classifier on three CyTOF
datasets. For each dataset, DeepCyTOF will (1) pick a reference sample and
train a feed-forward neural network classifier on the reference sample;
(2) train a MMD-ResNet for each pair of the remaining sample and the reference
sample; (3) predict the cell label information for each remaining sample.
Created on Jul 30, 2016
@author: urishaham, huaminli
'''

root = '/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data'

'''
packages to install:
- keras
- time
- fcsparser
- sklearn
- matplotlib
- numpy
- pandas
'''

import os

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = '-1'

from keras import backend as K
import numpy as np
import os.path
import time
import gzip, csv
from pathlib import Path

from glob import glob
import pandas as pd

os.chdir(root + "/src")

from Util import CostFunctions as cf
from Util import DataHandler as dh
from Util import denoisingAutoEncoder as dae
from Util import FileIO as io
from Util import feedforwadClassifier as net
from Util import MMDNet as mmd

os.chdir(root)


'''
Parameters.
dataSet           - A list of names of CyTOF data sets.
isCalibrate       - A indicator whether perform calibration or not.
denoise           - A indicator whether run de-noising auto-encoder or not.
loadModel         - A indicator whether load previous trained model or not.
hiddenLayersSizes - A list of 3 numbers to indicate the number of hidden
                    nodes in each hidden layer.
activation        - The activation function used to construct the classifier.
l2_penalty        - The regularization parameter to construct the classifier.
keepProb          - The keep probability for each single cell to be used as
                    the training set of de-noising autoencoder.
'''

isCalibrate = False
denoise = False
loadModel = False

hiddenLayersSizes = [12, 6, 3]
activation = 'softplus'
l2_penalty = 1e-4

'''
The user needs to specify the data set to run DeepCyTOF.
Make your choice here - an integer from 0 to 4.
0: MultiCenter study
1: 56 unnormalized
2: 56 normalized
3: 134 unnormalized
4: 134 normalized
'''

## load data set folders for nD and 2D data sets

# nD data folders
data_path_nD = 'data/nD/x'
data_dirs_nD = os.listdir(data_path_nD)
for i in range(len(data_dirs_nD)):
  data_dirs_nD[i] = data_path_nD + "/" + data_dirs_nD[i]

#2D data folders
data_path_2D = 'data/2D/x'
data_dirs_2D1 = os.listdir(data_path_2D)
data_dirs_2D = []
for i in range(len(data_dirs_2D1)):
  fold = data_path_2D + "/" + data_dirs_2D1[i]
  fold_ = os.listdir(fold)
  for fold_j in fold_:
    data_dirs_2D = data_dirs_2D + [fold + "/" + fold_j]

data_dirs = data_dirs_2D + data_dirs_nD
# for data_dir in data_dirs: ###################################################
data_dir = data_dirs_nD[0]
# data_dir = "src/MultiCenter_16sample"

data_paths = [y for x in os.walk(data_dir) for y in glob(os.path.join(x[0], '*.csv.gz'))]
actual_paths = [x.replace("/x/","/y/") for x in data_paths]

data = pd.read_csv(data_paths[0])
actual = pd.read_csv(actual_paths[0].replace("/x/", "/y/"))

# parameters
dataIndex = np.array(os.listdir(data_dir))
testIndex = dataIndex
trainNum = 10 # number of train samples
trainIndex = np.round(np.linspace(1, len(testIndex) - 1, trainNum)).astype(int)
trainIndex = dataIndex[trainIndex]
relevantMarkers = np.asarray(range(len(data.columns)))
mode = 'CSV.GZ'
numClasses = len(actual.columns)
keepProb = .8
'''
dataIndex = np.arange(1, 16 + 1)
trainIndex = dataIndex
testIndex = dataIndex
relevantMarkers = np.asarray([1, 2, 3, 4, 5, 6, 7, 8]) - 1
mode = 'CSV'
numClasses = 4
keepProb = .8
'''

'''
Choose the reference sample.
'''
print('Choose the reference sample')
refSampleInd = dh.chooseReferenceSample(data_dir, trainIndex, relevantMarkers, mode)
target = dh.loadDeepCyTOFData(data_dir, trainIndex[refSampleInd], relevantMarkers, mode)

# for data_path in data_paths:
#   data = pd.read_csv(data_path, compression='gzip', error_bad_lines=False)
#   actual = pd.read_csv(data_path.replace("/x/","/y/"), compression='gzip', error_bad_lines=False)

# Pre-process sample. Don't need to, my samples are cleaned and processed
target = dh.preProcessSamplesCyTOFData(target)

'''
Train the de-noising auto encoder.
'''
print('Train the de-noising auto encoder.')
res_dir = data_dir.replace("data/","results/").replace("/x/","/deepCyTOF_models/")
Path(res_dir).mkdir(parents=True, exist_ok=True)
DAE = dae.trainDAE(target, data_dir, refSampleInd, trainIndex,
                   relevantMarkers, mode, keepProb, denoise,
                   loadModel, res_dir)
denoiseTarget = dae.predictDAE(target, DAE, denoise)

'''
Train the feed-forward classifier on (de-noised) target.
'''
denoiseTarget, preprocessor = dh.standard_scale(denoiseTarget, preprocessor=None)

res_dir = data_dir.replace("data/", "results/").replace("/x/", "/deepCyTOF_models/")
Path(res_dir).mkdir(parents=True, exist_ok=True)
if loadModel:
  from keras.models import load_model
  cellClassifier = load_model(os.path.join(io.DeepLearningRoot(), res_dir + '/cellClassifier.h5'))
else:
  print('Train the classifier on de-noised Target')
  cellClassifier = net.trainClassifier(denoiseTarget, mode, refSampleInd,
                                       hiddenLayersSizes, activation, l2_penalty,
                                       res_dir)

'''
Test the performance with and without calibration.
'''
# Generate the output table.
dim = 2 if isCalibrate else 1
acc = np.zeros((testIndex.size, dim), np.float16)
F1 = np.zeros((testIndex.size, dim), np.float16)
mmd_before = np.zeros(testIndex.size)
mmd_after = np.zeros(testIndex.size)

for i in np.arange(testIndex.size): #################################################
  # i = 9
  # Load the source.
  sourceIndex = testIndex[i]
  source = dh.loadDeepCyTOFData(data_dir, sourceIndex, relevantMarkers, mode)
  source = dh.preProcessSamplesCyTOFData(source)
  
  # De-noising the source.
  denoiseSource = dae.predictDAE(source, DAE, denoise)
  denoiseSource, _ = dh.standard_scale(denoiseSource, preprocessor=preprocessor)
  
  # Predict the cell type of the source.
  print('Run the classifier on source ', str(sourceIndex),
        'without calibration')
  
  start = time.time()
  acc[i, 0], F1[i, 0], predLabel = net.prediction(denoiseSource, mode, i, cellClassifier)
  end = time.time()
  print(end - start)
  
  new_fold = data_dir.replace("data","results").replace("/x/","/deepCyTOF_labels/")
  Path(new_fold).mkdir(parents=True, exist_ok=True)
  np.savetxt(new_fold + "/" + str(testIndex[i]).replace(".gz",""), predLabel, delimiter=',', fmt='%i')
  
  sourceInds = np.random.randint(low=0, high=source.X.shape[0], size=1000)
  targetInds = np.random.randint(low=0, high=target.X.shape[0], size=1000)
  mmd_before[i] = K.eval(cf.MMD(denoiseSource.X, denoiseTarget.X).cost(
    K.variable(value=denoiseSource.X[sourceInds]),
    K.variable(value=denoiseTarget.X[targetInds])) )
  
  # f = open(dataPath + "/predlabel_nocal" + str(sourceIndex) + ".csv", 'w')
  # for item in predLabel:
  #     f.write(str(item.astype(int)) + '\n')
  # f.close()
  
  print('MMD before: ', str(mmd_before[i]))
  if isCalibrate:
    if loadModel:
      calibMMDNet = mmd.loadModel(denoiseTarget, denoiseSource,
                                  sourceIndex, predLabel, new_fold)
      calibrateSource = Sample(calibMMDNet.predict(denoiseSource.X),
                               denoiseSource.y)
      calibMMDNet = None
    else:
      calibrateSource = mmd.calibrate(denoiseTarget, denoiseSource,
                                      sourceIndex, predLabel, new_fold)
    
    print('Run the classifier on source ', str(sourceIndex),
          'with calibration')
    acc[i, 1], F1[i, 1], predLabell = net.prediction(calibrateSource,
                                                     mode, i,
                                                     cellClassifier)
    
    # f = open(dataPath + "/predlabel_cal" + str(sourceIndex) + ".csv", 'w')
    # for item in predLabell:
    #     f.write(str(item.astype(int)) + '\n')
    # f.close()
    
    mmd_after[i] = K.eval(cf.MMD(calibrateSource.X, denoiseTarget.X).cost(
      K.variable(value=calibrateSource.X[sourceInds]),
      K.variable(value=denoiseTarget.X[targetInds])))
    print('MMD after: ', str(mmd_after[i]))
    calibrateSource = None
  
  source = None
  denoiseSource = None
# end for loop
'''
Output the overall results.

CI = np.zeros(10000)
for i in range(10000):
  CI[i] = np.mean(np.random.choice(F1[:, 0], size=30))
CI = np.sort(CI)
print(mode, ', ', np.mean(CI), ' (', CI[250], CI[9750], ')')

if isCalibrate:
  CI = np.zeros(10000)
  for i in range(10000):
    CI[i] = np.mean(np.random.choice(F1[:, 1], size=30))
  CI = np.sort(CI)
  print(mode, ', ', np.mean(CI), ' (', CI[250], CI[9750], ')')
'''
#!/usr/bin/env python3
'''
Created on Oct 10, 2016

@author: huaminli
'''

import fcsparser
import numpy as np
from numpy import genfromtxt
from Util import FileIO as io
from Util import MMDNet as mmd
import os.path
import sklearn.preprocessing as prep
from sklearn.model_selection import train_test_split
import pandas as pd


class Sample:
    X = None
    y = None
    def __init__(self, X, y = None):
        self.X = X
        self.y = y

def preProcessSamplesCyTOFData(sample):
    sample.X = np.log(1 + np.abs(sample.X))    
    return sample

def standard_scale(sample, preprocessor = None):
    if preprocessor == None:
        preprocessor = prep.StandardScaler().fit(sample.X)
    sample.X = preprocessor.transform(sample.X)
    
    return sample, preprocessor

def loadDeepCyTOFData(dataPath, dataIndex, relevantMarkers, mode, skip_header = 0):
    X = []
    if mode == 'CSV.GZ':
        data_filename = dataPath + "/" + str(dataIndex) # I'm just going to give it the file name
        X = pd.read_csv(os.path.join(io.DeepLearningRoot(),data_filename)).values.tolist()
        print(np.shape(X))
        actual = pd.read_csv(os.path.join(io.DeepLearningRoot(),data_filename.replace("/x/","/y/")))
        labels = pd.DataFrame([0] * len(actual))
        for aci in range(len(actual.columns)):
            labels[actual[actual.columns[aci]] == 1] = aci + 1
        labels = [item for sublist in labels.values.tolist() for item in sublist]

    else:
        if mode == 'CSV':
            data_filename = dataPath + '/sample' + str(dataIndex)+'.csv'
            X = genfromtxt(os.path.join(io.DeepLearningRoot(),data_filename), delimiter=',', skip_header=skip_header)
        if mode == 'FCS':
            data_filename = dataPath + '/sample' + str(dataIndex)+'.fcs'
            _, X = fcsparser.parse(os.path.join(io.DeepLearningRoot() ,data_filename), reformat_meta=True)
            X = X.as_matrix()
        label_filename = dataPath + '/labels' + str(dataIndex) + '.csv'
        labels = genfromtxt(os.path.join(io.DeepLearningRoot(), label_filename), delimiter=',')
    labels = np.int_(labels)

    X = X[:, relevantMarkers]    
    sample = Sample(X, labels)
    
    return sample

def splitData(sample, test_size):
    data_train, data_test, label_train, label_test = train_test_split(sample.X, sample.y, test_size=test_size) 
    
    trainSample = Sample(data_train, label_train)
    testSample = Sample(data_test, label_test)
    return trainSample, testSample

def chooseReferenceSample(dataPath, dataIndex, relevantMarkers, mode):
    samples = []
    for i in dataIndex:
        sample = loadDeepCyTOFData(dataPath, i, relevantMarkers, mode)
        sample = preProcessSamplesCyTOFData(sample)
        samples.append(sample)
        
    numSamples = len(samples)
    norms = np.zeros(shape = [numSamples, numSamples])
    for i in range(numSamples):
        cov_i = np.cov(samples[i].X, rowvar = False)
        for j in range(numSamples):
            cov_j = np.cov(samples[j].X, rowvar = False)
            covDiff = cov_i - cov_j
            norms[i,j] = np.linalg.norm(covDiff, ord = 'fro')
            norms[j,i] = norms[i,j]
            avg = np.mean(norms, axis = 1)
            refSampleInd = np.argmin(avg)
        
    return refSampleInd
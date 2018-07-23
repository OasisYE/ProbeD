#!/usr/bin/env python
# -*-coding:utf-8 -*-
# ** author:Oasis
# *****************************
import re
from Bio.SeqUtils import MeltingTemp as mt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import numpy as np
from Bio.SeqUtils import GC

# LDA model information. 模型参数
tempList = [32, 37, 42, 47, 52, 57]
coefList = [[[-0.14494789, 0.18791679, 0.02588474]],
                  [[-0.13364364, 0.22510179, 0.05494031]],
                  [[-0.09006122, 0.25660706, 0.1078303]],
                  [[-0.01593182, 0.24498485, 0.15753649]],
                  [[0.01860365, 0.1750174, 0.17003374]],
                  [[0.03236755, 0.11624593, 0.24306498]]]
interList = [-1.17545204, -5.40436344, -12.45549846,
                   -19.32670233, -20.11992898, -23.98652919]
classList = [-1, 1]

# Convert lists to ndarrays.
coefArray = np.asarray(coefList)
interArray = np.asarray(interList)
classArray = np.asarray(classList)

# Determine which index to reference for model values.
# default
tempVal = 57
np_index = tempList.index(tempVal)

# Build model from encoded values.
clf = LinearDiscriminantAnalysis()
clf.coef_ = coefArray[np_index]
clf.intercept_ = interArray[np_index]
clf.classes_ = classArray

# Determine which classifier parameters to use.
clfT = tempList.index(tempVal)

# Make lists to hold data about candidates.
testList = []
testSet = set()
candsInfo = []


# the realignment result
samf = "/Users/yeweijian/Downloads/data/t.sam"

# Make a list to hold the output.
outList = []

# Open input file for reading
with open(samf, 'r') as f:
    file_read = [line.strip() for line in f]

# Process .sam file
for i in range(0, len(file_read), 1):
    if file_read[i][0] is not '@':
        chromField = file_read[i].split('\t')[2]
        chrom = file_read[i].split('\t')[0].split(':')[0]
        start = file_read[i].split('\t')[0].split(':')[1].split('-')[0]
        stop = file_read[i].split('\t')[0].split('-')[1].strip(' ')
        seq = file_read[i].split('\t')[9]
        Tm = mt.Tm_NN(seq)

        if re.match('\*', chromField) is None \
         and re.search('120M', file_read[i]) is not None:
            t = [float(len(seq)),
                 float(file_read[i].split('\t')[12].split(':')[2]),
                 GC(seq)]
            testList.append(t)
            print ('uniq alignment!')
            print ('%s\t%s\t%s\t%s\t%s' \
                   % (chrom, start, stop, seq, Tm))
        else:
            print ('aligned 0 time')


# Make ndarray for input into classifier.
testArray = np.asarray(testList)

# Create classifier
clf = LinearDiscriminantAnalysis()

# Load temperature-specific model information.
clf.coef_ = coefArray[clfT]
clf.intercept_ = interArray[clfT]
clf.classes_ = classArray

probs = clf.predict_proba(testArray)[:, 1]

for i in range(0, len(probs), 1):
    print ('The LDA score is: {}'.format(probs[i]))

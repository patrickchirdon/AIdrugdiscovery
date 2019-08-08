import numpy as np
from matplotlib import cm
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Draw import SimilarityMaps
from sklearn import svm
from sklearn.model_selection import train_test_split
import os
from sklearn.externals import joblib
# dense to sparse
from numpy import array
from scipy.sparse import csr_matrix

import keras
import tensorflow as tf
from sklearn import preprocessing
from keras import regularizers
abscaler= preprocessing.MaxAbsScaler()



def process(dropout, regularizers, lr, momentum, decay, nesterov, myloss, epochs, batchsize, filename, modelname):
     
    dataset = open(filename, "r" )
    dataset = [ line.rstrip().split(",") for line in dataset ][0:]
    mols = [ Chem.MolFromSmiles( line[0] ) for line in dataset ]


    u=0
    indexy=[]
    for y in mols:
        if y is not None:
            indexy.append(u)
            u+=1
        
    goodmols=[mols[k] for k in indexy]
    
    
    
   

    Y=[line[1] for line in dataset]
    goody=[Y[k] for k in indexy]
    trainmols, testmols, trainy, testy = train_test_split(goodmols, goody, test_size = 0.1, random_state=90  )
   
    #trainmols, testmols, trainy, testy = train_test_split(goodmols, Y, test_size = 0.1, random_state=90  )
    trainfps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in goodmols if m is not None ]
    testfps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in testmols if m is not None]
    
    
   
    u=0
    indexy=[]
    for y in trainfps:
        if y is not None:
            indexy.append(u)
            u+=1
    
    newy=array([int(goody[k]) for k in indexy])
    print(len(newy))
    
    np_fptrain=[]
    for fp in trainfps:
        arr=np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fptrain.append(arr)

    np_fptest=[]
    for fp in testfps:
        arr=np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fptest.append(arr)  
    
    
    a=csr_matrix(np_fptrain, dtype=np.int8).toarray()
    b=csr_matrix(np_fptest, dtype=np.int8).toarray()
    a=abscaler.fit_transform(a)
    b=abscaler.fit_transform(b)
    
    print(len(a))
    
    import matplotlib.pyplot as plt
    from sklearn.metrics import mean_squared_error, r2_score

    from sklearn import datasets

   

# In[ ]:


    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_validate
    from keras.models import Sequential
    from keras import layers
    from keras.optimizers import SGD




    from keras.models import Sequential
 
    from keras.callbacks import ReduceLROnPlateau
    from keras import regularizers as WeightRegularizer
    from keras.optimizers import SGD
    #SKlearn for metrics and datasplits
    from sklearn.model_selection import cross_validate
    from sklearn.metrics import roc_auc_score, roc_curve
    #Matplotlib for plotting
    from matplotlib import pyplot as plt
    
    mylist=[]



    model = Sequential()
    model.add(layers.Dense(236, input_dim=2048, activation='relu', kernel_regularizer=keras.regularizers.l2(l=regularizers)))
    model.add(layers.BatchNormalization())
    model.add(layers.Dropout(dropout))
    model.add(layers.Dense(236, activation='relu', kernel_regularizer=keras.regularizers.l2(l=regularizers)))
    model.add(layers.Dense(1,  activation='sigmoid'))
    
    
    
    sgd=SGD(lr, momentum, decay, nesterov)
    model.compile(loss=myloss, optimizer=sgd, metrics=['accuracy'])
    model.fit(a, newy, nb_epoch=epochs, batch_size=batchsize)
    
    
    model.save(modelname)
    
    
    
   
##########################################

    
    
    
    
 

    
    

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
abscaler= preprocessing.MaxAbsScaler()



def process(filename):
    
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
    




    model = Sequential()
    model.add(layers.Dense(236, input_dim=2048, activation='relu'))
    model.add(layers.BatchNormalization())
    model.add(layers.Dropout(.2))
    model.add(layers.Dense(236, activation='relu'))
    model.add(layers.Dense(1,  activation='sigmoid'))
    
    
    
    sgd=SGD(lr=.01, momentum=.8, decay=0.0, nesterov=False)
    model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])
    model.fit(a, newy, nb_epoch=7, batch_size=3)
    
    
    model.save('classifier.h5')
    from sklearn.metrics import confusion_matrix
    from keras.models import load_model
    model = load_model('classifier.h5')
   
    predy=model.predict(b)
    predy[predy > 0.5] = 1
    predy[predy <= 0.5] = 0
    testy=np.delete(testy, -1)
    testb=[int(h) for h in testy]
    print(confusion_matrix(testb, predy))
    
    
   
##########################################

process('gsk3.csv')

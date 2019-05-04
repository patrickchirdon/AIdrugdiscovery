#!/usr/bin/env python
# coding: utf-8

# In[3]:


from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import pandas as pd
import gzip

mols = pd.read_csv("solubility.csv", sep=" ")
mymols = [x for x in mols['smiles']  if x is not None]


moleculelist=[]
for u in mymols:
    if u is not None:
        moleculelist.append(Chem.MolFromSmiles(u))

moleculelist
fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in moleculelist if m is not None]
fps
len(fps)

myscores=[]
i=0
for y in moleculelist:
    if y is not None:
        myscores.append(i)    
    i+=1


thescores=[]
for g in myscores:
    thescores.append(mols['solubility'][g])




np_fp=[]
for fp in fps:
    arr=np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    np_fp.append(arr)

np_fp

from sklearn.linear_model import LinearRegression


np_scores=np.array(thescores)

np_scores


import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score

from sklearn import datasets

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(np_fp, np_scores, test_size=0.33, random_state=42)



# In[ ]:


from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_validate
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.optimizers import SGD




from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.callbacks import ReduceLROnPlateau
from keras import regularizers as WeightRegularizer
from keras.optimizers import SGD
#SKlearn for metrics and datasplits
from sklearn.model_selection import cross_validate
from sklearn.metrics import roc_auc_score, roc_curve
#Matplotlib for plotting
from matplotlib import pyplot as plt


# In[11]:


#Convert to Numpy arrays
from sklearn.preprocessing import StandardScaler
X = np.array(list(X_train))
y = y_train


# In[12]:


train_size=int(.7* X.shape[0])
X_train, X_test, y_train, y_test = X[0:train_size], X[train_size:],y[0:train_size],y[train_size:]



model = Sequential()
model.add(Dense(output_dim=5, input_dim=X.shape[1]))
model.add(Activation("sigmoid"))
model.add(Dense(output_dim=1))
model.add(Activation("linear"))

# In[34]:
model.compile(loss='mean_squared_error', optimizer=SGD(lr=0.01, momentum=0.9, nesterov=True))
history = model.fit(X_train, y_train, nb_epoch=500, batch_size=32)
y_pred = model.predict(X_test)
rms = (np.mean((y_test.reshape(-1,1) - y_pred)**2))**0.5
s = np.std(y_test -y_pred)
print("Neural Network RMS", rms)

print(s)

import matplotlib.pyplot as plt
plt.scatter(y_train,model.predict(X_train), label = 'Train', c='blue')
plt.title('Neural Network Predictor')
plt.xlabel('Measured Solubility')
plt.ylabel('Predicted Solubility')
plt.scatter(y_test,model.predict(X_test),c='lightgreen', label='Test', alpha = 0.8)
plt.legend(loc=4)
plt.show()

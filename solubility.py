from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
import os
from sklearn.externals import joblib
import glob
import numpy
import re
import os
import pandas as pd
import warnings; warnings.simplefilter('ignore')
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')

def solubilitypredictor(myinput):
    os.remove('output.txt')
    mylist=myinput.split('\n')
    for j in mylist:
        
        smile=j
        
        themol=Chem.MolFromSmiles(smile)
        estimator=joblib.load('solubility.pickle')
        fp = AllChem.GetMorganFingerprintAsBitVect(themol,2,nBits=2048)
        res = numpy.zeros(len(fp),numpy.int32)
        DataStructs.ConvertToNumpyArray(fp,res)
        probas = list(estimator.predict(res.reshape(1,-1))[0])
        f1=open('output.txt', 'a')
        print(smile, file=f1)
        print(probas, file=f1)



from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
import os
from sklearn.externals import joblib
import glob
import numpy
import re
import pandas as pd
import dropbox

thetargets=[]
theprobas=[]
theincorrect=Chem.SDMolSupplier("incorrect.sdf")

def getProba(fp, predictionFunction):
      return predictionFunction(fp)[0][1]
#themolecules=Chem.SDMolSupplier("validate.sdf")
smile="COc(cc1)cc(C(CS2)(N(Cc3cnccc3)C2=S)O)c1OC"

name="cob197"
mol=Chem.MolFromSmiles(smile)
info={}
def predictme(molecule):
    
    fp = AllChem.GetMorganFingerprintAsBitVect(molecule,2,nBits=2048, bitInfo=info)
    res = numpy.zeros(len(fp),numpy.int32)
    DataStructs.ConvertToNumpyArray(fp,res)
    
    for file in glob.glob('*.pkl'):
        
        line = re.sub('.pkl', '', file)
        estimator=joblib.load(file)
        probas = list(estimator.predict_proba(res.reshape(1,-1))[0])
        
        print(line, probas[1])
        
    

    for file in glob.glob('*.pkl'):
        line = re.sub('.pkl', '', file)
        estimator=joblib.load(file)
        fig, maxweight = SimilarityMaps.GetSimilarityMapForModel(mol, SimilarityMaps.GetMorganFingerprint, lambda x: getProba((x,), estimator.predict_proba))
        thename=name + '_' + line + '.png'
        probas = list(estimator.predict_proba(res.reshape(1,-1))[0])
        print(line, probas[1])
        print(line, probas[1], file=open("output/thepositives.txt", "a"))
        thetargets.append(line)
        theprobas.append(probas[1])
        #thelocation="output2/" +thename
        #fig.savefig(thelocation, bbox_inches='tight')


from rdkit.Chem import BRICS

from rdkit import Chem

from rdkit.Chem import DataStructs

from rdkit.Chem import AllChem

from rdkit.Chem import Descriptors

import numpy as np

import random

import os

import pandas as pd

 

def fragmenter(thefile): 
    os.remove('output.txt')
    id=[]

 

    for line in open(thefile):

        line=line.strip()

        id.append(line)

 

    

    df=pd.DataFrame()

 

    df=id

 

    count=0

    mylist=[]

 

    for y in df:

        base = Chem.MolFromSmiles(df[count])

    

        catalog = BRICS.BRICSDecompose(base)

        mcat = [Chem.MolFromSmiles(x) for x in catalog]

        ms = BRICS.BRICSBuild(mcat)

        for m in ms: 

            a=Chem.MolToSmiles(m)

            mylist.append(a)

        count=count +1

        

 

    df2 = pd.DataFrame({'smiles':mylist})
    f3=open('output.txt', 'w+')
    for j in mylist:
        print(j, file=f3)
    
    f3.close()
    return mylist
    


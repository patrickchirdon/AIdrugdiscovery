from rdkit.Chem import BRICS

from rdkit import Chem

from rdkit.Chem import DataStructs

from rdkit.Chem import AllChem

from rdkit.Chem import Descriptors

import numpy as np

import random

import pandas as pd

 

 

id=[]

 

for line in open("top500pocket3.smi"):

    line=line.strip()

    id.append(line)

 

    

df=pd.DataFrame()

 

df=id

 

count=1

mylist=[]

 

for y in df:

    base = Chem.MolFromSmiles(df[count])

    

    catalog = BRICS.BRICSDecompose(base)

    mcat = [Chem.MolFromSmiles(x) for x in catalog]

    ms = BRICS.BRICSBuild(mcat)

    for m in ms: 

        a=Chem.MolToSmiles(m)

        print(a)

    count=count +1

    mylist.append(a)

 

df2 = pd.DataFrame({'smiles':mylist})

df2.to_csv('mynewcompounds.csv')

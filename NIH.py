from rdkit.Chem.FilterCatalog import *
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, PandasTools, rdMolDescriptors, Lipinski
from rdkit.Chem.Draw import IPythonConsole
from rdkit import RDConfig
import pandas as pd
import numpy as np
import warnings; warnings.simplefilter('ignore')
mySMILESinput=pd.DataFrame(columns=['ID','my_smiles'])

params=FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
catalog=FilterCatalog(params)
suppl=Chem.SmilesMolSupplier("fda.smi")
inf=open('fda.smi','r')
sub_strct=[line.rstrip().split(" ") for line in inf]

ms = [x for x in suppl if x is not None]
i=0

for mol in ms:
    entry=catalog.GetFirstMatch(mol)
    sphybrid=Chem.rdMolDescriptors.CalcFractionCSP3(mol)
    if(entry is not None):
        print(i,sub_strct[i], "PAINS", entry.GetDescription(), "Fsp3", sphybrid)
    else:
        print(i,sub_strct[i], "PAINS OK", "Fsp3", sphybrid)
    i+=1

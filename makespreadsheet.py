import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

def makespreadsheet(thesmiles):
    mylist=thesmiles.split('\n')
    df = pd.DataFrame({'SMILES':mylist})

    df['Mol Image'] = [Chem.MolFromSmiles(s) for s in df['SMILES']]

    PandasTools.SaveXlsxFromFrame(df, 'test.xlsx', molCol='Mol Image')


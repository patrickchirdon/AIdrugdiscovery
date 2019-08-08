import numpy as np
from matplotlib import cm
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Draw import SimilarityMaps
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
import os
from sklearn.externals import joblib
def process(dropout, regularizers, lr, momentum, decay, nesterov, myloss, epochs, batchsize, filename, modelname):
    os.remove('output.txt')
    dataset = open(filename, "r" )
    dataset = [ line.rstrip().split("\t") for line in dataset ][0:]
    mols = [ Chem.MolFromSmiles( line[0][0] ) for line in dataset ]


    u=0
    indexy=[]
    for y in mols:
        if y is not None:
            indexy.append(u)
            u+=1
        
    goodmols=[mols[k] for k in indexy]


    Y= ([x[1] for x in dataset])
    myvalues=[Y[j] for j in indexy]
    Y=myvalues
   
    Y = [float(i) for i in Y]
    Y = np.asarray( Y )

    trainmols, testmols, trainy, testy = train_test_split(goodmols, Y, test_size = 0.1, random_state=90  )

# calc fingerprint

    trainfps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in trainmols if m is not None ]
    testfps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in testmols if m is not None]
    u=0
    indexy=[]
    for y in trainfps:
        if y is not None:
            indexy.append(u)
            u+=1
    goodY=[trainy[k] for k in indexy]

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

    cls = SVC( probability=True, C=500000)
    cls.fit( trainfps, goodY )
    pkl_name=filename + ".pkl"
    joblib.dump(cls, pkl_name)
    
        
    def getProba( fp, probabilityFunc ):
        return probabilityFunc( fp )[0][1]
 
    def mapperfunc( mol ):
        fig, weight = SimilarityMaps.GetSimilarityMapForModel( mol, SimilarityMaps.GetMorganFingerprint, lambda x: getProba( x, cls.predict_proba), colorMap=cm.bwr  )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 2 )
        print(fp)
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray( fp, arr )
        print(arr)
        res = cls.predict( arr )
        smi = Chem.MolToSmiles( mol )
        print(smi)
        
        if res[0] == 1:
            fig.savefig( "res/act_"+smi+"_.png", bbox_inches = "tight" )
        else:
            fig.savefig("res/nonact_"+smi+"_.png", bbox_inches = "tight" )
    m=Chem.MolFromSmiles("CN1C=CNC1=S")
    
    

    res= cls.predict(testfps)

    
    u=0
    indexy=[]
    for y in res:
        indexy.append(u)
        u+=1
    goodY=[testy[k] for k in indexy]
    from sklearn.metrics import confusion_matrix
    cm1=confusion_matrix(goodY,res)
    f=open("output.txt", "a")
    print(filename, file=open("predictions.txt", "a"))
    print(filename)
    print(confusion_matrix(goodY, res), file=f)
    print(confusion_matrix(goodY, res))
    from sklearn.metrics import accuracy_score
    print("fraction of correct predictions", accuracy_score(goodY, res), file=f)
    print("fraction of correct predictions", accuracy_score(goodY, res))

    sensitivity1 = cm1[0,0]/(cm1[0,0]+cm1[0,1])
    print('Sensitivity : ', sensitivity1, file=f)
    print('Sensitivity: ', sensitivity1)
    specificity1 = cm1[1,1]/(cm1[1,0]+cm1[1,1])
    print('Specificity : ', specificity1, file=f)
    f.close()
    print('specificity: ', specificity1)
############################################################


        
        

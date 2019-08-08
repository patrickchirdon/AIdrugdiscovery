import numpy as np
import torch
import os
import pandas as pd
from torch import nn
from torch.autograd import Variable
import torch.nn.functional as F
import torch.utils.data as D
import tables as tb
from sklearn.metrics import (matthews_corrcoef, 
                                 confusion_matrix, 
                                 f1_score, 
                                 roc_auc_score,
                                 accuracy_score,
                                 roc_auc_score)
import csv
i=[]
with open('chembltargets.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        i.append(row)
                          
def chemblpredict(theinput):
    
    f1=open('output.txt', 'a')
    #!/usr/bin/env python
# coding: utf-8

    # In[15]:


    
    

    # In[23]:

    
    mylist=theinput.split('\n')
    for j in mylist:

        # set the device to GPU if available
        device = torch.device('cpu')
        v=j

        # In[17]:


        MAIN_PATH = '.'
        DATA_FILE = 'mt_data.h5'
        MODEL_FILE = 'chembl_mt.model'
        N_WORKERS = 8 # Dataloader workers, prefetch data in parallel to have it ready for the model after each batch train
        BATCH_SIZE = 32 # 
        LR = 2 # Learning rate. Big value because of the way we are weighting the targets
        N_EPOCHS = 2


        # In[18]:


        class ChEMBLDataset(D.Dataset):

            def __init__(self, file_path):
                self.file_path = file_path
                with tb.open_file(self.file_path, mode='r') as t_file:
                    self.length = t_file.root.fps.shape[0]
                    self.n_targets = t_file.root.labels.shape[1]

            def __len__(self):
                return self.length

            def __getitem__(self, index):
                with tb.open_file(self.file_path, mode='r') as t_file:
                    structure = t_file.root.fps[index]
                    labels = t_file.root.labels[index]
                return structure, labels


        dataset = ChEMBLDataset(f"{MAIN_PATH}/{DATA_FILE}")
        validation_split = .2
        random_seed= 42

        dataset_size = len(dataset)
        indices = list(range(dataset_size))
        split = int(np.floor(validation_split * dataset_size))

        np.random.seed(random_seed)
        np.random.shuffle(indices)
        train_indices, test_indices = indices[split:], indices[:split]

        train_sampler = D.sampler.SubsetRandomSampler(train_indices)
        test_sampler = D.sampler.SubsetRandomSampler(test_indices)

        # dataloaders can prefetch the next batch if using n workers while
        # the model is tranining
        train_loader = torch.utils.data.DataLoader(dataset,
                                                   batch_size=BATCH_SIZE,
                                                   num_workers=N_WORKERS,
                                                   sampler=train_sampler)

        test_loader = torch.utils.data.DataLoader(dataset, 
                                                  batch_size=BATCH_SIZE,
                                                  num_workers=N_WORKERS,
                                                  sampler=test_sampler)


        # In[19]:


        class ChEMBLMultiTask(nn.Module):
            """
            Architecture borrowed from: https://arxiv.org/abs/1502.02072
            """
            def __init__(self, n_tasks):
                super(ChEMBLMultiTask, self).__init__()
                self.n_tasks = n_tasks
                self.fc1 = nn.Linear(1024, 2000)
                self.fc2 = nn.Linear(2000, 100)
                self.fc3 = nn.Linear(2000, 100)
                self.dropout = nn.Dropout(0.25)

                # add an independet output for each task int the output laer
                for n_m in range(self.n_tasks):
                    self.add_module(f"y{n_m}o", nn.Linear(100, 1))

            def forward(self, x):
                h1 = self.dropout(F.relu(self.fc1(x)))
                h2 = F.relu(self.fc2(h1))
                out = [torch.sigmoid(getattr(self, f"y{n_m}o")(h2)) for n_m in range(self.n_tasks)]
                return out

        # create the model, to GPU if available
        model = ChEMBLMultiTask(dataset.n_targets).to(device)

        # binary cross entropy
        # each task loss is weighted inversely proportional to its number of datapoints, borrowed from:
        # http://www.bioinf.at/publications/2014/NIPS2014a.pdf
        with tb.open_file(f"{MAIN_PATH}/{DATA_FILE}", mode='r') as t_file:
            weights = torch.tensor(t_file.root.weights[:])
            weights = weights.to(device)

        criterion = [nn.BCELoss(weight=w) for x, w in zip(range(dataset.n_targets), weights.float())]

        # stochastic gradient descend as an optimiser
        optimizer = torch.optim.SGD(model.parameters(), LR)


        # In[21]:





       


       


        model = ChEMBLMultiTask(560) # number of tasks
        model.load_state_dict(torch.load(f"./{MODEL_FILE}"))
        model.eval()


        # In[112]:


        from rdkit import Chem
        from rdkit.Chem import AllChem
        


        from rdkit import Chem, DataStructs
        


        # In[38]:


     
        # In[42]:


        def calc_fp(smiles, fp_size, radius):
            """
            calcs morgan fingerprints as a numpy array.
            """
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            mol.UpdatePropertyCache(False)
            Chem.GetSSSR(mol)
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_size)
            a = np.zeros((0,), dtype=np.float32)
            Chem.DataStructs.ConvertToNumpyArray(fp, a)
            return a


        # In[1]:


        FP_SIZE = 1024
        RADIUS = 2
        from rdkit.Chem import rdMolDescriptors


        # In[2]:

        
        u=calc_fp(j, 1024, 2)


        # In[130]:


        a=Variable(torch.from_numpy(u))


        # In[131]:

        b=model(a)
        
        print(str(v), file=f1)
        c=0
        thetargets =[]
        theprobas=[]
      


        for h in b:
            thetargets.append(str(i[c]).strip(',[]\'"'))
            g=np.array2string(Variable(h).data.cpu().numpy())[1:7]
            
            g=round(float(g), 0)
            
           
            theprobas.append(g)
            
            c+=1
            
        
        d = { 'Target':thetargets,'Probability':theprobas}
        df = pd.DataFrame(d)
        df2=df.sort_values(by="Probability", ascending=False)
        df2.to_csv(f1,  index=None, sep=' ')
        print('----------------------------------------------', file=f1)
         


   
    
   
    





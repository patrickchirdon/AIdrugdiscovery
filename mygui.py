import PySimpleGUI as sg
from NIH import painspredict
from solubility import solubilitypredictor
from fragmenter import fragmenter
from toxtest import toxtest
from chemblpredict import chemblpredict
from makespreadsheet import makespreadsheet

import warnings; warnings.simplefilter('ignore')
from rdkit import Chem, DataStructs
from rdkit.Chem.Draw import SimilarityMaps
import glob
from rdkit import rdBase
import os, io
from PIL import Image
from os import environ
from subprocess import call
from makeai import process
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import csv

def get_image_as_data(filename, width=None, height=None):
    # from PIL import Image         # use `pip install Pillow` to install PIL
    # import io
    im = Image.open(filename)
    if isinstance(width, int) and isinstance(height, int): # Resize if dimensions provided
        im = im.resize((width, height))
    im_bytes = io.BytesIO()
    im.save(im_bytes, format="PNG")
    return im_bytes.getvalue()
    
    
def get_image_files_list(folder):
    all_files = os.listdir(folder)
    image_files = []
    for file in all_files:
        extension = file.lower().split(".")[-1]
        print(file, extension)
        if extension in ["jpg", "png", "jpeg", "jpe"]:
            image_files.append(file)
    image_files.sort()
    return image_files
    
def demo_photo_picker3(default_folder, default_pic):
    folder = default_folder
    files_listing = get_image_files_list(folder)
    column1 = [
        [
            sg.Listbox(values=files_listing,
                change_submits=True, # trigger an event whenever an item is selected 
                size=(90, 30), 
                font=("Helvetica", 12),
                key="files_listbox")
        ]
    ]
    column2 = [
        [
            sg.Image( data=get_image_as_data(default_pic, 500, 500), 
                key="image", size=(500,500))
        ]
    ]
    layout = [
        [ 
            sg.Text("Select your photos folder"),
            sg.InputText(key="photo_folder", change_submits=True), # trigger an event whenever the item is changed 
            sg.FolderBrowse(target="photo_folder")
        ], [
            sg.Column( column1 ), 
            sg.Column( column2 )
        ], [
            sg.Button(button_text="Back")
        ]
    ]
    window = sg.Window('Pick a photo').Layout(layout)
    while True:
        event, values = window.Read()
        print(event)
        print(values)
        if event == "photo_folder":
            if values["photo_folder"] != "":
                if os.path.isdir(values["photo_folder"]):
                    folder = values["photo_folder"]
                    image_files = get_image_files_list(values["photo_folder"])
                    window.FindElement("files_listbox").Update(values=image_files)
                    if len(image_files) > 0:
                        full_filename = os.path.join(folder,image_files[0])
                        window.FindElement("image").Update(data=get_image_as_data(full_filename, 500, 500))
        if event == "files_listbox":
            full_filename = os.path.join(folder,values["files_listbox"][0])
            window.FindElement("image").Update(data=get_image_as_data(full_filename, 500, 500))
        if event is None or event == 'Back':
            return None


    
    


rdBase.DisableLog('rdApp.*')
menu_def=[['Help', ['About','Tutorial']]]
f2=open('output.txt', 'w+')
print("", file=f2)
f2.close()
form = sg.FlexForm('Biomedical Engineering Chemical Toolkit V 1.0', default_element_size=(80, 1))
menu_def=[['Help', ['About','Tutorial']]]
column1 = [[sg.Text('Column 1', background_color='#d3dfda', justification='center', size=(10,1))]]
folder="output"



layout = [
    
    [sg.Text('Create a Bioassay or Screen Your Compounds!', size=(80, 1), font=("Helvetica", 25))],
    [sg.Text('Enter the name of your molecule population here')],
    [sg.InputText('molecule name')],
    [sg.Multiline(default_text='Enter up to 100 smiles here. One smile only for toxicity map.', size=(70, 3))],
    [sg.Text('_'  * 80)],
    [sg.Listbox(values=('PAINS', 'Chembl Prediction', 'Toxicity Prediction', 'Toxicity Map', 'Build a SAR model', 'Predict solubility', 'Fragment Molecules', 'Make Spreadsheet'), size=(30, 8))],
    [sg.Output(size=(150,10), key='Output')], 
    [sg.Button('Submit')], [sg.Button('Clear')], [sg.Button('Images')], [sg.T('0%',key='output'), sg.ProgressBar(10000, orientation='h', size=(20,20), bar_color=('red', 'yellow green'), relief=sg.RELIEF_RAISED  ,key='progress')], [sg.Menu(menu_def, tearoff=True)]]
window=sg.Window('Biomedical Engineering Chemical Toolkit', layout)
progress_bar = window.FindElement('progress')
output = window.FindElement('output')
win2_active=False  
win3_active=False




  
    

while True:
    
    event, values = window.Read(timeout=100)
    
    if event is None:
        break
    if event == 'Images'  and not win2_active:
        win2_active = True
        window.Hide()
        layout2 = [[sg.Text('Go back to the toolkit?')],       
                   [sg.Button('Back')]]
                   
   
        default_pic = "IRF3.jpg"
        result = demo_photo_picker3(".", default_pic)
        print(result)

        win2 = sg.Window('Molecule Viewer').Layout(layout2)
        
        while True:
            ev2, vals2 = win2.Read()
            if ev2 is None or ev2 == 'Back':
                win2.Close()
                win2_active = False
                win3_active =False
                window.UnHide()
                window.FindElement('Output').Update('') 
                break
                
    if event == 'About':
        sg.Popup('About this program', 'Version 1.0', 'Created by Patrick Chirdon', 'at Ohio University', 'supported by Dr.Douglas Goetz and Dr.Sumit Sharma')
   
    if event == 'Tutorial':
        tutorial="Tutorial"
        spacer="======================================"
        targetprediction="Target Prediction"
        targetprediction1="46 models. Metrics for test data:"
        targetprediction2="accuracy 97.59 +/- 2.41"
        targetprediction3="sensitivity 91.9 +/- 8.1"
        targetprediction4="specificity 98.6 +/- 1.4"
        targetprediction5="215/247 87% correct mechanism on independent test set."
        chem="Chembl Target Prediction"
        chembl="Number of unique targets 560"
        chembl1="Ion channel 5"
        chembl2="kinase 96"
        chembl3="nuclear receptor 21"
        chembl4="GPCR 180"
        chembl5="Others 258"
        chembl6="accuracy .87"
        chembl7="auc .92"
        chembl8="sensitivity .76"
        chembl9="specificity .92"
        chembl10="precision .82"
        chembl11="225/225 100% correct mechanism on independent test set. Note-- 1 is considered positive and zero is negative for a given target."
        my_text="Interpreting output:  For the target predictions, the green represents a positive region for the molecule, the red represents a negative region of the molecule for a tested property, and gray represents no detection. For more on this method please read Similarity maps-- a visualization strategy for molecular fingerprints and machine learning methods."
        mytext2="Inside the target prediction folder, there should be .png images for each of the smiles in the output folder. Make sure to change the directory to the output directory of the targetprediction folder under the images menu. Since there are 46 models it is best to only use a few smiles at a time."
        mytext3="Creating your own models:  https://pubchem.ncbi.nlm.nih.gov/#query=interferon&tab=assay, Also see chembl bioassays. These assays must be saved as .txt files with two columns-- the first for the smiles and the next column for either 1 or zero (active and inactive respectively).  The text file with the smiles and 1's and 0's should be in the targetprediction folder.  The text file names should contain the name of the assay.  You want a model with both good sensitivity and specificity (as close to one as possible).  It is important to note that a model can appear highly accurate but if sensitivity is zero, then the model does not detect positives."
        mytext6="confusion matrix"
        mytext4="tn fp"
        mytext5="fn tp"
        mytext7="It is important to note that column 1, row 1 is NOT true positive as you might expect from stats class.  Sensitive models will not have 0 in the bottom right corner.  If you are not getting good sensitivity and specificity, then you may want to change the penalty C=500000 to some other value.  By default the SVC is set up to use a RBF best fit but this can be changed as per the scikit learn documentation.  The output files will be saved as .pkl files that can later be loaded for future use."
        mytext8="Pan Assay Interference"
        mytext9="See Seven Year Itch: Pan-Assay Interference Compounds (PAINS) in 2017—Utility and Limitations New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. Pan Assay Interference Compounds commonly result in false positives in biological screening assays."
        mytext10="Since they bind everything, they are not selective and therefore do not make good drug targets.  We found that the higher the drug score in Data Warrior http://www.openmolecules.org/datawarrior/ the lower the frequency of compounds containing PAINS.  Using data warrior’s evolutionary algorithm (be sure to use the wand tool if you want to fix the scaffold), evolve a few runs by taking the compounds with the top drug scores (macro → run macro → calculate properties) by taking the top 5 scoring compounds as starting points for evolution until you get drug scores greater than .9.  Select based on skelsphere similarity and the algorithm will generate a large number of compounds that have high drug scores, which are oftentimes painless." 
        mytext11="The program will tell you what functional groups for each compound were responsible for a positive PAINFUL test result.  The program also tells you the fraction of SP3 hybridized carbons.  Compounds with scores > .47 are more selective binders.  Note that double bonds reduce the fraction of sp3 hybridization, as they make the compound more flat.  See Escape from flatland: increasing saturation as an approach to improving clinical success. Pains are defined as follows:"
        mytext12="Doveston R, et al. A Unified Lead-oriented Synthesis of over Fifty Molecular Scaffolds. Org Biomol Chem 13 (2014) 859D65. doi:10.1039/C4OB02287D"
        mytext13="Jadhav A, et al, Quantitative Analyses of Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for Inhibitors of Thiol Protease.  J Med Chem 53 (2009) 37D51. doi:10.1021/jm901070c"
        mytext14="Fragmenter"
        mytext15="Input a list of smiles.  These will be recombined into new combinations.  When you take the lowest energy ligands from a docking program and recombine these there may be some compounds that bind with lower energy than the original."
        svm="Support Vector Machine"
        multiclass="Multiclass Classifier"
        spreadsheet="Make Spreadsheet"
        spreadsheet1="Input smiles.  The output will be a spreadsheet called test.xlsx in the target prediction folder that contains images of the molecules."
        solubility="Solubility"
        solubility1="Predicts log S.  Log S greater than -4 is soluble."
        solubility2="Root mean square error of 1.27 on a scale from -4 to 4."
        solubility3="linear regression"
        sar="Build a SAR model"
        sar1="cross entropy- default loss function for binary classification problems. Summarizes the average difference between the actual and predicted probability."
        sar2="hinge- alternative to cross entropy binary classification developed with SVM models used with support vector machine models"
        sar3="mse-default loss to use for regression problems. calculated as the average of the squared differences between the predicted and actual values"
        sar4="mae-for regression problems.  used in cases where there are outliers. average of the absolute difference between actual and predicted values"

        
        sg.PopupScrolled(tutorial,targetprediction, svm, targetprediction1, targetprediction2, targetprediction3, targetprediction4, targetprediction5, spacer,chem, multiclass, chembl, chembl1, chembl2, chembl3, chembl4, chembl5, chembl6, chembl7, chembl8, chembl9, chembl10, chembl11, spacer, my_text, spacer, mytext2, spacer, mytext3, spacer, mytext6, mytext4, mytext5, mytext7, spacer, mytext8, mytext9, mytext10, mytext11, mytext12, mytext13, spacer, mytext14, mytext15, spacer, spreadsheet, spreadsheet1, spacer, solubility, solubility1, solubility2, solubility3, spacer, sar, sar1, sar2, sar3, sar4)
        
        
        
    if event == 'Clear':
        
        window.FindElement('Output').Update('')    
       
    if event == 'Submit':
        
        smiles=values[1]
        f2=open('input.txt', 'w+')
        print(smiles, file=f2)
        f2.close()
        
        value=str(values[2])
        
        
        if value == "[\'PAINS\']":
            
            painspredict('input.txt', 'output.txt')
            f = open('output.txt', 'r+')
            lines = [line for line in f.readlines()]
            f.close()
            for i in lines:
                print(i)
            
        if value == "[\'Fragment Molecules\']": 
                  
            fragmenter('input.txt')
            f = open('output.txt', 'r+')
            lines = [line for line in f.readlines()]
            f.close()
            for i in lines:
                print(i)
                
        if value== "[\'Predict solubility\']":
            mysmiles=values[1] 
           
            solubilitypredictor(mysmiles)
            f = open('output.txt', 'r+')
            lines = [line for line in f.readlines()]
            f.close()
            for i in lines:
                print(i)   
                
        if value == "[\'Toxicity Prediction\']": 
            mysmiles=values[1] 
            mylist=mysmiles.split('\n')
            
           
         
                       
            for k in mylist:
                
                toxtest(k)
                f = open('output.txt', 'r+')
                lines = [line for line in f.readlines()]
                f.close()
                
                for i in lines:
                    
                    print(i)
        
        if value == "[\'Toxicity Map\']":
            environ['smile'] = values[1]
            l=0
             
            for file in glob.glob('*.pkl'):
                if(l <= 46):
                    for i in range(10000):
                        progress_bar.UpdateBar(i+1)
                        output.Update('{:2d}%'.format(int((i+1)/100)))
                        
                    environ['file'] = file
                
               
                    call(["runipy", "toxpredict.ipynb"])
                else:
                    break
                l+=1
                
        progress_bar.UpdateBar(0)         
        if value == "[\'Chembl Prediction\']": 
            mysmiles=values[1]
            mylist=mysmiles.split('\n')
            for k in mylist:
                
                chemblpredict(k)
                f = open('output.txt', 'r+')
                lines = [line for line in f.readlines()]
                f.close()
                for i in lines:
                    print(i)
                                       
        if value== "[\'Build a SAR model\']":
           
            mysmiles=values[1]
            
            sg.PopupOK('Must have smiles then a tab then 1 or 0 for active or inactive respectively.Recommended 100 active, 100 inactive')  
            
                 
            
              
            
            sg.SetOptions(text_justification='right')      
                                                   
            while True:
                
                layout = [[sg.Text('Machine Learning Command Line Parameters', font=('Helvetica', 16))],
                
                [sg.Text('Layers', font=('Helvetica', 15), justification='left')], 
                [sg.Text('relu activation', size=(15, 1)), sg.Text('number of nodes=', size=(18, 1)), sg.In(default_text='236', size=(10,1)), sg.Text('drop out', size=(18, 1)), sg.In(default_text='.2', size=(10,1)),      
                sg.Text('regularizers l=', size=(18, 1)), sg.In(default_text='.1', size=(10,1))],      
                 
                
                [sg.Text('Stochastic Gradient Descent', font=('Helvetica', 15), justification='left')],
                [sg.Text('learning rate', size=(18, 1)), sg.In(default_text='.01', size=(10,1))],
                [sg.Text('momentum', size=(18, 1)), sg.In(default_text='.8', size=(10,1))],      
                [sg.Text('decay', size=(18, 1)), sg.In(default_text='0', size=(10,1))],     
                [sg.Checkbox('nesterov', size=(20, 1), default=True)],      
                  
                   
                [sg.Text('Compile', font=('Helvetica', 15), justification='left')],      
                [sg.Radio('Cross-Entropy', 'loss', size=(12, 1), default=True)],       
                [sg.Radio('Hinge', 'loss', size=(12, 1)),      
                sg.Radio('MAE(L1)', 'loss', size=(12, 1))],      
                [sg.Radio('MSE(L2)', 'loss', size=(12, 1))],
                [sg.Text('number of epochs=', size=(18, 1)), sg.In(default_text='7', size=(10,1))],
                [sg.Text('batch size=', size=(18, 1)), sg.In(default_text='3', size=(10,1))],
                [sg.Text('Filename', size=(15, 1)), sg.In(default_text='test', size=(10, 1))],
                [sg.Text('Model name', size=(15, 1)), sg.In(default_text='my model', size=(10, 1))], 
                [sg.Button(button_text="Back")], [sg.Button(button_text="Submit")]]     
                  

                window3 = sg.Window('Machine Learning Front End', layout, font=("Helvetica", 12))
                
                event, values = window3.Read()
                nodes=values[0]
                dropout=values[1]
                regularizers=values[2]
                lr=values[3]
                momentum=values[4]
                decay=values[5]
                nesterov=values[6]
                crossentropy=values[7]
                hinge=values[8]
                mae=values[9]
                mse=values[10]
                epochs=values[11]
                batchsize=values[12]
                filename=values[13]
                modelname=values[14]
                
              
                
                       
                        
                if nesterov==True:
                    mynesterov="True"
                else:
                    mynesterov="False"
                
                if crossentropy==True:
                    mycrossentropy="True"
                  
                else:
                    mycrossentropy="False"
                
               
                if hinge == True:
                    myhinge="True"
                else:
                    myhinge="True"
                    
                lossfunctions=['binary_crossentropy', 'hinge', 'mean_absolute_error', 'mean_squared_error']
                

                
                mylosses=[crossentropy, hinge, mae, mse]
                
                for j in mylosses:
                    if j == True:
                        theloss=lossfunctions[j]
                    
              
                    
                args=[nodes, dropout, regularizers, lr, momentum, decay, mynesterov, theloss, epochs, batchsize, filename, modelname]   
                
                outF=open("csvfile.csv", "w")
                for line in args:
                    outF.write(line)
                    outF.write("\n")
                outF.close()
                      
                call(["runipy", "gsk.ipynb"])
                outN=open("samplefile.txt", "r")
                for line in outN:
                    print(line)
                outN.close()   
                
                
                
                
               
                
                if event == "Back":
                    
                    window3.Hide()
                    win3_active=False
                    window.UnHide()
                    break
                       
         
        
        if value== "[\'Make Spreadsheet\']":
            
            mysmiles=values[1] 
            
            
            makespreadsheet(mysmiles)
            sg.PopupOK('Your spreadsheet is named test.xlsx')            
                                            
     
window.Close()

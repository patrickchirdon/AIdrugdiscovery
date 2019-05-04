NIH.py	screen compounds for pan assay interference, fraction of sp3 hybridized carbons (>.47 is good), and total polar surface area
Rdockingscript	read rosetta output into r data frames
chirdonmontecarlo.py	my attempt at a monte carlo simulation of a NPT system of 256 molecules in a box
dockcopy.xml	required by rosetta for protein ligand docking
docking.sh	shell script for submitting rosetta docking jobs
fragmenter.py	fragment compounds into parts and recombine them.  useful when you know what binds well, take those compounds and fragment them and recombine to bet better ones.  also good for making new compounds from existing data bases.
moleculardynamicshw2.py molecular dynamics of NVT system of 256 molecules, lennard jones
moleculardynamicshw4.py	molecular dynamics, lennard jones
solubility.csv	required by solubilitypredict.py
solubilitypredict.py  can predict aqueous solubility of compounds in water.


---------------------------------------------------------------

# AIdrugdiscovery
Docking-- https://www.rosettacommons.org/
see wiki

https://www.rcsb.org/
crystal structures for docking

ADME-- http://www.swissadme.ch/index.php
see wiki

SEA search server-- http://sea.bkslab.org/
see wiki

Data Warrior-- http://www.openmolecules.org/datawarrior/
see wiki

SPACE-- http://www.bprc.ac.cn/space/


EPA TEST-- https://www.epa.gov/chemical-research/toxicity-estimation-software-tool-test
see wiki

Emoltox-- http://xundrug.cn/moltox


Chemicalize-- https://chemicalize.com/
Useful for obtaining the IUPAC names of SMILES


Open Babel--
https://openbabel.org/docs/dev/Installation/install.html

necessary for converting SMILES into 3D representations (note must use --gen3D).  Rosetta parameterizes .mol2 files


Also Useful--

https://www.emolecules.com/
library of virtual compounds

https://www.rdkit.org/
useful for computing tanimoto coefficients and molecular descriptors. 
---------------------------------

texts that might be useful-- 
https://www.amazon.com/Statistical-Mechanics-…/…/ref=sr_1_1…

https://www.amazon.com/Understanding-Molecula…/…/ref=sr_1_2…

https://www.amazon.com/Tutorials-Chemoinformatics-Alexandre-Varnek/dp/1119137969/ref=sr_1_1?ie=UTF8&qid=1547747839&sr=8-1&keywords=tutorials+in+chemoinformatics

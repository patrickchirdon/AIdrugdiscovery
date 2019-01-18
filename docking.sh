#!/bin/bash
echo -n "Enter the name of the receptor pdb file:"
read receptor
for files in ligand*.mol2
do
/users/PHS0297/ohu0515/rosetta/main/source/scripts/python/public/molfile_to_params.py -n "Out${files%.mol2}" -p "Out${files%.mol2}" $files
cat $receptor "Out${files%.mol2}_0001.pdb" > "Out${files%.mol2}_wReceptor.pdb"
~/rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease -in:path:database /users/PHS0297/ohu0515/rosetta/main/database/ -in:file:s "Out${files%.mol2}_wReceptor.pdb" -in:file:extra_res_fa "Out${files%.mol2}.params" -scorefile pocket1.sc -out:nstruct 3 -parser:protocol /users/PHS0297/ohu0515/emolecules/003/dockcopy.xml
mv "Out${files%.mol2}_wReceptor_0001.pdb" "${files%.mol2}_docked_0001.pdb"
mv "Out${files%.mol2}_wReceptor_0002.pdb" "${files%.mol2}_docked_0002.pdb"
mv "Out${files%.mol2}_wReceptor_0003.pdb" "${files%.mol2}_docked_0003.pdb"
awk '{print $2 "\t" $50 "\t" $NF}' pocket1.sc > reportcardpocket1.dat
done




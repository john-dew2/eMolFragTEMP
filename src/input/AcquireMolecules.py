from rdkit import Chem
from pathlib import Path
from Molecule import Molecule

import constants

#takes the contents of a file and puts it in a string for processing
def fileToString(file):

    contents = ""
    with open(file) as f:
        contents = f.read()
    f.close()
        
    return contents
   
#converts the contetns of a file to their respective molecule
def convertToRDkit(contents, extension):
    Chem.doKekule = False
    
    if (extension == constants.FASTA_FORMAT_EXT):
        return Chem.MolFromFASTA(contents) 

    if (extension == contants.YAML_FORMAT_EXT):
        return Chem.MolFromHELM(contents) 

    if (extension == contants.MOL2_FORMAT_EXT):
        return Chem.MolFromMol2Block(contents)

    if (extension == contants.MOL_FORMAT_EXT):
        return Chem.MolFromMolBlock(contents)

    if (extension == contants.PDB_FORMAT_EXT):
        return Chem.MolFromPDBBlock(contents) 

    if (extension == contants.SMARTS_FORMAT_EXT):
        return Chem.MolFromSmarts(contents) 

    if (extension == contants.SMILES_FORMAT_EXT):
        return Chem.MolFromSmiles(contents)

    if (extension == contants.TPL_FORMAT_EXT):
        return Chem.MolFromTPLBlock(contents)
        
    return None
    
def acquireMolecules(files):

    #parse each file to convert to RDKit molecule        
    for current_file in files:
        #get the contents of the file and the file type (extension) for processing
        file_contents = fileToSting(current_file)
        extension = current_file.suffix
            
        #get the molecule
        try:
            molecule = convertToRDkit(file_contents, extension)
        except:
            #print(f"RDKit failed to read {current_file.name}")
                
        #if the molecule didnt process let the user know
        if molecule == None:
            print(f'[Error] RDKit failed to convert {current_file.name} to a RDKit object.")

        #otherwise add it to our dataset and update the filenames we have
        else:
            molObject = Molecule(molecule, current_file.name)
            data.append(molObject)
            
    return data
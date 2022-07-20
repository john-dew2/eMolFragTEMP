from rdkit import Chem
from pathlib import Path
from eMolFragTEMP.src.representation import Molecule
from eMolFragTEMP.src.utilities import constants
from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.utilities import logging

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

    if (extension == constants.YAML_FORMAT_EXT):
        return Chem.MolFromHELM(contents) 

    if (extension == constants.MOL2_FORMAT_EXT):
        return Chem.MolFromMol2Block(contents)#, False)

    if (extension == constants.MOL_FORMAT_EXT):
        return Chem.MolFromMolBlock(contents)

    if (extension == constants.PDB_FORMAT_EXT):
        return Chem.MolFromPDBBlock(contents) 

    if (extension == constants.SMARTS_FORMAT_EXT):
        return Chem.MolFromSmarts(contents) 

    if (extension == constants.SMILES_FORMAT_EXT):
        return Chem.MolFromSmiles(contents)

    if (extension == constants.TPL_FORMAT_EXT):
        return Chem.MolFromTPLBlock(contents)
        
    return None
    
def acquireMolecules(files):
    data = []

    #parse each file to convert to RDKit molecule        
    for current_file in files:
      #get the contents of the file and the file type (extension) for processing
      file_contents = fileToString(current_file)
      extension = current_file.suffix
          
      #get the molecule
      try:
          molecule = convertToRDkit(file_contents, extension)
      except:
          logging.logger.error(f"RDKit failed to read {current_file.name}")

      #add it to our dataset and update the filenames we have
      else:
          molObject = Molecule.Molecule(molecule, current_file.name)
          data.append(molObject)
            
    return data
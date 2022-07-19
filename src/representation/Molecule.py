#That molecule class will contain the rdkit object, the name of the file it came from, as well as a list of 'equal other fragments'.
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from eMolFrag2.src.utilities import constants

class Molecule:
    def __init__(self, rdkit_object, file_name):
        self.rdkitObject = rdkit_object
        self.filename = file_name
        self.parent = None
        self.equalFragments = []
        
    def getParent(self):
        return self.parent
        
    def getRDKitObject(self):
        return self.rdkitObject
    
    def getFileName(self):
        return self.filename
    
    def getEqualFragments(self):
        return self.equalFragments
        
    def setEqualFragments(self, listOfFragments):
        self.equalFragments = listOfFragments

    def clearProperties(self):
        """
           Clean the rdkit molecule of all 'public' properties and 'private'
           Tripos ChargeType information

            @input: Rdkit.Mol
            @output: None
        """
        properties = self.rdkitObject.GetPropNames()
        properties.append(constants.ATOMTYPE_PROP) # Remove an errant property
        for property in properties:
            self.rdkitObject.ClearProp(property)

    def __eq__(self, molecule):
        # Acquire 
        fp1 = AllChem.GetMorganFingerprintAsBitVect(self.rdkitObject, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(molecule.rdkitObject, 3, nBits=2048)
        tc = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        
        return tc == 1
    
    def __str__(self):
        numAtoms = self.rdkitObject.GetNumAtoms()
        numBonds = self.rdkitObject.GetNumAtoms()
        
        return f"{self.filename} has {numAtoms} atoms and {numBonds} bonds"
        
        #do __eq__
        #do __str__ that meaningfully calls the rdkit molecule toString()
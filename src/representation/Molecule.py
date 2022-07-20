#
# That molecule class will contain the rdkit object, the name of the file it came from, as well as a list of 'equal other fragments'.
#
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from eMolFragTEMP.src.utilities import constants
from eMolFragTEMP.src.utilities import tc

class Molecule:
    def __init__(self, rdkit_mol,  file_name = None, parentMol = None):
        """
            @input: Rdkit.Mol object for this molecule
            @input: file_name -- file which this molecule originated
                       * In the case where our molecule is a fragment
                         (Brick / Linker) this will be None
            @input: parent -- Molecule that was fragmented to acquire this object
                       * In the case where our molecule is direct from a file, 
                         this will be None
        """
        self.rdkitObject = rdkit_mol
        self.filename = file_name
        self.parent = parentMol
        self.similar = []
        
    def getParent(self):
        return self.parent
        
    def getRDKitObject(self):
        return self.rdkitObject
    
    def getFileName(self):
        return self.filename
        
    def addTCSimilar(self, mol):
        """
            @input: mol -- a Molecule object that is checked for TC equivalence externally
        """
        return self.similar.append(mol)
    
    def clearProperties(self):
        """
           Clean the rdkit molecule of all 'public' properties and 'private'
           Tripos ChargeType information

            @input: Rdkit.Mol
            @output: None
        """
        properties = self.rdkitObject.GetPropNames()
        properties.append(constants.ATOMTYPE_PROP) # Remove an errant property
        properties.append(constants.ATOMTYPE_CHARGE_PROP) # Remove an errant property
        for property in properties:
            self.rdkitObject.ClearProp(property)

    def __hash__(self):       
        # hash(custom_object)
        return self.rdkitObject.__hash__()

    def __eq__(self, molecule):
        return tc.TCEquiv(self, molecule)
    
    def __str__(self):
        numAtoms = self.rdkitObject.GetNumAtoms()
        numBonds = self.rdkitObject.GetNumAtoms()
        
        return f"{self.filename} has {numAtoms} atoms and {numBonds} bonds"
        
    def makeFragmentFileName(self, \
                             file_name, \
                             prefix = "", \
                             numeric_suffix = 0, \
                             extension = constants.SDF_FORMAT_EXT):
        return f'{prefix}-{file_name}-{str(numeric_suffix).zfill(3)}{extension}'
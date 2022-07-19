import math 
from rdkit import DataStructs # For TC Computations

import sys
from eMolFrag2.src.utilities import constants 
from eMolFrag2.src.utilities import tc

from eMolFrag2.src.representation.Molecule import Molecule

#
# This class will simulate an equivalence class of molecules
#
#    We have a dictionary of the form <Molecule, [TC-'Equivalent' Molecules]>
#
class MoleculeDatabase(Molecule):

    def __init__(self, given_tc = constants.DEFAULT_TC_UNIQUENESS):
        self.database = {}
        
        if given_tc < 0 or given_tc > 1:
            print(f'Tanimoto coefficient constant {given_tc} is not in allowable range 0 <= tc <= 1.0')
            raise RuntimeError

        self.TC_THRESH = given_tc

    #
    # @input: 2 Molecule objects
    # @output: True if the molecules are TC-equivalent; False otherwise
    #          (we use math.isclose to assess floating-point-based equivalent)
    #
    def _TCEquiv(self, mol1, mol2):
        
        tanimoto = tc.TC(mol1, mol2)  
        # >   
        if (tanimoto > self.TC_THRESH):
          return True
        # =                                         # 4-decimal place equality 
        return math.isclose(tanimoto, self.TC_THRESH, rel_tol = 1e-5)

    
    #
    # Add one Molecule object to the database
    #
    # @output: if the given molecule is unique, return True (False if it is TC-redudant
    #
    def add(self, molecule):
      
        tc_equiv = list(filter(lambda db_mol : self._TCEquiv(molecule, db_mol), self.database.keys()))

        if len(tc_equiv) > 1:
            print(f'Internal MoleculeDatabase error; {len(tc_equiv)}-TC equivalent molecules')            

        # Empty ; we have a unique fragment; a new entry has { molecule, [] }
        if not tc_equiv:
            self.database[molecule] = []
            return True
            
        # Not a TC-unique fragment: { tc_equiv, [..., molecule] }
        self.database[tc_equiv[0]].append(molecule)
        return False
    
    #
    # Add a list of molecules
    # @output: the list of UNIQUE database additions
    #          (redundant count can be computed by user)
    #
    def addAll(self, molecules):
        return [mol for mol in molecules if self.add(mol)]
    
    def GetUniqueMolecules(self):
        return self.database.keys()

    #
    # Return all Molecule objects stored
    #
    def GetAllMolecules(self):
        all_mols = []
        for mol, tc_mols in self.database.items():
            all_mols.append(mol)
            all_mols += tc_mols

        return all_mols
        
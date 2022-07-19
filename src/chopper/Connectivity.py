from rdkit import Chem
from eMolFrag2.src.utilities import constants
from eMolFrag2.src.utilities import logging
import BRICS_custom

def getAtomTypeMap(mol):
  """ 
  Get the TriposAtomTypes of the molecule as map<atom_index, AtomType> 
  If input format molecule does not explicitly define AtomType (using
  Atom::GetProp(_TriposAtomType) we use Atom::GetSymbol as fallback.

      Parameters: 
                mol (Rdkit.Mol): Molecule (preferably from mol2 for TriposAtomType)

      Returns: 
                m (matrix): Molecule Adjacency Matrix
  """
  try:
    atom_atomType_map = [(id, atom.GetProp('_TriposAtomType')) for id, atom in enumerate(mol.GetAtoms())]
  except KeyError as e:
    logger.warning("Error: _TriposAtomType KeyError: Check input molecule type (.mol2 preferred).")
    atom_atomType_map = [(id, atom.GetSymbol()) for id, atom in enumerate(mol.GetAtoms())]

  return atom_atomType_map



def getBRICSBonds(snips, mol):
  """
      BRICS functionality requires bonds along with isotope information
          (bond, bond-isotopes) = ((int, int), (str, str))

      For our subset of BRICS bonds, acquire these bonds.

      Parameters:
                snips (set): Set of final BRICS snip list
                mol (Rdkit.Mol): Original Molecule

      Returns:
                bonds (set): Set of bonds formated for BreakBRICSBonds
  """
  return [(a, b) for a, b in list(Chem.BRICS.FindBRICSBonds(mol)) if tuple(sorted(a)) in snips]

def compute(parent_mol, bricks, linkers, snips):

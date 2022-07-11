from rdkit import Chem
import constants
    
#
# Return the number of atoms that are
#    (1) dummy atoms or
#    (2) radicals or
#    (3) hydrogens
#
# @input: rdkit fragment object
# @output: True if it meets Brick criteria (False indicates a Linker)
#
def CountAtomsToBeRemoved(rdkit_fragment):

    count = 0
    for atom in m.GetAtoms():

        sym = atom.GetSymbol()

        if sym == constants.DUMMY_ATOM_STR or \
           sym == constants.RADICAL_ATOM_STR or \
           sym == constants.HYDROGEN_ATOM_STR:
            count = count + 1

    return count

#
# Brick criteria:
#    Number of 'meaningful' atoms is 4 or more
#
# @input: rdkit fragment object
# @output: True if it meets Brick criteria (False indicates a Linker)
#
def isBrick(rdkit_frag):
    return rdkit_frag.GetNumAtoms() - CountAtomsToBeRemoved(rdkit_frag) >= constants.BRICK_MINIMUM_NUM_ATOMS


#
# @input: list of rdkit molecule fragment objects
# @output: a tuple of two lists (bricks, linkers) both lists are rdkit objects
#
def partition(rdkit_fragments):

    bricks = []
    linkers = []
    for rdkit_frag in rdkit_fragments:

        if isBrick(rdkit_frag):
            bricks.append(rdkit_frag)

        else:
            linkers.append(rdkit_frag)
            
    return bricks, linkers
from rdkit import Chem
import constants

#
# After chopping fragments, identify if any linker-linker neighbors have been chopped; reconnect them.
# i.e., if two linkers used to connect to each other, then connect them again to get larger linkers.
#
def combineLinkers(parent_rdkit_mol, rdkit_linker_list):


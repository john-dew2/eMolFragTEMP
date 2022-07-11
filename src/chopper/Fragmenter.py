from rdkit.Chem import BRICS
from rdkit import Chem

def BRICS_Fragment(molecule):

    fragmented_mol = Chem.FragmentOnBRICSBonds(molecule)

    return Chem.BRICS.GetMolFrags(fragemented_mol, asMols = True, sanitizeFrags = False)

#
# @input:  Molecule object
# @output: List of rdkit objects (fragments)
#
def fragment(mol):

    #
    # Fragment with the BRICS algorithm 
    #
    fragments = []
    try:
        fragments = BRICS_Fragment(molecule.getRDKitObject())
    except:
        error = f'molecule {molecule.GetFileName} failed to fragment with BRICS'
        print(error)

    return fragments
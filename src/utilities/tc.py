from rdkit import Chem
from rdkit import DataStructs # For TC Computations

#
# A utility function to compute the Tanimoto Coefficient
#
def TC(rdkit_mol1, rdkit_mol2)
    fp_1 = Chem.RDKFingerprint(rdkit_mol1)
    fp_2 = Chem.RDKFingerprint(rdkit_mol2)

    return = DataStructs.FingerprintSimilarity(fp_1, f_2)

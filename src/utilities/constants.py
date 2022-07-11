#Constants.py

#
# INPUT CONSTANTS
#

# RDKit can use mol2, mol, smiles, smarts, FASTA, HELM, PDB, PNG
FASTA_FORMAT_EXT = ".fasta"
YAML_FORMAT_EXT = ".yaml"
MOL2_FORMAT_EXT = ".mol2"
MOL_FORMAT_EXT = ".mol"
PDB_FORMAT_EXT = ".pdb"
SMARTS_FORMAT_EXT = ".sma"
SMILES_FORMAT_EXT = ".smi"
TPL_FORMAT_EXT = ".tpl"

ACCEPTED_FORMATS = [FASTA_FORMAT_EXT,\
                    YAML_FORMAT_EXT,\
                    MOL2_FORMAT_EXT,\
                    MOL_FORMAT_EXT,\
                    PDB_FORMAT_EXT,\
                    SMARTS_FORMAT_EXT,\
                    SMILES_FORMAT_EXT,\
                    TPL_FORMAT_EXT]

#
# ATOMIC CONSTANTS
#
HYGROGEN_ATOM_STR = "H"
CARBON_ATOM_STR = "C"
OXYGEN_ATOM_STR = "O"

RADICAL_ATOM_STR = "R"
DUMMY_ATOM_STR = "*"

#
# Linker / Brick CONSTANTS
#
LINKER_MAXIMUM_NUM_ATOMS = 3
BRICK_MINIMUM_NUM_ATOMS = LINKER_MAXIMUM_NUM_ATOMS + 1

DEFAULT_TC_UNIQUENESS = 1.0
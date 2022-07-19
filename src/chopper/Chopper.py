from rdkit import Chem
from eMolFrag2.src.utilities import logging
import Deconstructor
import Connectivity

def chop(rdkit_mol):
    """
        Chopping consists of the following algorithm:
            1. Build graph of molecule
            2. Find where BRICS would cleave (BRICS bonds) [snips in the code]
                 * We modify FindBRICSBonds code to ensure that L7 bonds are never
                   cleaved (according eMolFrag v1.0)
                 * Modified in BRICS_custom.py
            3. Break graph into fragments
                 Take the graph of molecule, substract proposed BRICS bonds
                 The resulting molecule is a disconnected set of subgraphs.
                 Each of these subgraphs is a fragment.
                 Those fragments with fewer than 4 atoms are the first phase of
                 linkers; all others are bricks
            4. Combine sequences of Linker-Linker fragments into single a Linker

            5. Sequentially BreakBRICSBonds over each snip
                 The resulting fragment is known to be either brick or linker.
                 When we break the bonds with BRICS, free radicals are added.

            6. Compute connectivy of freeradicals (from each snip calc atomtype of tuple pairs)

            7. Link atomtypes to atoms in our fragment
            
            
            Input: An rdkit molecule (ideally, with AtomType info from mol2 format)
            
            Output:
    """
    # Remove all hydrogens from our molecule for simplicity
    stripped_mol = parent_mol.GetRDKitObject().RemoveAllHs()

    #
    # Steps (1) - (5)
    # Deconstruct our molecule into the known fragments
    # These are sets of atom indices (bricks, linkers); snips are a set of bonds.
    #
    bricks, linkers, snips = Deconstructor.deconstruct(stripped_mol)

    #
    # Steps (6) - (7)
    # Compute connectivity and ensure atomtypes information is communicated to fragments
    #
    brick_db, linker_db = Connectivity.compute(rdkit_mol, bricks, linkers, snips)
    
    return brick_db, linker_db

#
# Create and add local Brick/Linker objects to the database; report results
#
def processFragments(bricks, linkers, parent):

    results = database.addAll([Brick(fragment, parent) for fragment in bricks])
    
    logging.logger.debug(f'Added {result.count(True)} TC-unique bricks; \
                         {result.count(False)} were TC-redundant')

    results = database.addAll([Linker(fragment, parent) for fragment in linkers])
    
    logging.logger.debug(f'Added {result.count(True)} TC-unique linkers; \
                         {result.count(False)} were TC-redundant')

#
# Chop many molecules
#
# @input: list of Molecule objects (containing rdkit objects)
# @output: MoleculeDatabase containing all fragments
#
def chopall(molList):

    database = MoleculeDatabase()

    for parent_mol in molList:

        logging.logger.debug(f'Processing molecule{parent_mol.GetFileName()}')

        # Actually chop and store the molecular fragments in the database
        bricks, linkers = chop(stripped_rdkit)

        processFragments(bricks, linkers, parent_mol)      

    return database
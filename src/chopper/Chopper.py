from rdkit import Chem
import utilities
import constants
import DoubleBondRecoonector
import Fragmenter
import FragmentPartitioner    

#
# Main chop routine:
#     (1) Fragment with BRICS
#     (2) Reconnect C.2 = C.2 double bonds broken by BRICS
#     (3) Identify bricks and linkers
#     (4) Reconnect adjacent linkers that share a summy bond
#     (5) Create internal representation of the fragments for downstream processing
#
# @input: A Molecule object (which contains an rdkit molecule object)
# @output: list of (local) Bricks and Linker objects
#
def chop(molecule):

    #
    # (1) Fragment
    #
    try: 
        rdkit_fragments = fragmenter.fragment(molecule.GetRDKitObject())

    except RDKitError:
        return []

    #
    # (2) Reconnect C.2 = C.2 Double Bonds
    #
    fragCountBefore = len(rdkit_fragments)
    rdkit_fragments = DoubleBondRecoonector.reconnect(molecule.GetRDKitObject(),
                                                      rdkit_fragments)

    if len(rdkit_fragments) != fragCountBefore:
        print(f'Reconnected {fragCountBefore - len(rdkit_fragments)} C.2-C.2 double bonds.')

    #
    # (3) Distribute the fragments into bricks / linkers
    #
    bricks, linkers = FragmentPartitioner.partition(fragments)

    #
    # (4) 'Adjacent' linkers will be combined into single linkers
    #
    linkers = LinkerCombiner.combineLinkers(molecule.GetRDKitObject(), bricks, linkers)

    #
    # (5) Create 'local' molecule representations
    #
    localBricks = [Brick(fragment, parent) for fragment in bricks]
    localLinkers = [Linker(fragment, parent) for fragment in linkers]
    
    return localBricks + localLinkers

#
# Chop many molecules
#
# @input: list of Molecule objects (containing rdkit objects)
# @output: MoleculeDatabase containing all fragments
#
def chopall(molList, level):

    database = MoleculeDatabase()
    for mol in molList:

        emit(level, f'Processing molecule{mol.GetFileName()}')

        result = database.addAll(chop(mol))
        
        emit(level, f'Added {result.count(True)} TC-unique fragments; \
             {result.count(False)} were TC-redundant')

    return database
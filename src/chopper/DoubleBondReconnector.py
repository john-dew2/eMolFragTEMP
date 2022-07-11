from itertools import filterfalse, tee
   
from rdkit import Chem
import fragmenter
import constants



#
# Does the bond consist of a carbon (C) and dummy created by BRICS (*)?
#
# @input: two RDKit Atom objects
#
def isCarbonDummyBond(rdkit_atom_begin, rdkit_atom_end):
    syms = [rdkit_atom_begin.GetSymbol(), rdkit_atom_end.GetSymbol()]

    return constants.CARBON_ATOM_STR in syms and
           constants.DUMMY_ATOM_STR in syms

#
# Does this rdkit molecule (fragment) have a C.2 = C.2 double bond
#
def isC2C2DoubleBond(rdkit_bond):

   if rdkit_bond != Chem.rdchem.BondType.DOUBLE return False
   
   return isCarbonDummyBond(rdkit_bond.GetBeginAtom(), \
                            rdkit_bond.GetEndAtom())

#
#  Check a single fragment for our desired bond-type
#
def hasC2C2DoubleBond(rdkit_frag):
   return next((True for bond in rdkit_frag.GetBonds() if isC2C2DoubleBond(bond)), False)

#
# Partition the fragments into two sublists;
#    (those with C.2 = C.2, and those without)
#
def partitionFragments(rdkit_frags):
   not_db = []
   db = []

   for frag in rdkit_frags:
       (db if hasC2C2DoubleBond(frag) else not_db).append(frag)

   return not_db, db
  
#
# @param: Parent molecule that has been fragmented
# @param: Set of fragments that contain C.2 = C.2 double bonds
# @output: A reconnected list of rdkit fragment objects
#
def reconnect_db(rdkit_parent, rdkit_db_fragments):
    
  
#
# The goal is to reconnect C.2 = C.2 double bonds that are broken by BRICS
#
def reconnect(parent_rdkit_mol, rdkit_fragment_list):

    # Partition in those having desired double-bonds (db) and those not (not_db)
    frags_not_db, frags_db = partitionFragments(rdkit_fragment_list)

    # No double-bonds to reconnect
    if not frags_not_db:
        return frags_db

    if len(frags_db) % 2 == 1:
        utilities.emitError(f'Odd number ({len(frags_db)}) of double-bonded fragments found.'")

    # Perform reconnection
    reconnected = reconnect_db(parent_rdkit_mol, frags_db)
    
    if len(reconnected) != len(frags_db) / 2:
        utilities.emitError(f'Expected to reconnect an even number fragments: \
                            {len(frags_db)} into {len(frags_db) / 2}')

    # Return the non-double-bonded fragments along with reconnected     
    return frags_not_db + reconnected
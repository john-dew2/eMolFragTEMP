from Molecule import Molecule
from eMolFrag2.src.utilities import constants

class Brick(Molecule):
    def __init__(self, rdkit_object, file_name):
        Molecule.__init__(self, rdkit_object, file_name)
        
    def getSDF():
        """
            Assuming a molecule with AtomType and connectivity
            information preseversed, output the molecule in SDF format
            with proper appendices

            @output: string in SDF format
        """
        clearProperties(mol)

        # Write to a string
        from rdkit.six import StringIO
        sio = StringIO()
        writer = Chem.SDWriter(sio)

        # ATOMTYPES appendix; take directly from Atom.Prop information
        mol.SetProp(SDF_OUTPUT_BRICK_ATOMTYPES, '\n'.join([atom.GetProp(ATOMTYPE_PROP) for atom in mol.GetAtoms()]))

        # Break connection string in separate cases...and join all together
        # Multiple connections may be stored as a comma-separated string for an atom
        # e.g., (1, C.3 C.3 N.4) -->
        #                2 C.3
        #                2 C.3
        #                2 N.4        
        atom_conn_map = {index : atom.GetProp(ATOM_CONNECTION_PROP) for index, atom in enumerate(mol.GetAtoms()) if atom.HasProp(ATOM_CONNECTION_PROP)}
        
        appendix = ''.join(''.join(f'{atom_index + 1} {atomtype}\n' for atomtype in atom_conn_map[atom_index].split(' ')) for atom_index in atom_conn_map)
        
        mol.SetProp(SDF_OUTPUT_BRICK_CONNECTIONS, appendix)

        # TODO: add similar fragments appendix

        writer.write(mol)
        writer.close()

        return sio.getvalue()
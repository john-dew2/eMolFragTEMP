from Molecule import Molecule

class Linker(Molecule):

    def __init__(self, rdkit_object, file_name):
        Molecule.__init__(self, rdkit_object, file_name):

    def __init__(self, rdkit_object, rdkit_parent):
        Molecule.__init__(self, rdkit_object, rdkit_parent):
        
    def getLinkerSDF():
        """
            Assuming a molecule with AtomType and connectivity
            information preseversed, output the molecule in SDF format
            with proper appendices

            @output: string in SDF format
        """
        clearProperties(self.rdkitObject)

        # Write to a string
        from rdkit.six import StringIO
        sio = StringIO()
        writer = Chem.SDWriter(sio)

        # A linker maintains a "maximum number of connections for each fragment"
        # All atoms must have an atomtype; not all atoms have connections
        appendix = '\n'.join(f'{str(len(atom.GetProp(ATOM_CONNECTION_PROP).split()) if atom.HasProp(ATOM_CONNECTION_PROP) else 0)} {atom.GetProp(ATOMTYPE_PROP)}'\
                             for atom in self.rdkitObject.GetAtoms())

        self.rdkitObject.SetProp(SDF_OUTPUT_LINKER_CONNECTIONS, appendix)

        # TODO: add similar fragments appendix

        writer.write(self.rdkitObject)
        writer.close()

        return sio.getvalue()
    
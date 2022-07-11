from Molecule import Molecule

class Brick(Molecule):
    def __init__(self, rdkit_object, file_name):
        Molecule.__init__(self, rdkit_object, file_name):
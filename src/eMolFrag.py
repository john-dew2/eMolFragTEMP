
import sys
from eMolFragTEMP.src.representation import Molecule
from pathlib import Path
from rdkit import Chem
from eMolFragTEMP.src.input import AcquireFiles, AcquireMolecules, Configuration, Options
#from eMolFragTEMP.src.input.Options import Options

def main():
    dataset = []
    
	#Verify Tools and Parse Command Line
    args = sys.argv
    initializer = Options.Options()
    Configuration.readConfigurationInput(initializer, args)
    
    #Input System
    files = AcquireFiles.acquireMoleculeFiles(initializer)
    dataset = AcquireMolecules.acquireMolecules(files)
    print(len(dataset))
    #Process
    
    #Post-Process
	
	





if __name__ == "__main__":
	main()
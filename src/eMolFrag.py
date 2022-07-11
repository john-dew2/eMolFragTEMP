
import sys
from representation import Molecule
from pathlib import Path
from rdkit import Chem
from input import AcquireFiles
from input import AcquireMolecules
from input import Configuration
from input.Options import Options

def main():
    dataset = []
    
	#Verify Tools and Parse Command Line
    args = sys.argv
    initializer = Options()
    Configuration.readConfigurationInput(initializer, args)
    
    #Input System
    files = AcquireFiles.acquireMoleculeFiles(initializer)
    dataset = AcquireMolecules.acquireMolecules(files)
    #Process
    
    #Post-Process
	
	





if __name__ == "__main__":
	main()

import sys
from Molecule import Molecule
from pathlib import Path
from rdkit import Chem
from input import AcquireFiles
from input import AcquireMolecules
from input import Configuration

def main():
    dataset = []
    
	#Verify Tools and Parse Command Line
    args = sys.argv
    Configuration.readConfigurationInput(args)
    
    #Input System
    files = AcquireFiles.acquireFiles(initializer)
    dataset = AcquireMolecules.acquireMolecules(files)
    
    #Process
    
    #Post-Process
	
	





if __name__ == "__main__":
	main()